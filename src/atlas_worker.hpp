#pragma once
/**
 * atlas_worker.hpp — ATLAS-GP basin-hopping worker
 *
 * Innovations:
 *   1. Stress-guided perturbation    — per-atom energy stress, CDF sampling
 *   2. Graph topology moves          — surface relocation, coordination swap
 *   3. Motif library                 — thread-safe, species-labelled fragments
 *   4. Taboo memory                  — full fingerprint, adaptive lambda
 *   5. GP surrogate                  — Matérn 5/2, shared across runs
 *   6. Acquisition-driven selection  — LCB scores N_cand moves, picks best
 *   7. Homotop search                — species-label swaps at fixed geometry
 *   8. Metropolis acceptance         — escapes local minima probabilistically
 *   9. Parallel tempering exchange   — hot/cold run communication
 *  10. Curvature-guided escape       — finite-difference Hessian soft modes
 */

#include "system.hpp"
#include "potential/potential_base.hpp"
#include "lbfgs_solver.hpp"
#include "fingerprint.hpp"
#include "gp_surrogate.hpp"
#include "acquisition.hpp"
#include "escape.hpp"
#include "active_learning.hpp"
#include "population.hpp"
#include "adaptive_descriptor.hpp"
#include "natural_gradient.hpp"
#include <vector>
#include <deque>
#include <algorithm>
#include <numeric>
#include <atomic>
#include <mutex>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <set>

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────
struct ATLASConfig {
    // Parallelism
    int    n_runs            = 8;
    int    n_workers         = 4;
    int    n_generations     = 300;

    // Local minimization
    int    local_maxiter     = 500;
    int    init_attempts     = 15;

    // Stress perturbation
    double stress_temp_init  = 1.0;
    double stress_temp_min   = 0.05;
    double stress_temp_decay = 0.990;
    double stress_n_frac     = 0.35;
    double stress_step_frac  = 0.25;
    double stress_step_decay = 0.998;

    // Graph moves
    double bond_cutoff       = 1.35;   // ×mean_sigma
    double graph_move_prob   = 0.30;
    double rotation_max_deg  = 60.0;

    // Motif library
    double motif_core_frac      = 0.40;
    int    motif_inject_interval= 25;
    int    motif_lib_maxsize    = 60;
    double motif_min_rmsd       = 0.15;

    // Taboo memory
    double taboo_lambda      = 0.05;
    double taboo_delta       = 0.25;
    int    taboo_max_size    = 200;
    double taboo_thresh      = 0.04;

    // Basin-hopping Metropolis
    double bh_temp           = 1.0;    // in reduced units
    double bh_temp_decay     = 0.9998;
    double bh_temp_min       = 0.20;

    // GP surrogate
    bool   use_gp            = true;
    double gp_kappa          = 2.0;
    int    gp_n_candidates   = 5;

    // Homotop search (heteroatomic)
    bool   homotop_search    = true;
    int    homotop_steps     = 50;

    // Isomer filtering
    double isomer_rmsd_thresh = 0.20;
    int    max_isomers        = 30;
    double energy_window      = 0.40;

    // Output
    std::string output_dir   = "atlas_results";
    bool        progress     = false;
    bool        verbose      = false;
    int         run_id       = 0;
    uint64_t    seed         = 0;

    // Session 3: parallel tempering
    bool   use_pt            = false;   // enable parallel tempering
    int    pt_n_replicas     = 4;       // replicas per run-group
    double pt_T_min          = 0.3;
    double pt_T_max          = 5.0;
    int    pt_exchange_interval = 20;

    // Session 3: eigenvector escape
    bool   use_eigenvec_escape  = true;    // Lanczos/dimer instead of diagonal
    int    lanczos_steps        = 20;

    // Session 3: active learning
    bool   use_active_learning  = false;
    int    oracle_budget        = 50;
    double oracle_uncertainty   = 2.0;

    // Session 4: population crossover
    bool   use_population       = true;   // evolutionary layer
    int    pop_size             = 20;
    int    pop_interval         = 15;     // crossover every N gens
    int    pop_n_children       = 3;      // children per crossover event
    double pop_diversity        = 0.25;   // min RMSD between pop members

    // Session 4: adaptive descriptor (ARD)
    bool   use_ard              = true;   // auto relevance determination
    int    ard_update_interval  = 80;

    // Session 4: natural gradient hyperparameter updates
    bool   use_ngd              = true;   // natural gradient vs grid search
};

// ─────────────────────────────────────────────────────────────────────────────
// Result
// ─────────────────────────────────────────────────────────────────────────────
struct RunResult {
    double best_energy = 1e30;
    System best_sys;
    std::vector<std::pair<double,System>> candidates;
    std::vector<double> energy_history;
    int  n_new_basins    = 0;
    int  moves_stress    = 0;
    int  moves_graph     = 0;
    int  moves_motif     = 0;
    int  moves_homotop   = 0;
    int  moves_gp_guided = 0;
    int  moves_gradient  = 0;
    int  moves_eigenvec  = 0;
    int  oracle_calls    = 0;
    int  moves_crossover = 0;   // Session 4: crossover children
    int  ard_updates     = 0;   // Session 4: ARD weight updates
    double revisit_rate  = 0.0;
};

// ─────────────────────────────────────────────────────────────────────────────
// Thread-safe RNG
// ─────────────────────────────────────────────────────────────────────────────
struct RNG {
    std::mt19937_64 eng;
    std::uniform_real_distribution<double> U{0.0,1.0};
    std::normal_distribution<double>       G{0.0,1.0};
    explicit RNG(uint64_t seed) : eng(seed) {}
    double rand()   { return U(eng); }
    double randn()  { return G(eng); }
    int randi(int lo, int hi) {
        if (lo >= hi) return lo;
        return std::uniform_int_distribution<int>(lo, hi-1)(eng);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Taboo memory — adaptive lambda, full fingerprint
// ─────────────────────────────────────────────────────────────────────────────
struct TabooMemory {
    double lam, delta, thresh;
    int    max_size;
    std::deque<Fingerprint> fps;
    int revisits=0, checks=0;

    TabooMemory(double l, double d, int ms, double th)
        : lam(l), delta(d), thresh(th), max_size(ms) {}

    void add(const Fingerprint& fp) {
        fps.push_back(fp);
        if ((int)fps.size() > max_size) fps.pop_front();
    }

    bool is_known(const Fingerprint& fp) const {
        for (auto& s : fps)
            if (fp.distance(s) < thresh) return true;
        return false;
    }

    double bias(const Fingerprint& fp) {
        ++checks;
        double b = 0.0, min_d = 1e30;
        for (auto& s : fps) {
            double d = fp.distance(s);
            min_d = std::min(min_d, d);
            double d2 = d*d / (2.0*delta*delta + 1e-20);
            b += std::exp(-d2);
        }
        if (min_d < thresh) ++revisits;
        return lam * b;
    }

    void adapt() {
        if (checks < 20) return;
        double r = (double)revisits / (checks + 1e-9);
        if      (r > 0.20) lam = std::min(lam * 1.05, 0.20);
        else if (r < 0.08) lam = std::max(lam * 0.95, 0.001);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Shared motif library (thread-safe, species-labelled)
// ─────────────────────────────────────────────────────────────────────────────
struct MotifLibrary {
    struct Motif { System fragment; double energy; };
    std::vector<Motif> lib;
    mutable std::mutex mtx;
    int max_size; double min_rmsd;

    MotifLibrary(int ms, double mr) : max_size(ms), min_rmsd(mr) {}

    void deposit(const System& sys, double core_frac, double cutoff_factor) {
        int M = std::max(2, (int)std::ceil(core_frac * sys.N));
        double cut2 = (cutoff_factor*sys.mean_sigma)*(cutoff_factor*sys.mean_sigma);

        // Coordination numbers
        std::vector<int> cn(sys.N, 0);
        for (int i=0;i<sys.N-1;++i) for (int j=i+1;j<sys.N;++j)
            if (sys.rij2(i,j) < cut2) { ++cn[i]; ++cn[j]; }

        std::vector<int> idx(sys.N); std::iota(idx.begin(),idx.end(),0);
        std::sort(idx.begin(),idx.end(),[&](int a,int b){return cn[a]>cn[b];});
        idx.resize(M);

        // Build fragment preserving species labels
        System frag;
        frag.is_molecular=false; frag.N=M;
        frag.species=sys.species; frag.mean_sigma=sys.mean_sigma; frag.R=sys.R;
        frag.pos.resize(3*M); frag.type_id.resize(M); frag.charges.resize(M,0);
        double cx=0,cy=0,cz=0;
        for (int k=0;k<M;++k) { cx+=sys.px(idx[k]); cy+=sys.py(idx[k]); cz+=sys.pz(idx[k]); }
        cx/=M; cy/=M; cz/=M;
        for (int k=0;k<M;++k) {
            frag.pos[3*k]  =sys.px(idx[k])-cx;
            frag.pos[3*k+1]=sys.py(idx[k])-cy;
            frag.pos[3*k+2]=sys.pz(idx[k])-cz;
            frag.type_id[k]=sys.type_id[idx[k]];
        }

        std::lock_guard<std::mutex> lk(mtx);
        for (auto& ex : lib) {
            if (ex.fragment.N != M) continue;
            if (kabsch_rmsd_species(frag, ex.fragment, sys.mean_sigma) < min_rmsd) return;
        }
        if ((int)lib.size() >= max_size) lib.erase(lib.begin());
        lib.push_back({frag, 0.0});
    }

    bool random_motif(Motif& out, RNG& rng) const {
        std::lock_guard<std::mutex> lk(mtx);
        if (lib.empty()) return false;
        out = lib[rng.randi(0,(int)lib.size())];
        return true;
    }

    int size() const { std::lock_guard<std::mutex> lk(mtx); return (int)lib.size(); }
};

// ─────────────────────────────────────────────────────────────────────────────
// Smart initialisation — species-aware random packing
// ─────────────────────────────────────────────────────────────────────────────
static System smart_init(const System& tmpl, const Potential& pot,
    int n_attempts, RNG& rng)
{
    const double min_dist = 0.82 * tmpl.mean_sigma;
    System best = tmpl; best.pos.assign(3*tmpl.N, 0.0);
    double best_e = 1e30;

    for (int att=0; att<n_attempts; ++att) {
        System s = tmpl; s.pos.assign(3*tmpl.N, 0.0);
        int placed = 0;
        for (int iter=0; iter<80000 && placed<tmpl.N; ++iter) {
            double x = (rng.rand()*2-1)*tmpl.R;
            double y = (rng.rand()*2-1)*tmpl.R;
            double z = (rng.rand()*2-1)*tmpl.R;
            if (x*x+y*y+z*z > tmpl.R*tmpl.R) continue;
            bool ok = true;
            for (int k=0;k<placed&&ok;++k) {
                double dx=x-s.pos[3*k],dy=y-s.pos[3*k+1],dz=z-s.pos[3*k+2];
                if (dx*dx+dy*dy+dz*dz < min_dist*min_dist) ok=false;
            }
            if (!ok) continue;
            s.pos[3*placed]=x; s.pos[3*placed+1]=y; s.pos[3*placed+2]=z;
            ++placed;
        }
        // Fallback: Gaussian packing
        if (placed < tmpl.N) {
            for (int i=placed; i<tmpl.N; ++i) {
                double x,y,z,r;
                do {
                    x=rng.randn(); y=rng.randn(); z=rng.randn();
                    r=std::sqrt(x*x+y*y+z*z);
                } while (r<1e-10);
                double scale = tmpl.R * std::cbrt(rng.rand()) / r;
                s.pos[3*i]=x*scale; s.pos[3*i+1]=y*scale; s.pos[3*i+2]=z*scale;
            }
        }
        double e = local_minimize(s, pot, 800);
        if (e < best_e && e < 1e29) { best_e=e; best=s; }
    }
    return best;
}

// ─────────────────────────────────────────────────────────────────────────────
// Stress-guided perturbation — uses per-atom energy stress
// ─────────────────────────────────────────────────────────────────────────────
static void stress_perturb(System& sys, const Potential& pot,
    double step, double T_stress, double n_frac, RNG& rng)
{
    const int N = sys.N;
    std::vector<double> ae(N, 0.0);
    double E0 = pot.energy_only(sys);
    const double h = 0.01 * sys.mean_sigma;

    for (int i=0; i<N; ++i) {
        double x0 = sys.pos[3*i];
        sys.pos[3*i] += h;
        double Ep = pot.energy_only(sys);
        sys.pos[3*i] = x0;
        ae[i] = std::abs(Ep - E0);
    }

    double smin = *std::min_element(ae.begin(),ae.end());
    double smax = *std::max_element(ae.begin(),ae.end());
    double srange = smax - smin + 1e-12;

    // Build CDF
    std::vector<double> cdf(N);
    double psum = 0;
    for (int i=0;i<N;++i) {
        double s = (ae[i]-smin)/srange;
        cdf[i] = std::exp(s/(T_stress+1e-12));
        psum += cdf[i];
    }
    for (int i=0;i<N;++i) cdf[i] /= psum;
    for (int i=1;i<N;++i) cdf[i] += cdf[i-1];

    int n_move = std::max(1,(int)std::round(n_frac*N));
    std::vector<bool> chosen(N, false);
    int moved=0, att=0;
    while (moved < n_move && att < N*4) {
        ++att;
        int idx = (int)(std::lower_bound(cdf.begin(),cdf.end(),rng.rand()) - cdf.begin());
        if (idx >= N) idx = N-1;
        if (chosen[idx]) continue;
        chosen[idx] = true;
        double amp = step * (0.5 + (ae[idx]-smin)/srange);
        sys.pos[3*idx]   += rng.randn()*amp;
        sys.pos[3*idx+1] += rng.randn()*amp;
        sys.pos[3*idx+2] += rng.randn()*amp;
        ++moved;
    }
    sys.project_sphere();
}

// ─────────────────────────────────────────────────────────────────────────────
// Graph topology move — surface relocation + coordination swap
// ─────────────────────────────────────────────────────────────────────────────
static void graph_move(System& sys, double cutoff_factor, RNG& rng) {
    const int N = sys.N;
    double cut2 = (cutoff_factor*sys.mean_sigma)*(cutoff_factor*sys.mean_sigma);
    std::vector<int> cn(N, 0);
    for (int i=0;i<N-1;++i) for (int j=i+1;j<N;++j)
        if (sys.rij2(i,j) < cut2) { ++cn[i]; ++cn[j]; }

    std::vector<int> idx(N); std::iota(idx.begin(),idx.end(),0);
    std::sort(idx.begin(),idx.end(),[&](int a,int b){return cn[a]<cn[b];});

    if (rng.rand() < 0.6) {
        // Surface relocation: move low-CN atom next to high-CN atom
        int pool = std::max(1, N/4);
        int src  = idx[rng.randi(0, pool)];
        int tgt  = idx[rng.randi(N-pool, N)];
        if (src == tgt) return;
        double r_eq = std::pow(2.0,1.0/6.0) * sys.mean_sigma;
        double dx=rng.randn(), dy=rng.randn(), dz=rng.randn();
        double n=std::sqrt(dx*dx+dy*dy+dz*dz);
        if (n<1e-10) return;
        double amp = r_eq*(0.85+0.30*rng.rand())/n;
        sys.pos[3*src]  = sys.pos[3*tgt]   + dx*amp;
        sys.pos[3*src+1]= sys.pos[3*tgt+1] + dy*amp;
        sys.pos[3*src+2]= sys.pos[3*tgt+2] + dz*amp;
    } else {
        // Position swap (preserves species labels — changes chemical ordering)
        int pool = std::max(1, N/4);
        int i = idx[rng.randi(0, pool)];
        int j = idx[rng.randi(N-pool, N)];
        for (int d=0;d<3;++d) std::swap(sys.pos[3*i+d], sys.pos[3*j+d]);
    }
    sys.project_sphere();
}

// ─────────────────────────────────────────────────────────────────────────────
// Homotop swap — swap species labels between two atoms of different species
// Geometry unchanged, chemical ordering changes
// ─────────────────────────────────────────────────────────────────────────────
static bool homotop_swap(System& sys, RNG& rng) {
    if (sys.n_species() < 2) return false;
    for (int att=0; att<200; ++att) {
        int i = rng.randi(0, sys.N);
        int j = rng.randi(0, sys.N);
        if (i==j || sys.type_id[i]==sys.type_id[j]) continue;
        std::swap(sys.type_id[i], sys.type_id[j]);
        return true;
    }
    return false;
}

// ─────────────────────────────────────────────────────────────────────────────
// Curvature-guided escape — finite-difference soft mode along Hessian
// eigenvector with lowest curvature → most efficient escape direction
// ─────────────────────────────────────────────────────────────────────────────
static void curvature_escape(System& sys, const Potential& pot,
    double step, RNG& rng)
{
    const int N3 = 3*sys.N;
    const double h = 0.02*sys.mean_sigma;

    // Finite-difference gradient at current point
    std::vector<double> g0(N3,0), gp(N3,0);
    pot.energy_gradient(sys, g0.data());

    // Build approximate diagonal Hessian (fast: N3 extra evaluations)
    std::vector<double> hess_diag(N3);
    for (int k=0; k<N3; ++k) {
        double x0 = sys.pos[k];
        sys.pos[k] += h;
        std::fill(gp.begin(),gp.end(),0);
        pot.energy_gradient(sys, gp.data());
        sys.pos[k] = x0;
        hess_diag[k] = (gp[k]-g0[k])/h;
    }

    // Soft mode = direction of smallest curvature
    // Perturb along N3-d vector weighted by 1/|H_kk|
    std::vector<double> soft(N3);
    double norm=0;
    for (int k=0;k<N3;++k) {
        soft[k] = rng.randn() / (std::abs(hess_diag[k]) + 0.1);
        norm += soft[k]*soft[k];
    }
    norm = std::sqrt(norm) + 1e-15;
    for (int k=0;k<N3;++k) sys.pos[k] += step * soft[k] / norm;
    sys.project_sphere();
}

// ─────────────────────────────────────────────────────────────────────────────
// Isomer filter — species-aware RMSD + energy window
// ─────────────────────────────────────────────────────────────────────────────
struct IsomerFilter {
    std::vector<std::pair<double,System>> isomers;
    double rmsd_thresh, energy_window;
    int    max_isomers;
    double best_energy = 1e30;

    IsomerFilter(double rt, double ew, int mi)
        : rmsd_thresh(rt), energy_window(ew), max_isomers(mi) {}

    void insert(double energy, const System& sys) {
        if (energy < best_energy) best_energy = energy;
        double cutoff = best_energy + std::abs(best_energy)*energy_window;
        if (energy > cutoff) return;

        // Check RMSD against existing isomers
        for (auto& [e,s] : isomers) {
            double rmsd = kabsch_rmsd_species(sys, s, sys.mean_sigma);
            if (rmsd < rmsd_thresh) {
                if (energy < e) e = energy; // update energy
                return;
            }
        }
        isomers.push_back({energy, sys});
        // Sort and trim
        std::sort(isomers.begin(),isomers.end(),
                  [](auto& a,auto& b){return a.first<b.first;});
        if ((int)isomers.size() > max_isomers) isomers.resize(max_isomers);
    }

    // Re-apply energy window after all insertions
    void trim() {
        if (isomers.empty()) return;
        best_energy = isomers[0].first;
        double cutoff = best_energy + std::abs(best_energy)*energy_window;
        while (!isomers.empty() && isomers.back().first > cutoff)
            isomers.pop_back();
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// ATLASWorker — single basin-hopping run
// ─────────────────────────────────────────────────────────────────────────────
class ATLASWorker {
public:
    ATLASWorker(const System& tmpl,
                const Potential& pot,
                const ATLASConfig& cfg,
                MotifLibrary& motif_lib,
                std::atomic<double>& shared_best,
                GPSurrogate& shared_gp,
                std::mutex& gp_mutex,
                ActiveLearner* active_learner = nullptr,
                Population* shared_pop = nullptr)
        : tmpl_(tmpl), pot_(pot), cfg_(cfg),
          rng_(cfg.seed + (uint64_t)cfg.run_id * 6364136223846793005ULL),
          motif_lib_(motif_lib),
          shared_best_(shared_best),
          shared_gp_(shared_gp),
          gp_mutex_(gp_mutex),
          taboo_(cfg.taboo_lambda, cfg.taboo_delta,
                 cfg.taboo_max_size, cfg.taboo_thresh),
          active_learner_(active_learner),
          shared_pop_(shared_pop)
    {
        step_     = cfg.stress_step_frac * tmpl.mean_sigma;
        T_bh_     = cfg.bh_temp * pot.get_mean_sigma(tmpl) * 10.0;
        T_stress_ = cfg.stress_temp_init;
        if (cfg.use_ard) {
            ard_.update_interval = cfg.ard_update_interval;
        }
        if (cfg.use_ngd) {
            ngd_opt_.use_ngd = true;
        }
    }
    RunResult run() {
        const int GENS = cfg_.n_generations;

        // ── Initialise ──────────────────────────────────────────────────────
        System cur = smart_init(tmpl_, pot_, cfg_.init_attempts, rng_);
        if (pot_.energy_only(cur) > 1e29) cur = tmpl_;
        double cur_e = pot_.energy_only(cur);
        Fingerprint cur_fp = Fingerprint::compute(cur);
        double cur_bias = taboo_.bias(cur_fp);
        double cur_eff  = cur_e + cur_bias;

        _register(cur, cur_e, cur_fp);
        result_.best_energy = cur_e;
        result_.best_sys    = cur;
        result_.candidates.push_back({cur_e, cur});
        result_.energy_history.reserve(GENS);

        int stag = 0;
        double T_bh = T_bh_;

        // ── Main generation loop ─────────────────────────────────────────────
        for (int gen=0; gen<GENS; ++gen) {

            // ── Generate candidate moves via gradient + random blend ──────────
            double best_score = 1e30;
            System best_trial = cur;
            int    best_mtype = 0;

            // Try gradient-informed proposal first (Session 2 innovation)
            bool used_gradient = false;
            if (cfg_.use_gp) {
                std::lock_guard<std::mutex> lk(gp_mutex_);
                if (shared_gp_.is_trained()) {
                    System grad_trial = cur;
                    bool ok = acq_.apply(grad_trial, cur_fp, shared_gp_,
                                         step_, cfg_.gp_kappa, rng_.eng);
                    if (ok) {
                        auto gfp = Fingerprint::compute(grad_trial);
                        double gscore = shared_gp_.lcb(gfp.data, cfg_.gp_kappa);
                        if (gscore < best_score) {
                            best_score = gscore;
                            best_trial = grad_trial;
                            best_mtype = 4;  // gradient move
                        }
                        ++result_.moves_gp_guided;
                        used_gradient = true;
                    }
                }
            }

            // Also generate N_cand random candidates and score them
            int n_cand = cfg_.gp_n_candidates;
            for (int c=0; c<n_cand; ++c) {
                System trial = cur;
                int mtype = 0;
                double r = rng_.rand();

                if (r < cfg_.graph_move_prob) {
                    graph_move(trial, cfg_.bond_cutoff, rng_);
                    mtype = 1;
                } else if (cfg_.homotop_search
                           && tmpl_.n_species() > 1
                           && r < cfg_.graph_move_prob + 0.15) {
                    if (homotop_swap(trial, rng_)) mtype = 3;
                    else { stress_perturb(trial,pot_,step_,T_stress_,
                                          cfg_.stress_n_frac,rng_); }
                } else {
                    stress_perturb(trial, pot_, step_, T_stress_,
                                   cfg_.stress_n_frac, rng_);
                    mtype = 0;
                }

                double score = 0.0;
                if (cfg_.use_gp) {
                    auto trial_fp = Fingerprint::compute(trial);
                    std::lock_guard<std::mutex> lk(gp_mutex_);
                    if (shared_gp_.is_trained()) {
                        // ARD-weighted LCB (Session 4)
                        auto wfp = cfg_.use_ard
                            ? ard_.apply(trial_fp.data)
                            : trial_fp.data;
                        score = shared_gp_.lcb(wfp, cfg_.gp_kappa);
                    } else {
                        score = rng_.rand();
                    }
                }

                if (score < best_score) {
                    best_score = score;
                    best_trial = trial;
                    best_mtype = mtype;
                }
            }

            // ── Local minimization ───────────────────────────────────────────
            double e_new = local_minimize(best_trial, pot_, cfg_.local_maxiter);
            if (e_new >= 1e29) {
                result_.energy_history.push_back(result_.best_energy);
                continue;
            }

            // ── Update GP (+ optional active learning oracle call) ────────────
            Fingerprint fp_new = Fingerprint::compute(best_trial);

            // Active learning: decide whether to call oracle on this structure
            if (cfg_.use_active_learning && active_learner_) {
                double e_corrected = active_learner_->call_oracle(
                    best_trial, fp_new, e_new, shared_gp_, gp_mutex_, gen);
                if (e_corrected < e_new) {
                    e_new = e_corrected;  // use oracle energy if available
                    ++result_.oracle_calls;
                }
            }

            {
                std::lock_guard<std::mutex> lk(gp_mutex_);
                shared_gp_.add(fp_new.data, e_new);
                if (shared_gp_.dirty && shared_gp_.n_obs() >= 3)
                    shared_gp_.rebuild();
            }

            // ── Taboo bias ───────────────────────────────────────────────────
            double bias_new = taboo_.bias(fp_new);
            double eff_new  = e_new + bias_new;

            // ── Metropolis acceptance ────────────────────────────────────────
            double dE = eff_new - cur_eff;
            bool accept = (dE <= 0.0)
                       || (rng_.rand() < std::exp(-dE / (T_bh + 1e-30)));

            if (accept) {
                cur = best_trial; cur_e = e_new;
                cur_fp = fp_new; cur_bias = bias_new; cur_eff = eff_new;
                _register(cur, cur_e, cur_fp);
                // Feed into population
                if (cfg_.use_population && shared_pop_)
                    shared_pop_->insert(cur, cur_e, cur_fp);
            }

            // ── Track move types ─────────────────────────────────────────────
            if (best_mtype==0) ++result_.moves_stress;
            else if (best_mtype==1) ++result_.moves_graph;
            else if (best_mtype==3) ++result_.moves_homotop;
            else if (best_mtype==4) ++result_.moves_gradient;

            // ── Update global best ───────────────────────────────────────────
            result_.candidates.push_back({cur_e, cur});
            if (cur_e < result_.best_energy) {
                result_.best_energy = cur_e;
                result_.best_sys    = cur;
                stag = 0;
                double old = shared_best_.load();
                while (cur_e < old &&
                       !shared_best_.compare_exchange_weak(old, cur_e)) {}
            } else {
                ++stag;
            }

            // ── Motif injection on stagnation ────────────────────────────────
            if (stag > 0 && stag % cfg_.motif_inject_interval == 0) {
                MotifLibrary::Motif m;
                if (motif_lib_.random_motif(m, rng_)) {
                    System inj = result_.best_sys;
                    // Perturb using best-known fragment as seed geometry
                    int M = m.fragment.N;
                    for (int k=0; k<std::min(M,inj.N); ++k) {
                        inj.pos[3*k]   = m.fragment.pos[3*k];
                        inj.pos[3*k+1] = m.fragment.pos[3*k+1];
                        inj.pos[3*k+2] = m.fragment.pos[3*k+2];
                    }
                    double e_inj = local_minimize(inj, pot_, cfg_.local_maxiter);
                    if (e_inj < 1e29) {
                        ++result_.moves_motif;
                        result_.candidates.push_back({e_inj, inj});
                        if (e_inj < result_.best_energy) {
                            result_.best_energy = e_inj;
                            result_.best_sys    = inj;
                            stag = 0;
                        }
                    }
                }
            }

            // ── Homotop search on stagnation ─────────────────────────────────
            if (cfg_.homotop_search && tmpl_.n_species() > 1
                && stag > 0 && stag % 50 == 0) {
                System ht = result_.best_sys;
                double ht_e = result_.best_energy;
                for (int hs=0; hs<cfg_.homotop_steps; ++hs) {
                    System trial = ht;
                    if (!homotop_swap(trial, rng_)) break;
                    double e_ht = local_minimize(trial, pot_, 300);
                    if (e_ht < 1e29) {
                        result_.candidates.push_back({e_ht, trial});
                        if (e_ht < ht_e) {
                            ht_e=e_ht; ht=trial; ++result_.moves_homotop;
                        }
                        if (e_ht < result_.best_energy) {
                            result_.best_energy=e_ht;
                            result_.best_sys=trial;
                            stag=0;
                        }
                    }
                }
            }

            // ── Population crossover (Session 4) ─────────────────────────────
            if (cfg_.use_population && shared_pop_ &&
                gen > 0 && gen % cx_engine_.interval == 0 &&
                shared_pop_->size() >= 2)
            {
                for (int ci=0; ci<cfg_.pop_n_children; ++ci) {
                    System child = cx_engine_.generate_child(*shared_pop_, rng_.eng);
                    double e_child = local_minimize(child, pot_, cfg_.local_maxiter);
                    if (e_child < 1e29) {
                        ++result_.moves_crossover;
                        Fingerprint fp_child = Fingerprint::compute(child);
                        result_.candidates.push_back({e_child, child});
                        shared_pop_->insert(child, e_child, fp_child);
                        {
                            std::lock_guard<std::mutex> lk(gp_mutex_);
                            shared_gp_.add(fp_child.data, e_child);
                        }
                        if (e_child < result_.best_energy) {
                            result_.best_energy = e_child;
                            result_.best_sys    = child;
                            stag = 0;
                        }
                    }
                }
                ++cx_engine_.n_crossovers;
            }

            // ── ARD weight update (Session 4) ─────────────────────────────────
            if (cfg_.use_ard && gen > 0 && gen % cfg_.ard_update_interval == 0) {
                std::lock_guard<std::mutex> lk(gp_mutex_);
                if (shared_gp_.is_trained()) {
                    ard_.update(shared_gp_);
                    ++result_.ard_updates;
                }
            }

            // ── Natural gradient hyperparam update (Session 4) ────────────────
            if (cfg_.use_ngd && gen > 0 && gen % 60 == 0) {
                std::lock_guard<std::mutex> lk(gp_mutex_);
                if (shared_gp_.n_obs() >= 12 && !shared_gp_.dirty) {
                    ngd_opt_.optimise(shared_gp_);
                }
            }


            int deep_stag_thresh = std::max(20, GENS/5);
            if (stag > 0 && stag % deep_stag_thresh == 0) {
                System esc = result_.best_sys;
                double esc_step = tmpl_.mean_sigma*(1.0+0.5*rng_.rand());

                if (cfg_.use_eigenvec_escape && rng_.rand() < 0.65) {
                    // Lanczos eigenvector following — guaranteed softest mode
                    eigenvector_escape(esc, pot_, esc_step,
                                       cfg_.seed ^ (uint64_t)gen ^ (uint64_t)stag);
                    ++result_.moves_eigenvec;
                } else {
                    // Random kick fallback
                    double kick = tmpl_.mean_sigma*(1.2+0.8*rng_.rand());
                    for (int i=0;i<esc.N;++i) {
                        esc.pos[3*i]  +=rng_.randn()*kick;
                        esc.pos[3*i+1]+=rng_.randn()*kick;
                        esc.pos[3*i+2]+=rng_.randn()*kick;
                    }
                }

                double e_esc = local_minimize(esc, pot_, cfg_.local_maxiter);
                if (e_esc < 1e29) {
                    Fingerprint fp_esc = Fingerprint::compute(esc);

                    // Active learning: possibly refine escape energy with oracle
                    if (cfg_.use_active_learning && active_learner_) {
                        e_esc = active_learner_->call_oracle(
                            esc, fp_esc, e_esc, shared_gp_, gp_mutex_, gen);
                        ++result_.oracle_calls;
                    }

                    {
                        std::lock_guard<std::mutex> lk(gp_mutex_);
                        shared_gp_.add(fp_esc.data, e_esc);
                        if (shared_gp_.dirty && shared_gp_.n_obs() >= 3)
                            shared_gp_.rebuild();
                    }
                    if (e_esc < result_.best_energy) {
                        result_.best_energy=e_esc; result_.best_sys=esc; stag=0;
                    }
                    cur=esc; cur_e=e_esc; cur_fp=fp_esc;
                    cur_bias=taboo_.bias(cur_fp); cur_eff=cur_e+cur_bias;
                }
            }

            // ── Cooling ──────────────────────────────────────────────────────
            T_stress_ = std::max(cfg_.stress_temp_min,
                                 T_stress_ * cfg_.stress_temp_decay);
            step_    *= cfg_.stress_step_decay;
            T_bh      = std::max(cfg_.bh_temp_min * 10.0,
                                 T_bh * cfg_.bh_temp_decay);
            taboo_.adapt();

            result_.energy_history.push_back(result_.best_energy);

            // ── JSON progress output ─────────────────────────────────────────
            if (cfg_.progress && (gen%5==0 || gen==GENS-1)) {
                #pragma omp critical
                {
                    std::cout << "{\"type\":\"progress\","
                              << "\"run\":"   << cfg_.run_id  << ","
                              << "\"gen\":"   << gen+1        << ","
                              << "\"total\":" << GENS         << ","
                              << "\"best\":"
                              << std::fixed<<std::setprecision(8)
                              << result_.best_energy           << ","
                              << "\"basins\":" << result_.n_new_basins << "}\n"
                              << std::flush;
                }
            }
        } // end generation loop

        // ── Final tight refinement of best structure ─────────────────────────
        double e_final = tight_minimize(result_.best_sys, pot_);
        if (e_final < result_.best_energy && e_final < 1e29) {
            result_.best_energy = e_final;
            result_.candidates.push_back({e_final, result_.best_sys});
        }

        result_.revisit_rate = (taboo_.checks > 0)
            ? (double)taboo_.revisits / taboo_.checks : 0.0;

        return result_;
    }

private:
    const System& tmpl_;
    const Potential& pot_;
    const ATLASConfig& cfg_;
    RNG rng_;
    MotifLibrary& motif_lib_;
    std::atomic<double>& shared_best_;
    GPSurrogate& shared_gp_;
    std::mutex& gp_mutex_;
    TabooMemory taboo_;
    RunResult result_{};
    double step_, T_bh_, T_stress_;
    AcquisitionEngine  acq_;
    ActiveLearner*     active_learner_ = nullptr;
    Population*        shared_pop_     = nullptr;   // Session 4
    CrossoverEngine    cx_engine_;                  // Session 4
    AdaptiveDescriptor ard_;                        // Session 4
    NaturalGradientOptimiser ngd_opt_;              // Session 4

    void _register(const System& sys, double e, const Fingerprint& fp) {
        if (!taboo_.is_known(fp)) {
            taboo_.add(fp);
            ++result_.n_new_basins;
            motif_lib_.deposit(sys, cfg_.motif_core_frac, cfg_.bond_cutoff);
        }
    }
};
