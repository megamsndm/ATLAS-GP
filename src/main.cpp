/**
 * main.cpp — ATLAS-GP
 * Usage: atlas <formula> [options]
 * Examples: atlas Ar13 --pot LJ --runs 8
 *           atlas Au32Ag6 --pot GUPTA --runs 16 --gens 500 --progress
 *           atlas Fe12O18 --pot MACE --runs 4 --gens 200 --progress
 */
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <omp.h>
#include <atomic>
#include <mutex>

#include "system.hpp"
#include "potential/potential_factory.hpp"
#include "fingerprint.hpp"
#include "gp_surrogate.hpp"
#include "atlas_worker.hpp"
#include "parallel_tempering.hpp"
#include "active_learning.hpp"
#include "population.hpp"
#include "adaptive_descriptor.hpp"
#include "natural_gradient.hpp"
#include "io.hpp"

static void print_banner() {
    std::cout <<
    "=================================================================\n"
    "  ATLAS-GP  |  Adaptive Topology-driven Landscape-Aware Search\n"
    "  GP-guided basin-hopping  |  Universal potential engine\n"
    "  All 103 elements  |  Atomic + Molecular  |  MACE-MP-0 ready\n"
    "=================================================================\n\n";
}

static void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " <formula> [options]\n\n"
    "Formula:  Ar13  Ni13  Au32Ag6  Fe12O18  Cu10Pt3  (any elements)\n\n"
    "  --pot LJ|MORSE|GUPTA|UFF|MACE|AUTO   potential (default AUTO)\n"
    "  --runs N         independent runs (default 8)\n"
    "  --gens N         generations per run (default 300)\n"
    "  --workers N      OpenMP threads (default all cores)\n"
    "  --seed N         RNG seed (0=random)\n"
    "  --no-gp          disable GP surrogate\n"
    "  --no-homotop     disable homotop search\n"
    "  --stress-step F  step size × sigma (default 0.25)\n"
    "  --stress-nfrac F fraction of atoms moved (default 0.35)\n"
    "  --graph-prob F   graph move probability (default 0.30)\n"
    "  --bond-cutoff F  bond cutoff × sigma (default 1.35)\n"
    "  --bh-temp F      BH temperature (default 1.0)\n"
    "  --taboo-lambda F taboo repulsion (default 0.05)\n"
    "  --gp-kappa F     LCB kappa (default 2.0)\n"
    "  --gp-cands N     GP candidate moves (default 5)\n"
    "  --rmsd-thresh F  isomer RMSD threshold (default 0.20)\n"
    "  --max-isomers N  max isomers (default 30)\n"
    "  --energy-window F energy window (default 0.40)\n"
    "  --output-dir D   output directory (default atlas_results)\n"
    "  --mace-script P  path to mace_server.py\n"
    "  --progress       emit JSON progress for GUI\n"
    "  --verbose        verbose statistics\n"
    "  --pt             enable parallel tempering\n"
    "  --pt-replicas N  replicas per PT group (default 4)\n"
    "  --pt-tmin F      PT coldest temperature (default 0.3)\n"
    "  --pt-tmax F      PT hottest temperature (default 5.0)\n"
    "  --no-eigenvec    disable eigenvector escape (use diagonal approx)\n"
    "  --active         enable active learning (MACE oracle refinement)\n"
    "  --oracle-budget N max MACE calls when using --active (default 50)\n"
    "  --no-pop         disable population crossover\n"
    "  --pop-size N     population size (default 20)\n"
    "  --no-ard         disable adaptive descriptor weighting\n"
    "  --no-ngd         disable natural gradient hyperparameter updates\n\n";
}

int main(int argc, char* argv[]) {
    print_banner();

    if (argc < 2 || std::string(argv[1])=="-h" || std::string(argv[1])=="--help") {
        print_usage(argv[0]); return 0;
    }

    std::vector<std::pair<std::string,int>> composition;
    try { composition = parse_formula(argv[1]); }
    catch (const std::exception& e) {
        std::cerr << "ERROR parsing formula '" << argv[1] << "': " << e.what() << "\n";
        return 1;
    }

    ATLASConfig cfg;
    std::string pot_name = "AUTO";
    std::string mace_script = "mace_server.py";

    for (int i=2; i<argc; ++i) {
        std::string a = argv[i];
        auto nx  = [&]()->std::string{ if(i+1>=argc){std::cerr<<"Missing value for "<<a<<"\n";exit(1);}return argv[++i]; };
        auto nxd = [&](){ return std::stod(nx()); };
        auto nxi = [&](){ return std::stoi(nx()); };

        if      (a=="--pot")           pot_name              = nx();
        else if (a=="--runs")          cfg.n_runs            = nxi();
        else if (a=="--gens")          cfg.n_generations     = nxi();
        else if (a=="--workers")       cfg.n_workers         = nxi();
        else if (a=="--seed")          cfg.seed              = (uint64_t)nxi();
        else if (a=="--no-gp")         cfg.use_gp            = false;
        else if (a=="--no-homotop")    cfg.homotop_search    = false;
        else if (a=="--homotop")       cfg.homotop_search    = true;
        else if (a=="--stress-step")   cfg.stress_step_frac  = nxd();
        else if (a=="--stress-nfrac")  cfg.stress_n_frac     = nxd();
        else if (a=="--graph-prob")    cfg.graph_move_prob   = nxd();
        else if (a=="--bond-cutoff")   cfg.bond_cutoff       = nxd();
        else if (a=="--bh-temp")       cfg.bh_temp           = nxd();
        else if (a=="--bh-decay")      cfg.bh_temp_decay     = nxd();
        else if (a=="--taboo-lambda")  cfg.taboo_lambda      = nxd();
        else if (a=="--taboo-delta")   cfg.taboo_delta       = nxd();
        else if (a=="--gp-kappa")      cfg.gp_kappa          = nxd();
        else if (a=="--gp-cands")      cfg.gp_n_candidates   = nxi();
        else if (a=="--rmsd-thresh")   cfg.isomer_rmsd_thresh= nxd();
        else if (a=="--max-isomers")   cfg.max_isomers       = nxi();
        else if (a=="--energy-window") cfg.energy_window     = nxd();
        else if (a=="--output-dir")    cfg.output_dir        = nx();
        else if (a=="--mace-script")   mace_script           = nx();
        else if (a=="--progress")      cfg.progress          = true;
        else if (a=="--verbose")       cfg.verbose           = true;
        // Session 3 flags
        else if (a=="--pt")            cfg.use_pt            = true;
        else if (a=="--pt-replicas")   cfg.pt_n_replicas     = nxi();
        else if (a=="--pt-tmin")       cfg.pt_T_min          = nxd();
        else if (a=="--pt-tmax")       cfg.pt_T_max          = nxd();
        else if (a=="--no-eigenvec")   cfg.use_eigenvec_escape = false;
        else if (a=="--active")        cfg.use_active_learning = true;
        else if (a=="--oracle-budget") cfg.oracle_budget     = nxi();
        // Session 4 flags
        else if (a=="--no-pop")        cfg.use_population    = false;
        else if (a=="--pop-size")      cfg.pop_size          = nxi();
        else if (a=="--no-ard")        cfg.use_ard           = false;
        else if (a=="--no-ngd")        cfg.use_ngd           = false;
        else { std::cerr<<"Unknown arg: "<<a<<"\nRun atlas --help\n"; return 1; }
    }

    if (cfg.seed==0)
        cfg.seed=(uint64_t)std::chrono::high_resolution_clock::now()
                 .time_since_epoch().count();

    // Build system template
    System tmpl = System::build_atomic(composition);
    double sigma_sum=0;
    for (auto& sp:tmpl.species) sigma_sum+=get_uff(sp).sigma;
    tmpl.mean_sigma = sigma_sum/tmpl.species.size();
    tmpl.R = 1.3*tmpl.mean_sigma*std::cbrt((double)tmpl.N);

    // Build potential
    std::unique_ptr<Potential> pot;
    try {
        pot = make_potential(pot_name, tmpl, mace_script);
        tmpl.mean_sigma = pot->get_mean_sigma(tmpl);
        tmpl.R = 1.3*tmpl.mean_sigma*std::cbrt((double)tmpl.N);
    } catch (const std::exception& e) {
        std::cerr<<"ERROR building potential: "<<e.what()<<"\n"; return 1;
    }

    std::cout << "[ATLAS] Formula   : " << tmpl.formula() << "\n"
              << "[ATLAS] Atoms     : " << tmpl.N << "\n"
              << "[ATLAS] Species   : ";
    for (auto& s:tmpl.species) std::cout<<s<<" ";
    std::cout << "\n"
              << "[ATLAS] Potential : " << pot->name() << "\n"
              << "[ATLAS] Runs      : " << cfg.n_runs << "\n"
              << "[ATLAS] Gens/run  : " << cfg.n_generations << "\n"
              << "[ATLAS] GP        : " << (cfg.use_gp?"ON":"OFF") << "\n"
              << "[ATLAS] Homotop   : "
              << (cfg.homotop_search&&tmpl.n_species()>1?"ON":"OFF") << "\n\n";

    if (cfg.progress) {
        std::cout<<"{\"type\":\"header\",\"formula\":\""<<tmpl.formula()
                 <<"\",\"N\":"<<tmpl.N<<",\"potential\":\""<<pot->name()
                 <<"\",\"runs\":"<<cfg.n_runs<<",\"gens\":"<<cfg.n_generations
                 <<"}\n"<<std::flush;
    }

    // Shared state
    std::atomic<double> shared_best{1e30};
    GPSurrogate         shared_gp;
    std::mutex          gp_mutex;
    MotifLibrary        motif_lib(cfg.motif_lib_maxsize, cfg.motif_min_rmsd);
    IsomerFilter        filter(cfg.isomer_rmsd_thresh * tmpl.mean_sigma,
                               cfg.energy_window, cfg.max_isomers);
    std::mutex          filter_mutex;
    // Session 4: shared population for crossover
    Population          shared_pop(cfg.pop_size, cfg.pop_diversity);

    // Session 3: Active learner (optional oracle)
    ActiveLearner active_learner;
    if (cfg.use_active_learning) {
        // Build fast + oracle potentials
        // Oracle = MACE, fast = AUTO (will be UFF/Gupta/LJ)
        std::unique_ptr<Potential> oracle_pot;
        try {
            oracle_pot = make_potential("MACE", tmpl, mace_script);
            active_learner.set_oracle(oracle_pot.release(),
                                      cfg.oracle_budget,
                                      cfg.oracle_uncertainty);
            std::cout << "[AL] Active learning enabled. Oracle budget: "
                      << cfg.oracle_budget << " MACE calls\n";
        } catch (...) {
            std::cout << "[AL] Warning: MACE not available, disabling active learning\n";
            cfg.use_active_learning = false;
        }
    }

    omp_set_num_threads(cfg.n_workers);
    auto t0 = std::chrono::high_resolution_clock::now();

    std::vector<RunResult> all_results(cfg.n_runs);

    if (!cfg.use_pt) {
        // ── Standard parallel independent runs ───────────────────────────
        #pragma omp parallel for schedule(dynamic,1)
        for (int r=0; r<cfg.n_runs; ++r) {
            ATLASConfig rc=cfg; rc.run_id=r;
            ActiveLearner* al = cfg.use_active_learning ? &active_learner : nullptr;
            Population*    pp = cfg.use_population      ? &shared_pop     : nullptr;
            ATLASWorker worker(tmpl, *pot, rc, motif_lib,
                               shared_best, shared_gp, gp_mutex, al, pp);
            all_results[r]=worker.run();
            {
                std::lock_guard<std::mutex> lk(filter_mutex);
                for (auto& [e,s]:all_results[r].candidates) filter.insert(e,s);
            }
            if (!cfg.progress) {
                #pragma omp critical
                std::cout<<"[Run "<<std::setw(2)<<r+1<<"/"<<cfg.n_runs<<"]  best="
                         <<std::fixed<<std::setprecision(8)<<all_results[r].best_energy
                         <<"  basins="<<all_results[r].n_new_basins<<"\n"<<std::flush;
            }
        }
    } else {
        // ── Parallel tempering: groups of replicas ────────────────────────
        // Divide n_runs into PT groups of pt_n_replicas each
        int n_groups = std::max(1, cfg.n_runs / cfg.pt_n_replicas);
        int R        = cfg.pt_n_replicas;

        PTConfig pt_cfg = make_pt_config(cfg, R);
        PTEngine pt_engine(pt_cfg, cfg.seed);

        std::cout << "[PT] " << n_groups << " group(s) of "
                  << R << " replicas each\n"
                  << "[PT] T ladder: ";
        for (double t : pt_engine.temps) std::cout<<std::fixed<<std::setprecision(2)<<t<<" ";
        std::cout << "\n";

        // Run groups sequentially, replicas within a group in parallel
        for (int g=0; g<n_groups; ++g) {
            // Within this group: run R replicas with different temperatures
            std::vector<RunResult> group_results(R);

            #pragma omp parallel for schedule(static,1) num_threads(std::min(R, cfg.n_workers))
            for (int r=0; r<R; ++r) {
                ATLASConfig rc     = cfg;
                rc.run_id          = g*R + r;
                rc.bh_temp         = pt_engine.temps[r];
                rc.bh_temp_min     = pt_engine.temps[0] * 0.5;
                rc.use_pt          = false;  // prevent recursion

                ActiveLearner* al = cfg.use_active_learning ? &active_learner : nullptr;
                Population*    pp = cfg.use_population      ? &shared_pop     : nullptr;
                ATLASWorker worker(tmpl, *pot, rc, motif_lib,
                                   shared_best, shared_gp, gp_mutex, al, pp);
                group_results[r] = worker.run();
            }

            // Collect group results
            for (int r=0; r<R; ++r) {
                int global_r = g*R + r;
                if (global_r < (int)all_results.size())
                    all_results[global_r] = group_results[r];
                std::lock_guard<std::mutex> lk(filter_mutex);
                for (auto& [e,s]:group_results[r].candidates) filter.insert(e,s);
            }

            // Attempt PT swaps between adjacent replicas
            // Populate replica states from group results
            std::vector<System>  rep_sys(R);
            std::vector<double>  rep_e(R);
            for (int r=0;r<R;++r) {
                rep_sys[r] = group_results[r].best_sys;
                rep_e[r]   = group_results[r].best_energy;
            }
            pt_engine.init(rep_sys, rep_e);

            int swaps = pt_engine.attempt_swaps(g % 2);
            if (cfg.verbose)
                std::cout << "[PT] Group " << g+1 << "/" << n_groups
                          << "  swaps=" << swaps << "\n" << std::flush;
        }

        if (cfg.verbose) pt_engine.print_stats();
    }

    filter.trim();
    if (!filter.isomers.empty()) {
        std::cout << "  Isomers found: " << filter.isomers.size()
                  << "  Best E: " << filter.isomers[0].first << "\n";
    }
    double elapsed=std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now()-t0).count();

    std::cout << "\n=== RESULTS ===\n"
              << "  Global min  : "<<std::fixed<<std::setprecision(8)
              <<(filter.isomers.empty()?1e30:filter.isomers[0].first)<<"\n"
              << "  Isomers     : " << filter.isomers.size() << "\n"
              << "  Time        : " << std::setprecision(1) << elapsed << " s\n\n";

    if (cfg.verbose) {
        int ts=0,tg=0,tm=0,th=0,tgp=0,tgrad=0,tev=0,tor=0,tcx=0,tard=0,tb=0;
        for (auto& r:all_results){ts+=r.moves_stress;tg+=r.moves_graph;
            tm+=r.moves_motif;th+=r.moves_homotop;tgp+=r.moves_gp_guided;
            tgrad+=r.moves_gradient;tev+=r.moves_eigenvec;
            tor+=r.oracle_calls;tcx+=r.moves_crossover;tard+=r.ard_updates;
            tb+=r.n_new_basins;}
        std::cout<<"  Stress moves   : "<<ts<<"\n  Graph moves    : "<<tg
                 <<"\n  Motif inj.     : "<<tm<<"\n  Homotop swaps  : "<<th
                 <<"\n  GP-guided      : "<<tgp<<"\n  Gradient moves : "<<tgrad
                 <<"\n  Eigenvec esc.  : "<<tev<<"\n  Crossovers     : "<<tcx
                 <<"\n  ARD updates    : "<<tard<<"\n  Oracle calls   : "<<tor
                 <<"\n  Total basins   : "<<tb<<"\n\n";
        if (cfg.use_population) shared_pop.print_summary();
    }

    if (cfg.use_active_learning) active_learner.print_summary();

    // Convert isomers to pairs for writer
    // Compatible with both old (pair<double,System>) and new (Isomer struct) filter
    std::vector<std::pair<double,System>> iso_pairs;
    for (auto& iso : filter.isomers)
        iso_pairs.push_back(std::make_pair(iso.first, iso.second));
    write_isomers(cfg.output_dir, tmpl.formula(), iso_pairs, pot->name());

    if (cfg.progress) {
        double be=filter.isomers.empty()?1e30:filter.isomers[0].first;
        std::cout<<"{\"type\":\"result\",\"formula\":\""<<tmpl.formula()
                 <<"\",\"potential\":\""<<pot->name()
                 <<"\",\"best_energy\":"<<std::setprecision(8)<<be
                 <<",\"n_isomers\":"<<filter.isomers.size()
                 <<",\"elapsed\":"<<std::setprecision(2)<<elapsed<<"}\n"<<std::flush;
    }
    return 0;
}
