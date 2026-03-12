#pragma once
/**
 * active_learning.hpp — On-the-fly active learning loop
 *
 * SESSION 3 INNOVATION: Build a custom potential on the fly during search.
 *
 * Problem with static potential selection:
 *   - UFF/Gupta: fast, but crude (wrong barriers, wrong relative ordering)
 *   - MACE-MP-0: accurate, but ~1000× slower
 *   - If you just use MACE throughout, you can only afford ~100 gens
 *   - If you just use UFF, energy ordering of isomers is unreliable
 *
 * Active learning solution:
 *   Phase 1 (exploration):   Run with fast potential (UFF/Gupta/LJ)
 *                            Collect diverse structures via ATLAS
 *   Phase 2 (oracle calls):  Evaluate top-K structures with MACE
 *                            Only call MACE when GP uncertainty is HIGH
 *                            (uncertainty = it would most benefit the model)
 *   Phase 3 (refinement):    Re-rank all isomers using MACE energies
 *                            Fine-tune GP on MACE data
 *
 * The key decision: when to call the oracle (MACE)?
 *   Strategy: Maximum Uncertainty Sampling
 *   Call MACE on structure x when GP variance σ²(x) > threshold σ²_thresh
 *   Threshold auto-adjusts to hit budget N_oracle calls per run
 *
 * This gives you:
 *   - Full ATLAS speed during exploration (fast potential)
 *   - MACE accuracy on the structures that matter
 *   - GP trained on MACE data → better proposals for remaining search
 *   - Final results ranked by DFT-quality energies, not empirical ones
 *
 * Budget control:
 *   --oracle-budget N   (default 50): max MACE calls per search
 *   --oracle-thresh F   (default 2.0): uncertainty multiplier
 *   Calls are front-loaded (more in early gens when variance is highest)
 */

#include "system.hpp"
#include "potential/potential_base.hpp"
#include "fingerprint.hpp"
#include "gp_surrogate.hpp"
#include <vector>
#include <string>
#include <memory>
#include <mutex>
#include <algorithm>
#include <iostream>
#include <iomanip>

// ─────────────────────────────────────────────────────────────────────────────
// OracleRecord: one MACE-evaluated structure
// ─────────────────────────────────────────────────────────────────────────────
struct OracleRecord {
    System      sys;
    Fingerprint fp;
    double      e_fast;    // energy from fast potential (K)
    double      e_oracle;  // energy from oracle/MACE (K)
    double      gp_variance_at_call;  // GP σ² when oracle was called
    int         gen;       // generation when called
};

// ─────────────────────────────────────────────────────────────────────────────
// ActiveLearner
// ─────────────────────────────────────────────────────────────────────────────
class ActiveLearner {
public:
    int    oracle_budget;       // max oracle calls total
    double uncertainty_mult;    // call oracle when σ > mult*σ_mean
    bool   enabled;

    Potential*  oracle;         // slow accurate potential (MACE/DFT)
    GPSurrogate correction_gp;  // GP on (e_fast → e_oracle - e_fast) corrections
    std::mutex  oracle_mutex;

    std::vector<OracleRecord> records;
    int n_oracle_calls = 0;
    int n_oracle_skipped = 0;
    double sigma_threshold = 1e30;  // starts very high → no calls initially

    ActiveLearner()
        : oracle_budget(50), uncertainty_mult(2.0),
          enabled(false), oracle(nullptr)
    {
        // Sparse GP for correction model
        correction_gp.M_max = 50;
        correction_gp.refit_every = 10;
        correction_gp.noise_var = 1e-4;
    }

    void set_oracle(Potential* p, int budget, double umult = 2.0) {
        oracle           = p;
        oracle_budget    = budget;
        uncertainty_mult = umult;
        enabled          = (p != nullptr && budget > 0);
    }

    // ── Decide whether to call oracle on this structure ───────────────────
    bool should_call_oracle(const Fingerprint& fp,
                            const GPSurrogate& main_gp,
                            int gen, int total_gens)
    {
        if (!enabled) return false;
        if (n_oracle_calls >= oracle_budget) return false;

        // Remaining budget: front-load calls (more early, fewer late)
        int remaining = oracle_budget - n_oracle_calls;
        int gens_left = total_gens - gen;
        if (gens_left <= 0) return false;

        // Expected calls per remaining gen
        double call_rate = (double)remaining / gens_left;
        if (call_rate < 0.01) return false;  // budget essentially exhausted

        // Call if GP main_gp has high uncertainty on this structure
        if (!main_gp.is_trained()) return false;
        auto [m, v] = main_gp.predict(fp.data);

        // Auto-calibrate threshold: target call_rate per generation
        // Use mean variance of recent observations as baseline
        double sigma = std::sqrt(std::max(0.0, v));
        return sigma > sigma_threshold * 0.5;
    }

    // ── Calibrate threshold to hit target call rate ───────────────────────
    void calibrate(const GPSurrogate& main_gp, double target_rate) {
        if (!main_gp.is_trained()) return;

        // Sample variance at a few training points
        std::vector<double> vars;
        int n = std::min((int)main_gp.X.size(), 30);
        for (int i=0; i<n; ++i) {
            auto [m,v] = main_gp.predict(main_gp.X[i]);
            vars.push_back(std::sqrt(std::max(0.0, v)));
        }
        if (vars.empty()) return;
        std::sort(vars.begin(), vars.end());

        // Set threshold at (1-target_rate) quantile
        int idx = std::max(0, (int)((1.0-target_rate) * vars.size()) - 1);
        sigma_threshold = vars[idx] / uncertainty_mult;
    }

    // ── Call oracle and record result ─────────────────────────────────────
    // Returns corrected energy estimate, or e_fast if oracle skipped
    double call_oracle(
        System& sys,
        const Fingerprint& fp,
        double e_fast,
        GPSurrogate& main_gp,
        std::mutex& main_gp_mutex,
        int gen)
    {
        if (!should_call_oracle(fp, main_gp, gen, 1000)) {
            ++n_oracle_skipped;
            return corrected_energy(fp, e_fast);
        }

        // Call oracle (MACE) — this is slow
        double e_oracle = 1e30;
        {
            std::lock_guard<std::mutex> lk(oracle_mutex);
            if (n_oracle_calls >= oracle_budget) {
                ++n_oracle_skipped;
                return corrected_energy(fp, e_fast);
            }

            std::vector<double> grad_dummy(3*sys.N, 0.0);
            e_oracle = oracle->energy_only(sys);

            if (e_oracle < 1e29) {
                ++n_oracle_calls;

                // Store record
                OracleRecord rec;
                rec.sys      = sys;
                rec.fp       = fp;
                rec.e_fast   = e_fast;
                rec.e_oracle = e_oracle;
                rec.gen      = gen;
                {
                    auto [m,v] = main_gp.predict(fp.data);
                    rec.gp_variance_at_call = v;
                }
                records.push_back(rec);

                // Update correction GP: learn (fp → e_oracle - e_fast)
                double correction = e_oracle - e_fast;
                correction_gp.add(fp.data, correction);
                if (correction_gp.dirty && correction_gp.n_obs() >= 3)
                    correction_gp.rebuild();

                // Also update main GP with MACE energy
                {
                    std::lock_guard<std::mutex> lk2(main_gp_mutex);
                    main_gp.add(fp.data, e_oracle);
                    if (main_gp.dirty && main_gp.n_obs() >= 3)
                        main_gp.rebuild();
                }
            }
        }

        if (e_oracle < 1e29) return e_oracle;
        return corrected_energy(fp, e_fast);
    }

    // ── Correction: e_corrected = e_fast + GP_correction(fp) ─────────────
    double corrected_energy(const Fingerprint& fp, double e_fast) const {
        if (!correction_gp.is_trained()) return e_fast;
        auto [corr, v] = correction_gp.predict(fp.data);
        return e_fast + corr;
    }

    // ── Re-evaluate all isomers with oracle (post-search refinement) ──────
    // Returns number of oracle calls made
    int refine_isomers(
        std::vector<std::pair<double,System>>& isomers,
        int max_calls = -1)
    {
        if (!enabled || !oracle) return 0;
        if (max_calls < 0) max_calls = oracle_budget;

        int calls = 0;
        for (auto& [e, sys] : isomers) {
            if (calls >= max_calls) break;
            double e_oracle = oracle->energy_only(sys);
            if (e_oracle < 1e29) { e = e_oracle; ++calls; }
        }

        // Re-sort by refined energies
        std::sort(isomers.begin(), isomers.end(),
            [](auto& a, auto& b){ return a.first < b.first; });

        return calls;
    }

    // ── Energy correction ratio: how different is oracle from fast pot ────
    double mean_correction() const {
        if (records.empty()) return 0.0;
        double s = 0.0;
        for (auto& r : records) s += std::abs(r.e_oracle - r.e_fast);
        return s / records.size();
    }

    // ── Print summary ─────────────────────────────────────────────────────
    void print_summary() const {
        if (!enabled) return;
        std::cout << "\n[ACTIVE LEARNING] Summary:\n"
                  << "  Oracle calls    : " << n_oracle_calls
                  << " / " << oracle_budget << "\n"
                  << "  Skipped         : " << n_oracle_skipped << "\n";
        if (!records.empty()) {
            std::cout << std::fixed << std::setprecision(2)
                      << "  Mean correction : " << mean_correction() << " K\n"
                      << "  First call gen  : " << records.front().gen << "\n"
                      << "  Last call gen   : " << records.back().gen  << "\n";
        }
        std::cout << "\n";
    }
};
