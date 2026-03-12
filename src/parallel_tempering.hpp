#pragma once
/**
 * parallel_tempering.hpp — Replica Exchange Monte Carlo
 *
 * SESSION 3 INNOVATION: Parallel tempering (PT) is the single most
 * impactful improvement for rugged energy landscapes with deep funnels.
 *
 * Problem with independent basin-hopping runs:
 *   Each run has ONE temperature. Hot runs explore but don't converge.
 *   Cold runs converge fast but get permanently trapped in funnels.
 *   They cannot learn from each other.
 *
 * Parallel tempering solution:
 *   Run N_replica walkers simultaneously at different temperatures:
 *     T_1 < T_2 < T_3 < ... < T_R   (geometric ladder)
 *   Periodically attempt to SWAP configurations between adjacent replicas:
 *     Accept swap with probability min(1, exp(-(β_i - β_j)(E_j - E_i)))
 *   Result:
 *     - Hot replicas explore widely, find new funnels
 *     - Cold replicas refine within funnels precisely
 *     - Swaps allow cold replicas to escape funnels via hot replicas
 *     - Provably converges to Boltzmann distribution
 *
 * Implementation:
 *   - PTEngine manages R replicas sharing one GPSurrogate and MotifLibrary
 *   - Each replica is an ATLASWorker with its own temperature
 *   - Exchange attempts every PT_INTERVAL generations
 *   - Uses even/odd swap scheme to avoid deadlocks in parallel execution
 *   - Exchange statistics tracked to auto-tune temperature ladder
 *     (target acceptance rate: 20-40% → adjust T spacing automatically)
 *
 * Temperature ladder:
 *   T_i = T_min * (T_max/T_min)^(i/(R-1))   geometric spacing
 *   Auto-tuned: if acceptance rate < 15% → compress ladder
 *               if acceptance rate > 50% → expand ladder
 */

#include "system.hpp"
#include "potential/potential_base.hpp"
#include "fingerprint.hpp"
#include "gp_surrogate.hpp"
#include "atlas_worker.hpp"
#include <vector>
#include <mutex>
#include <atomic>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

// ─────────────────────────────────────────────────────────────────────────────
// PTConfig: parallel tempering parameters
// ─────────────────────────────────────────────────────────────────────────────
struct PTConfig {
    int    n_replicas        = 4;      // number of temperature replicas
    double T_min             = 0.3;    // coldest replica temperature (reduced)
    double T_max             = 5.0;    // hottest replica temperature
    int    exchange_interval = 20;     // generations between swap attempts
    bool   auto_tune_temps   = true;   // auto-adjust T ladder
    double target_acc_rate   = 0.25;   // target swap acceptance rate
    int    tune_interval     = 100;    // generations between T adjustments
};

// ─────────────────────────────────────────────────────────────────────────────
// Replica state: current system, energy, temperature
// ─────────────────────────────────────────────────────────────────────────────
struct Replica {
    System sys;
    double energy;
    double T;           // current temperature (reduced units)
    int    id;
    int    n_swaps_attempted = 0;
    int    n_swaps_accepted  = 0;

    double swap_acceptance() const {
        return n_swaps_attempted > 0
            ? (double)n_swaps_accepted / n_swaps_attempted
            : 0.0;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// PTEngine: manages replica exchange across all walkers
// ─────────────────────────────────────────────────────────────────────────────
class PTEngine {
public:
    PTConfig               cfg;
    std::vector<double>    temps;       // temperature ladder
    std::vector<Replica>   replicas;    // current replica states
    std::mutex             swap_mutex;
    std::mt19937_64        rng;
    int                    gen_count = 0;

    PTEngine() : rng(42) {}
    explicit PTEngine(const PTConfig& c, uint64_t seed = 42)
        : cfg(c), rng(seed)
    {
        _build_temp_ladder();
    }

    // ── Initialise replicas from a set of systems ─────────────────────────
    void init(const std::vector<System>& systems,
              const std::vector<double>& energies)
    {
        int R = cfg.n_replicas;
        replicas.resize(R);
        for (int i=0;i<R;++i) {
            replicas[i].id  = i;
            replicas[i].T   = temps[i];
            int src = i % (int)systems.size();
            replicas[i].sys    = systems[src];
            replicas[i].energy = energies[src];
        }
    }

    // ── Attempt swap between adjacent replicas (even/odd scheme) ─────────
    // Call after every PT_INTERVAL generations
    // Returns number of accepted swaps
    int attempt_swaps(int parity) {
        // parity=0: swap (0,1),(2,3),...
        // parity=1: swap (1,2),(3,4),...
        int R = (int)replicas.size();
        int n_accepted = 0;
        std::uniform_real_distribution<double> U(0.0, 1.0);

        std::lock_guard<std::mutex> lk(swap_mutex);
        for (int i = parity; i < R-1; i += 2) {
            Replica& ri = replicas[i];
            Replica& rj = replicas[i+1];

            ++ri.n_swaps_attempted;
            ++rj.n_swaps_attempted;

            // Metropolis-Hastings swap criterion:
            // P(swap) = min(1, exp(-(β_i - β_j)(E_j - E_i)))
            double beta_i = 1.0 / (ri.T + 1e-15);
            double beta_j = 1.0 / (rj.T + 1e-15);
            double dE     = rj.energy - ri.energy;
            double log_p  = -(beta_i - beta_j) * dE;

            bool accept = (log_p >= 0.0) || (std::log(U(rng)) < log_p);

            if (accept) {
                // Swap configurations but keep temperatures
                std::swap(ri.sys,    rj.sys);
                std::swap(ri.energy, rj.energy);
                ++ri.n_swaps_accepted;
                ++rj.n_swaps_accepted;
                ++n_accepted;
            }
        }
        ++gen_count;
        return n_accepted;
    }

    // ── Auto-tune temperature ladder ──────────────────────────────────────
    // Adjust T spacing so swap acceptance stays near target (20-40%)
    void auto_tune() {
        if (!cfg.auto_tune_temps) return;
        int R = (int)replicas.size();
        if (R < 2) return;

        bool any_change = false;
        for (int i=0; i<R-1; ++i) {
            double acc = replicas[i].swap_acceptance();
            if (replicas[i].n_swaps_attempted < 20) continue;

            if (acc < 0.12) {
                // Too few swaps: compress ladder (bring temps closer)
                double mid = 0.5*(temps[i] + temps[i+1]);
                temps[i]   = mid * 0.9 + temps[i]   * 0.1;
                temps[i+1] = mid * 0.9 + temps[i+1] * 0.1;
                any_change = true;
            } else if (acc > 0.55) {
                // Too many swaps: expand ladder (increase T spacing)
                double gap = temps[i+1] - temps[i];
                temps[i]   -= gap * 0.1;
                temps[i+1] += gap * 0.1;
                temps[i]   = std::max(cfg.T_min * 0.5, temps[i]);
                temps[i+1] = std::min(cfg.T_max * 2.0, temps[i+1]);
                any_change = true;
            }
        }

        if (any_change) {
            std::sort(temps.begin(), temps.end());
            for (int i=0;i<R;++i) replicas[i].T = temps[i];
            // Reset acceptance counters
            for (auto& r : replicas) {
                r.n_swaps_attempted = 0;
                r.n_swaps_accepted  = 0;
            }
        }
    }

    // ── Get coldest replica ───────────────────────────────────────────────
    const Replica& coldest() const {
        return *std::min_element(replicas.begin(), replicas.end(),
            [](const Replica& a, const Replica& b){ return a.T < b.T; });
    }

    // ── Get best (lowest energy) replica ─────────────────────────────────
    const Replica& best() const {
        return *std::min_element(replicas.begin(), replicas.end(),
            [](const Replica& a, const Replica& b){
                return a.energy < b.energy; });
    }

    // ── Print PT statistics ───────────────────────────────────────────────
    void print_stats() const {
        std::cout << "[PT] Replica exchange statistics:\n";
        std::cout << std::fixed << std::setprecision(3);
        for (auto& r : replicas) {
            std::cout << "  Replica " << r.id
                      << "  T=" << std::setw(6) << r.T
                      << "  E=" << std::setw(12) << std::setprecision(6) << r.energy
                      << "  swap_acc=" << std::setw(5) << std::setprecision(3)
                      << r.swap_acceptance() << "\n";
        }
        std::cout << "\n";
    }

    int n_replicas() const { return (int)replicas.size(); }

private:
    void _build_temp_ladder() {
        int R = cfg.n_replicas;
        temps.resize(R);
        if (R == 1) { temps[0]=cfg.T_min; return; }
        for (int i=0;i<R;++i) {
            double frac = (double)i/(R-1);
            temps[i] = cfg.T_min * std::pow(cfg.T_max/cfg.T_min, frac);
        }
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// integrate_pt_into_config: add PT parameters to ATLASConfig
// Called by main() when --pt flag is given
// ─────────────────────────────────────────────────────────────────────────────
inline PTConfig make_pt_config(const ATLASConfig& cfg, int n_replicas=4) {
    PTConfig pt;
    pt.n_replicas       = n_replicas;
    pt.T_min            = cfg.bh_temp_min;
    pt.T_max            = cfg.bh_temp * 8.0;
    pt.exchange_interval= 20;
    pt.auto_tune_temps  = true;
    return pt;
}
