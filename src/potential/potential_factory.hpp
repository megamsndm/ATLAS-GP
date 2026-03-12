#pragma once
/**
 * potential_factory.hpp
 * Creates the right Potential object from a name string.
 * Also implements the cascade strategy: cheap potential for search,
 * MACE for final refinement.
 */

#include "potential_base.hpp"
#include "lj_potential.hpp"
#include "morse_potential.hpp"
#include "gupta_potential.hpp"
#include "uff_potential.hpp"
#include "mace_potential.hpp"
#include <memory>
#include <string>
#include <algorithm>
#include <set>
#include <iostream>

// Normalise potential name to uppercase
static std::string norm_pot(const std::string& s) {
    std::string r=s;
    for (char& c:r) c=(char)std::toupper((unsigned char)c);
    return r;
}

inline std::unique_ptr<Potential> make_potential(
    const std::string& name,
    const System& sys,
    const std::string& mace_script = "")
{
    std::string n = norm_pot(name);

    if (n=="LJ" || n=="LJ126" || n=="LJ-12-6")
        return std::make_unique<LJPotential>(sys);

    if (n=="MORSE")
        return std::make_unique<MorsePotential>(sys);

    if (n=="GUPTA" || n=="RGL" || n=="GUPTA-RGL")
        return std::make_unique<GuptaPotential>(sys);

    if (n=="UFF" || n=="UNIVERSAL")
        return std::make_unique<UFFPotential>(sys);

    if (n=="MACE" || n=="MACE-MP-0" || n=="MACE-MP") {
        std::string script = mace_script.empty() ? "mace_server.py" : mace_script;
        auto pot = std::make_unique<MACEPotential>(sys, script);
        if (!pot->is_available()) {
            std::cerr << "[ATLAS] MACE not available: " << pot->get_fallback_reason()
                      << "\n[ATLAS] Falling back to UFF.\n";
            return std::make_unique<UFFPotential>(sys);
        }
        std::cout << "[ATLAS] Using MACE-MP-0 (near-DFT accuracy, 89 elements)\n";
        return pot;
    }

    // AUTO: pick best available for this system
    if (n=="AUTO") {
        // Try MACE first
        if (!mace_script.empty() || [](){
            FILE*f=fopen("mace_server.py","r");
            if(f){fclose(f);return true;}return false;
        }()) {
            auto pot = std::make_unique<MACEPotential>(sys, mace_script.empty()?"mace_server.py":mace_script);
            if (pot->is_available()) {
                std::cout << "[ATLAS] AUTO selected: MACE-MP-0\n";
                return pot;
            }
        }
        // Pick based on elements
        bool all_metals = true;
        std::set<std::string> gupta_supported = {
            "Ni","Cu","Rh","Pd","Ag","Ir","Pt","Au","Al","Pb","Fe","Co","Ti","Mo","W","Re"
        };
        for (auto& el : sys.species)
            if (gupta_supported.find(el)==gupta_supported.end()) { all_metals=false; break; }
        if (all_metals && sys.n_species()<=2) {
            std::cout << "[ATLAS] AUTO selected: Gupta-RGL\n";
            return std::make_unique<GuptaPotential>(sys);
        }
        // Noble gases → LJ
        std::set<std::string> noble = {"He","Ne","Ar","Kr","Xe","Rn"};
        bool all_noble = true;
        for (auto& el : sys.species)
            if (noble.find(el)==noble.end()) { all_noble=false; break; }
        if (all_noble) {
            std::cout << "[ATLAS] AUTO selected: LJ-12-6\n";
            return std::make_unique<LJPotential>(sys);
        }
        // Default: UFF (covers all 103 elements)
        std::cout << "[ATLAS] AUTO selected: UFF\n";
        return std::make_unique<UFFPotential>(sys);
    }

    throw std::runtime_error("Unknown potential: '"+name+"'. Options: LJ, Morse, Gupta, UFF, MACE, AUTO");
}

// ─────────────────────────────────────────────────────────────────────────────
// Cascade potential: use cheap_pot for search, mace_pot for final refinement
// ─────────────────────────────────────────────────────────────────────────────
struct CascadePotential {
    std::unique_ptr<Potential> cheap;
    std::unique_ptr<Potential> accurate;
    bool has_accurate;

    CascadePotential(std::unique_ptr<Potential> c, std::unique_ptr<Potential> a, bool ha)
        : cheap(std::move(c)), accurate(std::move(a)), has_accurate(ha) {}

    // Use cheap potential during search
    double search_energy_gradient(const System& sys, double* grad) const {
        return cheap->energy_gradient(sys, grad);
    }
    double search_energy(const System& sys) const {
        return cheap->energy_only(sys);
    }

    // Use accurate potential for final refinement
    double refine_energy_gradient(const System& sys, double* grad) const {
        if (has_accurate) return accurate->energy_gradient(sys, grad);
        return cheap->energy_gradient(sys, grad);
    }
    double refine_energy(const System& sys) const {
        if (has_accurate) return accurate->energy_only(sys);
        return cheap->energy_only(sys);
    }

    std::string search_name() const { return cheap->name(); }
    std::string refine_name() const { return has_accurate ? accurate->name() : cheap->name(); }
    double get_mean_sigma(const System& sys) const { return cheap->get_mean_sigma(sys); }
};

inline CascadePotential make_cascade(
    const std::string& search_pot,
    const System& sys,
    const std::string& mace_script = "")
{
    auto cheap = make_potential(search_pot, sys, mace_script);

    // Try to make MACE for refinement
    std::string script = mace_script.empty() ? "mace_server.py" : mace_script;
    FILE* f = fopen(script.c_str(), "r");
    bool mace_available = (f != nullptr);
    if (f) fclose(f);

    if (mace_available && norm_pot(search_pot) != "MACE") {
        auto mace = std::make_unique<MACEPotential>(sys, script);
        if (mace->is_available()) {
            std::cout << "[ATLAS] Cascade: search=" << cheap->name()
                      << " + refine=MACE-MP-0\n";
            return CascadePotential(std::move(cheap), std::move(mace), true);
        }
    }
    std::cout << "[ATLAS] Single potential: " << cheap->name() << "\n";
    return CascadePotential(std::move(cheap),
                            make_potential(search_pot, sys, mace_script), false);
}
