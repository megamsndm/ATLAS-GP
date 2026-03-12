#pragma once
/**
 * gupta_potential.hpp — Gupta / RGL many-body embedded-atom potential
 *
 * E = Σᵢ { −[Σⱼ≠ᵢ ξ²ᵢⱼ exp(−2qᵢⱼ(rᵢⱼ/r₀ᵢⱼ − 1))]^½  +  Σⱼ>ᵢ Aᵢⱼ exp(−pᵢⱼ(rᵢⱼ/r₀ᵢⱼ − 1)) }
 *
 * Band energy: Eb = −√(ρᵢ),  ρᵢ = Σⱼ ξ²exp(−2q(r/r₀−1))
 * Repulsion:   Er = Σ A exp(−p(r/r₀−1))
 *
 * Parameters: Cleri & Rosato, PRB 48:22 (1993) — 16 FCC/HCP metals
 * Heteroatomic mixing: geometric mean for ξ,A; arithmetic mean for p,q,r₀
 *
 * Analytical gradient derived from band energy via chain rule.
 */

#include "potential_base.hpp"

struct GuptaElemParams {
    double A, xi, p, q, r0;
};

inline GuptaElemParams get_gupta_params(const std::string& el) {
    // Cleri-Rosato (1993) parameters — FCC metals
    // A[eV], xi[eV], p, q, r0[Å]  → convert A,xi to K: 1 eV = 11604.5 K
    static const std::map<std::string,GuptaElemParams> db = {
    //  el       A[eV]    xi[eV]   p        q       r0[Å]
    {"Ni",  {0.0376,  1.0700,  16.999,  1.189,  2.4911}},
    {"Cu",  {0.0855,  1.2240,  10.960,  2.278,  2.5560}},
    {"Rh",  {0.0629,  1.6600,  18.450,  1.867,  2.6900}},
    {"Pd",  {0.1746,  1.7180,  10.867,  3.742,  2.7485}},
    {"Ag",  {0.1028,  1.1780,  10.928,  3.139,  2.8885}},
    {"Ir",  {0.1156,  2.2890,  16.980,  2.691,  2.7150}},
    {"Pt",  {0.2975,  2.6950,  10.612,  4.004,  2.7747}},
    {"Au",  {0.2061,  1.7900,  10.229,  4.036,  2.8840}},
    {"Al",  {0.1221,  1.3160,   8.612,  2.516,  2.8637}},
    {"Pb",  {0.0980,  0.9140,   9.576,  3.648,  3.5012}},
    {"Fe",  {0.1300,  1.6000,  10.500,  2.600,  2.553 }},  // BCC approx
    {"Co",  {0.0950,  1.4880,  14.300,  1.970,  2.507 }},
    {"Ti",  {0.1530,  1.5920,   8.620,  2.700,  2.900 }},
    {"Mo",  {0.2100,  1.8700,  15.050,  2.500,  2.726 }},
    {"W",   {0.3060,  2.2600,  15.350,  2.290,  2.741 }},
    {"Re",  {0.2200,  2.1000,  16.000,  2.400,  2.760 }},
    };
    auto it = db.find(el);
    if (it != db.end()) return it->second;
    // Fallback: construct approximate Gupta from UFF
    auto& u = get_uff(el);
    return { u.epsilon/11604.5 * 0.5, u.epsilon/11604.5 * 2.0,
             10.0, 2.5, u.sigma * std::pow(2.0,1.0/6.0) };
}

struct GuptaPotential : public Potential {
    struct PairP { double A, xi2, p, q, r0; };
    std::vector<std::vector<PairP>> pp;
    std::vector<double> sigmas;
    static constexpr double eV_to_K = 11604.5221;
    int n_types=0;

    explicit GuptaPotential(const System& sys) {
        n_types = sys.n_species();
        pp.resize(n_types, std::vector<PairP>(n_types));
        sigmas.resize(n_types);
        for (int i=0;i<n_types;++i) {
            auto pi = get_gupta_params(sys.species[i]);
            sigmas[i] = pi.r0;
            for (int j=i;j<n_types;++j) {
                auto pj = get_gupta_params(sys.species[j]);
                PairP p;
                p.A   = std::sqrt(pi.A  * pj.A)  * eV_to_K;
                p.xi2 = std::sqrt(pi.xi * pj.xi) * eV_to_K;
                p.xi2 *= p.xi2;   // store xi² directly
                p.p   = 0.5*(pi.p + pj.p);
                p.q   = 0.5*(pi.q + pj.q);
                p.r0  = 0.5*(pi.r0+ pj.r0);
                pp[i][j] = pp[j][i] = p;
            }
        }
    }

    double get_mean_sigma(const System& sys) const override {
        double s=0; for (double si:sigmas) s+=si; return s/sigmas.size();
    }
    std::string name() const override { return "Gupta-RGL"; }

    double energy_gradient(const System& sys, double* grad) const override {
        const int N=sys.N;
        // Band energy density per atom
        std::vector<double> rho(N, 0.0);
        // Repulsive energy (accumulate directly)
        double Er=0.0;

        for (int i=0;i<N-1;++i) {
            for (int j=i+1;j<N;++j) {
                double dx=sys.pos[3*i]-sys.pos[3*j];
                double dy=sys.pos[3*i+1]-sys.pos[3*j+1];
                double dz=sys.pos[3*i+2]-sys.pos[3*j+2];
                double r=std::sqrt(dx*dx+dy*dy+dz*dz);
                if (r<1e-10) continue;
                const auto& p=pp[sys.type_id[i]][sys.type_id[j]];
                double x=r/p.r0-1.0;
                double ep = std::exp(-p.p*x);
                double eq = std::exp(-2.0*p.q*x);
                Er     += p.A*ep;
                rho[i] += p.xi2*eq;
                rho[j] += p.xi2*eq;
            }
        }

        // Band energy
        double Eb=0.0;
        std::vector<double> inv_sqrt_rho(N,0.0);
        for (int i=0;i<N;++i){
            if (rho[i]>1e-30){
                double sr=std::sqrt(rho[i]);
                Eb -= sr;
                inv_sqrt_rho[i] = 0.5/sr;
            }
        }

        // Gradient
        for (int i=0;i<N-1;++i){
            for (int j=i+1;j<N;++j){
                double dx=sys.pos[3*i]-sys.pos[3*j];
                double dy=sys.pos[3*i+1]-sys.pos[3*j+1];
                double dz=sys.pos[3*i+2]-sys.pos[3*j+2];
                double r2=dx*dx+dy*dy+dz*dz;
                if (r2<1e-20) continue;
                double r=std::sqrt(r2);
                const auto& p=pp[sys.type_id[i]][sys.type_id[j]];
                double x=r/p.r0-1.0;
                // dEr/dr = -p.A*p.p/r0 * exp(-p*x)
                double dep = -p.A*p.p/p.r0 * std::exp(-p.p*x);
                // dEb/drij via chain rule: -0.5/sqrt(rho_i)*drhoi/drij + i<->j
                double deq = -2.0*p.q/p.r0 * p.xi2 * std::exp(-2.0*p.q*x);
                double dEb_dr = -(inv_sqrt_rho[i]+inv_sqrt_rho[j]) * deq;
                double fc = (dep + dEb_dr) / r;
                grad[3*i]  +=fc*dx; grad[3*j]  -=fc*dx;
                grad[3*i+1]+=fc*dy; grad[3*j+1]-=fc*dy;
                grad[3*i+2]+=fc*dz; grad[3*j+2]-=fc*dz;
            }
        }
        return Er + Eb;
    }

    double energy_only(const System& sys) const override {
        const int N=sys.N;
        std::vector<double> rho(N,0.0);
        double Er=0.0;
        for (int i=0;i<N-1;++i){
            for (int j=i+1;j<N;++j){
                double dx=sys.pos[3*i]-sys.pos[3*j];
                double dy=sys.pos[3*i+1]-sys.pos[3*j+1];
                double dz=sys.pos[3*i+2]-sys.pos[3*j+2];
                double r=std::sqrt(dx*dx+dy*dy+dz*dz);
                if (r<1e-10) continue;
                const auto& p=pp[sys.type_id[i]][sys.type_id[j]];
                double x=r/p.r0-1.0;
                Er     += p.A*std::exp(-p.p*x);
                double eq=std::exp(-2.0*p.q*x);
                rho[i] += p.xi2*eq;
                rho[j] += p.xi2*eq;
            }
        }
        double Eb=0.0;
        for (int i=0;i<N;++i) if (rho[i]>1e-30) Eb-=std::sqrt(rho[i]);
        return Er+Eb;
    }
};
