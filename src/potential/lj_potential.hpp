#pragma once
/**
 * lj_potential.hpp — Lennard-Jones 12-6
 *
 * Homoatomic:   E = 4ε Σᵢ<ⱼ [(σ/rᵢⱼ)¹² − (σ/rᵢⱼ)⁶]
 * Heteroatomic: εᵢⱼ = √(εᵢεⱼ),  σᵢⱼ = (σᵢ+σⱼ)/2  (Lorentz-Berthelot)
 *
 * Parameters from UFF database (all 103 elements supported).
 * Works for noble gases, simple molecular pairs, approximate metals.
 */

#include "potential_base.hpp"

struct LJPotential : public Potential {
    // Pre-computed pair parameters (upper triangle only)
    struct PairParams { double eps4, sig2; };
    std::vector<std::vector<PairParams>> pair_params; // [type_i][type_j]
    std::vector<double> sigmas;  // per type, for mean_sigma
    int n_types = 0;

    explicit LJPotential(const System& sys) {
        n_types = sys.n_species();
        sigmas.resize(n_types);
        pair_params.resize(n_types, std::vector<PairParams>(n_types));
        for (int i=0;i<n_types;++i) {
            auto& pi = get_uff(sys.species[i]);
            sigmas[i] = pi.sigma;
            for (int j=i;j<n_types;++j) {
                auto& pj = get_uff(sys.species[j]);
                double sij = lb_sigma  (pi.sigma,   pj.sigma);
                double eij = lb_epsilon(pi.epsilon, pj.epsilon);
                PairParams p;
                p.eps4 = 4.0 * eij;
                p.sig2 = sij * sij;
                pair_params[i][j] = pair_params[j][i] = p;
            }
        }
    }

    double get_mean_sigma(const System& sys) const override {
        double s=0; for (double si : sigmas) s+=si; return s/sigmas.size();
    }

    std::string name() const override { return "LJ-12-6"; }

    double energy_gradient(const System& sys, double* grad) const override {
        double E = 0.0;
        const int N = sys.N;
        const auto& pos = sys.pos;
        const auto& tid = sys.type_id;
        for (int i=0;i<N-1;++i) {
            const double xi=pos[3*i], yi=pos[3*i+1], zi=pos[3*i+2];
            const int ti = tid[i];
            for (int j=i+1;j<N;++j) {
                const double dx=xi-pos[3*j], dy=yi-pos[3*j+1], dz=zi-pos[3*j+2];
                const double r2 = dx*dx+dy*dy+dz*dz;
                if (r2<1e-20) continue;
                const auto& p = pair_params[ti][tid[j]];
                const double ir2  = p.sig2/r2;
                const double ir6  = ir2*ir2*ir2;
                const double ir12 = ir6*ir6;
                E += p.eps4*(ir12-ir6);
                const double fc = p.eps4*(-12.0*ir12+6.0*ir6)/r2;
                grad[3*i]  +=fc*dx; grad[3*j]  -=fc*dx;
                grad[3*i+1]+=fc*dy; grad[3*j+1]-=fc*dy;
                grad[3*i+2]+=fc*dz; grad[3*j+2]-=fc*dz;
            }
        }
        return E;
    }

    double energy_only(const System& sys) const override {
        double E=0.0;
        const int N=sys.N;
        const auto& pos=sys.pos;
        const auto& tid=sys.type_id;
        for (int i=0;i<N-1;++i){
            const double xi=pos[3*i],yi=pos[3*i+1],zi=pos[3*i+2];
            const int ti=tid[i];
            for (int j=i+1;j<N;++j){
                const double dx=xi-pos[3*j],dy=yi-pos[3*j+1],dz=zi-pos[3*j+2];
                const double r2=dx*dx+dy*dy+dz*dz;
                if (r2<1e-20) continue;
                const auto& p=pair_params[ti][tid[j]];
                const double ir2=p.sig2/r2, ir6=ir2*ir2*ir2;
                E+=p.eps4*(ir6*ir6-ir6);
            }
        }
        return E;
    }
};
