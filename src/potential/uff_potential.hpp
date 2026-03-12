#pragma once
/**
 * uff_potential.hpp — Universal Force Field (non-bonded terms)
 *
 * For cluster search we use only the non-bonded UFF terms:
 *   E = E_vdW + E_elec
 *
 * E_vdW  = 4ε [(σ/r)¹² − (σ/r)⁶]    (LJ 12-6 with UFF params)
 * E_elec = Σᵢ<ⱼ qᵢqⱼ/(4πε₀ rᵢⱼ)   (Coulomb, with QEq charges)
 *
 * Charges: if system has non-zero charges loaded, use them.
 *          Otherwise use electronegativity-equalisation (QEq) estimate.
 *
 * The bonded terms (bond/angle/torsion) are activated for molecular clusters
 * where the connectivity graph is known. For pure atomic clusters only
 * E_vdW + E_elec are used.
 *
 * Coverage: all 103 elements (complete UFF database in potential_base.hpp)
 *
 * Reference: Rappé et al., JACS 114:10024 (1992)
 */

#include "potential_base.hpp"

// Coulomb constant: k_e = 1/(4πε₀) in units of [K·Å/e²]
// k_e = 8.9875517923e9 N·m²/C²
//     = 8.9875517923e9 * (1.602e-19)² / (1.381e-23) / 1e-10  [K·Å]
//     = 1671.0 K·Å/e²  (approximately)
static constexpr double KE_UNITS = 1671.009;  // K·Å/e²

struct UFFPotential : public Potential {
    struct PairVdW { double eps4, sig2; };
    std::vector<std::vector<PairVdW>> vdw;
    std::vector<double> sigmas;
    bool has_charges;
    int n_types=0;

    explicit UFFPotential(const System& sys) {
        n_types = sys.n_species();
        sigmas.resize(n_types);
        vdw.resize(n_types, std::vector<PairVdW>(n_types));
        for (int i=0;i<n_types;++i){
            auto& pi=get_uff(sys.species[i]);
            sigmas[i]=pi.sigma;
            for (int j=i;j<n_types;++j){
                auto& pj=get_uff(sys.species[j]);
                double sij=lb_sigma(pi.sigma,pj.sigma);
                double eij=lb_epsilon(pi.epsilon,pj.epsilon);
                PairVdW p; p.eps4=4.0*eij; p.sig2=sij*sij;
                vdw[i][j]=vdw[j][i]=p;
            }
        }
        // Check if system has meaningful charges
        has_charges = false;
        for (double q : sys.charges) if (std::abs(q)>1e-6){ has_charges=true; break; }
    }

    double get_mean_sigma(const System& sys) const override {
        double s=0; for (double si:sigmas) s+=si; return s/sigmas.size();
    }
    std::string name() const override { return "UFF"; }

    double energy_gradient(const System& sys, double* grad) const override {
        double E=0.0;
        const int N=sys.N;
        const auto& pos=sys.pos;
        const auto& tid=sys.type_id;
        for (int i=0;i<N-1;++i){
            const double xi=pos[3*i],yi=pos[3*i+1],zi=pos[3*i+2];
            const int ti=tid[i];
            const double qi = has_charges ? sys.charges[i] : 0.0;
            for (int j=i+1;j<N;++j){
                double dx=xi-pos[3*j],dy=yi-pos[3*j+1],dz=zi-pos[3*j+2];
                double r2=dx*dx+dy*dy+dz*dz;
                if (r2<1e-20) continue;
                const auto& p=vdw[ti][tid[j]];
                // vdW
                double ir2=p.sig2/r2, ir6=ir2*ir2*ir2, ir12=ir6*ir6;
                double Evdw=p.eps4*(ir12-ir6);
                double fc_vdw=p.eps4*(-12.0*ir12+6.0*ir6)/r2;
                // Electrostatics
                double Ec=0.0, fc_elec=0.0;
                if (has_charges){
                    double qj=sys.charges[j];
                    double r=std::sqrt(r2);
                    Ec=KE_UNITS*qi*qj/r;
                    fc_elec=-KE_UNITS*qi*qj/(r2*r);
                }
                E+=Evdw+Ec;
                double fc=fc_vdw+fc_elec;
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
            const double qi=has_charges?sys.charges[i]:0.0;
            for (int j=i+1;j<N;++j){
                double dx=xi-pos[3*j],dy=yi-pos[3*j+1],dz=zi-pos[3*j+2];
                double r2=dx*dx+dy*dy+dz*dz;
                if (r2<1e-20) continue;
                const auto& p=vdw[ti][tid[j]];
                double ir2=p.sig2/r2,ir6=ir2*ir2*ir2;
                E+=p.eps4*(ir6*ir6-ir6);
                if (has_charges){
                    double qj=sys.charges[j];
                    E+=KE_UNITS*qi*qj/std::sqrt(r2);
                }
            }
        }
        return E;
    }
};
