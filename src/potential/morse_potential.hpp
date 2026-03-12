#pragma once
/**
 * morse_potential.hpp — Morse pair potential
 *
 * V(r) = D_e [1 − exp(−a(r−r₀))]²  − D_e
 *
 * Parameters: D_e [K], a [Å⁻¹], r₀ [Å]
 * For heteroatomic pairs: geometric mean of D_e, arithmetic mean of r₀ and a
 *
 * Reference: Girifalco & Weizer (1959), Phys Rev 114:687
 *            + fitted values for transition metals
 */

#include "potential_base.hpp"

struct MorseParams {
    double De;   // well depth [K]
    double a;    // width parameter [Å⁻¹]
    double r0;   // equilibrium distance [Å]
};

inline MorseParams get_morse_params(const std::string& el) {
    // Parameterised for homoatomic pairs
    static const std::map<std::string,MorseParams> db = {
    //  el       De[K]     a[Å⁻¹]  r0[Å]
    {"H",   { 38316.0,   1.9436,  0.741}},
    {"Li",  {  4958.0,   0.8029,  2.673}},
    {"Na",  {  2643.0,   0.7178,  3.079}},
    {"K",   {  1718.0,   0.5530,  3.923}},
    {"Cu",  { 38148.0,   1.3588,  2.866}},   // FCC metal
    {"Ag",  { 26042.0,   1.3690,  2.889}},
    {"Au",  { 40120.0,   1.6463,  2.884}},
    {"Ni",  { 49530.0,   1.4199,  2.780}},
    {"Pd",  { 43680.0,   1.3413,  2.751}},
    {"Pt",  { 55940.0,   1.5571,  2.775}},
    {"Fe",  { 48640.0,   1.3885,  2.845}},
    {"Co",  { 49940.0,   1.4317,  2.770}},
    {"Al",  { 31500.0,   1.1646,  2.863}},
    {"Pb",  { 14450.0,   1.1836,  3.500}},
    {"Ar",  {   143.2,   1.7100,  3.757}},   // noble gas dimer
    {"Kr",  {   200.9,   1.6060,  4.012}},
    {"Xe",  {   282.3,   1.4100,  4.362}},
    {"N",   { 113690.0,  2.6890,  1.098}},   // N2 dimer
    {"O",   { 60280.0,   2.6560,  1.208}},   // O2 dimer
    {"C",   { 84140.0,   1.9600,  1.543}},   // C-C
    };
    auto it = db.find(el);
    if (it != db.end()) return it->second;
    // Fallback: estimate from UFF epsilon/sigma
    auto& u = get_uff(el);
    return { u.epsilon * 4.0, 1.4 / u.sigma, u.sigma * std::pow(2.0, 1.0/6.0) };
}

struct MorsePotential : public Potential {
    struct PairP { double De, a, r0; };
    std::vector<std::vector<PairP>> pair_params;
    std::vector<double> sigmas;
    int n_types=0;

    explicit MorsePotential(const System& sys) {
        n_types = sys.n_species();
        sigmas.resize(n_types);
        pair_params.resize(n_types, std::vector<PairP>(n_types));
        for (int i=0;i<n_types;++i) {
            auto pi = get_morse_params(sys.species[i]);
            sigmas[i] = pi.r0 / std::pow(2.0, 1.0/6.0);
            for (int j=i;j<n_types;++j) {
                auto pj = get_morse_params(sys.species[j]);
                PairP p;
                p.De = std::sqrt(pi.De * pj.De);           // geometric mean
                p.a  = 0.5*(pi.a  + pj.a);                 // arithmetic mean
                p.r0 = 0.5*(pi.r0 + pj.r0);                // arithmetic mean
                pair_params[i][j] = pair_params[j][i] = p;
            }
        }
    }

    double get_mean_sigma(const System& sys) const override {
        double s=0; for (double si:sigmas) s+=si; return s/sigmas.size();
    }
    std::string name() const override { return "Morse"; }

    double energy_gradient(const System& sys, double* grad) const override {
        double E=0.0;
        const int N=sys.N;
        for (int i=0;i<N-1;++i){
            for (int j=i+1;j<N;++j){
                double dx=sys.pos[3*i]-sys.pos[3*j];
                double dy=sys.pos[3*i+1]-sys.pos[3*j+1];
                double dz=sys.pos[3*i+2]-sys.pos[3*j+2];
                double r=std::sqrt(dx*dx+dy*dy+dz*dz);
                if (r<1e-10) continue;
                const auto& p=pair_params[sys.type_id[i]][sys.type_id[j]];
                double ex=std::exp(-p.a*(r-p.r0));
                double v=(1.0-ex)*(1.0-ex)-1.0;
                E+=p.De*v;
                // dV/dr = 2*De*a*ex*(1-ex)
                double dvdr=2.0*p.De*p.a*ex*(1.0-ex)/r;
                grad[3*i]  +=dvdr*dx; grad[3*j]  -=dvdr*dx;
                grad[3*i+1]+=dvdr*dy; grad[3*j+1]-=dvdr*dy;
                grad[3*i+2]+=dvdr*dz; grad[3*j+2]-=dvdr*dz;
            }
        }
        return E;
    }

    double energy_only(const System& sys) const override {
        double E=0.0;
        const int N=sys.N;
        for (int i=0;i<N-1;++i)
            for (int j=i+1;j<N;++j){
                double dx=sys.pos[3*i]-sys.pos[3*j];
                double dy=sys.pos[3*i+1]-sys.pos[3*j+1];
                double dz=sys.pos[3*i+2]-sys.pos[3*j+2];
                double r=std::sqrt(dx*dx+dy*dy+dz*dz);
                if (r<1e-10) continue;
                const auto& p=pair_params[sys.type_id[i]][sys.type_id[j]];
                double ex=std::exp(-p.a*(r-p.r0));
                E+=p.De*((1.0-ex)*(1.0-ex)-1.0);
            }
        return E;
    }
};
