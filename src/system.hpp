#pragma once
/**
 * system.hpp — Universal cluster descriptor
 *
 * Supports:
 *   - Homoatomic clusters:    Ar13, Ni55, Au38
 *   - Heteroatomic clusters:  Au32Ag6, Ni12Cu1, Fe12O18
 *   - Molecular clusters:     (H2O)13, (NH3)6, benzene·water
 *   - Mixed metal/molecular:  Au13·(H2O)6
 *
 * Design:
 *   All algorithms (search, GP, filter) operate on System&.
 *   The potential layer reads type_id[] to get per-pair parameters.
 *   Molecular mode adds orientation quaternions and a bond graph.
 */

#include <vector>
#include <string>
#include <array>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <map>
#include <iostream>

// ─────────────────────────────────────────────────────────────────────────────
// Quaternion for molecular orientation (w, x, y, z)
// ─────────────────────────────────────────────────────────────────────────────
struct Quat {
    double w=1, x=0, y=0, z=0;
    Quat() = default;
    Quat(double w, double x, double y, double z) : w(w),x(x),y(y),z(z) {}

    Quat operator*(const Quat& q) const {
        return { w*q.w - x*q.x - y*q.y - z*q.z,
                 w*q.x + x*q.w + y*q.z - z*q.y,
                 w*q.y - x*q.z + y*q.w + z*q.x,
                 w*q.z + x*q.y - y*q.x + z*q.w };
    }
    void normalise() {
        double n = std::sqrt(w*w+x*x+y*y+z*z);
        if (n > 1e-15) { w/=n; x/=n; y/=n; z/=n; }
    }
    // Rotate vector v by this quaternion
    std::array<double,3> rotate(double vx, double vy, double vz) const {
        // v' = q v q*
        double t2=x*vy-y*vx, t1=w*vz+y*vx-x*vy;
        // Rodrigues form (faster)
        double rx = vx + 2.0*(w*(vy*z - vz*y) + x*(x*vx+y*vy+z*vz) - vx*(x*x+y*y+z*z));
        double ry = vy + 2.0*(w*(vz*x - vx*z) + y*(x*vx+y*vy+z*vz) - vy*(x*x+y*y+z*z));
        double rz = vz + 2.0*(w*(vx*y - vy*x) + z*(x*vx+y*vy+z*vz) - vz*(x*x+y*y+z*z));
        (void)t1; (void)t2;
        return {rx, ry, rz};
    }
    static Quat from_axis_angle(double ax, double ay, double az, double angle) {
        double s = std::sin(angle*0.5);
        double n = std::sqrt(ax*ax+ay*ay+az*az);
        if (n < 1e-15) return Quat{};
        return { std::cos(angle*0.5), ax/n*s, ay/n*s, az/n*s };
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Atom site within a molecule template
// ─────────────────────────────────────────────────────────────────────────────
struct MolSite {
    std::string element;
    double rx, ry, rz;   // position relative to molecule CoM, in Å
    double charge;        // partial charge (e)
};

// ─────────────────────────────────────────────────────────────────────────────
// Molecule template (rigid body)
// ─────────────────────────────────────────────────────────────────────────────
struct MolTemplate {
    std::string name;              // e.g. "H2O", "NH3", "benzene"
    std::vector<MolSite> sites;    // atom positions in body frame

    // Total atoms in one molecule
    int n_atoms() const { return (int)sites.size(); }

    // Moment of inertia (for orientation moves)
    std::array<double,3> inertia_eigenvalues() const {
        double Ixx=0,Iyy=0,Izz=0;
        for (auto& s : sites) {
            Ixx += s.ry*s.ry + s.rz*s.rz;
            Iyy += s.rx*s.rx + s.rz*s.rz;
            Izz += s.rx*s.rx + s.ry*s.ry;
        }
        return {Ixx, Iyy, Izz};
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// System — the universal cluster representation
// ─────────────────────────────────────────────────────────────────────────────
struct System {
    // ── Atomic degrees of freedom ──────────────────────────────────────────
    int                  N;           // total atoms (all species combined)
    std::vector<double>  pos;         // flat [3N]: x0,y0,z0, x1,y1,z1, ...
    std::vector<int>     type_id;     // [N]: index into species[] per atom
    std::vector<std::string> species; // unique element names, indexed by type_id
    std::vector<double>  charges;     // [N]: partial charges (0 if not set)

    // ── Molecular degrees of freedom (empty for pure atomic) ──────────────
    bool                 is_molecular = false;
    std::vector<double>  com;         // [3*N_mol]: centre of mass per molecule
    std::vector<Quat>    orient;      // [N_mol]:   orientation quaternion
    std::vector<int>     mol_id;      // [N]: which molecule does atom i belong to
    std::vector<int>     mol_type;    // [N_mol]: molecule template index
    std::vector<MolTemplate> mol_templates; // molecule type definitions

    // ── Geometry ────────────────────────────────────────────────────────────
    double               R;           // confinement sphere radius [Å]
    double               mean_sigma;  // mean σ for sphere scaling / bond cutoffs

    // ── Constructors ────────────────────────────────────────────────────────
    System() : N(0), R(0), mean_sigma(1.0) {}

    // Build atomic system from species list + counts
    // e.g. build_atomic({{"Au",32},{"Ag",6}})
    static System build_atomic(
        const std::vector<std::pair<std::string,int>>& composition,
        double sphere_factor = 1.3,
        double mean_sig = 2.5)
    {
        System s;
        s.is_molecular = false;
        s.N = 0;
        for (auto& [el,n] : composition) {
            int tid = (int)s.species.size();
            s.species.push_back(el);
            for (int i=0;i<n;++i) {
                s.type_id.push_back(tid);
                s.N++;
            }
        }
        s.pos.assign(3*s.N, 0.0);
        s.charges.assign(s.N, 0.0);
        s.mol_id.assign(s.N, -1);
        s.mean_sigma = mean_sig;
        s.R = sphere_factor * mean_sig * std::cbrt((double)s.N);
        return s;
    }

    // ── Accessors ───────────────────────────────────────────────────────────
    double  px(int i) const { return pos[3*i];   }
    double  py(int i) const { return pos[3*i+1]; }
    double  pz(int i) const { return pos[3*i+2]; }
    double& px(int i)       { return pos[3*i];   }
    double& py(int i)       { return pos[3*i+1]; }
    double& pz(int i)       { return pos[3*i+2]; }

    int n_species() const { return (int)species.size(); }
    int count(int tid) const {
        return (int)std::count(type_id.begin(), type_id.end(), tid);
    }

    // ── Geometry helpers ────────────────────────────────────────────────────
    void centre() {
        double cx=0,cy=0,cz=0;
        for (int i=0;i<N;++i){cx+=pos[3*i];cy+=pos[3*i+1];cz+=pos[3*i+2];}
        cx/=N; cy/=N; cz/=N;
        for (int i=0;i<N;++i){pos[3*i]-=cx;pos[3*i+1]-=cy;pos[3*i+2]-=cz;}
    }

    void project_sphere() {
        for (int i=0;i<N;++i){
            double r2=pos[3*i]*pos[3*i]+pos[3*i+1]*pos[3*i+1]+pos[3*i+2]*pos[3*i+2];
            if (r2>R*R){ double s=R/std::sqrt(r2); pos[3*i]*=s;pos[3*i+1]*=s;pos[3*i+2]*=s; }
        }
    }

    double rij2(int i, int j) const {
        double dx=pos[3*i]-pos[3*j], dy=pos[3*i+1]-pos[3*j+1], dz=pos[3*i+2]-pos[3*j+2];
        return dx*dx+dy*dy+dz*dz;
    }
    double rij(int i, int j) const { return std::sqrt(rij2(i,j)); }

    // ── Deep copy ───────────────────────────────────────────────────────────
    System clone() const { return *this; }

    // ── Composition string for display ──────────────────────────────────────
    std::string formula() const {
        std::string f;
        for (int t=0;t<(int)species.size();++t){
            f += species[t];
            int n = count(t);
            if (n>1) f += std::to_string(n);
        }
        return f;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Formula parser: "Au32Ag6Cu5" → [{"Au",32},{"Ag",6},{"Cu",5}]
//                 "Ni13"       → [{"Ni",13}]
//                 "Ar"         → [{"Ar",1}]   (implicit count 1)
// ─────────────────────────────────────────────────────────────────────────────
inline std::vector<std::pair<std::string,int>> parse_formula(const std::string& s) {
    if (s.empty()) throw std::runtime_error("Empty formula");
    std::vector<std::pair<std::string,int>> result;
    size_t i = 0;
    while (i < s.size()) {
        if (!std::isupper((unsigned char)s[i]))
            throw std::runtime_error("Expected uppercase letter at position "+std::to_string(i)+" in '"+s+"'");
        std::string elem;
        elem += s[i++];
        while (i < s.size() && std::islower((unsigned char)s[i]))
            elem += s[i++];
        std::string numstr;
        while (i < s.size() && std::isdigit((unsigned char)s[i]))
            numstr += s[i++];
        int n = numstr.empty() ? 1 : std::stoi(numstr);
        if (n < 1) throw std::runtime_error("Zero count for element "+elem);
        result.push_back({elem, n});
    }
    if (result.empty()) throw std::runtime_error("No elements found in '"+s+"'");
    int total = 0; for (auto& [e,n] : result) total += n;
    if (total < 2) throw std::runtime_error("Need at least 2 atoms total");
    if (total > 1000) throw std::runtime_error("Too many atoms (max 1000): "+std::to_string(total));
    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// Molecule template database
// ─────────────────────────────────────────────────────────────────────────────
inline MolTemplate get_molecule_template(const std::string& name) {
    // TIP4P-like geometry for water, standard geometries for others
    static std::map<std::string,MolTemplate> db;
    if (db.empty()) {
        // Water H2O — TIP4P geometry
        {
            MolTemplate m; m.name = "H2O";
            // O at origin, H atoms at ±52.26° from z axis, bond length 0.9572 Å
            double rOH = 0.9572, angle = 104.52*M_PI/180.0;
            double hx = rOH*std::sin(angle*0.5);
            double hz = rOH*std::cos(angle*0.5);
            m.sites = {
                {"O",  0.0,  0.0,  0.0,   -0.8476},
                {"H",  hx,   0.0,  hz,     0.4238},
                {"H", -hx,   0.0,  hz,     0.4238},
            };
            db["H2O"] = m; db["water"] = m;
        }
        // Ammonia NH3
        {
            MolTemplate m; m.name = "NH3";
            double rNH=1.012, angle=106.67*M_PI/180.0;
            double hx=rNH*std::sin(angle*0.5), hz=-rNH*std::cos(angle*0.5);
            m.sites = {
                {"N",  0.0,    0.0,    0.0,   -1.02},
                {"H",  hx,     0.0,    hz,     0.34},
                {"H", -hx*0.5, hx*0.866, hz,   0.34},
                {"H", -hx*0.5,-hx*0.866, hz,   0.34},
            };
            db["NH3"] = m; db["ammonia"] = m;
        }
        // CO2
        {
            MolTemplate m; m.name = "CO2";
            m.sites = {
                {"C",  0.0,  0.0,  0.0,    0.6512},
                {"O",  1.163,0.0,  0.0,   -0.3256},
                {"O", -1.163,0.0,  0.0,   -0.3256},
            };
            db["CO2"] = m;
        }
        // N2
        {
            MolTemplate m; m.name = "N2";
            m.sites = {{"N", 0.55, 0.0, 0.0, 0.0}, {"N", -0.55, 0.0, 0.0, 0.0}};
            db["N2"] = m;
        }
        // HF
        {
            MolTemplate m; m.name = "HF";
            m.sites = {{"F", 0.0, 0.0, 0.0, -0.36}, {"H", 0.917, 0.0, 0.0, 0.36}};
            db["HF"] = m;
        }
        // HCl
        {
            MolTemplate m; m.name = "HCl";
            m.sites = {{"Cl",0.0,0.0,0.0,-0.2}, {"H",1.275,0.0,0.0,0.2}};
            db["HCl"] = m;
        }
        // Methane CH4
        {
            MolTemplate m; m.name = "CH4";
            double r=1.089, s=r/std::sqrt(3.0);
            m.sites = {
                {"C", 0.0, 0.0, 0.0, -0.24},
                {"H",  s,   s,   s,   0.06},
                {"H", -s,  -s,   s,   0.06},
                {"H", -s,   s,  -s,   0.06},
                {"H",  s,  -s,  -s,   0.06},
            };
            db["CH4"] = m; db["methane"] = m;
        }
    }
    auto it = db.find(name);
    if (it == db.end())
        throw std::runtime_error("Molecule template '"+name+"' not found. Available: H2O, NH3, CO2, N2, HF, HCl, CH4");
    return it->second;
}
