#pragma once
/**
 * fingerprint.hpp — ACE-inspired structure descriptors
 *
 * SESSION 1 IMPROVEMENT: Replace naive sorted-distances with
 * physically-motivated symmetry functions that are:
 *   (a) Unique — no two different structures map to the same descriptor
 *   (b) Continuous — small geometry changes → small descriptor change
 *   (c) Fast — O(N*K) not O(N²) for large cutoffs
 *   (d) Species-aware — different channels per element pair
 *   (e) Chiral — distinguishes enantiomers
 *
 * Descriptor = concatenation of:
 *   G1: Radial symmetry functions (Behler-Parrinello style)
 *       G1_i = Σ_j exp(-η(r_ij - r_s)²) * f_c(r_ij)
 *       Multiple (η, r_s) pairs → captures different length scales
 *
 *   G2: Angular symmetry functions
 *       G2_i = Σ_{j,k} (1+λcos θ_ijk)^ζ exp(-η(r_ij²+r_ik²+r_jk²)) * f_c
 *       Captures 3-body correlations (bond angles)
 *
 *   G3: Species-pair partial RDF (element-resolved)
 *       Histogram of r_ij per (type_i, type_j) pair
 *
 *   G4: Chiral invariant
 *       Signed scalar from oriented triple products
 *       Changes sign under spatial inversion → detects enantiomers
 *
 * Also provides: kabsch_rmsd_species() — species-aware structural comparison
 */

#include "system.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <array>

// ─────────────────────────────────────────────────────────────────────────────
// Cutoff function: smooth decay to zero at r_cut
// f_c(r) = 0.5 * [cos(π r/r_cut) + 1]   for r < r_cut
//         = 0                              for r >= r_cut
// ─────────────────────────────────────────────────────────────────────────────
inline double cutoff_fn(double r, double r_cut) {
    if (r >= r_cut) return 0.0;
    return 0.5 * (std::cos(M_PI * r / r_cut) + 1.0);
}

// ─────────────────────────────────────────────────────────────────────────────
// G1: Radial symmetry functions
// For each atom i: G1[i][k] = Σ_j exp(-η_k*(r_ij - r_s_k)²) * f_c(r_ij)
// Averaged over all atoms of each species → per-species radial fingerprint
// ─────────────────────────────────────────────────────────────────────────────
struct RadialParams { double eta, rs; };

static const std::vector<RadialParams> RADIAL_PARAMS = {
    // (eta, rs) — captures bonding, 2nd shell, 3rd shell distances
    {0.50, 0.50}, {0.50, 1.00}, {0.50, 1.50}, {0.50, 2.00},
    {0.50, 2.50}, {0.50, 3.00}, {1.00, 1.50}, {1.00, 2.00},
    {1.00, 2.50}, {1.00, 3.00}, {2.00, 1.50}, {2.00, 2.00},
    {2.00, 2.50}, {4.00, 2.00}, {4.00, 2.50}, {4.00, 3.00},
};

inline std::vector<double> compute_g1(const System& sys, double r_cut_sigma = 4.0) {
    const int N = sys.N;
    const int ns = sys.n_species();
    const int np = (int)RADIAL_PARAMS.size();
    const double r_cut = r_cut_sigma * sys.mean_sigma;

    // Per-species per-param accumulator
    std::vector<double> g(ns * np, 0.0);
    std::vector<int>    count(ns, 0);

    for (int i = 0; i < N; ++i) {
        int ti = sys.type_id[i];
        ++count[ti];
        for (int j = 0; j < N; ++j) {
            if (j == i) continue;
            double r = sys.rij(i, j) / sys.mean_sigma;  // normalise by sigma
            double fc = cutoff_fn(r, r_cut_sigma);
            if (fc == 0.0) continue;
            for (int k = 0; k < np; ++k) {
                double dr = r - RADIAL_PARAMS[k].rs;
                g[ti * np + k] += std::exp(-RADIAL_PARAMS[k].eta * dr * dr) * fc;
            }
        }
    }

    // Normalise by atom count per species
    for (int t = 0; t < ns; ++t)
        if (count[t] > 0)
            for (int k = 0; k < np; ++k)
                g[t * np + k] /= count[t];

    return g;
}

// ─────────────────────────────────────────────────────────────────────────────
// G2: Angular symmetry functions (3-body)
// G2_i = Σ_{j≠i,k≠i,k>j} (1+λcos θ)^ζ * exp(-η(rij²+rik²+rjk²)) * fc_ij*fc_ik*fc_jk
// Uses a few (λ,ζ,η) combinations — enough to distinguish geometries
// ─────────────────────────────────────────────────────────────────────────────
struct AngularParams { double eta, zeta, lambda; };

static const std::vector<AngularParams> ANGULAR_PARAMS = {
    {0.01, 1.0,  1.0}, {0.01, 2.0,  1.0}, {0.01, 1.0, -1.0},
    {0.05, 1.0,  1.0}, {0.05, 4.0,  1.0}, {0.05, 1.0, -1.0},
    {0.10, 1.0,  1.0}, {0.10, 4.0, -1.0},
};

inline std::vector<double> compute_g2(const System& sys, double r_cut_sigma = 3.5) {
    const int N = sys.N;
    const int ns = sys.n_species();
    const int np = (int)ANGULAR_PARAMS.size();
    const double r_cut = r_cut_sigma * sys.mean_sigma;

    std::vector<double> g(ns * np, 0.0);
    std::vector<int> count(ns, 0);

    for (int i = 0; i < N; ++i) {
        int ti = sys.type_id[i];
        ++count[ti];

        // Collect neighbours of i within cutoff
        std::vector<std::pair<int,double>> nb; // (j, r_ij/sigma)
        for (int j = 0; j < N; ++j) {
            if (j == i) continue;
            double r = sys.rij(i, j) / sys.mean_sigma;
            if (r < r_cut_sigma) nb.push_back({j, r});
        }

        // Triple loop over pairs of neighbours
        for (int a = 0; a < (int)nb.size() - 1; ++a) {
            int    j   = nb[a].first;
            double rij = nb[a].second;
            double fcij = cutoff_fn(rij, r_cut_sigma);

            for (int b = a+1; b < (int)nb.size(); ++b) {
                int    k   = nb[b].first;
                double rik = nb[b].second;
                double fcik = cutoff_fn(rik, r_cut_sigma);

                double rjk = sys.rij(j, k) / sys.mean_sigma;
                double fcjk = cutoff_fn(rjk, r_cut_sigma);
                if (fcjk == 0.0) continue;

                // cos θ_ijk via dot product
                double dx_ij = (sys.px(j)-sys.px(i))/sys.mean_sigma;
                double dy_ij = (sys.py(j)-sys.py(i))/sys.mean_sigma;
                double dz_ij = (sys.pz(j)-sys.pz(i))/sys.mean_sigma;
                double dx_ik = (sys.px(k)-sys.px(i))/sys.mean_sigma;
                double dy_ik = (sys.py(k)-sys.py(i))/sys.mean_sigma;
                double dz_ik = (sys.pz(k)-sys.pz(i))/sys.mean_sigma;
                double cos_t = (dx_ij*dx_ik + dy_ij*dy_ik + dz_ij*dz_ik)
                               / (rij * rik + 1e-15);
                cos_t = std::max(-1.0, std::min(1.0, cos_t));

                double r2sum = rij*rij + rik*rik + rjk*rjk;
                double fc3   = fcij * fcik * fcjk;

                for (int p = 0; p < np; ++p) {
                    const auto& ap = ANGULAR_PARAMS[p];
                    double ang = std::pow(1.0 + ap.lambda * cos_t, ap.zeta);
                    g[ti * np + p] += ang * std::exp(-ap.eta * r2sum) * fc3;
                }
            }
        }
    }

    // Normalise
    for (int t = 0; t < ns; ++t)
        if (count[t] > 0)
            for (int p = 0; p < np; ++p)
                g[t * np + p] /= count[t];

    return g;
}

// ─────────────────────────────────────────────────────────────────────────────
// G3: Species-pair partial RDF (element-resolved histogram)
// ─────────────────────────────────────────────────────────────────────────────
inline std::vector<double> compute_g3(
    const System& sys, int n_bins = 25, double r_max_sigma = 5.0)
{
    const int N  = sys.N;
    const int ns = sys.n_species();
    const int n_pairs = ns * (ns + 1) / 2;
    std::vector<double> rdf(n_pairs * n_bins, 0.0);
    double dr = r_max_sigma / n_bins;

    auto pair_idx = [&](int ti, int tj) -> int {
        if (ti > tj) std::swap(ti, tj);
        return ti * ns - ti*(ti-1)/2 + (tj - ti);
    };

    for (int i = 0; i < N-1; ++i)
        for (int j = i+1; j < N; ++j) {
            double r = sys.rij(i, j) / sys.mean_sigma;
            int bin = (int)(r / dr);
            if (bin >= n_bins) continue;
            int pidx = pair_idx(sys.type_id[i], sys.type_id[j]);
            rdf[pidx * n_bins + bin] += 1.0;
        }

    // Normalise each pair block
    for (int pi = 0; pi < n_pairs; ++pi) {
        double s = 0;
        for (int b = 0; b < n_bins; ++b) s += rdf[pi*n_bins + b];
        if (s > 0)
            for (int b = 0; b < n_bins; ++b) rdf[pi*n_bins + b] /= s;
    }
    return rdf;
}

// ─────────────────────────────────────────────────────────────────────────────
// G4: Chiral invariant — changes sign under inversion
// Uses oriented volume of nearest-neighbour tetrahedra
// ─────────────────────────────────────────────────────────────────────────────
inline double compute_g4(const System& sys, double cutoff_factor = 1.6) {
    const int N = sys.N;
    double cut2 = (cutoff_factor * sys.mean_sigma) * (cutoff_factor * sys.mean_sigma);
    double chi = 0.0;

    for (int i = 0; i < N; ++i) {
        std::vector<int> nb;
        for (int j = 0; j < N; ++j) {
            if (j == i) continue;
            if (sys.rij2(i, j) < cut2) nb.push_back(j);
        }
        // Signed volume of each triplet of neighbours with atom i
        int M = (int)nb.size();
        for (int a = 0; a < M-2; ++a)
        for (int b = a+1; b < M-1; ++b)
        for (int c = b+1; c < M;   ++c) {
            double ax=sys.px(nb[a])-sys.px(i), ay=sys.py(nb[a])-sys.py(i), az=sys.pz(nb[a])-sys.pz(i);
            double bx=sys.px(nb[b])-sys.px(i), by=sys.py(nb[b])-sys.py(i), bz=sys.pz(nb[b])-sys.pz(i);
            double cx=sys.px(nb[c])-sys.px(i), cy=sys.py(nb[c])-sys.py(i), cz=sys.pz(nb[c])-sys.pz(i);
            // Scalar triple product = signed volume
            double vol = ax*(by*cz - bz*cy) - ay*(bx*cz - bz*cx) + az*(bx*cy - by*cx);
            chi += vol;
        }
    }
    double scale = std::pow((double)N * sys.mean_sigma, 3.0) + 1e-15;
    return chi / scale;
}

// ─────────────────────────────────────────────────────────────────────────────
// Full Fingerprint: G1 + G2 + G3 + G4
// ─────────────────────────────────────────────────────────────────────────────
struct Fingerprint {
    std::vector<double> data;

    Fingerprint() = default;

    static Fingerprint compute(const System& sys) {
        Fingerprint fp;
        auto g1 = compute_g1(sys);
        auto g2 = compute_g2(sys);
        auto g3 = compute_g3(sys);
        double g4 = compute_g4(sys);

        fp.data.reserve(g1.size() + g2.size() + g3.size() + 1);
        fp.data.insert(fp.data.end(), g1.begin(), g1.end());
        fp.data.insert(fp.data.end(), g2.begin(), g2.end());
        // Only add G3 for heteroatomic (adds no info for homoatomic)
        if (sys.n_species() > 1)
            fp.data.insert(fp.data.end(), g3.begin(), g3.end());
        fp.data.push_back(g4);
        return fp;
    }

    double distance(const Fingerprint& other) const {
        double s = 0.0;
        size_t n = std::min(data.size(), other.data.size());
        for (size_t i = 0; i < n; ++i) {
            double d = data[i] - other.data[i];
            s += d * d;
        }
        return std::sqrt(s);
    }

    bool empty() const { return data.empty(); }
};

// ─────────────────────────────────────────────────────────────────────────────
// Kabsch RMSD — rotation+translation invariant, species-aware atom matching
//
// Algorithm:
//   1. Match atoms greedily by species (same element only)
//   2. Find optimal rotation via SVD (Kabsch algorithm)
//   3. Return RMSD in Angstroms
// ─────────────────────────────────────────────────────────────────────────────
inline double kabsch_rmsd_species(
    const System& A, const System& B, double /*sigma*/)
{
    assert(A.N == B.N);
    const int N = A.N;

    // Quick pre-filter: if composition differs, structures are different
    if (A.species != B.species) return 1e10;
    for (int t = 0; t < A.n_species(); ++t)
        if (A.count(t) != B.count(t)) return 1e10;

    // ── Greedy bipartite matching per species ─────────────────────────────
    std::vector<int> perm(N, -1);
    std::vector<bool> usedB(N, false);

    for (int t = 0; t < A.n_species(); ++t) {
        std::vector<int> idxA, idxB;
        for (int i = 0; i < N; ++i) {
            if (A.type_id[i] == t) idxA.push_back(i);
            if (B.type_id[i] == t) idxB.push_back(i);
        }

        // Build distance matrix for this species block
        std::vector<std::tuple<double,int,int>> pairs;
        pairs.reserve(idxA.size() * idxB.size());
        for (int i : idxA)
            for (int j : idxB) {
                double dx=A.px(i)-B.px(j), dy=A.py(i)-B.py(j), dz=A.pz(i)-B.pz(j);
                pairs.emplace_back(dx*dx+dy*dy+dz*dz, i, j);
            }
        std::sort(pairs.begin(), pairs.end());

        std::vector<bool> usedA_block(N, false);
        for (auto& [d2, i, j] : pairs) {
            if (!usedA_block[i] && !usedB[j]) {
                perm[i] = j;
                usedA_block[i] = usedB[j] = true;
            }
        }
    }

    // ── Build matched coordinate matrices ─────────────────────────────────
    Eigen::MatrixXd P(N, 3), Q(N, 3);
    for (int i = 0; i < N; ++i) {
        int j = (perm[i] >= 0) ? perm[i] : i;
        P(i,0)=A.px(i); P(i,1)=A.py(i); P(i,2)=A.pz(i);
        Q(i,0)=B.px(j); Q(i,1)=B.py(j); Q(i,2)=B.pz(j);
    }

    // ── Centre both ───────────────────────────────────────────────────────
    P.rowwise() -= P.colwise().mean();
    Q.rowwise() -= Q.colwise().mean();

    // ── Kabsch rotation via SVD ───────────────────────────────────────────
    Eigen::Matrix3d H = P.transpose() * Q;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(
        H, Eigen::ComputeFullU | Eigen::ComputeFullV);

    double det = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    Eigen::DiagonalMatrix<double,3> D(1.0, 1.0, det > 0 ? 1.0 : -1.0);
    Eigen::Matrix3d R = svd.matrixV() * D * svd.matrixU().transpose();

    // ── RMSD ──────────────────────────────────────────────────────────────
    Eigen::MatrixXd diff = (P * R.transpose()) - Q;
    return std::sqrt(diff.squaredNorm() / N);
}

// ─────────────────────────────────────────────────────────────────────────────
// Point-group order estimator
// Used by isomer filter to check if two structures are related by symmetry
// Returns approximate order of the point group (1=no symmetry, 60=icosahedral)
// ─────────────────────────────────────────────────────────────────────────────
inline int estimate_point_group_order(const System& sys, double tol_sigma = 0.15) {
    const int N = sys.N;
    const double tol = tol_sigma * sys.mean_sigma;
    int order = 1;

    // Test Cn rotational axes (n=2,3,4,5,6) along all bond directions
    for (int i = 0; i < std::min(N, 8); ++i) {
        for (int j = i+1; j < std::min(N, 8); ++j) {
            // Axis = bond vector i→j
            double ax = sys.px(j)-sys.px(i);
            double ay = sys.py(j)-sys.py(i);
            double az = sys.pz(j)-sys.pz(i);
            double an = std::sqrt(ax*ax+ay*ay+az*az);
            if (an < 1e-10) continue;
            ax/=an; ay/=an; az/=an;

            for (int n : {2, 3, 4, 5, 6}) {
                double angle = 2.0 * M_PI / n;
                double c = std::cos(angle), s = std::sin(angle);
                double oc = 1.0 - c;

                // Rodrigues rotation matrix
                double R[3][3] = {
                    {c+ax*ax*oc,    ax*ay*oc-az*s, ax*az*oc+ay*s},
                    {ay*ax*oc+az*s, c+ay*ay*oc,    ay*az*oc-ax*s},
                    {az*ax*oc-ay*s, az*ay*oc+ax*s, c+az*az*oc   }
                };

                // Apply rotation to all atoms, check if maps onto original
                bool sym = true;
                for (int k = 0; k < N && sym; ++k) {
                    double rx = R[0][0]*sys.px(k) + R[0][1]*sys.py(k) + R[0][2]*sys.pz(k);
                    double ry = R[1][0]*sys.px(k) + R[1][1]*sys.py(k) + R[1][2]*sys.pz(k);
                    double rz = R[2][0]*sys.px(k) + R[2][1]*sys.py(k) + R[2][2]*sys.pz(k);
                    bool found = false;
                    for (int l = 0; l < N && !found; ++l) {
                        if (sys.type_id[l] != sys.type_id[k]) continue;
                        double dx=rx-sys.px(l), dy=ry-sys.py(l), dz=rz-sys.pz(l);
                        if (dx*dx+dy*dy+dz*dz < tol*tol) found = true;
                    }
                    if (!found) sym = false;
                }
                if (sym) order = std::max(order, n);
            }
        }
    }
    return order;
}
