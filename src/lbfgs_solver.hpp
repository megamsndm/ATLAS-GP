#pragma once
/**
 * lbfgs_solver.hpp — Local minimizer using L-BFGS-B
 *
 * Potential-agnostic: works with ANY Potential implementation.
 * Two modes:
 *   local_minimize()  — standard convergence (search phase)
 *   tight_minimize()  — very tight convergence (final refinement)
 */
#include <LBFGSB.h>
#include "system.hpp"
#include "potential/potential_base.hpp"
#include <iomanip>

using namespace LBFGSpp;

struct PotentialFunctor {
    const Potential& pot;
    System& sys;
    PotentialFunctor(const Potential& p, System& s) : pot(p), sys(s) {}
    double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) const {
        for (int i = 0; i < 3*sys.N; ++i) sys.pos[i] = x(i);
        grad.setZero();
        return pot.energy_gradient(sys, grad.data());
    }
};

inline double local_minimize(
    System& sys,
    const Potential& pot,
    int    max_iter = 500,
    double g_tol    = 1e-8)
{
    sys.centre();
    sys.project_sphere();

    LBFGSBParam<double> param;
    param.max_iterations = max_iter;
    param.epsilon        = g_tol;
    param.epsilon_rel    = 1e-14;
    param.max_linesearch = 200;
    param.m              = 10;

    LBFGSBSolver<double> solver(param);
    PotentialFunctor f(pot, sys);

    Eigen::VectorXd lb = Eigen::VectorXd::Constant(3*sys.N, -1.5*sys.R);
    Eigen::VectorXd ub = Eigen::VectorXd::Constant(3*sys.N,  1.5*sys.R);
    Eigen::VectorXd x  = Eigen::Map<Eigen::VectorXd>(sys.pos.data(), 3*sys.N);

    double fx = pot.energy_only(sys);
    try { solver.minimize(f, x, fx, lb, ub); } catch (...) {}

    for (int i = 0; i < 3*sys.N; ++i) sys.pos[i] = x(i);
    sys.centre();

    // Reject unphysical overlap
    const double min_r2 = (0.6*sys.mean_sigma)*(0.6*sys.mean_sigma);
    for (int i = 0; i < sys.N-1; ++i)
        for (int j = i+1; j < sys.N; ++j)
            if (sys.rij2(i,j) < min_r2) return 1e30;

    return pot.energy_only(sys);
}

inline double tight_minimize(System& sys, const Potential& pot) {
    return local_minimize(sys, pot, 2000, 1e-12);
}
