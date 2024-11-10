#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include "LineSearch.hh"

//== NAMESPACES ===============================================================

namespace AOPT
{
    /**
     * @brief NewtonMethods is a collection of functions implementing several variations of
     * Newton's method.
     */
    class NewtonMethods
    {
    public:
        typedef FunctionBaseSparse::Vec Vec;   // dense vector arbitrary size
        typedef FunctionBaseSparse::Mat Mat;   // dense matrix arbitrary size
        typedef FunctionBaseSparse::T T;       // Triplets
        typedef FunctionBaseSparse::SMat SMat; // sparse matrix arbitrary size

        /**
         * @brief solve
         * \param _problem pointer to any function/problem inheriting from FunctionBaseSparse
         *        on which the basic Newton Method will be applied
         * \param _initial_x starting point of the method
         * \param _eps epsilon under which the method stops
         * \param _max_iters maximum iteration of the method
         */
        static Vec solve(FunctionBaseSparse *_problem, const Vec &_initial_x, const double _eps = 1e-4, const int _max_iters = 1000000)
        {
            std::cout << "******** Newton Method ********" << std::endl;

            // Squared epsilon for stopping criterion
            double e2 = 2 * _eps * _eps;

            int n = _problem->n_unknowns();

            // Initial point
            Vec x = _initial_x;

            // Allocate gradient and Hessian storage
            Vec g(n);
            SMat H(n, n);

            // Allocate storage for search direction
            Vec delta_x(n);
            int iter = 0;

            // Initialize the sparse Cholesky solver
            Eigen::SimplicialLLT<SMat> solver;

            //------------------------------------------------------//
            // TODO: Implement Newton method

            while (iter < _max_iters)
            {
                // Evaluate gradient and Hessian at current point x
                _problem->eval_gradient(x, g);
                _problem->eval_hessian(x, H);

                // Check for convergence: if gradient norm squared is less than epsilon
                if (g.squaredNorm() < e2)
                {
                    std::cout << "Convergence achieved after " << iter << " iterations." << std::endl;
                    break;
                }

                // Solve H * delta_x = -g for delta_x
                solver.compute(H);
                if (solver.info() != Eigen::Success)
                {
                    std::cerr << "Failed to decompose Hessian matrix; it may not be positive definite." << std::endl;
                    break;
                }
                delta_x = solver.solve(-g);

                // Perform backtracking line search to find a suitable step size
                double t = LineSearch::backtracking_line_search(_problem, x, g, delta_x, 1.0);

                // Update x
                x += t * delta_x;

                ++iter;
            }

            if (iter == _max_iters)
            {
                std::cout << "Reached maximum iterations without convergence." << std::endl;
            }

            //------------------------------------------------------//

            return x;
        }

        /**
         * @brief solve with the Projected Hessian method
         * \param _problem pointer to any function/problem inheriting from FunctionBaseSparse.
         *        This problem MUST provide a working eval_hessian() function for this method to work.
         *
         * \param _initial_x starting point of the method
         * \param _gamma scaling factor for the increase of delta
         * \param _eps epsilon under which the method stops
         * \param _max_iters maximum iteration of the method
         */
        static Vec solve_with_projected_hessian(FunctionBaseSparse *_problem, const Vec &_initial_x, const double _gamma = 10.0,
                                                const double _eps = 1e-4, const int _max_iters = 1000000)
        {
            bool converged = false;
            return solve_with_projected_hessian(_problem, converged, _initial_x, _gamma, _eps, _max_iters);
        }

        static Vec solve_with_projected_hessian(FunctionBaseSparse *_problem, bool &_converged, const Vec &_initial_x, const double _gamma = 10.0,
                                                const double _eps = 1e-4, const int _max_iters = 1000000)
        {
            std::cout << "******** Newton Method with Projected Hessian ********" << std::endl;

            // Squared epsilon for stopping criterion
            double e2 = 2 * _eps * _eps;

            int n = _problem->n_unknowns();

            // Initial point
            Vec x = _initial_x;

            // Allocate gradient and Hessian storage
            Vec g(n);
            SMat H(n, n);

            // Allocate search direction vector storage
            Vec delta_x(n);
            int iter = 0;

            // Identity matrix and scaling factor for adding positive values to the diagonal
            SMat I(n, n);
            I.setIdentity();

            _converged = false;

            // Initialize sparse Cholesky solver
            Eigen::SimplicialLLT<SMat> solver;

            //------------------------------------------------------//

            // Compute initial gradient and Hessian
            _problem->eval_gradient(x, g);
            _problem->eval_hessian(x, H);

            //------------------------------------------------------//
            // TODO: Implement Newton with Projected Hessian method
            // Hint: if the factorization fails, add delta * I to the Hessian and retry.
            //       Increase delta by multiplying by _gamma (recommended gamma > 1).
            //
            // Set initial delta based on trace(H) as:
            // delta = 1e-3 * |trace(H)| / n
            // Initial delta based on the trace of H
            double trace_H = H.diagonal().sum();         // Sum of the diagonal of H (approx. trace)
            double delta = 1e-3 * std::abs(trace_H) / n; // Initial value for delta
            std::cout << "Initial delta calculated: " << delta << std::endl;

            while (iter < _max_iters)
            {
                // Attempt Cholesky factorization of H + delta * I
                SMat H_mod = H + delta * I;
                solver.compute(H_mod);

                // Check if factorization succeeded
                if (solver.info() == Eigen::Success)
                {
                    // Compute search direction: delta_x = -H^-1 * g
                    delta_x = solver.solve(-g);

                    // Update the current solution
                    x += delta_x;

                    // Check for convergence
                    if (g.norm() < _eps)
                    {
                        _converged = true;
                        std::cout << "Convergence reached in " << iter << " iterations." << std::endl;
                        break;
                    }

                    // Recompute gradient and Hessian for next iteration
                    _problem->eval_gradient(x, g);
                    _problem->eval_hessian(x, H);

                    // Force delta to 1e-1 after initialization
                    delta = 1e-1;
                }
                else
                {
                    // Factorization failed, increase delta by gamma
                    delta *= _gamma;
                    std::cout << "Factorization failed, increasing delta to " << delta << std::endl;
                }

                // Increment iteration count
                iter++;

                if (iter >= _max_iters)
                {
                    std::cout << "Max iterations reached without convergence." << std::endl;
                }
            }

            //------------------------------------------------------//

            return x;
        }
        // end of solve_with_projected_hessian
    };
} // namespace AOPT