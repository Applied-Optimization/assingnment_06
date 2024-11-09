#pragma once

#include "SpringElement2DWithLength.hh"

#include <Eigen/Eigenvalues>
//== NAMESPACES ===============================================================

namespace AOPT {

//== CLASS DEFINITION =========================================================

/**
*   Class that overrides the hessian of the non-convex energy of the spring element
 * by fixing the negative eigen values of the hessian matrix
*/

class SpringElement2DWithLengthPSDHess : public SpringElement2DWithLength {
public:

    SpringElement2DWithLengthPSDHess(): SpringElement2DWithLength() {}

    inline virtual void eval_hessian(const Vec &_x, const Vec &_coeffs, Mat &_H) override {
        //------------------------------------------------------//
        //TODO: compute the hessian matrix and project it to a positve definite matrix
        //Hint: 1. to compute the eigen decomposition, use
        //          Eigen::SelfAdjointEigenSolver<Mat> solver(A);
        //          Mat evecs = solver.eigenvectors();  //this matrix contains the eigenvectors in its columns
        //          Vec evals = solver.eigenvalues();
        //      2. to convert a vector d to a (dense) diagonal matrix D, use
        //          D = d.asDiagonal()

        // 1.Compute the Hessian matrix
        SpringElement2DWithLength::eval_hessian(_x, _coeffs, _H);
        // 2.Compute the eigen decomposition of the Hessian matrix
        Eigen::SelfAdjointEigenSolver<Mat> solver(_H);
        // 3.Get the eigen values
        Vec evals = solver.eigenvalues();
        // 4.Get the eigen vectors
        Mat evecs = solver.eigenvectors();
        // 5. Project the negative eigen values to a small positive value
        for (int i = 0; i < evals.size(); ++i) {
            if (evals[i] < 0) {
                evals[i] = m_eps;
            }
        }
        // 6. Reconstruct the Hessian matrix
        _H = evecs * evals.asDiagonal() * evecs.transpose();


        //------------------------------------------------------//
    }

    static constexpr double m_eps = 1e-7;
};

//=============================================================================

}
