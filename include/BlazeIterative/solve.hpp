// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_SOLVE_HPP
#define BLAZE_ITERATIVE_SOLVE_HPP

#include "IterativeCommon.hpp"
#include "IterativeTag.hpp"
#include "solvers/solvers.hpp"
#include <type_traits>
#include <cstring>
#include <utility>

#include <iostream>

BLAZE_NAMESPACE_OPEN
ITERATIVE_NAMESPACE_OPEN

/**
 * Solve a linear system using a preallocated buffer "x".
 * The values in "x" are used as the initial guess for
 * the iterative algorithm.
 */
template<typename MatrixType, typename T, typename TagType>
void solve_inplace(DynamicVector<T> &x,
                   const MatrixType &A,
                   const DynamicVector<T> &b,
                   TagType &tag)
{
    //Compile-time assertions
    BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE(MatrixType);
    static_assert(std::is_same<T, typename MatrixType::ElementType>::value,
                  "Matrix and vector data types must be the same");

    //Run-time assertions checking conditions that would be problems later anyway
    assert(A.columns() == b.size() && "A and b must have consistent dimensions");
    assert(x.size() == b.size() && "x and b must be the same length");
    assert(A.rows() == A.columns() && "A must be a square matrix");

    // Call specific solver
    detail::solve_impl(x, A, b, tag);
};


template<typename MatrixType, typename T, typename TagType>
void solve_inplace(DynamicVector<T> &x,
                   const MatrixType &A,
                   const DynamicVector<T> &b,
                   TagType &tag,
                   std::string Preconditioner)
{
    //Compile-time assertions
    BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE(MatrixType);
    static_assert(std::is_same<T, typename MatrixType::ElementType>::value,
                  "Matrix and vector data types must be the same");

    //Run-time assertions checking conditions that would be problems later anyway
    assert(A.columns() == b.size() && "A and b must have consistent dimensions");
    assert(x.size() == b.size() && "x and b must be the same length");
    assert(A.rows() == A.columns() && "A must be a square matrix");

    // Call specific solver
    detail::solve_impl(x, A, b, tag,Preconditioner);
};

    // For Arnoldi
    // Lanczos
    // GMRES
    template<typename MatrixType, typename T, typename TagType>
    void solve_inplace(DynamicVector<T> &x,
                       const MatrixType &A,
                       const DynamicVector<T> &b,
                       TagType &tag,
                       const std::size_t &n)
    {
        //Compile-time assertions
        BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE(MatrixType);
        static_assert(std::is_same<T, typename MatrixType::ElementType>::value,
                      "Matrix and vector data types must be the same");

        //Run-time assertions checking conditions that would be problems later anyway
        assert(A.columns() == b.size() && "A and b must have consistent dimensions");
        assert(n >= 1 && "n must be larger than or equal to 1");

        // Call specific solver
        detail::solve_impl(x, A, b, tag, n);
    };


/**
 * \brief Solver the linear system \f$ Ax = b \f$ using an iterative solver.
 *
 * Top-level entry point to the iterative solver collection.
 * This method returns a newly-allocated blaze::DynamicVector<T>
 * by allocating the return vector and calling the solve_inplace
 * method. All data integrity checking (that can be done independent
 * of the specific iterative method) is done in the call to solve_inplace.
 * The initial guess is zeros(length(b),1) (in matlab parlance).
 * To provide a different initial guess, use the solve_inplace method.
 *
 * The specific iterative solver is called via tag dispatch. The tag
 * holds methods to control logging and convergence criteria as well.
 * See blaze::iterative::IterativeTag for data members available to
 * all iterative methods. See the method-specific tag types for
 * members/methods available only to a particular iterative method.
 */
template<typename MatrixType, typename T, typename TagType>
DynamicVector<T> solve(const MatrixType &A, const DynamicVector<T> &b, TagType &tag)
{
    DynamicVector<T> x(b.size(), 0.0);
    solve_inplace(x, A, b, tag);

    return x;
};

/**
 * \brief Solver the linear system \f$ Ax = b \f$ using an iterative solver using a preconditioner.
 *
 * Top-level entry point to the iterative solver collection.
 * This method returns a newly-allocated blaze::DynamicVector<T>
 * by allocating the return vector and calling the solve_inplace
 * method. All data integrity checking (that can be done independent
 * of the specific iterative method) is done in the call to solve_inplace.
 * The initial guess is zeros(length(b),1) (in matlab parlance).
 * To provide a different initial guess, use the solve_inplace method.
 *
 * The specific iterative solver is called via tag dispatch. The tag
 * holds methods to control logging and convergence criteria as well.
 * See blaze::iterative::IterativeTag for data members available to
 * all iterative methods. See the method-specific tag types for
 * members/methods available only to a particular iterative method.
 */
template<typename MatrixType, typename T, typename TagType>
DynamicVector<T> solve(const MatrixType &A, 
                       const DynamicVector<T> &b, 
                       TagType &tag, 
                       std::string Preconditioner)
{
    DynamicVector<T> x(b.size(), 0.0);
    solve_inplace(x, A, b, tag, Preconditioner);

    return x;
};


// For Arnoldi
// For Lanczos
// For GMRES
    template<typename MatrixType, typename T, typename TagType>
    DynamicVector<T> solve(const MatrixType &A, const DynamicVector<T> &b, TagType &tag, const std::size_t &n)
    {
        if(typeid(tag).name() == typeid(ArnoldiTag).name()) {
            DynamicVector<T> x(n, 0.0);
            solve_inplace(x, A, b, tag, n);
            return x;

        } else if(typeid(tag).name() == typeid(LanczosTag).name()) {
            DynamicVector<T> x(n, 0.0);
            solve_inplace(x, A, b, tag, n);
            return x;
        } else {
            DynamicVector<T> x(n, 0.0);
            solve_inplace(x, A, b, tag, n);
            return x;
        }

    };


ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_SOLVE_HPP
