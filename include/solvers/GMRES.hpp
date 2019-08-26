// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_GMRES_HPP
#define BLAZE_ITERATIVE_GMRES_HPP

#include <iostream>
#include "IterativeCommon.hpp"
#include "GMRESTag.hpp"
#include <utility>
#include "Arnoldi.hpp"

BLAZE_NAMESPACE_OPEN
    ITERATIVE_NAMESPACE_OPEN

        namespace detail {

            template< typename MatrixType, typename T, typename VectorType>
            std::pair< VectorType,  VectorType> arnoldi( const MatrixType &A, MatrixType &Q, int &k)
            {
                std::size_t m = A.columns();
                DynamicVector<T> q(m);
                DynamicVector<T> h(k+2);

                q = A * column(Q,k);
                for(int i = 0; i <= k; ++i){
                    h[i] = ctrans(q) * column(Q,i);
                    q -= h[i] * column(Q,i);
                }
                h[k+1] = norm(q);
                q = q / h[k+1];
                return std::make_pair(h,q);
            }

            template<typename T>
            std::pair< T, T> givens_rotation( T &v1, T &v2)
            {
                T cs_k, sn_k;
                if (v1 == 0){
                    cs_k = 0;
                    sn_k = 0;

                } else {
                    auto t = sqrt(v1 * v1 + v2 * v2);
                    cs_k = abs(v1) / t;
                    sn_k = cs_k * v2 / v1;
                }

                return std::make_pair(cs_k,sn_k);
            }

            template<typename T, typename VectorType>
            std::tuple< VectorType,  T, T> apply_givens_rotation( VectorType &h, VectorType &cs, VectorType &sn, int &k)
            {
                for(int i = 0; i < k; ++i){
                    auto temp = cs[i] * h[i] + sn[i] * h[i+1];
                    h[i+1] = -sn[i] * h[i] + cs[i] * h[i+1];
                    h[i] = temp;
                }

                T cs_k, sn_k;
                auto res = givens_rotation(h[k], h[k+1]);
                cs_k = res.first;
                sn_k = res.second;

                h[k] = cs_k * h[k] + sn_k * h[k+1];
                h[k+1] = 0.0;

                return std::make_tuple(h, cs_k, sn_k);
            }

            template<typename MatrixType, typename T>
            void  solve_impl(
             //       MatrixType &h,
              //      MatrixType &Q,
                    DynamicVector<T> &x,
                    const MatrixType &A,
                    const DynamicVector<T> &b,
                    const DynamicVector<T> &x0,
                    GMRESTag &tag,
                    const std::size_t &n
            )
            {

                BLAZE_INTERNAL_ASSERT(A.isSymmetric(), "A must be a symmetric matrix")
                
                BLAZE_INTERNAL_ASSERT(n >= 1, "n must larger than or equal to 1")

                // A: m * m matrix; n is max_iteration;


                std::size_t m = A.columns();
                DynamicVector<T> r(m);
                DynamicMatrix<T> H(n+1,n);
                DynamicMatrix<T> Q(m,n+1);
                DynamicVector<T> q(m);
                DynamicVector<T> sn(n,0);
                DynamicVector<T> cs(n,0);
                DynamicVector<T> e1(n+1,0);
                DynamicVector<T> beta(n+1);


                r = b - A * x0;
                auto err = norm(r) / norm(b);

                e1[0] = 1;
                column(Q,0) = r / norm(r);
                beta = norm(r) * e1;

                for(int k = 0; k < n; ++k){
                    auto res_1 = arnoldi(A, Q, k);
                    column(H,k) = res_1.first;
                    column(Q,k+1) = res_1.second;

                    auto res_2 = apply_givens_rotation(column(H,k), cs, sn, k);
                    column(H,k) = res_2.first;

                    beta[k+1] = -sn[k] * beta[k];
                    beta[k] = cs[k] * beta[k];
                    err = abs(beta[k+1]) / norm(b);

                    double eps = 1e-12;
                    if (err <= eps)
                        break;

                }
                auto H_sub = submatrix(H, 0, 0, n, n);
                auto beta_sub = subvector(beta, 0, n);
                auto y = inv(H_sub) *  beta_sub;
                auto Q_sub = submatrix(Q, 0, 0, n, n);
                x += Q_sub * y;



            }; // end solve_imple function

        } //end namespace detail

    ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_GMRES_HPP
