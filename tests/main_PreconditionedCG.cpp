// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "BlazeIterative.hpp"
#include <iostream>
#include <cstdlib>

using namespace blaze;
using namespace blaze::iterative;

int main() {

    DynamicMatrix<double,false> A{{2,-1,0},{-1,2,-1},{0,-1,1}};
    DynamicVector<double> b{0, 0, 1};
    DynamicVector<double> x1{1, 2, 3};


    PreconditionCGTag tag;
    tag.do_log() = true;
    auto x2 = solve(A,b,tag, "Jacobi");

    auto error = norm(x1 - x2);

    if (error < EPSILON){
        std::cout << " Pass test of Preconditioned CG" << std::endl;
        return EXIT_SUCCESS;
    } else{
        std::cout << "Fail test of Preconditioned CG" << std::endl;
        return EXIT_FAILURE;
    }

}

