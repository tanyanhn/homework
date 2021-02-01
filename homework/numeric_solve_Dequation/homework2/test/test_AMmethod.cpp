#include "catch.hpp"
#include "../src/ThreeBody.h"
#include <vector>
#include "equalrangematcher.h"
#include <fstream>
#include <memory>

using namespace std;
using namespace Catch;
using namespace solver_ty;
using namespace integrator_ty;

SCENARIO("fdf caculate", "[fdf]"){
    Real m = 1/81.45;
    fdf test(m);
    int size = 36;
    Real u[] = {1, 2, 3, 4, 5, 6};
    Real A[36],
            aim[] = {0, 0, 0, 1.000000000000000, 0, 0,
        0 ,0, 0, 0, 1.000000000000000, 0,
        0, 0, 0, 0, 0, 1.000000000000000,
        0.985054292842976, 0.008145739754421, 0.012218609631632, 0, 2.000000000000000, 0,
        0.008145739754421, 0.997264894715356, 0.024499101271573, -2.000000000000000, 0, 0,
        0.012218609631632, 0.024499101271573, 0.017680812441667, 0, 0, 0};
    test.df(u, A);
    CHECK_THAT(vector<Real>(A, A+size), EqualsRange(vector<Real>(aim, aim+size), 1e-5));
    m = 1/179.9;
    fdf testm(m);
    Real um[] = {10.1,-0.4, -8.3, 11.1, 9.99, 7.3};
    Real aimm[] = {0, 0, 0, 1.000000000000000, 0, 0,
        0, 0, 0, 0, 1.000000000000000, 0,
        0, 0, 0, 0, 0, 1.000000000000000,
        1.000352605058543, -0.000031671787999, -0.000657189600981, 0, 2.000000000000000, 0,
        -0.000031671787999, 0.999554240261516, 0.000026032240203, -2.000000000000000, 0, 0,
        -0.000657189600981, 0.000026032240203, 0.000093154679941, 0, 0, 0};
    test.df(um, A);
    CHECK_THAT(vector<Real>(A, A+size), EqualsRange(vector<Real>(aimm, aimm+size), 1e-5));
}

SCENARIO("ThreeBody AM method for initial 1","[AMmethod][initial1]"){
    int size = 6, step = 24000;
    Real m = 1/81.45, T = 17.06521656015796, timeSt = T/step;
    Real u[] ={0.994, 0, 0, 0, -2.0015851063790825224, 0};
    ThreeBody<1> tb1(m, timeSt, u, AM);
    ThreeBody<2> tb2(m, timeSt, u, AM);
    ThreeBody<3> tb3(m, timeSt, u, AM);
    ThreeBody<4> tb4(m, timeSt, u, AM);
    std::unique_ptr<Real[]> out1(tb1.solver(T));
    std::unique_ptr<Real[]> out2(tb2.solver(T));
    std::unique_ptr<Real[]> out3(tb3.solver(T));
    std::unique_ptr<Real[]> out4(tb4.solver(T));
    fstream of1("../output/A-Moulton/p1init1.txt",fstream::out | fstream::trunc);
    fstream of2("../output/A-Moulton/p2init1.txt",fstream::out | fstream::trunc);
    fstream of3("../output/A-Moulton/p3init1.txt",fstream::out | fstream::trunc);
    fstream of4("../output/A-Moulton/p4init1.txt",fstream::out | fstream::trunc);
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of1 << out1[i * size + k] << "\t";
        }
        of1 << "\n";
    }
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of2 << out2[i * size + k] << "\t";
        }
        of2 << "\n";
    }
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of3 << out3[i * size + k] << "\t";
        }
        of3 << "\n";
    }
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of4 << out4[i * size + k] << "\t";
        }
        of4 << "\n";
    }
}

SCENARIO("ThreeBody AM method for initial 2","[AMmethod][initial2]"){
    int size = 6, step = 24000;
    Real m = 1/81.45, T = 19.14045706162071, timeSt = T/step;
    Real u[] ={0.87978, 0, 0, 0, -0.3797, 0};
    ThreeBody<1> tb1(m, timeSt, u, AM);
    ThreeBody<2> tb2(m, timeSt, u, AM);
    ThreeBody<3> tb3(m, timeSt, u, AM);
    ThreeBody<4> tb4(m, timeSt, u, AM);
    std::unique_ptr<Real[]> out1(tb1.solver(T));
    std::unique_ptr<Real[]> out2(tb2.solver(T));
    std::unique_ptr<Real[]> out3(tb3.solver(T));
    std::unique_ptr<Real[]> out4(tb4.solver(T));
    fstream of1("../output/A-Moulton/p1init2.txt",fstream::out | fstream::trunc);
    fstream of2("../output/A-Moulton/p2init2.txt",fstream::out | fstream::trunc);
    fstream of3("../output/A-Moulton/p3init2.txt",fstream::out | fstream::trunc);
    fstream of4("../output/A-Moulton/p4init2.txt",fstream::out | fstream::trunc);
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of1 << out1[i * size + k] << "\t";
        }
        of1 << "\n";
    }
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of2 << out2[i * size + k] << "\t";
        }
        of2 << "\n";
    }
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of3 << out3[i * size + k] << "\t";
        }
        of3 << "\n";
    }
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of4 << out4[i * size + k] << "\t";
        }
        of4 << "\n";
    }
}