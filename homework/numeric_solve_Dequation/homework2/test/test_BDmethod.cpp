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

SCENARIO("ThreeBody BDF method for initial 1","[BDmethod][initial1]"){
    int size = 6, step = 170000;
    Real m = 1/81.45, T = 17.06521656015796, timeSt = T/step;
    Real u[] ={0.994, 0, 0, 0, -2.0015851063790825224, 0};
    ThreeBody<1> tb1(m, timeSt, u, BDF);
    ThreeBody<2> tb2(m, timeSt, u, BDF);
    ThreeBody<3> tb3(m, timeSt, u, BDF);
    ThreeBody<4> tb4(m, timeSt, u, BDF);
    std::unique_ptr<Real[]> out1(tb1.solver(T));
    std::unique_ptr<Real[]> out2(tb2.solver(T));
    std::unique_ptr<Real[]> out3(tb3.solver(T));
    std::unique_ptr<Real[]> out4(tb4.solver(T));
    fstream of1("../output/BDF/p1init1.txt",fstream::out | fstream::trunc);
    fstream of2("../output/BDF/p2init1.txt",fstream::out | fstream::trunc);
    fstream of3("../output/BDF/p3init1.txt",fstream::out | fstream::trunc);
    fstream of4("../output/BDF/p4init1.txt",fstream::out | fstream::trunc);
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

SCENARIO("ThreeBody BDF method for initial 2","[BDmethod][initial2]"){
    int size = 6, step = 24000;
    Real m = 1/81.45, T = 19.14045706162071, timeSt = T/step;
    Real u[] ={0.87978, 0, 0, 0, -0.3797, 0};
    ThreeBody<1> tb1(m, timeSt, u, BDF);
    ThreeBody<2> tb2(m, timeSt, u, BDF);
    ThreeBody<3> tb3(m, timeSt, u, BDF);
    ThreeBody<4> tb4(m, timeSt, u, BDF);
    std::unique_ptr<Real[]> out1(tb1.solver(T));
    std::unique_ptr<Real[]> out2(tb2.solver(T));
    std::unique_ptr<Real[]> out3(tb3.solver(T));
    std::unique_ptr<Real[]> out4(tb4.solver(T));
    fstream of1("../output/BDF/p1init2.txt",fstream::out | fstream::trunc);
    fstream of2("../output/BDF/p2init2.txt",fstream::out | fstream::trunc);
    fstream of3("../output/BDF/p3init2.txt",fstream::out | fstream::trunc);
    fstream of4("../output/BDF/p4init2.txt",fstream::out | fstream::trunc);
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