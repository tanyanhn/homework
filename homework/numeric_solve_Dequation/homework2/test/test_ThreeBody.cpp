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



SCENARIO("ThreeBody::f","[f]"){
    int size = 6;
    Real m = 1/81.45;
    Real u[] = {1, 2, 3, 4, 5, 6},
            aim[] = { 4.000000000000000 ,  5.000000000000000  , 6.000000000000000 ,
        10.980959986346534 , -6.038135678931384 , -0.057203518397077};
    Real* v = new Real[size];
    ThreeBody<2> tb(m, 1, u, RK);
    tb.f(u,v);
    CHECK_THAT(vector<Real>(v, v+size), EqualsRange(vector<Real>(aim, aim+size), 1e-5));
}


SCENARIO("ThreeBody RK method for initial 1","[RKmethod][initial1]"){
    int size = 6, step = 24000;
    Real m = 1/81.45, T = 17.06521656015796, timeSt = T/step;
    Real u[] ={0.994, 0, 0, 0, -2.0015851063790825224, 0};
    ThreeBody<1> tb1(m, timeSt, u, RK);
    ThreeBody<2> tb2(m, timeSt, u, RK);
    std::unique_ptr<Real[]> out1(tb1.solver(T));
    fstream of1("../output/Runge-Kutta/modifyEuler1.txt",fstream::out | fstream::trunc);
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of1 << out1[i * size + k] << "\t";
        }
        of1 << "\n";
    }
    std::unique_ptr<Real[]> out2(tb2.solver(T));
    fstream of2("../output/Runge-Kutta/improveEuler1.txt",fstream::out | fstream::trunc);
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of2 << out2[i * size + k] << "\t";
        }
        of2 << "\n";
    }
    step = 24000; timeSt = T/step;
    ThreeBody<4> tb4(m, timeSt, u, RK);
    std::unique_ptr<Real[]> out4(tb4.solver(T));
    fstream of4("../output/Runge-Kutta/classicalRK1.txt",fstream::out | fstream::trunc);
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of4 << out4[i * size + k] << "\t";
        }
        of4 << "\n";
    }
    // step = 2400000; timeSt = T/step;
    // ThreeBody<4> tb4r(m, timeSt, u, RK);
    // std::unique_ptr<Real[]> out4r(tb4r.solver(T));
    // fstream of4r("../output/1right.txt",fstream::out | fstream::trunc);
    // for(auto i = 0; i <= step; ++i){
    //     for(auto k = 0; k < size; ++k){
    //         of4r << out4r[i * size + k] << "\t";
    //     }
    //     of4r << "\n";
    // }
}

SCENARIO("ThreeBody RK method for initial 2","[RKmethod][initial2]"){
    int size = 6, step = 24000;
    Real m = 1/81.45, T = 19.14045706162071, timeSt = T/step;
    Real u[] ={0.87978, 0, 0, 0, -0.3797, 0};
    ThreeBody<1> tb1(m, timeSt, u, RK);
    ThreeBody<2> tb2(m, timeSt, u, RK);
    std::unique_ptr<Real[]> out1(tb1.solver(T));
    fstream of1("../output/Runge-Kutta/modifyEuler2.txt",fstream::out | fstream::trunc);
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of1 << out1[i * size + k] << "\t";
        }
        of1 << "\n";
    }
    std::unique_ptr<Real[]> out2(tb2.solver(T));
    fstream of2("../output/Runge-Kutta/improveEuler2.txt",fstream::out | fstream::trunc);
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of2 << out2[i * size + k] << "\t";
        }
        of2 << "\n";
    }
    step = 24000; timeSt = T/step;
    ThreeBody<4> tb4(m, timeSt, u, RK);
    std::unique_ptr<Real[]> out4(tb4.solver(T));
    fstream of4("../output/Runge-Kutta/classicalRK2.txt",fstream::out | fstream::trunc);
    for(auto i = 0; i <= step; ++i){
        for(auto k = 0; k < size; ++k){
            of4 << out4[i * size + k] << "\t";
        }
        of4 << "\n";
    }
    // step = 2400000; timeSt = T/step;
    // ThreeBody<4> tb4r(m, timeSt, u, RK);
    // std::unique_ptr<Real[]> out4r(tb4r.solver(T));
    // fstream of4r("../output/2right.txt",fstream::out | fstream::trunc);
    // for(auto i = 0; i <= step; ++i){
    //     for(auto k = 0; k < size; ++k){
    //         of4r << out4r[i * size + k] << "\t";
    //     }
    //     of4r << "\n";
    // }
}

SCENARIO("point of point operaor","[point]"){
    Real y0[3] = {1, 2, 3},
            y1[3] = {11, 12, 13},
            y2[4] = {21, 22, 23, 24},
            y3[2] = {31, 32};
    Real* y[] = {y0, y1, y2, y3};
    CHECK(y[2][3] == 24);
    CHECK(y[3][1] == 32);
    y[2] = y[0];
    CHECK(y[2][1] == 2);
}



