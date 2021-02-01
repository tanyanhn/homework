#include "catch.hpp"
#include "../src/possionsolver.h"
#include <cmath>
#include "equalrangematcher.h"
#include "../src/is_odd.h"

using namespace std;
using namespace solver_ty;



SCENARIO("possionsolver::updatarvalue","[possionsolver][updatarvalue]"){
    WHEN("Dim = 1"){
        int n[]{64};
        Real h[]{M_PI / n[0]};
        int size = n[0] + 1;
        unique_ptr<Real[]> v(new Real[n[0] + 1]),
                f(new Real[n[0] + 1]);
        for(auto k = 0; k < size; ++k){
            f[k] = sin(k * h[0]);
            v[k] = sin(k * h[0]);
        }
        possionsolver<1> psolver1(n, h, move(f), move(v));
        psolver1.updatervalue();
        auto rvalue = psolver1.test_getvalue("rvalue");
        CHECK_THAT(vector<Real>(rvalue.get(), rvalue.get()+size), EqualsRange(vector<Real>(size, 0.0), 1e-3));
    }
    WHEN("Dim = 2"){
        int n[]{64, 64};
        Real h[]{M_PI / n[0], M_PI / n[1]};
        int size = (n[0] + 1) * (n[1] + 1);
        unique_ptr<Real[]> v(new Real[size]),
                f(new Real[size]);
        for(auto i = 0; i < n[0] + 1; ++i){
            for(auto j = 0; j< n[1] + 1; ++j){
                int id = i*(n[1] + 1) + j;
                f[id] = 2 * sin(i * h[0]) * sin(j * h[1]);
                v[id] = sin(i * h[0]) * sin(j * h[1]);
            }
        }
        possionsolver<2> psolver1(n, h, move(f), move(v));
        psolver1.updatervalue();
        auto rvalue = psolver1.test_getvalue("rvalue");
        CHECK_THAT(vector<Real>(rvalue.get(), rvalue.get()+size), EqualsRange(vector<Real>(size, 0.0), 1e-3));
    }
}


SCENARIO("possionsolver::Relax", "[possionsolver][relax]"){
    WHEN("Dim = 1"){
        int n[]{64};
        Real h[]{M_PI / n[0]};
        int size = n[0] + 1;
        THEN("GaussSeidel relax"){
            unique_ptr<Real[]> v(new Real[size]),
                    f(new Real[size]),
                    aimv1(new Real[size]),
                    aimf1(new Real[size]);
            for(auto k = 0; k < size; ++k){
                f[k] = sin(k * h[0]);
                v[k] = sin(k * h[0]);
                aimv1[k] = v[k];
                aimf1[k] = f[k];
            }
            possionsolver<1> psolver1(n, h, move(f), move(v));
            psolver1.chooseoption("relax_type", GaussSeidel);
            GIVEN("relax 1 times"){
                for(auto k = 1; k < size-1; ++k){
                    aimv1[k] = (aimv1[k-1] + aimv1[k+1] + h[0] * h[0] * aimf1[k]) / 2;
                }
                psolver1.Relax();
                auto outv1 = psolver1.test_getvalue("value");
                REQUIRE_THAT(vector<Real>(aimv1.get(), aimv1.get()+size), EqualsRange(vector<Real>(outv1.get(), outv1.get()+size), 1e-15));
            }
            GIVEN("relax 5 times"){
                for(auto t = 0; t < 5; ++t){
                    for(auto k = 1; k < size-1; ++k){
                        aimv1[k] = (aimv1[k-1] + aimv1[k+1] + h[0] * h[0] * aimf1[k]) / 2;
                    }
                    psolver1.Relax();
                }
                auto outv1 = psolver1.test_getvalue("value");
                REQUIRE_THAT(vector<Real>(aimv1.get(), aimv1.get()+size), EqualsRange(vector<Real>(outv1.get(), outv1.get()+size), 1e-15));
            }
        }
    }
    WHEN("Dim = 2"){
        int n[]{4, 4};
        Real h[]{4.0 / n[0], 4.0 / n[1]};
        int size = (n[0] + 1) * (n[1] + 1);
        THEN("GuassSediel relax"){
            unique_ptr<Real[]> v(new Real[size]),
                    f(new Real[size]),
                    aimv(new Real[size]),
                    aimf(new Real[size]);
            for(auto k = 0; k < size; ++k){
                f[k] = 0;
                v[k] = 0;
                aimv[k] = v[k];
                aimf[k] = f[k];
            }
            int count = 1;
            for(auto i = 1; i < n[0]; ++i){
                for(auto j = 1; j < n[1]; ++j){
                    int id = i * (n[1] + 1) + j;
                    f[id] = count++;
                    v[id] = f[id];
                    aimv[id] = f[id];
                    aimf[id] = f[id];
                }
            }
            possionsolver<2> psolver1(n, h, move(f), move(v));
            psolver1.chooseoption("relax_type", GaussSeidel);
            GIVEN("relax 1 times"){
                aimv[6] = (aimv[1] + aimv[5] + aimv[7] + aimv[11] + h[0]*h[0]*aimf[6]) / 4;
                aimv[7] = (aimv[2] + aimv[6] + aimv[8] + aimv[12] + h[0]*h[0]*aimf[7]) / 4;
                aimv[8] = (aimv[3] + aimv[7] + aimv[9] + aimv[13] + h[0]*h[0]*aimf[8]) / 4;
                aimv[11] = (aimv[6] + aimv[10] + aimv[12] + aimv[16] + h[0]*h[0]*aimf[11]) / 4;
                aimv[12] = (aimv[7] + aimv[11] + aimv[13] + aimv[17] + h[0]*h[0]*aimf[12]) / 4;
                aimv[13] = (aimv[8] + aimv[12] + aimv[14] + aimv[18] + h[0]*h[0]*aimf[13]) / 4;
                aimv[16] = (aimv[11] + aimv[15] + aimv[17] + aimv[21] + h[0]*h[0]*aimf[16]) / 4;
                aimv[17] = (aimv[12] + aimv[16] + aimv[18] + aimv[22] + h[0]*h[0]*aimf[17]) / 4;
                aimv[18] = (aimv[13] + aimv[17] + aimv[19] + aimv[23] + h[0]*h[0]*aimf[18]) / 4;
                psolver1.Relax();
                auto outv = psolver1.test_getvalue("value");
                REQUIRE_THAT(vector<Real>(aimv.get(), aimv.get()+size), EqualsRange(vector<Real>(outv.get(), outv.get()+size), 1e-15));
            }
        }
    }
}


SCENARIO("possionsolver::Restrict","[possionsolver][restrict]"){
    WHEN("Dim = 1"){
        int n[]{64};
        Real h[]{M_PI / n[0]};
        int finesize = n[0] + 1,
                coarsize = n[0] / 2 + 1;
        THEN("full_restrict"){
            unique_ptr<Real[]> v(new Real[finesize]),
                    f(new Real[finesize]),
                    r(new Real[finesize]),
                    aimv1(new Real[coarsize]),
                    aimf1(new Real[coarsize]);
            for(auto k = 0; k < finesize; ++k){
                f[k] = sin(k * h[0]);
                v[k] = sin(k * h[0]);
                r[k] = 0;
            }
            for(auto k = 0; k < coarsize; ++k){
                aimv1[k] = 0;
                aimf1[k] = 0;
            }
            for(auto k = 1; k < finesize - 1; ++k){
                r[k] = f[k] + (v[k - 1] + v[k + 1] - 2 * v[k]) / (h[0] * h[0]);
            }
            GIVEN("restrict 1 times"){
                for(auto k = 1; k < coarsize-1; ++k){
                    //aimv1[k] = (v[2 * k - 1] + v[2 * k + 1] + 2 * v[2 * k]) / 4;
                    aimf1[k] = (r[2 * k - 1] + r[2 * k + 1] + 2 * r[2 * k]) / 4;
                }
                possionsolver<1> psolver1(n, h, move(f), move(v));
                auto coarsolver1 = psolver1.Restrict();
                auto outv1 = coarsolver1->test_getvalue("value"),
                        outf1 = coarsolver1->test_getvalue("fvalue");
                //REQUIRE_THAT(vector<Real>(aimv1.get(), aimv1.get()+coarsize), EqualsRange(vector<Real>(outv1.get(), outv1.get()+coarsize), 1e-10));
                REQUIRE_THAT(vector<Real>(aimf1.get(), aimf1.get()+coarsize), EqualsRange(vector<Real>(outf1.get(), outf1.get()+coarsize), 1e-10));
            }
        }
    }
    WHEN("Dim = 2"){
        int n[]{128, 128};
        Real h[]{M_PI / n[0], M_PI /n[1]};
        int finesize = (n[0] + 1) * (n[1] + 1),
                coarsize = (n[0] / 2 + 1) * (n[1] / 2 + 1);
        THEN("full_restrict"){
            unique_ptr<Real[]> v(new Real[finesize]),
                    f(new Real[finesize]),
                    aimv2(new Real[coarsize]),
                    aimf2(new Real[coarsize]);
            for(auto i = 0; i < n[0] + 1; ++i){
                for(auto j = 0; j< n[1] + 1; ++j){
                    int id = i*(n[1] + 1) + j;
                    f[id] = 2 * sin(i * h[0]) * sin(j * h[1]);
                    v[id] = sin(i * h[0]) * sin(j * h[1]);
                }
            }
            for(auto k = 0; k < coarsize; ++k){
                aimv2[k] = 0;
                aimf2[k] = 0;
            }
            GIVEN("restrict 1 times"){
                // for(auto i = 1; i < n[0] / 2; ++i){
                //     for(auto j = 1; j < n[1] / 2; ++j){
                //         int coarseid = i * (n[1] / 2 + 1) + j,
                //                 id = 2 * i * (n[1] + 1) + 2 * j;
                //         aimv2[coarseid] = ((v[id - n[1] - 2] + v[id - n[1]] +
                //                             v[id + n[1]] + v[id + n[1] + 2]) +
                //                            (v[id - n[1] - 1] + v[id - 1] +
                //                             v[id + 1] + v[id + n[1] + 1]) * 2 +
                //                            v[id] * 4) / 16;
                //     }
                // }
                possionsolver<2> psolver2(n, h, move(f), move(v));
                psolver2.updatervalue();
                auto r = psolver2.test_getvalue("rvalue");
                for(auto i = 1; i < n[0] / 2; ++i){
                    for(auto j = 1; j < n[1] / 2; ++j){
                        int coarseid = i * (n[1] / 2 + 1) + j,
                                id = 2 * i * (n[1] + 1) + 2 * j;
                        aimf2[coarseid] = ((r[id - n[1] - 2] + r[id - n[1]] +
                                            r[id + n[1]] + r[id + n[1] + 2]) +
                                           (r[id - n[1] - 1] + r[id - 1] +
                                            r[id + 1] + r[id + n[1] + 1]) * 2 +
                                           r[id] * 4) / 16;
                    }
                }
                auto coarsolver2 = psolver2.Restrict();
                auto outv1 = coarsolver2->test_getvalue("value"),
                        outf1 = coarsolver2->test_getvalue("fvalue");
                //REQUIRE_THAT(vector<Real>(aimv2.get(), aimv2.get()+coarsize), EqualsRange(vector<Real>(outv1.get(), outv1.get()+coarsize), 1e-10));
                REQUIRE_THAT(vector<Real>(aimf2.get(), aimf2.get()+coarsize), EqualsRange(vector<Real>(outf1.get(), outf1.get()+coarsize), 1e-10));
            }
        }
    }
}


SCENARIO("possionsolver::Interpolate","[possionsolver][.interpolate]"){
    WHEN("Dim = 1"){
        int n[]{64};
        Real h[]{M_PI / n[0]};
        int coarsize = n[0] + 1,
                finesize = n[0] * 2 + 1;
        THEN("Linearinterpolate"){
            unique_ptr<Real[]> v(new Real[coarsize]),
                    f(new Real[coarsize]),
                    aimv1(new Real[finesize]);
            for(auto k = 0; k < coarsize; ++k){
                f[k] = sin(k * h[0]);
                v[k] = sin(k * h[0]);
            }
            GIVEN("interpolate 1 times"){
                aimv1[2 * coarsize - 2] = v[coarsize - 1]; 
                for(auto k = 0; k < coarsize - 1; ++k){
                    aimv1[2 * k] = v[k];
                    aimv1[2 * k + 1] = (v[k] + v[k + 1]) / 2;
                }
                possionsolver<1> psolver1(n, h, move(f), move(v));
                auto coarsolver1 = psolver1.Interpolate();
                        // coarsolver2 = psolver1.Linearinterpolate(),
                        // coarsolver3 = psolver1.Linearinterpolatet();
                auto outv1 = coarsolver1->test_getvalue("value");
                        // outv2 = coarsolver2->test_getvalue("value"),
                        // outv3 = coarsolver3->test_getvalue("value");
                REQUIRE_THAT(vector<Real>(aimv1.get(), aimv1.get()+coarsize), EqualsRange(vector<Real>(outv1.get(), outv1.get()+coarsize), 1e-10));
                // REQUIRE_THAT(vector<Real>(aimv1.get(), aimv1.get()+coarsize), EqualsRange(vector<Real>(outv2.get(), outv2.get()+coarsize), 1e-10));
                // REQUIRE_THAT(vector<Real>(aimv1.get(), aimv1.get()+coarsize), EqualsRange(vector<Real>(outv3.get(), outv3.get()+coarsize), 1e-10));
            }
        }
    }
    WHEN("Dim = 2"){
        int n[]{2, 2};
        Real h[]{M_PI / n[0], M_PI / n[1]};
        int coarsize = (n[0] + 1) * (n[1] + 1),
                finesize = (n[0] * 2 + 1) * (n[1] * 2 + 1);
        THEN("Linearinterpolate"){
            unique_ptr<Real[]> v(new Real[coarsize]),
                    f(new Real[coarsize]),
                    aimv2(new Real[finesize]);
            for(auto i = 0; i < n[0] + 1; ++i){
                for(auto j = 0; j < n[1] + 1; ++j){
                    int id = i * (n[1] + 1) + j;
                    f[id] = 2 * sin(i * h[0]) * sin(j * h[1]);
                    v[id] = sin(i * h[0]) * sin(j * h[1]);
                }
            }
            for(auto k = 0; k < finesize; ++k){
                aimv2[k] = 0;
            }
            GIVEN("interpolate 1 times"){
                REQUIRE(is_odd(99));
                REQUIRE(is_odd(-99));
                REQUIRE(!is_odd(0));
                REQUIRE(!is_odd(-2));
                REQUIRE(!is_odd(100));
                for(auto i = 1; i < n[0] * 2; ++i){
                    for(auto j = 1; j < n[1] * 2; ++j){
                        int fineid = i * (n[1] * 2 + 1) + j,
                                id = (i / 2) * (n[1] + 1) + (j / 2);
                        if(!is_odd(i) && !is_odd(j)){
                            aimv2[fineid] =  v[id];
                        }
                        else if(is_odd(i) && !is_odd(j)){
                            aimv2[fineid] = (v[id] + v[id + n[1] + 1]) / 2;
                        }
                        else if(!is_odd(i) && is_odd(j)){
                            aimv2[fineid] = (v[id] + v[id + 1]) / 2;
                        }
                        else if(is_odd(i) && is_odd(j)){
                            aimv2[fineid] = (v[id] + v[id + n[1] + 1] +
                                             v[id + 1] + v[id + n[1] + 2]) / 4;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                }
                possionsolver<2> psolver1(n, h, move(f), move(v));
                auto coarsolver1 = psolver1.Interpolate();
                        // coarsolver2 = psolver1.Linearinterpolate(),
                        // coarsolver3 = psolver1.Linearinterpolatet();
                auto outv1 = coarsolver1->test_getvalue("value");
                        // outv2 = coarsolver2->test_getvalue("value"),
                        // outv3 = coarsolver3->test_getvalue("value");
                REQUIRE_THAT(vector<Real>(aimv2.get(), aimv2.get()+finesize), EqualsRange(vector<Real>(outv1.get(), outv1.get()+finesize), 1e-10));
                // REQUIRE_THAT(vector<Real>(aimv2.get(), aimv2.get()+finesize), EqualsRange(vector<Real>(outv2.get(), outv2.get()+finesize), 1e-10));
                // REQUIRE_THAT(vector<Real>(aimv2.get(), aimv2.get()+finesize), EqualsRange(vector<Real>(outv3.get(), outv3.get()+finesize), 1e-10));
            }
        }
    }
}
