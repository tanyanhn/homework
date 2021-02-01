#include "../src/cycle.h"
#include "../src/possionsolver.h"
#include "equalrangematcher.h"
#include "catch.hpp"
#include "../src/norm.h"

using namespace std;
using namespace solver_ty;


SCENARIO("norm","[norm]"){
    Real v[]{-1,7};
    REQUIRE(norm<2>(v, 2) == 5);
    REQUIRE(norm<1>(v, 2) == 4);
}

SCENARIO("Vcycle, sin(x)","[cycle][Vcycle][sin(x)]"){
    int v02[]{0, 2}, v12[]{1, 2}, v20[]{2, 0}, v21[]{2, 1}, v22[]{2, 2};
    int loops = 15,
            num = 256;
    GIVEN(" -(sin)\" = sin"){
        WHEN("Dim = 1"){
            int n[]{num};
            Real h[]{M_PI / n[0]};
            int size = n[0] + 1;
            unique_ptr<Real[]> f(new Real[size]);
            vector<Real> aim(size);
            for(auto k = 0; k < size; ++k){
                f[k] = sin(k * h[0]);
                aim[k] = sin(k * h[0]);
            }
            WHEN("using GuassSediel"){
                possionsolver<1> psolver21(n, h, move(f), nullptr);
                psolver21.chooseoption("relax_type", GaussSeidel);
                psolver21.solver(v21, loops);
                auto result = psolver21.getresult();
                auto e = move(result.first), r = move(result.second);
                vector<Real> ve(size), vr(size);
                for(auto k = 0; k < size; ++k){
                    e[k] -= aim[k];
                    ve[k] = e[k];
                    vr[k] = r[k];
                }
                Real norme = norm<2>(e.get(), size),
                        normr = norm<2>(r.get(), size);
                CHECK(norme < 1e-5);
                CHECK(normr < 1e-7);
            }
            WHEN("using Jacobirelax"){
                possionsolver<1> psolver21(n, h, move(f), nullptr);
                psolver21.chooseoption("relax_type", Jacobi);
                psolver21.solver(v21, loops);
                auto result = psolver21.getresult();
                auto e = move(result.first), r = move(result.second);
                vector<Real> ve(size), vr(size);
                for(auto k = 0; k < size; ++k){
                    e[k] -= aim[k];
                    ve[k] = e[k];
                    vr[k] = r[k];
                }
                Real norme = norm<2>(e.get(), size),
                        normr = norm<2>(r.get(), size);
                CHECK(norme < 1e-4);
                CHECK(normr < 1);
                //CHECK_THAT(vector<Real>(get<0>(result).get(), get<0>(result).get()+size), EqualsRange(aim, 1e-1));
            }
        }
        WHEN("Dim = 2"){
            int n[]{num, num};
            Real h[]{M_PI / n[0], M_PI / n[1]};
            int size = (n[0] + 1) * (n[1] + 1);
            unique_ptr<Real[]> aim(new Real[size]),
                    f(new Real[size]);
            for(auto i = 0; i < n[0] + 1; ++i){
                for(auto j = 0; j < n[1] + 1; ++j){
                    int id = i*(n[1] + 1) + j;
                    f[id] = 2 * sin(i * h[0]) * sin(j * h[1]);
                    aim[id] = sin(i * h[0]) * sin(j * h[1]);
                }
            }
            WHEN("using GaussSeidel"){
                possionsolver<2> psolver21(n, h, move(f), nullptr);
                psolver21.chooseoption("relax_type", GaussSeidel);
                psolver21.solver(v21, loops);
                auto result = psolver21.getresult();
                auto e = move(result.first), r = move(result.second);
                vector<Real> ve(size), vr(size);
                for(auto k = 0; k < size; ++k){
                    e[k] -= aim[k];
                    ve[k] = e[k];
                    vr[k] = r[k];
                }
                Real norme = norm<2>(e.get(), size),
                        normr = norm<2>(r.get(), size);
                CHECK(norme < 1e-5);
                CHECK(normr < 1e-7);
            }
            WHEN("using Jacobirelax"){
                possionsolver<2> psolver21(n, h, move(f), nullptr);
                psolver21.chooseoption("relax_type", Jacobi);
                psolver21.solver(v21, loops);
                auto result = psolver21.getresult();
                auto e = move(result.first), r = move(result.second);
                vector<Real> ve(size), vr(size);
                for(auto k = 0; k < size; ++k){
                    e[k] -= aim[k];
                    ve[k] = e[k];
                    vr[k] = r[k];
                }
                Real norme = norm<2>(e.get(), size),
                        normr = norm<2>(r.get(), size);
                CHECK(norme < 1e-5);
                CHECK(normr < 1e-4);
            }
        }
    }
}

SCENARIO("Vcycle , exp(sin(x))","[Vcycle][exp(sin(x))]"){
    int v02[]{0, 2}, v12[]{1, 2}, v20[]{2, 0}, v21[]{2, 1}, v22[]{2, 2};
    int loops = 15,
            num =1024;
    GIVEN(" -exp(sin)\" = -sin * sin * exp(sin) + cos * exp(sin)"){
        WHEN("Dim = 1"){
            int n[]{num};
            Real h[]{1.0 / n[0]};
            int size = n[0] + 1;
            unique_ptr<Real[]> f(new Real[size]),
                    v(new Real[size]);
            vector<Real> aim(size);
            for(auto k = 0; k < size; ++k){
                f[k] = - cos(k * h[0]) * cos(k * h[0]) * exp(sin(k * h[0])) +
                        sin(k * h[0]) * exp(sin(k * h[0]));
                aim[k] = exp(sin(k * h[0]));
                v[k] =  0; //(expm1( sin(1) )) * k / (size - 1) + 1;
            }
            v[0] = 1; v[size-1] = exp(sin(1));
            WHEN("using GuassSediel"){
                possionsolver<1> psolver21(n, h, move(f), move(v));
                psolver21.chooseoption("relax_type", GaussSeidel);
                psolver21.solver(v21, loops);
                auto result = psolver21.getresult();
                auto e = move(result.first), r = move(result.second);
                vector<Real> ve(size), vr(size);
                for(auto k = 0; k < size; ++k){
                    e[k] -= aim[k];
                    ve[k] = e[k];
                    vr[k] = r[k];
                }
                Real norme = norm<2>(e.get(), size),
                        normr = norm<2>(r.get(), size);
                CHECK(norme < 1e-7);
                CHECK(normr < 1e-6);
            }
            WHEN("using Jacobirelax"){
                possionsolver<1> psolver21(n, h, move(f), move(v));
                psolver21.chooseoption("relax_type", weightJacobi);
                psolver21.solver(v21, loops);
                auto result = psolver21.getresult();
                auto e = move(result.first), r = move(result.second);
                vector<Real> ve(size), vr(size);
                for(auto k = 0; k < size; ++k){
                    e[k] -= aim[k];
                    ve[k] = e[k];
                    vr[k] = r[k];
                }
                Real norme = norm<2>(e.get(), size),
                        normr = norm<2>(r.get(), size);
                CHECK(norme < 1e-7);
                CHECK(normr < 1e-6);
            }
        }
        WHEN("Dim = 2"){
            int n[]{num, num};
            Real h[]{1.0 / n[0], 1.0 / n[1]};
            int size = (n[0] + 1) * (n[1] + 1);
            unique_ptr<Real[]> aim(new Real[size]),
                    f(new Real[size]),
                    v(new Real[size]);
            for(auto i = 0; i < n[0] + 1; ++i){
                for(auto j = 0; j < n[1] + 1; ++j){
                    int id = i*(n[1] + 1) + j;
                    f[id] = - cos(i * h[0]) * cos(i * h[0]) * sin(j * h[1]) * sin(j * h[1]) *
                            exp(sin(i * h[0]) * sin(j * h[1]))
                            - cos(j * h[1]) * cos(j * h[1]) * sin(i * h[0]) * sin(i * h[0]) *
                            exp(sin(i * h[0]) * sin(j * h[1])) +
                            sin(i * h[0]) * sin(j * h[1]) * exp(sin(i * h[0]) * sin(j * h[1])) +
                            sin(j * h[1]) * sin(i * h[0]) * exp(sin(i * h[0]) * sin(j * h[1]));
                    aim[id] = exp(sin(i * h[0]) * sin(j * h[1]));
                    if(i == 0 || i == n[0] || j == 0 || j == n[1]){
                        v[id] = aim[id];
                    }
                    else{
                        v[id] = 0;
                    }
                }
            }
            WHEN("using GaussSeidel"){
                possionsolver<2> psolver21(n, h, move(f), move(v));
                psolver21.chooseoption("relax_type", RB_GaussSeidel);
                //for(auto i = 0; i < loops; ++i){
                //   psolver21.solver(v22, 1);
                    psolver21.solver(v22, loops);
                    auto result = psolver21.getresult();
                    auto e = move(result.first), r = move(result.second);
                    vector<Real> ve(size), vr(size);
                    for(auto k = 0; k < size; ++k){
                        e[k] -= aim[k];
                        ve[k] = e[k];
                        vr[k] = r[k];
                    }
                    Real norme = norm<2>(e.get(), size),
                            normr = norm<2>(r.get(), size);
                    //   INFO("i is " << i);
                    CHECK(norme < 1e-6);
                    CHECK(normr < 1e-6);
                    // }
            }
            WHEN("using Jacobirelax"){
                possionsolver<2> psolver21(n, h, move(f), move(v));
                psolver21.chooseoption("relax_type", weightJacobi);
                //for(auto i = 0; i < loops; ++i){
                //psolver21.solver(v22, 1);
                psolver21.solver(v22, loops);
                auto result = psolver21.getresult();
                auto e = move(result.first), r = move(result.second);
                vector<Real> ve(size), vr(size);
                for(auto k = 0; k < size; ++k){
                    e[k] -= aim[k];
                    ve[k] = e[k];
                    vr[k] = r[k];
                }
                Real norme = norm<2>(e.get(), size),
                        normr = norm<2>(r.get(), size);
                // INFO("i is " << i);
                CHECK(norme < 1e-6);
                CHECK(normr < 1e-5);
                //}
            }
        }
    }
}

SCENARIO("FMGcycle , exp(sin(x))","[FMGcycle][exp(sin(x))]"){
    int v21[]{1, 2, 1}, v22[]{1, 2, 2};
    int loops = 3,
            num = 256;
    GIVEN(" -exp(sin)\" = -sin * sin * exp(sin) + cos * exp(sin)"){
        WHEN("Dim = 1"){
            int n[]{num};
            Real h[]{1.0 / n[0]};
            int size = n[0] + 1;
            unique_ptr<Real[]> f(new Real[size]),
                    v(new Real[size]);
            vector<Real> aim(size);
            for(auto k = 0; k < size; ++k){
                f[k] = - cos(k * h[0]) * cos(k * h[0]) * exp(sin(k * h[0])) +
                        sin(k * h[0]) * exp(sin(k * h[0]));
                aim[k] = exp(sin(k * h[0]));
                v[k] =  0; //(expm1( sin(1) )) * k / (size - 1) + 1;
            }
            v[0] = 1; v[size-1] = exp(sin(1));
            WHEN("using GuassSediel"){
                possionsolver<1> psolver21(n, h, move(f), move(v));
                psolver21.chooseoption("relax_type", GaussSeidel);
                psolver21.chooseoption("cycle_type", full_multi_cycle);
                psolver21.solver(v21, loops);
                auto result = psolver21.getresult();
                auto e = move(result.first), r = move(result.second);
                vector<Real> ve(size), vr(size);
                for(auto k = 0; k < size; ++k){
                    e[k] -= aim[k];
                    ve[k] = e[k];
                    vr[k] = r[k];
                }
                Real norme = norm<2>(e.get(), size),
                        normr = norm<2>(r.get(), size);
                CHECK(norme < 1e-5);
                //CHECK(normr < 1e-10);
            }
        }
        WHEN("Dim = 2"){
            int n[]{num, num};
            Real h[]{1.0 / n[0], 1.0 / n[1]};
            int size = (n[0] + 1) * (n[1] + 1);
            unique_ptr<Real[]> aim(new Real[size]),
                    f(new Real[size]),
                    v(new Real[size]);
            for(auto i = 0; i < n[0] + 1; ++i){
                for(auto j = 0; j < n[1] + 1; ++j){
                    int id = i*(n[1] + 1) + j;
                    f[id] = - cos(i * h[0]) * cos(i * h[0]) * sin(j * h[1]) * sin(j * h[1]) *
                            exp(sin(i * h[0]) * sin(j * h[1]))
                            - cos(j * h[1]) * cos(j * h[1]) * sin(i * h[0]) * sin(i * h[0]) *
                            exp(sin(i * h[0]) * sin(j * h[1])) +
                            sin(i * h[0]) * sin(j * h[1]) * exp(sin(i * h[0]) * sin(j * h[1])) +
                            sin(j * h[1]) * sin(i * h[0]) * exp(sin(i * h[0]) * sin(j * h[1]));
                    aim[id] = exp(sin(i * h[0]) * sin(j * h[1]));
                    if(i == 0 || i == n[0] || j == 0 || j == n[1]){
                        v[id] = aim[id];
                    }
                    else{
                        v[id] = 0;
                    }
                }
            }
            WHEN("using GaussSeidel"){
                possionsolver<2> psolver21(n, h, move(f), move(v));
                psolver21.chooseoption("relax_type", GaussSeidel);
                psolver21.chooseoption("cycle_type", full_multi_cycle);
                psolver21.solver(v21, loops);
                auto result = psolver21.getresult();
                auto e = move(result.first), r = move(result.second);
                vector<Real> ve(size), vr(size);
                for(auto k = 0; k < size; ++k){
                    e[k] -= aim[k];
                    ve[k] = e[k];
                    vr[k] = r[k];
                }
                Real norme = norm<2>(e.get(), size),
                        normr = norm<2>(r.get(), size);
                CHECK(norme < 1e-5);
                //CHECK(normr < 1e-11);
            }
        }
    }
}