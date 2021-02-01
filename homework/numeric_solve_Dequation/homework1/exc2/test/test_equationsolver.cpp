#include "../src/equationsolver.h"
#include "catch.hpp"
#include "equalrangematcher.h"

using namespace solver_ty;
using namespace std;

SCENARIO("equationsolver test", "[equationsolver][initiallize][chooseoption]"){
    WHEN("Dim = 1"){
        int n[]{16};
        Real h[]{0.35};
        unique_ptr<Real[]> f(new Real[n[0] + 1]),
                v(new Real[n[0] + 1]);
        equationsolver<1> solver1(n, h, move(f), move(v));
        THEN("chooseoption"){
            vector<string> vs{"boundary_type", "restrict_operator", "interpolat_operator",
                "cycle_type", "stop_criteria", "relax_type"};
            for(auto k = 0; k < 5; ++k){
                for(auto i = 0; i < 2; ++i){
                    REQUIRE(solver1.chooseoption(vs[k], i));
                }
            }
            for(auto i = 0; i < 4; ++i){
                REQUIRE(solver1.chooseoption(vs[5], i));
            }
        }
    }
    WHEN("Dim = 2 without initialize"){
        int n[]{16, 64};
        Real h[]{0.15, 0.05};
        int size = (n[0]+1) * (n[1]+1);
        unique_ptr<Real[]> f(new Real[(n[0]+1) * (n[1]+1)]);
                //v(new Real[(n[0]+1) * (n[1]+1)]);
        equationsolver<2> solver2(n, h, move(f), (nullptr));
        THEN("chooseoption"){
            vector<string> vs{"boundary_type", "restrict_operator", "interpolat_operator",
                "cycle_type", "stop_criteria", "relax_type"};
            for(auto k = 0; k < 5; ++k){
                for(auto i = 0; i < 2; ++i){
                    REQUIRE(solver2.chooseoption(vs[k], i));
                }
            }
            for(auto i = 0; i < 4; ++i){
                REQUIRE(solver2.chooseoption(vs[5], i));
            }
        }
        THEN("Verify inner value without initial"){
            auto sn = solver2.test_getvalue("n"),
                    sh = solver2.test_getvalue("h"),
                    svalue = solver2.test_getvalue("value");
            REQUIRE(sn[0] == n[0]);
            REQUIRE(sn[1] == n[1]);
            REQUIRE_THAT(sh[0], Catch::Matchers::WithinRel(h[0]));
            REQUIRE_THAT(sh[1], Catch::Matchers::WithinRel(h[1]));
            CHECK_THAT(vector<Real>(svalue.get(), svalue.get()+size), EqualsRange(vector<Real>(size, 0.0), 1e-5));
        }
    }
    WHEN("Dim = 2 with initialize"){
        int n[]{16, 64};
        Real h[]{0.15, 0.05};
        int size = (n[0]+1) * (n[1]+1);
        unique_ptr<Real[]> f(new Real[(n[0]+1) * (n[1]+1)]),
                v(new Real[(n[0]+1) * (n[1]+1)]),
                tv(new Real[(n[0]+1) * (n[1]+1)]),
                tf(new Real[(n[0]+1) * (n[1]+1)]);
        for(auto i = 0; i < size; ++i){
            v[i] = i; f[i] = i * i;
        }
        std::copy(v.get(),v.get()+size,tv.get());
        std::copy(f.get(),f.get()+size,tf.get());
        equationsolver<2> solver2(n, h, move(f), move(v));
        THEN("Verify inner value"){
            auto sn = solver2.test_getvalue("n"),
                    sh = solver2.test_getvalue("h"),
                    svalue = solver2.test_getvalue("value"),
                    sfvalue = solver2.test_getvalue("fvalue");
            REQUIRE(sn[0] == n[0]);
            REQUIRE(sn[1] == n[1]);
            REQUIRE_THAT(sh[0], Catch::Matchers::WithinRel(h[0]));
            REQUIRE_THAT(sh[1], Catch::Matchers::WithinRel(h[1]));
            CHECK_THAT(vector<Real>(svalue.get(), svalue.get()+size), EqualsRange(vector<Real>(tv.get(), tv.get()+size), 1e-5));
            CHECK_THAT(vector<Real>(sfvalue.get(), sfvalue.get()+size), EqualsRange(vector<Real>(tf.get(), tf.get()+size), 1e-5));
        }
    }
}


SCENARIO("equationsolver::mathoperator","[equationsolver][mathoperator][operator=]"){
    WHEN("Dim = 1"){
        int n[]{64};
        Real h[]{M_PI / n[0]};
        int size = n[0] + 1;
        unique_ptr<Real[]> vsin(new Real[size]),
                v1(new Real[size]),
                v7(new Real[size]);
        vector<Real> aim(size);
        for(auto k = 0; k < size; ++k){
            vsin[k] = sin(k * h[0]);
            v1[k] = 1;
            v7[k] = 7;
            aim[k] = vsin[k];
        }
        equationsolver<1> psolversin(n, h, nullptr, move(vsin)),
                psolver1(n, h, nullptr, move(v1)),
                psolver7(n, h, nullptr, move(v7)),
                psolver0(n, h, nullptr);
        psolversin.plusvalue(psolver0);
        auto valuesin = psolversin.test_getvalue("value");
        CHECK_THAT(vector<Real>(valuesin.get(), valuesin.get()+size), EqualsRange(aim, 1e-10));
        auto copypsolversin(psolversin);
        auto copyvaluesin = copypsolversin.test_getvalue("value");
        CHECK_THAT(vector<Real>(copyvaluesin.get(), copyvaluesin.get()+size), EqualsRange(aim, 1e-10));
        psolver1.assignvalue(psolver7);
        auto value1 = psolver1.test_getvalue("value");
        CHECK_THAT(vector<Real>(value1.get(), value1.get()+size), EqualsRange(vector<Real>(size, 7.0), 1e-10));
        psolversin.minusvalue(psolver1);
        for(auto k = 0; k < size; ++k){
            aim[k] -= 7.0;
        }
        valuesin.reset(psolversin.test_getvalue("value").release());
        CHECK_THAT(vector<Real>(valuesin.get(), valuesin.get()+size), EqualsRange(aim, 1e-10));
        copypsolversin = psolversin;
        auto copyvaluesin2 = copypsolversin.test_getvalue("value");
        CHECK_THAT(vector<Real>(copyvaluesin2.get(), copyvaluesin2.get()+size), EqualsRange(aim, 1e-10));
    }
    WHEN("Dim = 2"){
        int n[]{64, 64};
        Real h[]{M_PI / n[0], M_PI / n[1]};
        int size = (n[0] + 1) * (n[1] + 1);
        unique_ptr<Real[]> vsin(new Real[size]),
                v1(new Real[size]),
                v7(new Real[size]);
        vector<Real> aim(size);
        for(auto i = 0; i < n[0] + 1; ++i){
            for(auto j = 0; j< n[1] + 1; ++j){
                int id = i*(n[1] + 1) + j;
                vsin[id] = sin(i * h[0]) * sin(j * h[1]);
                v1[id] = 1;
                v7[id] = 7;
                aim[id] = vsin[id];
            }
        }
        equationsolver<2> psolversin(n, h, nullptr, move(vsin)),
                psolver1(n, h, nullptr, move(v1)),
                psolver7(n, h, nullptr, move(v7)),
                psolver0(n, h, nullptr);
        psolversin.plusvalue(psolver0);
        auto valuesin = psolversin.test_getvalue("value");
        CHECK_THAT(vector<Real>(valuesin.get(), valuesin.get()+size), EqualsRange(aim, 1e-10));
        psolver1.assignvalue(psolver7);
        auto value1 = psolver1.test_getvalue("value");
        CHECK_THAT(vector<Real>(value1.get(), value1.get()+size), EqualsRange(vector<Real>(size, 7.0), 1e-10));
        psolversin.minusvalue(psolver1);
        for(auto k = 0; k < size; ++k){
            aim[k] -= 7.0;
        }
        valuesin.reset(psolversin.test_getvalue("value").release());
        CHECK_THAT(vector<Real>(valuesin.get(), valuesin.get()+size), EqualsRange(aim, 1e-10));
    }
}