#ifndef TIMEINTEGRATOR
#define TIMEINTEGRATOR

#include "real.h"
#include <vector>
#include <string>

namespace integrator_ty{

enum method{AB = 0, AM = 1, BDF = 2, RK = 3};

class TimeIntegrator  {
  protected:
    Real timeSt;
    Real* value = nullptr;
    int size;
    method met = RK;
    Real aimT;
  public:
    TimeIntegrator(){};
    TimeIntegrator(Real st, Real* v, int sz, method m = RK, Real TT = -1)
            : timeSt(st), size(sz), met(m), aimT(TT){
        value = new Real[size];
        for(auto i = 0; i < size; ++i){
            value[i] = v[i];
        }
    };
    TimeIntegrator(const TimeIntegrator& T)
            : timeSt(T.timeSt), size(T.size), met(T.met), aimT(T.aimT){
        value = new Real[size];
        for(auto i = 0; i < size; ++i){
            value[i] = T.value[i];
        }
    };
    TimeIntegrator& operator=(const TimeIntegrator& T){
        timeSt = T.timeSt; size = T.size; met = T.met; aimT = T.aimT;
        if(value != nullptr)
        delete[] value;
        value = new Real[size];
        for(auto i = 0; i < size; ++i){
            value[i] = T.value[i];
        }
        return *this;
    };
    template<typename T>
    T setvalue(std::string, T);
    virtual Real* solver(Real T = -1);
    virtual ~TimeIntegrator(){ delete[] value;}
};

template<typename T>
T TimeIntegrator::setvalue(std::string name, T v){
    if(name == "timeSt"){
        auto temp = timeSt;
        timeSt = v;
        return temp;
    }
    if(name == "size"){
        auto temp = size;
        size = int(v);
        return temp;
    }
    if(name == "met"){
        auto temp = met;
        met = method(v);
        return temp;
    }
    if(name == "aimT"){
        auto temp = aimT;
        aimT = v;
        return temp;
    }
    return 0;
}



}

#endif