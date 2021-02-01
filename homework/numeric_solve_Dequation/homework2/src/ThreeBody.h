#ifndef THREEBODY
#define THREEBODY


#include "TimeIntegrator.h"
#include <vector>

namespace integrator_ty{

class fdf {
  public:
    Real mu;
    fdf(){mu = 0;}
    explicit fdf(Real m) : mu(m) {};
    virtual void f(const Real* u, Real* result);
    virtual void df(const Real* u, Real* result);
    ~fdf(){}
};

template<int p>
class ThreeBody : public TimeIntegrator {
  public:
    Real mu;
    fdf fop;
    ThreeBody(){}
    explicit ThreeBody(Real m, Real ts, Real* v, method me, Real TT = -1)
            : mu(m), TimeIntegrator(ts, v, 6, me, TT), fop(mu)
    {}
    template<int op>
    ThreeBody(const ThreeBody<op>& T) {
        *this = T;
    }
    template<int op>
    ThreeBody& operator=(const ThreeBody<op>& T){
        TimeIntegrator::operator=(T);
        mu = T.mu;
        fop.mu = mu;
        return *this;
    }
    Real* solver(Real T = -1) override ;
    Real* LMMsolver(Real T);
    void RKmethod(Real* u, Real ts, Real* result);
    void ABmethod(Real** u, Real ts, Real** y, Real* result);
    void AMmethod(Real** u, Real ts, Real** y, Real* result);
    void BDmethod(Real** u, Real ts, Real* result);
    void f(const Real* U, Real* result);
    void df(const Real* U, Real* result);
    ~ThreeBody(){}
};

}


#endif
