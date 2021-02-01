#ifndef STEPOPERATOR
#define STEPOPERATOR

#include "real.h"

namespace integrator_ty{
class Thrbodyoperator  {
  public:
    Thrbodyoperator(){}
    void RKmethod(Real* u, Real ts, Real* result);
    void ABmethod(Real** u, Real ts, Real* result);
    void f(Real* U, Real* result);
    ~Thrbodyoperator(){}
};
}

#endif