#ifndef INVERSE
#define INVERSE
#include "real.h"
#include <lapack.h>
#include <lapacke.h>

namespace integrator_ty{


class Inverse  {
  public:
    int dim;
    explicit Inverse(int i) : dim(i){}
    void operator()(Real*);
    void Solverequations(Real* A, Real* r);
    ~Inverse(){}
};



}

#endif