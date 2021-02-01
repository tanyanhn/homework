#ifndef CYCLE_TY
#define CYCLE_TY

#include "equationsolver.h"

namespace solver_ty{

template<int dim>
class cycle {
  public:
    cycle(){}
    cycle(const cycle&){}
    cycle& operator=(const cycle&) { return *this;}
    virtual void operator()(equationsolver<dim>&,const int * const){
        assert(false && "cycle::operator() need override.");
    };
    virtual ~cycle() = 0;
};

template<int dim>
cycle<dim>::~cycle(){}

template<int dim>
class Vcycle : public cycle<dim> {
  public:
    Vcycle(){}
    Vcycle(const Vcycle&){}
    Vcycle& operator=(const Vcycle&){ return *this;}
    ~Vcycle(){}
    void operator()(equationsolver<dim>&, const int * const) override;
};

template<int dim>
class FMGcycle : public cycle<dim> {
  public:
    FMGcycle(){}
    FMGcycle(const FMGcycle&){}
    FMGcycle& operator=(const FMGcycle&){ return *this;}
    ~FMGcycle(){}
    void operator()(equationsolver<dim>&, const int * const) override;
};
}

#endif