#include "cycle.h"

using namespace std;

namespace solver_ty{

template<int dim>
//equationsolver<dim>
void Vcycle<dim>::operator()(equationsolver<dim>& solver, const int* const v) {
    for(auto k = 0; k < v[0]; ++k){
        solver.Relax();
    }
    if(solver.is_coarsest()){
        solver.setvaluezero();
    }
    else{
        std::unique_ptr<equationsolver<dim>> csolver(solver.Restrict());
        Vcycle<dim> Vc;
        csolver->setvaluezero();
        Vc(*(csolver.get()), v);
        std::unique_ptr<equationsolver<dim>> fsolver = csolver->Interpolate();
        solver.plusvalue(*(fsolver.get()));
    }
    for(auto k = 0; k < v[1]; ++k){
        solver.Relax();     
    }
}

template<int dim>
void FMGcycle<dim>::operator()(equationsolver<dim>& solver, const int* const v){
    Vcycle<dim> Vc;
    if(solver.is_coarsest()){
        solver.setvaluezero(false);
    }
    else{
        std::unique_ptr<equationsolver<dim>> csolver(solver.Restrict());
        FMGcycle<dim> FMGc;
        FMGc(*(csolver.get()), v);
        std::unique_ptr<equationsolver<dim>> fsolver(csolver->Interpolate());
        solver.assignvalue(*(fsolver.get()), false);
    }
    for(auto k = 0; k < v[0]; ++k){
        Vc(solver, v+1);
    }
}

template class Vcycle<1>;
template class Vcycle<2>;
template class FMGcycle<1>;
template class FMGcycle<2>;

}