#include "TimeIntegrator.h"

using namespace std;

namespace integrator_ty{

Real* TimeIntegrator::solver(Real T){
    return nullptr;
}



template<>
Real TimeIntegrator::setvalue<Real>(string, Real);

template<>
Real* TimeIntegrator::setvalue<Real*>(string, Real*);

template<>
method TimeIntegrator::setvalue<method>(string, method);


}