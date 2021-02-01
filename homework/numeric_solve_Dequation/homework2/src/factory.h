#ifndef FACTORY
#define FACTORY

#include "ThreeBody.h"
#include <string>

namespace integrator_ty{

class Factory  {
  public:
    TimeIntegrator* operator()(std::string, int p = 0);
    ~Factory(){}
};
}
#endif