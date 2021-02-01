#ifndef NORM_TY
#define NORM_TY

#include "real.h"
#include <cmath>

namespace solver_ty{

template<int dim, typename T>
Real norm(T v, int n){
    Real sum = 0;
    for(auto i = 0; i < n; ++i){
        sum += std::pow(fabs(v[i]), dim);
    }
    return std::pow(sum / n, 1.0 / dim);
}

 
}

#endif