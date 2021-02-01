#include "norm.h"

template<>
Real solver_ty::norm<0, Real*>(Real* v, int n){
    Real sum = 0;
    for(auto i = 0; i < n; ++i){
        if(sum < fabs(v[i])){
            sum = fabs(v[i]);
        }
    }
    return sum;
}