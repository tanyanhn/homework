#ifndef IS_ODD_TY
#define IS_ODD_TY

namespace solver_ty{

bool is_odd(int a){
    return a - a / 2 * 2;
}
}

#endif