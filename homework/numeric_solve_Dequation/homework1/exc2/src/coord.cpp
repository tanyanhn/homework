#include "coord.h"

using namespace solver_ty;


template<>
std::unique_ptr<int[]> coordinate<1>::getnearcoord(int co) const {
    assert(co > 0 && co < n[0]);
    return std::unique_ptr<int[]>(new int[2]{co-1, co+1});
}

template<>
std::unique_ptr<int[]> coordinate<1>::getdiagnolcoord(int) const {
    assert(false && "Shouldn't use getdiagnolcoord in dim = 1.");
    return nullptr;
}

template<>
std::unique_ptr<int[]> coordinate<2>::getnearcoord(int co) const {
    assert(co > n[1] && co < n[0]*(n[1]+1));
    assert(co - std::floor(co/double(n[1]+1))*(n[1]+1) > 0);
    assert(std::floor(-co/double(n[1]+1))*(-n[1]-1) - co > 0);
    return std::unique_ptr<int[]>(new int[4]{co-n[1]-1 ,co-1, co+1, co+n[1]+1});
}


template<>
std::unique_ptr<int[]> coordinate<2>::getdiagnolcoord(int co) const {
    assert(co > n[1] && co < n[0]*(n[1]+1));
    assert(co - std::floor(co/double(n[1]+1))*(n[1]+1) > 0);
    assert(std::floor(-co/double(n[1]+1))*(-n[1]-1) - co > 0);
    return std::unique_ptr<int[]>(new int[4]{co-n[1]-2 ,co-n[1], co+n[1], co+n[1]+2});
}




