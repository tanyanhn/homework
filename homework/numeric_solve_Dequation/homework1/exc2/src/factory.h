#ifndef FACTORY_TY
#define FACTORY_TY

#include "possionsolver.h"
#include <memory.h>
#include <iostream>
#include <fstream>

namespace solver_ty{

template<int dim>
equationsolver<dim>* solverFactory(std::string file){
    std::fstream is(file, std::ios_base::in);
    equationsolver<dim>* psolver;
    std::string t;
    Real re;
    int in;
    int size = 1;
    int n[dim];
    Real h[dim];
    std::string type;
    while(is >> t){
        if(t == "solvertype:"){
            is >> type;
            if(type == "possionsolver"){
                for(auto i = 0; i < dim; ++i){
                    is >> in;
                    n[i] = in;
                    size *= (in+1);
                    is >> re;
                    h[i] = re;
                }
            }
        }
        else if(t == "solver_argument:"){
            while(true){
                char c;
                is.get(c);
                if(c == '%'){
                    char comment[256];
                    is.getline(comment, 256);
                }
                else{
                    is.unget();
                    is >> t;
                    if(t != "end"){
                        is >> in;
                        psolver->chooseoption(t, in);
                    }
                    else{
                        break;
                    }
                }
            }
        }
        else if(t == "initial:"){
            if(type == "possionsolver"){
                std::fstream is1, is2;
                is >> t;
                is1.open(t,std::ios_base::in);
                is >> t;
                is2.open(t,std::ios_base::in);
                Real* v(new Real[size]), * f(new Real[size]);
                for(auto k = 0; k < size; ++k){
                    is1 >> re;
                    v[k] = re;
                    is2 >> re;
                    f[k] = re;
                }
                psolver = new possionsolver<dim>(n, h,
                                                 std::unique_ptr<Real[]>(f),
                                                 std::unique_ptr<Real[]>(v));
                is1.close(); is2.close();
            }   
        }
    }
    is.close();
    return psolver;
}

}

#endif