#include "equationsolver.h"
#include <assert.h>
#include <iostream>

namespace solver_ty{

template<int dim>
equationsolver<dim>::equationsolver(int nn[dim], Real hh[dim],
                                    std::unique_ptr<Real[]> f,
                                    std::unique_ptr<Real[]> initialvalue)
        : pcoord(new coordinate<dim>(nn)) {
    int size = 1;
    for(auto i = 0; i < dim; ++i){
        size = size * (nn[i] + 1);
    }
    rvalue.reset(new Real[size]);
    if(f != nullptr){
        fvalue.swap(f);
    }
    else{
        fvalue.reset(new Real[size]);
        for(auto i = 0; i < size; ++i){
            fvalue[i] = 0;
        }
    }
    for(auto i = 0; i < dim; ++i){
        n[i] = nn[i];
        h[i] = hh[i];
    }
    if(initialvalue != nullptr){
        value.swap(initialvalue);
    }
    else{
        value.reset(new Real[size]);
        for(auto i = 0; i < size; ++i){
            value[i] = 0;
        }
    }
}

template<int dim>
equationsolver<dim>::equationsolver(const equationsolver& solver){
    *this = solver;
}

template<int dim>
equationsolver<dim>& equationsolver<dim>::operator=(const equationsolver& solver){
    btype = solver.btype;
    restrictop = solver.restrictop;
    interpolatop = solver.interpolatop;
    ctype = solver.ctype;
    stopcondi = solver.stopcondi;
    rtype = solver.rtype;
    int size = 1;
    for(auto i = 0; i < dim; ++i){
        size *= (solver.n[i] + 1);
        n[i] = solver.n[i];
        h[i] = solver.h[i];
    }
    pcoord.reset(new coordinate<dim>(n));
    std::unique_ptr<Real[]> tempfvalue(new Real[size]),
            tempvalue(new Real[size]),
            temprvalue(new Real[size]),
            tempdfvalue(new Real[size]);
    for(auto k = 0; k < size; ++k){
        tempfvalue[k] = solver.fvalue[k];
        tempvalue[k] = solver.value[k];
        temprvalue[k] = solver.rvalue[k];
        if(solver.dfvalue != nullptr){
            tempdfvalue[k] = solver.dfvalue[k];
        }
    }
    fvalue.swap(tempfvalue);
    value.swap(tempvalue);
    rvalue.swap(temprvalue);
    dfvalue.swap(tempdfvalue);
    return *this;
}

template<int dim>
bool equationsolver<dim>::chooseoption(std::string s, int type){
    if(s == "boundary_type"){
        if(type > 1){
            return false;
        }
        btype = boundary_type(type);
    }
    else if(s == "restrict_operator"){
        if(type > 1){
            return false;
        }
        restrictop = restrict_operator(type);
    }
    else if(s == "interpolat_operator"){
        if(type > 1){
            return false;
        }
        interpolatop = interpolat_operator(type);
    }
    else if(s == "cycle_type"){
        if(type > 1){
            return false;
        }
        ctype = cycle_type(type);
    }
    else if(s == "stop_criteria"){
        if(type > 1){
            return false;
        }
        stopcondi = stop_criteria(type);
    }
    else if(s == "relax_type"){
        if(type > 3){
            return false;
        }
        rtype = relax_type(type);
    }
    else {
        return false;
    }
    return true;
}

template<int dim>
std::unique_ptr<Real[]> equationsolver<dim>::test_getvalue(std::string s) const{
    int size = 1;
    for(auto i = 0 ; i < dim; ++i){
        size = size * (n[i] + 1);
    }
    if(s == "n"){
        std::unique_ptr<Real[]> output(new Real[dim]);
        for(auto i = 0; i < dim; ++i){
            output[i] = n[i];
        }
        return output;
    }
    else if(s == "h"){
        std::unique_ptr<Real[]> output(new Real[dim]);
        for(auto i = 0; i < dim; ++i){
            output[i] = h[i];
        }
        return output;
    }
    else if(s == "fvalue"){
        std::unique_ptr<Real[]> output(new Real[size]);
        for(auto i = 0; i < size; ++i){
            output[i] = fvalue[i];
        }
        return output;
    }
    else if(s == "value"){
        std::unique_ptr<Real[]> output(new Real[size]);
        for(auto i = 0; i < size; ++i){
            output[i] = value[i];
        }
        return output;
    }
    else if(s == "rvalue"){
        std::unique_ptr<Real[]> output(new Real[size]);
        for(auto i = 0; i < size; ++i){
            output[i] = rvalue[i];
        }
        return output;
    }
    else if(s == "dfvalue"){
        std::unique_ptr<Real[]> output(new Real[size]);
        for(auto i = 0; i < size; ++i){
            output[i] = dfvalue[i];
        }
        return output;
    }
    else{
        return nullptr;
    }
}
template<int dim>
bool equationsolver<dim>::is_coarsest() const{
    bool output = false;
    for(auto i = 0; i < dim; ++i){
        output = (output || (n[i] <= 4));
    }
    return output;
}

template<int dim>
void equationsolver<dim>::plusvalue(const equationsolver& addv, bool boundary){
    int size = 1;
    for(auto k = 0; k < dim; ++k){
        size *= (n[k] + 1);
    }
    for(auto i = 0; i < size; ++i){
        if(boundary){
            value[i] += addv.value[i];
        }
        else{
            if(!is_boundary(i)){
                value[i] += addv.value[i];
            }
        }
    }
}

template<int dim>
void equationsolver<dim>::minusvalue(const equationsolver& addv, bool boundary){
    int size = 1;
    for(auto k = 0; k < dim; ++k){
        size *= (n[k] + 1);
    }
    for(auto i = 0; i < size; ++i){
        if(boundary){
            value[i] -= addv.value[i];
        }
        else{
            if(!is_boundary(i)){
                value[i] -= addv.value[i];
            }
        }
    }
}

template<int dim>
void equationsolver<dim>::assignvalue(const equationsolver& addv, bool boundary){
    int size = 1;
    for(auto k = 0; k < dim; ++k){
        size *= (n[k] + 1);
    }
    for(auto i = 0; i < size; ++i){
        if(boundary){
            value[i] = addv.value[i];
        }
        else{
            if(!is_boundary(i)){
                value[i] = addv.value[i];
            }
        }
    }
}

template<int dim>
void equationsolver<dim>::solver(int*, Real){
    assert(false && "solver need overrider!");
}

template<int dim>
void equationsolver<dim>::updatervalue(){
    assert(false && "updatarvalue need overrider!");
}

template<int dim>
void equationsolver<dim>::Relax(){
    assert(false && "Relax need overrider!");
}

template<int dim>
std::unique_ptr<equationsolver<dim>> equationsolver<dim>::Restrict(){
    assert(false && "Restrict need overrider!");
}

template<int dim>
std::unique_ptr<equationsolver<dim>> equationsolver<dim>::Interpolate(){
    assert(false && "Interpolate need overrider!");
}

template<int dim>
void equationsolver<dim>::copyflag(equationsolver<dim>& e){
    e.btype = btype;
    e.restrictop = restrictop;
    e.interpolatop = interpolatop;
    e.ctype = ctype;
    e.stopcondi = stopcondi;
    e.rtype = rtype;
}

template<int dim>
std::pair<std::unique_ptr<Real[]>, std::unique_ptr<Real[]>> equationsolver<dim>::getresult(){
    updatervalue();
    auto v = test_getvalue("value"),
            rv = test_getvalue("rvalue");
    return std::make_pair(std::move(v), std::move(rv));
}

template<>
bool equationsolver<1>::is_boundary(int id) const{
    return id == 0 || id == n[0];
}
template<>
bool equationsolver<2>::is_boundary(int id) const{
    if(id <= n[1] || id >= n[0] * (n[1] + 1)){
        return true;
    }
    int a = id - round(double(id) / double(n[1] + 1)) * (n[1] + 1);
    return a == 0 || a == -1;
}

template<int dim>
bool equationsolver<dim>::is_boundary(int id) const{
    Real facesize = 1;
    int newn[dim-1];
    Real newh[dim-1];
    for(auto i = 1; i < dim; ++i){
        facesize *= (this->n[i] + 1);
        newn[i-1] = n[i];
        newh[i-1] = h[i];
    }
    int i0 = std::floor(double(id) / facesize);
    if(i0 == 0 || i0 == n[0]){
        return true;
    }
    int newid = id - i0 * facesize;
    equationsolver<dim - 1> face(newn, newh, nullptr, nullptr);
    return face.is_boundary(newid);
}

template<int dim>
void equationsolver<dim>::setvaluezero(bool boundary){
    int size = 1;
    for(auto i = 0 ; i < dim; ++i){
        size = size * (n[i] + 1);
    }
    for(auto k = 0; k < size; ++k){
        if(boundary){
            value[k] = 0;
        }
        else{
            if(!is_boundary(k)){
                value[k] = 0;
            }
        }
    }
}

template<int dim>
bool equationsolver<dim>::is_homogeneous() const{
    int size = 1;
    for(auto i = 0 ; i < dim; ++i){
        size = size * (n[i] + 1);
    }
    for(auto k = 0; k < size; ++k){
        if(is_boundary(k)){
            if(value[k] != 0){
                return false;
            }
        }
    }
    return true;
}

template class equationsolver<1>;
template class equationsolver<2>;
}