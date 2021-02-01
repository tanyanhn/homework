#include "possionsolver.h"
#include <iostream>
#include "cycle.h"
#include "norm.h"


namespace solver_ty{

template<int dim>
possionsolver<dim>::possionsolver(int nn[dim], Real hh[dim],
                                  std::unique_ptr<Real[]> f,
                                  std::unique_ptr<Real[]> v) :
        equationsolver<dim>(nn, hh, move(f), move(v)) {
    if(dim > 1){
        for(auto i = 0; i < dim-1; ++i){
            if(hh[i] != hh[i+1]){
                std::cout << "Haven't implement situation that hx != hy.";
            }
        }
    }
}

template<>
void possionsolver<1>::updatervalue() {
    if(btype == Diri){
        Real nearv[2];
        rvalue[0] = 0;
        rvalue[n[0]] = 0;
        for(auto i = 1; i < n[0]; ++i){
            int id = i;
            auto co = pcoord->getnearcoord(id);
            nearv[0] = value[co[0]];
            nearv[1] = value[co[1]];
            rvalue[id] = (h[0] * h[0] * fvalue[id] + nearv[0] + nearv[1] - 2 * value[id])
                    / (h[0] * h[0]);
        }
    }
    else{
        assert(false && "Haven't define relax type.");
    }
}

template<>
void possionsolver<2>::updatervalue(){
    if(btype == Diri){
        Real nearv[4];
        for(auto i = 0; i < n[0] + 1; ++i){
            for(auto j = 0; j < n[1] + 1; ++j){
                int id = i * (n[1] + 1) + j;
                if(i == 0 || i == n[0] || j == 0 || j == n[1]){
                    rvalue[id] = 0;
                } else {
                    auto co = pcoord->getnearcoord(id);
                    nearv[0] = value[co[0]];
                    nearv[1] = value[co[1]];
                    nearv[2] = value[co[2]];
                    nearv[3] = value[co[3]];
                    rvalue[id] = ((h[0] * h[0] *fvalue[id] +
                                   (nearv[0] + nearv[3] - 2 * value[id]) +
                                   (nearv[1] + nearv[2] - 2 * value[id])) /
                                  (h[0] * h[0]));
                }
            }
        }
    }
    else{
        assert(false && "Haven't define boundary type function.");
    }
}

template<>
void possionsolver<1>::Fullrestrict(Real* v, Real* result){
    *result = (v[0] + v[1] + 2 * v[2]) / 4;
}

template<>
void possionsolver<2>::Fullrestrict(Real* v, Real* result){
    *result = ((v[0] + v[1] + v[2] + v[3]) + 2 * (v[4] + v[5] + v[6] + v[7]) + 4 * v[8]) / 16;
}

template<>
void possionsolver<1>::Jacobirelax(Real* v, Real f, Real* result){
    *result = (v[0] + v[1] + h[0] * h[0] * (f)) / 2;
}

template<>
void possionsolver<2>::Jacobirelax(Real* v, Real f, Real* result){
    *result = (v[0] + v[1] + v[2] + v[3] + h[0] * h[1] * (f)) / 4;
}

template<int dim>
void possionsolver<dim>::Relax() {
    if(this->rtype == GaussSeidel) {
        Relax(this->value);
    }
    else if(this->rtype == RB_GaussSeidel){
        RB_GaussSeidel_Relax();
    }
    else if(this->rtype == Jacobi || this->rtype == weightJacobi){
        int size = 1;
        for(auto i = 0; i < dim; ++i){
            size *= (this->n[i] + 1);
        }
        std::unique_ptr<Real[]> v(new Real[size]);
        Relax(v);
        this->value.swap(v);
    }
}

template<>
void possionsolver<1>::RB_GaussSeidel_Relax(){
    if(this->btype == Diri){
        Real nearv[2];
        for(auto m = 0; m < 2; ++m){
            for(auto i = 1 + m; i < n[0]; i += 2){
                int id = i;
                auto co = this->pcoord->getnearcoord(id);
                for(auto k = 0; k < 2; ++k){
                    nearv[k] = this->value[co[k]];
                }
                Jacobirelax(nearv, this->fvalue[id], this->value.get()+id);
            }
        }
    }
    else{
        assert(false && "Haven't define boundary type.");
    }
}

template<>
void possionsolver<2>::RB_GaussSeidel_Relax(){
    if(this->btype == Diri){
        Real nearv[4];
        for(auto m = 0; m < 2; ++m){
            for(auto i = 1; i < n[0]; i += 2){
                for(auto j = 1 + m; j < n[1]; j += 2){
                    int id = i*(n[1] + 1) + j;
                    auto co = this->pcoord->getnearcoord(id);
                    for(auto k = 0; k < 4; ++k){
                        nearv[k] = this->value[co[k]];
                    }
                    Jacobirelax(nearv, this->fvalue[id], this->value.get()+id);
                }
            }
            for(auto i = 2; i < n[0]; i += 2){
                for(auto j = 2 - m; j < n[1]; j += 2){
                    int id = i*(n[1] + 1) + j;
                    auto co = this->pcoord->getnearcoord(id);
                    for(auto k = 0; k < 4; ++k){
                        nearv[k] = this->value[co[k]];
                    }
                    Jacobirelax(nearv, this->fvalue[id], this->value.get()+id);
                }
            }
        }
    }
    else{
        assert(false && "Haven't define boundary type.");
    }
}

template<>
void possionsolver<1>::Relax(std::unique_ptr<Real[]>& v){
    if(btype == Diri){
        Real nearv[2];
        for(auto i = 0; i < n[0] + 1; ++i){
            int id = i;
            if(i == 0 || i == n[0]){
                v[id] = this->value[id];
            }
            else {
                auto co = pcoord->getnearcoord(id);
                if(rtype == Jacobi || rtype == GaussSeidel){
                    nearv[0] = this->value[co[0]];
                    nearv[1] = this->value[co[1]];
                    Jacobirelax(nearv, fvalue[id], v.get()+id);
                }
                else if(rtype == weightJacobi){
                    nearv[0] = this->value[co[0]];
                    nearv[1] = this->value[co[1]];
                    Jacobirelax(nearv, fvalue[id], v.get()+id);
                    v[id] = (2 *v[id] + this->value[id]) / 3.0; 
                }
                else{
                    assert(false && "Haven't define relax type.");
                }
            }
        }
    }
    else{
        assert(false && "Haven't define boundary type.");
    }
}

template<>
void possionsolver<2>::Relax(std::unique_ptr<Real[]>& v){
    if(btype == Diri){
        Real nearv[4];
        Real tempv;
        for(auto i = 0; i < n[0] + 1; ++i){
            for(auto j = 0; j < n[1] + 1; ++j){
                int id = i*(n[1] + 1) + j;
                if(i == 0 || i == n[0] || j == 0 || j == n[1]){
                    v[id] = this->value[id];
                }
                else{
                    auto co = pcoord->getnearcoord(id);
                    if(rtype == Jacobi || rtype == GaussSeidel){
                        nearv[0] = this->value[co[0]];
                        nearv[1] = this->value[co[1]];
                        nearv[2] = this->value[co[2]];
                        nearv[3] = this->value[co[3]];
                        Jacobirelax(nearv, fvalue[id], v.get()+id);
                    }
                    else if(rtype == weightJacobi){
                        nearv[0] = this->value[co[0]];
                        nearv[1] = this->value[co[1]];
                        nearv[2] = this->value[co[2]];
                        nearv[3] = this->value[co[3]];
                        Jacobirelax(nearv, fvalue[id], &tempv);
                        v[id] = (2 * tempv + this->value[id]) / 3.0; 
                    }
                    else{
                        assert(false && "Haven't define relax type.");
                    }
                }
            }
        }
    }
    else{
        assert(false && "Haven't define boundary type function.");
    }
}

template<>
std::unique_ptr<equationsolver<1>> possionsolver<1>::Restrict(){
    updatervalue();
    int coarsen[1];
    Real coarseh[1];
    for(auto i = 0; i < 1; ++i){
        coarsen[i] = equationsolver<1>::n[i] / 2;
        coarseh[i] = equationsolver<1>::h[i] * 2;
    }
    int size = 1;
    for(auto i = 0; i < 1; ++i){
        size = size * (coarsen[i] + 1);
    }
    std::unique_ptr<Real[]> initialv = std::unique_ptr<Real[]>(new Real[size]),
            initialf = std::unique_ptr<Real[]>(new Real[size]);
    Real nearv[3];
    if(btype == Diri){
        initialv[0] = value[0];
        initialv[coarsen[0]] = value[n[0]];
        initialf[0] = 0;
        initialf[coarsen[0]] = 0;
        for(auto i = 1; i < coarsen[0]; ++i){
            int coarseid = i,
                    id =  2 * i;
            auto co = pcoord->getnearcoord(id);
            if(restrictop == full){
                // nearv[0] = value[co[0]];
                // nearv[1] = value[co[1]];
                // nearv[2] = value[id];
                // Fullrestrict(nearv, initialv.get()+coarseid);
                initialv[coarseid] = 0;
                nearv[0] = rvalue[co[0]];
                nearv[1] = rvalue[co[1]];
                nearv[2] = rvalue[id];
                Fullrestrict(nearv, initialf.get()+coarseid);
            }
            else if(restrictop == inject){
                initialf[coarseid] = rvalue[id];
            }
            else{
                assert(false && "Haven't define restrict type.");
            }
        }
    } else{
        assert(false && "Haven't define boundary type.");
    }
    
    return  std::unique_ptr<equationsolver<1>>(
        new possionsolver<1>(coarsen, coarseh,std::move(initialf),std::move(initialv)));
}

template<>
std::unique_ptr<equationsolver<2>> possionsolver<2>::Restrict(){
    //possionsolver<2>* possionsolver<2>::Restrict(){
    updatervalue();
    int coarsen[2];
    Real coarseh[2];
    for(auto i = 0; i < 2; ++i){
        coarsen[i] = equationsolver<2>::n[i] / 2;
        coarseh[i] = equationsolver<2>::h[i] * 2;
    }
    int size = 1;
    for(auto i = 0; i < 2; ++i){
        size = size * (coarsen[i] + 1);
    }
    std::unique_ptr<Real[]> initialv = std::unique_ptr<Real[]>(new Real[size]),
            initialf = std::unique_ptr<Real[]>(new Real[size]);
    Real nearv[9];
    if(btype == Diri){
        for(auto k = 0; k < coarsen[1] + 1; ++k){
            initialv[k] = value[2 * k];
            initialv[k + coarsen[0] * (coarsen[1] + 1)] = value[2 * k + n[0] * (n[1] + 1)];
            initialf[k] = 0;
            initialf[k + coarsen[0] * (coarsen[1] + 1)] = 0;
        }
        for(auto k = 0; k < coarsen[0] + 1; ++k){
            initialv[k * (coarsen[1] + 1)] = value[2 * k * (n[1] + 1)];
            initialv[k * (coarsen[1] + 1) + coarsen[1]] = value[2 * k * (n[1] + 1) + n[1]];
            initialf[k * (coarsen[1] + 1)] = 0;
            initialf[k * (coarsen[1] + 1) + coarsen[1]] = 0;
        }
        for(auto i = 1; i < coarsen[0]; ++i){
            for(auto j = 1; j < coarsen[1]; ++j){
                int coarseid = i * (coarsen[1] + 1) + j,
                        id = 2 * i * (n[1] + 1) + 2 * j;
                if(restrictop == full){
                    auto nearco = pcoord->getnearcoord(id),
                            diagco = pcoord->getdiagnolcoord(id);
                    // nearv[0] = value[diagco[0]];
                    // nearv[1] = value[diagco[1]];
                    // nearv[2] = value[diagco[2]];
                    // nearv[3] = value[diagco[3]];
                    // nearv[4] = value[nearco[0]];
                    // nearv[5] = value[nearco[1]];
                    // nearv[6] = value[nearco[2]];
                    // nearv[7] = value[nearco[3]];
                    // nearv[8] = value[id];
                    // Fullrestrict(nearv, initialv.get()+coarseid);
                    initialv[coarseid] = 0;
                    nearv[0] = rvalue[diagco[0]];
                    nearv[1] = rvalue[diagco[1]];
                    nearv[2] = rvalue[diagco[2]];
                    nearv[3] = rvalue[diagco[3]];
                    nearv[4] = rvalue[nearco[0]];
                    nearv[5] = rvalue[nearco[1]];
                    nearv[6] = rvalue[nearco[2]];
                    nearv[7] = rvalue[nearco[3]];
                    nearv[8] = rvalue[id];
                    Fullrestrict(nearv, initialf.get()+coarseid);
                }
                else if(restrictop == inject){
                    initialf[coarseid] = rvalue[id];
                }
                else{
                    assert(false && "Haven't define restrict type.");
                }
            }
        }
    }
    else{
        assert(false && "Haven't define boundary type.");
    }
    return std::unique_ptr<equationsolver<2>>(
        new possionsolver<2>(coarsen, coarseh,move(initialf),move(initialv)));
}

template<int dim>
std::unique_ptr<equationsolver<dim>> possionsolver<dim>::Interpolate(){
    //possionsolver<dim>* possionsolver<dim>::Interpolate(){
    if(this->btype == Diri){
        if(this->interpolatop == linear){
            return move(Linearinterpolate());
        } else{
            assert(false && "Haven't define interpolat type.");
        }
    } else{
        assert(false && "Haven't define boundary type.");
    }
}

template<>
std::unique_ptr<equationsolver<1>> possionsolver<1>::Linearinterpolatet(){
    int finen[1] = {n[0] * 2};
    Real fineh[1] = {h[0] / 2};
    auto initialv = std::unique_ptr<Real[]>(new Real[finen[0] + 1]),
            initialf = std::unique_ptr<Real[]>(nullptr);
    initialv[0] = 0;
    initialv[finen[0]] = 0;
    initialv[finen[0] - 1] = (value[n[0] - 1] + value[n[0]]) / 2;
    for(auto i = 1; i < n[0]; ++i){
        initialv[2 * i - 1] = (value[i-1] + value[i]) / 2;
        initialv[2 * i] = value[i];
    }
    return std::unique_ptr<equationsolver<1>>(
        new possionsolver<1>(finen, fineh,move(initialf),move(initialv)));
}

template<>
std::unique_ptr<equationsolver<1>> possionsolver<1>::Linearinterpolate(){
    int finen[1] = {n[0] * 2};
    Real fineh[1] = {h[0] / 2};
    coordinate<1> finecoord(finen);
    auto initialv = std::unique_ptr<Real[]>(new Real[finen[0] + 1]),
            initialf = std::unique_ptr<Real[]>(nullptr);
    for(auto i = 0; i < finen[0] + 1; ++i){
        initialv[i] = 0;
    }
    initialv[0] = 0;
    initialv[finen[0]] = 0;
    for(auto i = 1; i < n[0]; ++i){
        int id = 2 * i;
        initialv[id] = value[i];
        auto co = finecoord.getnearcoord(id);
        initialv[co[0]] += value[i] / 2;
        initialv[co[1]] += value[i] / 2;
    }
    return std::unique_ptr<equationsolver<1>>(
        new possionsolver<1>(finen, fineh,move(initialf),move(initialv)));
}

template<>
std::unique_ptr<equationsolver<2>> possionsolver<2>::Linearinterpolatet(){
    int finen[2] = {n[0] * 2, n[1] * 2};
    Real fineh[2] = {h[0] / 2, h[1] / 2};
    auto initialv = std::unique_ptr<Real[]>(new Real[(finen[0] + 1) * (finen[1] + 1)]),
            initialf = std::unique_ptr<Real[]>(nullptr);
    for(auto k = 0; k < finen[1] + 1; ++k){
        initialv[k] = 0;
        initialv[k + finen[0] * (finen[1] + 1)] = 0;
    }
    for(auto k = 0; k < finen[0] + 1; ++k){
        initialv[k * (finen[1] + 1)] = 0;
        initialv[k * (finen[1] + 1) + finen[1]] = 0;
    }
    for(auto k = 1; k < n[1]; ++k){
        initialv[finen[1] + k * 2] = (value[k - 1] + value[k] +
                                      value[n[1] + k] + value[n[1] + 1 + k]) / 4;
        initialv[finen[1] + k * 2 + 1] = (value[k] + value[n[1] + 1 + k]) / 2;
    }
    for(auto k = 1; k < n[0]; ++k){
        initialv[(2 * k - 1) * (finen[1] + 1) + finen[1] - 1] =
                (value[k * (n[1] + 1) - 2] + value[k * (n[1] + 1) - 1] +
                 value[(k + 1) * (n[1] + 1) - 2] + value[(k + 1) * (n[1] + 1) - 1]) / 4;
        initialv[(2 * k) * (finen[1] + 1) + finen[1] - 1] =
                (value[(k + 1) * (n[1] + 1) - 2] + value[(k + 1) * (n[1] + 1) - 1]) / 2;
    }
    initialv[finen[0] * (finen[1] + 1) - 2] =
            (value[n[0] * (n[1] + 1) - 2] + value[n[0] * (n[1] + 1) - 1] +
             value[(n[0] + 1) * (n[1] + 1) - 2] + value[(n[0] + 1) * (n[1] + 1) - 1]) / 4;
    for(auto i = 1; i < n[0]; ++i){
        for(auto j = 1; j < n[1]; ++j){
            initialv[2 * i * (finen[1] + 1) + j * 2 - 1] =
                    (value[i * (n[1] + 1) + j - 1] +
                     value[i * (n[1] + 1) + j]) / 2;
            initialv[2 * i * (finen[1] + 1) + j * 2] =
                    (value[i * (n[1] + 1) + j]);
            initialv[(2 * i + 1) * (finen[1] + 1) + j * 2 - 1] =
                    (value[i * (n[1] + 1) + j - 1] +
                     value[i * (n[1] + 1) + j] +
                     value[(i + 1) * (n[1] + 1) + j - 1] +
                     value[(i + 1) * (n[1] + 1) + j]) / 4;
            initialv[(2 * i + 1) * (finen[1] + 1) + j * 2] =
                    (value[i * (n[1] + 1) + j] +
                     value[(i + 1) * (n[1] + 1) + j]) / 2;
        }
    }
    return std::unique_ptr<equationsolver<2>>(
        new possionsolver<2>(finen, fineh,move(initialf),move(initialv)));
}

template<>
std::unique_ptr<equationsolver<2>> possionsolver<2>::Linearinterpolate(){
    int finen[2] = {n[0] * 2, n[1] * 2};
    Real fineh[2] = {h[0] / 2, h[1] / 2};
    int finesize = (finen[0] + 1) * (finen[1] + 1);
    auto initialv = std::unique_ptr<Real[]>(new Real[finesize]),
            initialf = std::unique_ptr<Real[]>(nullptr);
    coordinate<2> finecoord(finen);
    for(auto k = 0; k < finesize; ++k){
        initialv[k] = 0;
    }
    for(auto i = 1; i < n[0]; ++i){
        for(auto j = 1; j < n[1]; ++j){
            int fineid = (2 * i) * (finen[1] + 1) + 2 * j,
                    coarid = i * (n[1] + 1) + j;
            initialv[fineid] = value[coarid];
            auto nearco = finecoord.getnearcoord(fineid),
                    diagco = finecoord.getdiagnolcoord(fineid);
            for(auto k = 0; k < 4; ++k){
                initialv[nearco[k]] += value[coarid] / 2;
                initialv[diagco[k]] += value[coarid] / 4;
            }
        }
    }
    return std::unique_ptr<equationsolver<2>>(
        new possionsolver<2>(finen, fineh,move(initialf),move(initialv)));
}


template<int dim>
void possionsolver<dim>::solver(int* v, Real st){
    if(this->stopcondi == maxstep){
        int step = round(st);
        for(auto istep = 0; istep < step; ++istep){
            if(this->ctype == V_cycle){
                Vcycle<dim> Vc;
                Vc(*this, v);
            }
            else if(this->ctype == full_multi_cycle){
                FMGcycle<dim> FMGc;
                if(this->is_homogeneous()){
                    FMGc(*(this), v); 
                }
                else{
                    Vcycle<dim> Vc;
                    Vc(*this, v+1);
                    auto csolver = this->Restrict();
                    csolver->setvaluezero();
                    FMGc(*(csolver.get()), v);
                    auto fsolver = csolver->Interpolate();
                    this->plusvalue(*(fsolver.get()));
                    for(auto k = 0; k < v[2]; ++k){
                        this->Relax();     
                    }
                }
            }
        }
    }
    else if(this->stopcondi == relate_error){
        int size = 1;
        for(auto i = 0; i < 1; ++i){
            size = size * (this->n[i] + 1);
        }
        this->updatervalue();
        auto rnorm = norm<0>(this->rvalue.get(), size);
        int step = 0;
        bool stop = false;
        while(true){
            auto oldvalue = this->getresult();
            if(this->ctype == V_cycle){
                Vcycle<dim> Vc;
                Vc(*this, v);
            }
            else if(this->ctype == full_multi_cycle){
                FMGcycle<dim> FMGc;
                if(step == 0){
                    FMGc((*this), v);
                }
                else{
                    std::unique_ptr<Real[]> tempfvalue(new Real[size]);
                    tempfvalue.swap(this->rvalue);
                    possionsolver<dim> temp(this->n, this->h, nullptr, std::move(tempfvalue));
                    FMGc((temp), v);
                    this->plusvalue(temp);
                }
            }
            else {
                assert(false && "Haven't define stop type with given cycle type.");
            }
            for(auto k = 0; k < size; ++k){
                this->updatervalue();
                auto newrnorm = norm<0>(this->rvalue.get(), size);
                if(newrnorm / rnorm < st){
                    stop = true;
                }
            }
            if(stop || ++step > 100){
                break;
            }
        }
    }
}


template class possionsolver<1>;
template class possionsolver<2>;

}
