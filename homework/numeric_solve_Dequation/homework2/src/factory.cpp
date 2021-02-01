#include "factory.h"
#include <fstream>
#include <iostream>
#include <memory>

using namespace std;

namespace integrator_ty{

TimeIntegrator* Factory::operator()(string file , int p){
    ifstream ifs(file);
    Real mu, timestep, aimT;
    char c, comment[1024];
    string str, type;
    int met;
    unique_ptr<Real[]> v;
    c = ifs.get();
    while(c == '%'){
        ifs.unget();
        ifs.getline(comment, 1024);
        c = ifs.get();
    }
    ifs.unget();
    while(ifs >> str){
        if(str == "Integratortype:"){
            ifs >> type;
            if(type == "ThreeBody"){
                ifs >> mu;
                ifs >> timestep;
                ifs >> met;
            }
        }
        else if(str == "Initial:"){
            if(type == "ThreeBody"){
                v.reset(new Real[6]);
                Real temp;
                for(auto i = 0; i < 6; ++i){
                    ifs >> temp;
                    v[i] = temp;
                }
            }
        }
        else if(str == "AimTime:"){
            ifs >> aimT;
        }
        else if(str == "Orderp:"){
            if(p == 0){
                ifs >> p;
            }
            else{
                ifs.get();
            }
        }
    }
    if(type == "ThreeBody"){
        if(p == 1){
            return new ThreeBody<1>(mu, timestep, v.get(), method(met), aimT);
        }
        if(p == 2){
            return new ThreeBody<2>(mu, timestep, v.get(), method(met), aimT);
        }
        if(p == 3){
            return new ThreeBody<3>(mu, timestep, v.get(), method(met), aimT);
        }
        if(p == 4){
            return new ThreeBody<4>(mu, timestep, v.get(), method(met), aimT);
        }
    }
    cout << "Input file haven't choose p";
    return nullptr;
}
}