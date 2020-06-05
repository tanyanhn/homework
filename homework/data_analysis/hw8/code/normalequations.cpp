#include<iostream>
#include<algorithm>
#include "normalequations.h"

using namespace std;


vector<double> Normalequation::work() const{
    int n = a;
    double G[(n+1) * (n+1)], c[n+1];
    vector<double> anwser;
    Generatematrix(G, n);
    Generateconstant(c, n);
    Equationsolver(G, c, n+1);
    for_each(c, c + n+1, [&anwser](double i){ anwser.push_back(i);});
    return anwser;
}


void Normalequation::Generatematrix(double* G, const int n) const{
    //vector<double> x = DLS::x;
    for(auto i = 0; i < n+1 ; i++){
        double sum = 0;
         for(auto k = x.cbegin(); k != x.cend(); k++){
            sum += pow(*k, i);
        }
        for(auto j =0; j <= i; ++j)
            G[i + n*j] = sum;
    }
    for(auto i = 1; i < n+1 ; i++){
        double sum = 0;
        for(auto k = x.cbegin(); k != x.cend(); k++){
            sum += pow(*k,n + i);
        }
        for(auto j =0; j <= n-i; ++j)
            G[n + i*(n+1) + j*n] = sum;
    }
    //for(auto i = 0; i < n+1; ++i)
}

void Normalequation::Generateconstant(double* c, const int n) const{
    vector<double> x = DLS::x, y = DLS::y;
    for(auto i = 0; i < n+1; ++i){
        double sum = 0;
        auto j = y.begin();
        for(auto k = x.cbegin(); k!= x.cend(); ++k){
            sum += pow(*k, i) * (*j);
            ++j;
        }
        c[i] = sum;
    }
}