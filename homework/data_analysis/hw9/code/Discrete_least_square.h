#ifndef DLS_H
#define DLS_H

#include<vector>
#include<iostream>
#include<lapacke.h>

//using namespace std;

class DLS {
protected:
    std::vector<double> x,y;
    int a;
public:
    DLS(){};
    DLS(const std::vector<double>& xx, const std::vector<double>& yy, const int aa) : x(xx), y(yy), a(aa){}
    //DLS(const DLS& temp) : x(temp.x), y(temp.y) {}
    virtual ~DLS(){}
    void initialize(std::istream& is){
        int aa = 0;
        is >> aa;
        a = aa;
        std::vector<double> xx,yy;
        double tempx,tempy;
        while(is >> tempx >> tempy){
            xx.push_back(tempx);
            yy.push_back(tempy);
        }
        x = xx;
        y = yy;
    };
    virtual std::vector<double> work(double* ) const = 0;
    void Equationsolver(double* A, double* B, const int N) const{
        int ipiv[N];
        LAPACKE_dgesv(LAPACK_COL_MAJOR, N, 1, A, N, ipiv, B, N);
    }
};



#endif