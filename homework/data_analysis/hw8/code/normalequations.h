#ifndef NORMALEQAUTION_H
#define NORMALEQAUTION_H




#include"Discrete_least_square.h"




class Normalequation : public DLS {
public:
    Normalequation(){}
    Normalequation(const std::vector<double>& xx, const std::vector<double>& yy, const int aa) : DLS(xx, yy, aa) {}
    ~Normalequation(){}
    std::vector<double> work() const;
    void Generatematrix(double* G, const int N) const;
    void Generateconstant(double* c, const int N) const;
};


#endif