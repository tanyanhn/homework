#ifndef QR_FACTOR_H
#define QR_FACTOR_H

#include"Discrete_least_square.h"


class QRfactor: public DLS{
public:
    QRfactor(){}
    QRfactor(const std::vector<double>& xx, const std::vector<double>& yy, const int aa) : DLS(xx, yy, aa) {}
    ~QRfactor(){}
    std::vector<double> work(double* ) const;
    void Generatematrix(double* G, double* c, const int m,  const int n) const;
    void Maqr(double** A, double** Q, const int m, const int n) const;
};



#endif