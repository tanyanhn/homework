#include "ThreeBody.h"
#include <cmath>
#include "norm.h"
#include "Inverse.h"
#include <iostream>

using namespace std;
using namespace solver_ty;

namespace integrator_ty{

void fdf::f(const Real* u, Real* y1){
    y1[0] = u[3];
    y1[1] = u[4];
    y1[2] = u[5];
    y1[3] = 2 * u[4] + u[0] -
            mu * (u[0] + mu - 1) /
            pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu - 1, 2), 1.5) -
            (1 - mu) * (u[0] + mu) /
            pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu, 2), 1.5);
    y1[4] = - 2 * u[3] + u[1] -
            mu * u[1] /
            pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu - 1, 2), 1.5) -
            (1 - mu) * (u[1]) /
            pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu, 2), 1.5);
    y1[5] = - mu * u[2] /
            pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu - 1, 2), 1.5) -
            (1 - mu) * (u[2]) /
            pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu, 2), 1.5);
}

void fdf::df(const Real* u, Real* output){
    output[0] = 0; output[1] = 0; output[2] = 0;
    output[3] = 1; output[4] = 0; output[5] = 0;
    output[6] = 0; output[7] = 0; output[8] = 0;
    output[9] = 0; output[10] = 1; output[11] = 0;
    output[12] = 0; output[13] = 0; output[14] = 0;
    output[15] = 0; output[16] = 0; output[17] = 1;
    
    output[18] = (mu - 1)/pow(pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2), (3.0/2)) - mu/pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0/2))
            + (3*mu*(2*mu + 2*u[0] - 2)*(mu + u[0] - 1))/(2*pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0/2))) - (3*(2*mu + 2*u[0])*(mu + u[0])*(mu - 1))/(2*pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0/2))) + 1;
    
    output[19] = (mu * (u[0] + mu - 1) * 3 * u[1]) / (pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu - 1, 2), 2.5)) +
            ((1 - mu) * (u[0] + mu) * 3 * u[1]) / (pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu, 2), 2.5));
    
    output[20] = (mu * (u[0] + mu - 1) * 3 * u[2]) / (pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu - 1, 2), 2.5)) +
            ((1 - mu) * (u[0] + mu) * 3 * u[2]) / (pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu, 2), 2.5));

    output[21] = 0;
    output[22] = 2;
    output[23] = 0;

    output[24] = (mu * u[1] * 3 * (u[0] + mu - 1)) / (pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu - 1, 2), 2.5)) + ((1 - mu) * u[1] * 3 * (u[0] + mu)) / (pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu, 2), 2.5));

    output[25] = (mu - 1)/pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0/2)) - mu/pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0/2)) - (3*pow(u[1], 2)*(mu - 1))/pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0/2)) + (3*mu*pow(u[1],2))/pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0/2)) + 1;

    output[26] = (mu * u[1] * 3 * u[2]) / (pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu - 1, 2), 2.5)) +
            ((1 - mu) * u[1] * 3 * u[2]) / (pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu, 2), 2.5));

    output[27] = -2;
    output[28] = 0;
    output[29] = 0;

    output[30] = (mu * u[2] * (u[0] + mu - 1) * 3) / (pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu - 1, 2), 2.5)) +
            ((1 - mu) * u[2] * (u[0] + mu) * 3) / (pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu, 2), 2.5));

    output[31] = (mu * u[2] * u[1] * 3) / (pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu - 1, 2), 2.5)) +
            ((1 - mu) * u[2] * u[1] * 3) / (pow(pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu, 2), 2.5));

    output[32] = (mu - 1)/pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0/2)) - mu/pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0/2)) - (3*pow(u[2], 2)*(mu - 1))/pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0/2)) + (3*mu*pow(u[2], 2))/pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0/2));
    
    output[33] = 0;
    output[34] = 0;
    output[35] = 0;
}

template<int p>
void ThreeBody<p>::f(const Real* u, Real* y1){
    return fop.f(u, y1);
}

template<int p>
void ThreeBody<p>::df(const Real* u, Real* y1){
    return fop.df(u, y1);
}


template<>
void ThreeBody<1>::RKmethod(Real* u, Real ts, Real* out){
    Real y1[size], y2[size], u1[size];
    f(u, y1);
    for(auto i = 0; i < size; ++i){
        u1[i] = u[i] + ts / 2 * y1[i];
    }
    f(u1, y2);
    for(auto i = 0; i < size; ++i){
        out[i] = u[i] + ts * (y2[i]);
    }
}

template<>
void ThreeBody<2>::RKmethod(Real* u, Real ts, Real* out){
    Real y1[size], y2[size], u1[size];
    f(u, y1);
    for(auto i = 0; i < size; ++i){
        u1[i] = u[i] + ts * y1[i];
    }
    f(u1, y2);
    for(auto i = 0; i < size; ++i){
        out[i] = u[i] + ts / 2 * (y1[i] + y2[i]);
    }
}
template<>
void ThreeBody<3>::RKmethod(Real*, Real, Real*){
}

template<>
void ThreeBody<4>::RKmethod(Real* u, Real ts, Real* out){
    Real y1[size], y2[size], y3[size], y4[size], u1[size], u2[size], u3[size];
    f(u, y1);
    for(auto i = 0; i < size; ++i){
        u1[i] = u[i] + ts / 2 * y1[i];
    }
    f(u1, y2);
    for(auto i = 0; i < size; ++i){
        u2[i] = u[i] + ts / 2 * y2[i];
    }
    f(u2, y3);
    for(auto i = 0; i < size; ++i){
        u3[i] = u[i] + ts * y3[i];
    }
    f(u3, y4);
    for(auto i = 0; i < size; ++i){
        out[i] = u[i] + ts / 6 * (y1[i] + 2 * y2[i] + 2 * y3[i] + y4[i]);
    }
}

template<>
void ThreeBody<1>::ABmethod(Real** u, Real ts, Real** y, Real* out){
    // Real y1[size];
    //f(u[0], y1);
    for(auto i = 0; i < size; ++i){
        out[i] = u[0][i] + ts * y[0][i];
    }
}

template<>
void ThreeBody<2>::ABmethod(Real** u, Real ts, Real** y, Real* out){
    // Real y1[size], y2[size];
    // f(u[0], y1);
    // f(u[1], y2);
    for(auto i = 0; i < size; ++i){
        out[i] = u[0][i] + ts * (3.0 / 2 * y[0][i] - 1.0 / 2 * y[1][i]);
    }
}

template<>
void ThreeBody<3>::ABmethod(Real** u, Real ts, Real** y, Real* out){
    // Real y1[size], y2[size], y3[size];
    // f(u[0], y1);
    // f(u[1], y2);
    // f(u[2], y3);
    for(auto i = 0; i < size; ++i){
        out[i] = u[0][i] + ts * (23.0 / 12 * y[0][i] - 16.0 / 12 * y[1][i] + 5.0 / 12 * y[2][i]);
    }
}

template<>
void ThreeBody<4>::ABmethod(Real** u, Real ts, Real** y, Real* out){
    // Real y1[size], y2[size], y3[size], y4[size];
    // f(u[0], y1);
    // f(u[1], y2);
    // f(u[2], y3);
    // f(u[3] ,y4);
    for(auto i = 0; i < size; ++i){
        out[i] = u[0][i] + ts * (55.0 / 24 * y[0][i] - 59.0 / 24 * y[1][i]
                                 + 37.0 / 24 * y[2][i] - 9.0 / 24 * y[3][i]);
    }
}

template<>
void ThreeBody<1>::AMmethod(Real** u, Real ts, Real** y, Real* out){
    for(auto k = 0; k < size; ++k){
        out[k] = u[0][k];
    }
    Real A[size * size], c[size], r[size];
    int NTst = 5;
    Inverse inv(size);
    for(auto i = 0; i < size; ++i){
        c[i] = u[0][i] + ts * (1.0 / 2 * y[0][i]);
    }
    int count = 1;
    while(true){
        fop.f(out, r);
        fop.df(out, A);
        for(auto i = 0; i < size; ++i){
            r[i] = -(ts * 1.0 / 2 * r[i] + c[i] - out[i]);
            for(auto j = 0; j < size; ++j){
                A[i * size + j] = ts * 1.0 / 2.0 * A[i * size + j];
            }
            A[i * size + i] -= 1;
        }
        for(auto i = 0; i < size; ++i){
            for(auto j = 0; j < size; ++j){
                r[i] += A[i * size + j] * out[j];
            }
        }
        inv.Solverequations(A, r);
        for(auto i = 0; i < size; ++i){
            out[i] = r[i];
        }
        if(count++ > NTst){
            break;
        }
    }
}

template<>
void ThreeBody<2>::AMmethod(Real** u, Real ts, Real** y, Real* out){
    for(auto k = 0; k < size; ++k){
        out[k] = u[0][k];
    }
    Real A[size * size], c[size], r[size];
    int NTst = 5;
    Inverse inv(size);
    for(auto i = 0; i < size; ++i){
        c[i] = u[0][i] + ts * (8.0 / 12 * y[0][i] - 1.0 / 12 * y[1][i]);
    }
    int count = 1;
    while(true){
        fop.f(out, r);
        fop.df(out, A);
        for(auto i = 0; i < size; ++i){
            r[i] = -(ts * 5.0 / 12 * r[i] + c[i] - out[i]);
            for(auto j = 0; j < size; ++j){
                A[i * size + j] = ts * 5.0 / 12.0 * A[i * size + j];
            }
            A[i * size + i] -= 1;
        }
        for(auto i = 0; i < size; ++i){
            for(auto j = 0; j < size; ++j){
                r[i] += A[i * size + j] * out[j];
            }
        }
        inv.Solverequations(A, r);
        for(auto i = 0; i < size; ++i){
            out[i] = r[i];
        }
        if(count++ > NTst){
            break;
        }
    }
}

template<>
void ThreeBody<3>::AMmethod(Real** u, Real ts, Real** y, Real* out){
    for(auto k = 0; k < size; ++k){
        out[k] = u[0][k];
    }
    Real A[size * size], c[size], r[size];
    int NTst = 5;
    Inverse inv(size);
    for(auto i = 0; i < size; ++i){
        c[i] = u[0][i] + ts * (19.0 / 24 * y[0][i] - 5.0 / 24 * y[1][i] + 1.0 / 24 * y[2][i]);
    }
    int count = 1;
    while(true){
        fop.f(out, r);
        fop.df(out, A);
        for(auto i = 0; i < size; ++i){
            r[i] = -(ts * 9.0 / 24 * r[i] + c[i] - out[i]);
            for(auto j = 0; j < size; ++j){
                A[i * size + j] = ts * 9.0 / 24.0 * A[i * size + j];
            }
            A[i * size + i] -= 1;
        }
        for(auto i = 0; i < size; ++i){
            for(auto j = 0; j < size; ++j){
                r[i] += A[i * size + j] * out[j];
            }
        }
        inv.Solverequations(A, r);
        for(auto i = 0; i < size; ++i){
            out[i] = r[i];
        }
        if(count++ > NTst){
            break;
        }
    }
}


template<>
void ThreeBody<4>::AMmethod(Real** u, Real ts, Real** y, Real* out){
    for(auto k = 0; k < size; ++k){
        out[k] = u[0][k];
    }
    Real A[size * size], c[size], r[size];
    int NTst = 5;
    Inverse inv(size);
    for(auto i = 0; i < size; ++i){
        c[i] = u[0][i] + ts * (646.0 / 720 * y[0][i] - 264.0 / 720 * y[1][i] + 106.0 / 720 * y[2][i] - 19.0 / 720 * y[3][i]);
    }
    int count = 1;
    while(true){
        fop.f(out, r);
        fop.df(out, A);
        for(auto i = 0; i < size; ++i){
            r[i] = -(ts * 251.0 / 720 * r[i] + c[i] - out[i]);
            for(auto j = 0; j < size; ++j){
                A[i * size + j] = ts * 251.0 / 720.0 * A[i * size + j];
            }
            A[i * size + i] -= 1;
        }
        for(auto i = 0; i < size; ++i){
            for(auto j = 0; j < size; ++j){
                r[i] += A[i * size + j] * out[j];
            }
        }
        inv.Solverequations(A, r);
        for(auto i = 0; i < size; ++i){
            out[i] = r[i];
        }
        if(count++ > NTst){
            break;
        }
    }
}

template<>
void ThreeBody<1>::BDmethod(Real** u, Real ts, Real* out){
    for(auto k = 0; k < size; ++k){
        out[k] = u[0][k];
    }
    Real A[size * size], c[size], r[size];
    int NTst = 5;
    Inverse inv(size);
    for(auto i = 0; i < size; ++i){
        c[i] = u[0][i];
    }
    int count = 1;
    while(true){
        fop.f(out, r);
        fop.df(out, A);
        for(auto i = 0; i < size; ++i){
            r[i] = -(ts * r[i] + c[i] - out[i]);
            for(auto j = 0; j < size; ++j){
                A[i * size + j] = ts * A[i * size + j];
            }
            A[i * size + i] -= 1;
        }
        for(auto i = 0; i < size; ++i){
            for(auto j = 0; j < size; ++j){
                r[i] += A[i * size + j] * out[j];
            }
        }
        inv.Solverequations(A, r);
        for(auto i = 0; i < size; ++i){
            out[i] = r[i];
        }
        if(count++ > NTst){
            break;
        }
    }
}

template<>
void ThreeBody<2>::BDmethod(Real** u, Real ts, Real* out){
    for(auto k = 0; k < size; ++k){
        out[k] = u[0][k];
    }
    Real A[size * size], c[size], r[size];
    int NTst = 5;
    Inverse inv(size);
    for(auto i = 0; i < size; ++i){
        c[i] = 4.0 / 3 * u[0][i] - 1.0 / 3 * u[1][i];
    }
    int count = 1;
    while(true){
        fop.f(out, r);
        fop.df(out, A);
        for(auto i = 0; i < size; ++i){
            r[i] = -(ts * 2.0 / 3 * r[i] + c[i] - out[i]);
            for(auto j = 0; j < size; ++j){
                A[i * size + j] = ts * 2.0 / 3 * A[i * size + j];
            }
            A[i * size + i] -= 1;
        }
        for(auto i = 0; i < size; ++i){
            for(auto j = 0; j < size; ++j){
                r[i] += A[i * size + j] * out[j];
            }
        }
        inv.Solverequations(A, r);
        for(auto i = 0; i < size; ++i){
            out[i] = r[i];
        }
        if(count++ > NTst){
            break;
        }
    }
}

template<>
void ThreeBody<3>::BDmethod(Real** u, Real ts, Real* out){
    for(auto k = 0; k < size; ++k){
        out[k] = u[0][k];
    }
    Real A[size * size], c[size], r[size];
    int NTst = 5;
    Inverse inv(size);
    for(auto i = 0; i < size; ++i){
        c[i] = 18.0 / 11 * u[0][i] - 9.0 / 11 * u[1][i] + 2.0 / 11 * u[2][i];
    }
    int count = 1;
    while(true){
        fop.f(out, r);
        fop.df(out, A);
        for(auto i = 0; i < size; ++i){
            r[i] = -(ts * 6.0 / 11 * r[i] + c[i] - out[i]);
            for(auto j = 0; j < size; ++j){
                A[i * size + j] = ts * 6.0 / 11 * A[i * size + j];
            }
            A[i * size + i] -= 1;
        }
        for(auto i = 0; i < size; ++i){
            for(auto j = 0; j < size; ++j){
                r[i] += A[i * size + j] * out[j];
            }
        }
        inv.Solverequations(A, r);
        for(auto i = 0; i < size; ++i){
            out[i] = r[i];
        }
        if(count++ > NTst){
            break;
        }
    }
}

template<>
void ThreeBody<4>::BDmethod(Real** u, Real ts, Real* out){
    for(auto k = 0; k < size; ++k){
        out[k] = u[0][k];
    }
    Real A[size * size], c[size], r[size];
    int NTst = 5;
    Inverse inv(size);
    for(auto i = 0; i < size; ++i){
        c[i] = 48.0 / 25 * u[0][i] - 36.0 / 25 * u[1][i] + 16.0 / 25 * u[2][i]
                - 3.0 / 25 * u[3][i];
    }
    int count = 1;
    while(true){
        fop.f(out, r);
        fop.df(out, A);
        for(auto i = 0; i < size; ++i){
            r[i] = -(ts * 12.0 / 25 * r[i] + c[i] - out[i]);
            for(auto j = 0; j < size; ++j){
                A[i * size + j] = ts * 12.0 / 25 * A[i * size + j];
            }
            A[i * size + i] -= 1;
        }
        for(auto i = 0; i < size; ++i){
            for(auto j = 0; j < size; ++j){
                r[i] += A[i * size + j] * out[j];
            }
        }
        inv.Solverequations(A, r);
        for(auto i = 0; i < size; ++i){
            out[i] = r[i];
        }
        if(count++ > NTst){
            break;
        }
    }
}

template<int p>
Real* ThreeBody<p>::LMMsolver(Real T){
    int step = floor(T / timeSt - 1);
    Real lastts = T - timeSt * step;
    Real** U = new Real*[p];
    Real** y = new Real*[p];
    Real* output = new Real[size * (step+2)],
            *ally = new Real[size * (step+2)];
    ThreeBody<4> initgenerator(*this);
    for(auto i = 0; i < p; ++i){
        if(i == 0){
            for(auto k = 0; k < size; ++k){
                output[k] = value[k];
            }
            U[p-1] = output;
            fop.f(U[p-1], ally);
            y[p-1] = ally;
        }
        else{
            initgenerator.RKmethod(U[p-i], timeSt, output + i * size);
            U[p-i-1] = output + i * size;
            fop.f(U[p-i-1], ally + i * size);
            y[p-i-1] = ally + i * size;
        }
    }
    int istep = p-1;
    while(true){
        ++istep;
        if(istep >= (step  + 1)){
            if(met == AB){
                ABmethod(U, lastts, y, output + istep * size);
            }
            else if(met == AM){
                AMmethod(U, lastts, y, output + istep * size);
            }
            else if(met == BDF){
                BDmethod(U, lastts, output + istep * size);
            }
            break;
        }
        else{
            if(met == AB){
                ABmethod(U, timeSt, y, output + istep * size);
            }
            else if(met == AM){
                AMmethod(U, timeSt, y, output + istep * size);
            }
            else if(met == BDF){
                BDmethod(U, timeSt, output + istep * size);
            }
            for(auto i = p-1; i >0; --i){
                U[i] = U[i-1];
                y[i] = y[i-1];
            }
            U[0] = output + istep * size;
            fop.f(U[0], ally + istep * size);
            y[0] = ally + istep * size;
        }
    }
    delete[] U;
    delete[] y;
    delete[] ally;
    return output;
}

template<int p>
Real* ThreeBody<p>::solver(Real T){
    if(T < 0){
        T = aimT;
        if(T < 0){
            cout << "AimTime T < 0";
            return nullptr;
        }
    }
    if(timeSt > 1){
        timeSt = T / timeSt;
    }
    if(met != RK){
        return LMMsolver(T);
    }
    int step = floor(T / timeSt - 1);
    Real lastts = T - timeSt * step;
    Real U0[size], U1[size];
    Real* output = new Real[size * (step+2)];
    for(auto k = 0; k < size; ++k){
        output[k] = value[k];
        U0[k] = value[k];
    }
    int istep = 0;
    while(true){
        ++istep;
        if(istep == (step  + 1)){
            RKmethod(U0, lastts, U1);
            for(auto k = 0; k < size; ++k){
                output[istep * size + k] = U1[k];
                U0[k] = U1[k];
            }
            break;
        }
        else{
            RKmethod(U0, timeSt, U1);
            for(auto k = 0; k < size; ++k){
                output[istep * size + k] = U1[k];
                U0[k] = U1[k];
            }
        }
    }
    return output;
}

template class ThreeBody<1>;
template class ThreeBody<2>;
template class ThreeBody<3>;
template class ThreeBody<4>;

}