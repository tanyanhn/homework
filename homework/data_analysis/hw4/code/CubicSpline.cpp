#include"CubicSpline.h"
#include"gsl_function.h"
#include<lapacke.h>
#include<algorithm>
#include<iterator>
using namespace std;


template<class T>
CubicSpline<T>::CubicSpline(T ff, void* params){
    f = f_factory(ff, params);
}

template<class T>
vector<my_f_params> CubicSpline<T>::operator()(const int N) const {
    double X[N], B[N-2], A[(N-2)*(N-2)];
    for_each(A[0], A[(N-2)*(N-2)-1], [](double& d){ d = 0;});
    for(auto i = 0; i < N; i++){
        X[i] = CubicSpline<T>::a + (CubicSpline<T>::b - CubicSpline<T>::a) * i / (N-1);
    }
    vector<pair<double, double>> interpolateValue;
    interpolateValue =  Equationgenerator(A, B, X, N);
    Equationsolver(A, B, N-2);
    for(auto i = 0; i < N-2; i++){
        interpolateValue[i+1].second = B[i];
    }
    vector<my_f_params> anwser;
    for(auto i = 0; i != N-1; i++){
        vector<double> x, f, fp;
        x.push_back(X[i]);
        x.push_back(X[i+1]);
        f.push_back(interpolateValue[i].first);
        f.push_back(interpolateValue[i+1].first);
        fp.push_back(interpolateValue[i].second);
        fp.push_back(interpolateValue[i+1].second);
        my_f_params params = CubicSpline<T>::Hermiteintepolation(x, f, fp);
        anwser.push_back(params);
    }
    return anwser;
}

template<class T>
my_f_params CubicSpline<T>::Hermiteintepolation(vector<double> x, vector<double> f, vector<double> fp) const {
    vector<vector<double>> a;
    vector<double> params;
    double K = (f[0] - f[1]) / (x[0] - x[1]);
    params.push_back((fp[0] + fp[1] - 2 * K) /
                     ((x[1] - x[0]) * (x[1] - x[0])));
    params.push_back((3 * K - 2 * fp[0] - fp[1]) /
                     (x[1] - x[0]));
    params.push_back(fp[0]);
    params.push_back(f[0]);
    a.push_back(params);
    my_f_params anwser(a);
    anwser.getmid() = x[0];
    return anwser;
}


template<class T>
vector<pair<double, double>> CompleteCubic<T>::Equationgenerator(double* A, double* B, double* X, const int N) const {
    double mu[N-2], lam[N-2];
    for(auto i = 0; i < N-2; i++){
        mu[i] = (X[i] - X[i-1]) / (X[i+1] - X[i-1]);
        lam[i] = (X[i+1] - X[i]) / (X[i+1] - X[i-1]);
    }
    A[0] = 2;
    A[1] = lam[1];
    A[(N-2)*(N-2) - 1] = 2;
    A[(N-2)*(N-2) - 2] = mu[N-4];
    for(auto i = 1; i < N-3; i++){
        int k = i * (N-2) + i;
        A[k] = 2;
        A[k-1] = mu[i-1];
        A[k+1] = lam[i+1];
    }
    double fprime[N];
    double f[N];
    gsl_deriv_forward(&CubicSpline<T>::f, CubicSpline<T>::a, 0.1, &fprime[0], 0.00001);
    gsl_deriv_backward(&CubicSpline<T>::f, CubicSpline<T>::b, 0.1, &fprime[N-1], 0.00001);
    for(auto i = 0; i < N; i++){
        f[i] = GSL_FN_EVAL(&CubicSpline<T>::f , X[i]);
    }
    for(auto i = 0; i < N-2; i++){
        B[i] = 3 * mu[i] * (f[i] - f[i+1]) / (X[i] - X[i+1]) +
            3 * lam[i] * (f[i-1] - f[i]) / ( X[i-1] - X[i]);
    }
    B[0] = B[0] - lam[0] * fprime[0];
    B[N-3] = B[N-3] - mu[N-3] * fprime[N-1];
    vector<pair<double, double>> anwser;
    anwser.push_back(make_pair(f[0], fprime[0]));
    for(auto i = 1; i < N-1; i++){
        anwser.push_back(make_pair(f[i], 0));
    }
    anwser.push_back(make_pair(f[N-1], fprime[N-1]));
    return anwser;
}

template<class T>
void CompleteCubic<T>::Equationsolver(double* A, double* B, const int N) const{
    int ipiv[N];
    LAPACKE_dgesv(LAPACK_COL_MAJOR, N, 1, A, N, ipiv, B, N);
    cout << "The solution is : \n";
    copy(B[0], B[N-1],
             ostream_iterator<double>(cout, " "));
    cout << endl;
}