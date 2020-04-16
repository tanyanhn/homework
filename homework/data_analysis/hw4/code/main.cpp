#include"CubicSpline.h"
#include"gsl_function.h"
#include<fstream>
using namespace std;

//template<class T>
//double my_f(double x, void* p);
double computeMaxError(gsl_function f1, gsl_function f2, const double min, const double max, double h){
    double anwser = abs(GSL_FN_EVAL(&f1, min) - GSL_FN_EVAL(&f2, min));
    for( auto i = min + h; i < max; i+=h){
        double temp = abs(GSL_FN_EVAL(&f1, i) - GSL_FN_EVAL(&f2, i));
        if(anwser < temp){
            anwser = temp;
        }
    }
    return anwser;
}

const double h = 1e-5;

int main(int argc, char *argv[])
{
    //ofstream os(argv[1]);
    vector<double> deno;
    deno.push_back(25);
    deno.push_back(0);
    deno.push_back(1);
    vector<vector<double>> params;
    params.push_back(deno);
    my_f_params fparams(params);
    fparams.getmid() = 0;
    double a = -1, b = 1;
    gsl_function g = f_factory<Fractional>(&my_f<Fractional>, &fparams);
    CompleteCubic<gslfunction> F(g);
    //double f, fprime, error = 1e-5;
    //f = GSL_FN_EVAL(&g, a);
    //gsl_deriv_forward(&(g), a, 0.1, &fprime, &error);
    // CompleteCubic<gslfunction> F(&my_f<Fractional>, &fparams);
    F.setInterval(a, b);
    int K[5] = {6, 11, 21, 41, 81};
    vector<vector<my_f_params>> anwser;
    for(auto i = 0; i < 5; i ++){
        vector<my_f_params> p = F.CubicSpline::operator()(K[i]);
        anwser.push_back(p);
    }
    int m = 1;
    for(auto i = anwser.begin(); i != anwser.end(); i++){
        string s("output/");
        s += argv[m];
        ofstream os(s);
        m++;
        int k = i->size();
        double x[k+1];
        double error = 0;
        for(auto i = 0; i < k+1; i++){
            x[i] = a + (b - a) * i / (k);
        }
        for(auto j = 0; j != k; j++){
            os << x[j] << " " << x[j+1] << endl;
            gsl_function interpolate = f_factory<Polynomial>(&my_f<Polynomial>, &(*i)[j]);
            double temp = computeMaxError(g, interpolate, x[j], x[j+1], h);
            if( error < temp){
                error = temp;
            }
            os << (*i)[j];
        }
        cout << "N = " << k + 1 << ": max-norm error is " << error << endl;
    }
    return 0;
}
