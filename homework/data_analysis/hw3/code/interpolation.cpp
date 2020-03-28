#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include"function.h"

using namespace std;

double FUN0(const vector<double>& parameter){
    double x = parameter[0];
    return 1 / (1 + x * x);
}

double FUN2(const vector<double>& parameter){
    double x = parameter[0];
    return 1 / (1 + 25 * x * x);
}


int main(int argc, char *argv[])
{
    F f(FUN0);
    function classfun(f);
    classfun.setab(-5, 5);
    string outputnewton("Data/newton");
    for(int n = 2; n < 10; n += 2){
        classfun.clear();
        classfun.newtoninterpolation(n);
        vector<pair<double,double>> coefficient = classfun.getcoefficient();
        string s = outputnewton + " " + to_string(n);
        ofstream stream(s);
        for(auto i = coefficient.begin(); i != coefficient.end(); i++){
            stream << i->first << " ";
        }
        stream << endl;
        for(auto i = coefficient.begin(); i != coefficient.end(); i++){
            stream << i->second << " ";
        }
    }
    F f2(FUN2);
    function classfun2(f2);
    classfun2.setab(-1, 1);
    function classfun3 = classfun2;
    string outputnewton2("Data/newton2"),
        outputchebyshevi3("Data/chebyshev3");
    for(int n = 5; n < 25; n += 5){
        classfun2.clear();
        classfun3.clear();
        classfun2.newtoninterpolation(n);
        classfun3.chebyshevinterpolation(n);
        vector<pair<double,double>> coefficient2 = classfun2.getcoefficient(),
            coefficient3 = classfun3.getcoefficient();
        string s2 = outputnewton2 + " " + to_string(n),
            s3 = outputchebyshevi3 + " " + to_string(n);
        ofstream stream2(s2),
            stream3(s3);
        for(auto i = coefficient2.begin(); i != coefficient2.end(); i++){
            stream2 << i->first << " ";
        }
        stream2 << endl;
        for(auto i = coefficient2.begin(); i != coefficient2.end(); i++){
            stream2 << i->second << " ";
        }
        for(auto i = coefficient3.begin(); i != coefficient3.end(); i++){
            stream3 << i->first << " ";
        }
        stream3 << endl;
        for(auto i = coefficient3.begin(); i != coefficient3.end(); i++){
            stream3 << i->second << " ";
        }
    }
    return 0;
}
