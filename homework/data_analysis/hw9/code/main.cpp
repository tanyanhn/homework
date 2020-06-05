#include"Discrete_least_square.h"
#include"normalequations.h"
#include"QR_factorization.h"
#include<fstream>
#include<iostream>

using namespace std;


int main(int argc, char *argv[])
{
    ifstream is1("data/file1");
    ifstream is2("data/file1");
    ofstream os1("data/file1.output");
    ofstream os2("data/file2.output");
    ofstream os3("data/matrix1"), os4("data/matrix2");
    Normalequation solver;
    QRfactor QRsolver;
    solver.initialize(is1);
    QRsolver.initialize(is2);
    double A1[3*3], A2[3*3];
    vector<double> anwser = solver.work(A1);
    vector<double> QRanwser = QRsolver.work(A2);
    for(auto i = 0; i < 9; ++i){
        os3 << A1[i] << " ";
        os4 << A2[i] << " ";
    }
    int k = 0;
    cout << "normalequations anwser:" << endl;
    for(auto i = anwser.begin(); i!= anwser.end(); ++i){
        // cout << "normalequations anwser:" << endl;
        os1 << *i << " ";
        cout << *i << " x^" << k++ << " + ";
    }
    cout << endl;
    k = 0;
    cout << "QRfactor anwser:" << endl;
    for(auto i = QRanwser.begin(); i!= QRanwser.end(); ++i){
        //  cout << "QRfactor anwser:" << endl;
        os2 << *i << " ";
        cout << *i << " x^" << k++ << " + ";
    }
    cout << endl;
    return 0;
}
