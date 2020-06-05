#include"Discrete_least_square.h"
#include"normalequations.h"
#include<fstream>
#include<iostream>

using namespace std;


int main(int argc, char *argv[])
{
    ifstream is1("data/file1");
    ofstream os1("data/file1.output");
    Normalequation solver;
    solver.initialize(is1);
    vector<double> anwser = solver.work();
    int k = 0;
    for(auto i = anwser.begin(); i!= anwser.end(); ++i){
        os1 << *i << " ";
        cout << *i << " x^" << k++ << " + ";
    }
    cout << endl;
    return 0;
}
