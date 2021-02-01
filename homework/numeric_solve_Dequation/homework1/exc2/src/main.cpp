#include "coord.h"
#include <iostream>
#include "equationsolver.h"
#include "factory.h"
#include <fstream>
#include "norm.h"

using namespace std;
using namespace solver_ty;


int main(int argc, char *argv[])
{
    fstream fs(argv[1]);
    while(true){
        char c;
        fs.get(c);
        if(c == '%'){
            char comment[256];
            fs.getline(comment, 256);
        }
        else{
            fs.unget();
            break;
        }
    }
    int dim;
    int loops;
    int v[3];
    Real error;
    Real norme = 0, normr = 0, radio;
    fs >> dim;
    fs >> loops;
    fs >> error;
    int count = 0;
    while(true){
        int temp;
        fs >> temp;
        if(temp < 0){
            break;
        }
        v[count++] = temp;
    }
    std::string file, aim;
    if(dim == 1){
        while(fs >> aim){
            fs >> file;
            fstream fa(aim);
            std::unique_ptr<equationsolver<1>> p(solverFactory<1>(file));
            int size = p->size();
            std::unique_ptr<Real[]> m(new Real[size]);
            cout << "size of grid : " << size-1 << ", dim is: " << dim << endl << endl;
            for(auto k = 0; k < size; ++k){
                Real tmp;
                fa >> tmp;
                m[k] = tmp;
            }
            for(auto k = 1; k < loops+1; ++k){
                p->solver(v, 1);
                auto result = p->getresult();
                auto e = move(result.first), r = move(result.second);
                for(auto k = 0; k < size; ++k){
                    e[k] -= m[k];
                }
                Real tmpnorme = norm<0>(e.get(), size),
                        tmpnormr = norm<0>(r.get(), size);
                if(k == 1){
                    cout << "After " << k << "times solver, norm(e) is: " << tmpnorme
                         << " ,norm(r) is: " << tmpnormr << endl;
                }
                else{
                    radio = tmpnormr / normr;
                    cout << "After " << k << "times solver, norm(e) is: " << tmpnorme
                         << " ,norm(r) is: " << tmpnormr << " , rate is: " << radio << endl;
                }
                normr = tmpnormr; norme = tmpnorme;
                if(normr < error){
                    cout << "\n Already reach stop criteria.\n";
                }
                else if(k == loops){
                    cout << "\n Haven't' reach stop criteria.\n";
                }
            }
            cout << endl << "========================================================================" << endl << endl;
        }
    }
    return 0;
}
