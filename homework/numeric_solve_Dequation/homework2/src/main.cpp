#include "factory.h"
#include "ThreeBody.h"
#include "norm.h"
#include <fstream>
#include <iostream>
#include <memory>

using namespace std;
using namespace integrator_ty;
using namespace solver_ty;

int main(int argc, char *argv[])
{
    fstream fs(argv[1]);
    char comment[1024];
    while(true){
        char c;
        fs.get(c);
        if(c == '%'){
            fs.getline(comment, 1024);
        }
        else{
            fs.unget();
            break;
        }
    }
    string file, aimf;
    Factory fac;
    unique_ptr<Real[]> aim(new Real[2400001 * 6]),
            error1(new Real[2400001 * 6]),
            error2(new Real[2400001 * 6]);
    double lastn[6];   double p1, p2; double start[4], end[4], time[4];
    auto count = 0;
    while(fs >> file){
        fs >> aimf;
        ++count;
        cout << "\n\n\n In initial value " << count << endl;
        int step[] = {6000, 12000, 24000, 48000, 96000, 120000, 240000, 480000};
        int fstep= 2400000, size = 6;
        fstream afs(aimf);
        for(auto i = 0; i < 2400001; ++i){
            for(auto j = 0; j < 6; ++j){
                Real temp;
                afs >> temp;
                aim[i * 6 + j]  = temp;
            }
        }
        unique_ptr<TimeIntegrator> pRKsolver(fac(file)),
                pEusolver(fac(file, 1));
        for(auto stepk = 0; stepk < 8; ++stepk){
            int cstep = step[stepk];
            pRKsolver->setvalue<Real>(string("timeSt"), cstep);
            pEusolver->setvalue<Real>(string("timeSt"), cstep);
            auto start1 = clock();
            std::unique_ptr<Real[]> RKout(pRKsolver->solver());
            auto end1 = clock();
            auto start2 = clock();
            unique_ptr<Real[]> Euout(pEusolver->solver());
            auto end2 = clock();
            double time1 = (double)(end1-start1)/CLOCKS_PER_SEC,
                    time2 = (double)(end2-start2)/CLOCKS_PER_SEC;
            int multi = round(double(fstep) / cstep);
            for(auto istep = 0; istep <= cstep; ++istep){
                for(auto i = 0; i < 6; i++){
                    error1[istep * size + i] = RKout[istep * size + i] -  aim[istep * multi * size + i];
                    error2[istep * size + i] = Euout[istep * size + i] -  aim[istep * multi * size + i];
                }
            }
            auto norme1 = norm<0>(error1.get(), (cstep+1) * size),
                    norme2 = norm<0>(error2.get(), (cstep+1) * size);
            if(stepk == 0){
                lastn[4] = norme1; lastn[5] = norme2;
            }
            else {
                p1 = log(norme1 / lastn[4]) / log(double(step[stepk-1]) / step[stepk]);
                p2 = log(norme2 / lastn[5]) / log(double(step[stepk-1]) / step[stepk]);
                lastn[4] = norme1; lastn[5] = norme2;
            }
            cout << "\n Steps is : " << cstep << endl << endl
                 << "norm of error in RK method : " << norme1 << endl
                 << "CPU calculate time :" << time1 << " ms. \n";
            if(stepk !=0){
                cout << "convergence rate :" << p1 << endl <<  endl;
            }
            cout << "norm of error in Eu method : " << norme2 << endl
                 << "CPU calculate time :" << time2 << " ms. \n";
            if(stepk !=0){
                cout << "convergence rate :" << p2 << endl << endl;
            }
            for(int p = 0; p <4; ++p){
                unique_ptr<TimeIntegrator> psolver(fac(file, p+1));
                psolver->setvalue("met", 0);
                psolver->setvalue<Real>(string("timeSt"), cstep);
                start[p] = clock();
                unique_ptr<Real[]> pout(psolver->solver());
                end[p] = clock();
                time[p] = (double)(end[p]-start[p])/CLOCKS_PER_SEC;
                for(auto istep = 0; istep <= cstep; ++istep){
                    for(auto i = 0; i < 6; i++){
                        error1[istep * size + i] =
                                pout[istep * size + i] -
                                aim[istep * multi * size + i];
                    }
                }
                auto norme1 = norm<0>(error1.get(), (cstep+1) * size);
                if(stepk == 0){
                    lastn[p] = norme1;
                }
                else {
                    p1 = log(norme1 / lastn[p]) / log(double(step[stepk-1]) / step[stepk]);
                    lastn[p] = norme1;
                }
                cout << "norm of error in A-Bashforth method p = " << p+1
                     << ", error : " << norme1 << endl
                     << "CPU calculate time :" << time[p] << " ms. \n";
                if(stepk !=0){
                    cout << "convergence rate :" << p1 << endl << endl;
                }
            }
            // if(norme1 < 1e-3 && norme2 < 1e-3){
            //     cout << "\n Arealdy reach accuratly " << 1e-3 << endl;
            //     break;
            // }
        }
        cout << "============================================================" << endl;
    }
    return 0;
}
