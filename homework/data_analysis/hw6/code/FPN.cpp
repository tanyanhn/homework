#include<fstream>
#include<cmath>
#include <iostream>
#include<vector>
#include<algorithm>
using namespace std;

double printBinary(const vector<int>& m, const int beta, const int e){
    double anwser = 0;
    int k = e;
    for_each(m.begin(), m.end(), [&](int a){
                                     anwser += a * pow(beta, k--);
                                 });
    return anwser;
}


class FPN{
    int beta, p, L, U;
public:
    FPN(const int bb, const int pp, const int ll, const int uu)
        : beta(bb), p(pp), L(ll), U(uu) {}
    ~FPN(){}
    double UFL() const {
        return pow(beta, L);
    }
    double OFL() const {
        return pow(beta, U) * (beta - pow(beta, 1 - p));
    }
    void enumerate(ostream& os) const {
        vector<int> m(p);
        int num = 0;
        for(auto i = m.begin(); i != m.end(); i++){
            *i = 0;
        }
        os << 0 << " ";
        ++num;
        m[0] = 1;
        for(auto e = L; !(e > U); e++){
            bool pluse = false;
            while(pluse == false){
                double y = printBinary(m, beta, e);
                os << y << " " << -y << " ";
                ++num;
                num++;
                bool plusm = false;
                while(plusm == false){
                    for(auto i = p-1; i > -1; i--){
                        if(m[i] != beta -1){
                            m[i] += 1;
                            plusm = true;
                            break;
                        }
                        else if(i != 0){
                            m[i] = 0;
                        }
                        else{
                            plusm = true;
                            pluse = true;
                            m[0] = 1;
                            break;
                        }
                    }
                }
            }
        }
        os << endl << num << endl;
    }
    void enumeratesubnormal(ostream& os) const {
        vector<int> m(p);
        int num = 0;
        for(auto i = m.begin(); i != m.end(); i++){
            *i = 0;
        }
        m[p-1] = 1;
        bool pluse = false;
        while(pluse == false){
            double y = printBinary(m, beta, L);
            os << y << " " << -y << " ";
            ++num;
            num++;
            bool plusm = false;
            while(plusm == false){
                for(auto i = p-1; i > -1; i--){
                    if(m[i] != beta -1){
                        m[i] += 1;
                        plusm = true;
                        break;
                    }
                    else if(i != 1){
                        m[i] = 0;
                    }
                    else{
                        plusm = true;
                        pluse = true;
                        //m[0] = 1;
                        break;
                    }
                }
            }
        }
        os << endl << num << endl;
    }
};





int main(int argc, char *argv[])
{
    int beta = 2, p = 3, L = -1, U = 1;
    ofstream osf("output/F1"),
        osexf("output/F2");
    //double UFL = pow(beta, L),
    //    OFL = pow(beta, U) * (beta - pow(beta, 1 - p));
    FPN F(beta, p, L, U);
    cout << "UFL: " << F.UFL() << endl
         << "OFL: " << F.OFL() << endl;
    F.enumerate(cout);
    F.enumeratesubnormal(cout);
    F.enumerate(osf);
    F.enumerate(osexf);
    F.enumeratesubnormal(osexf);
    return 0;
}
