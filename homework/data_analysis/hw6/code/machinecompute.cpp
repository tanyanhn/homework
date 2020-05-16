#include<fstream>
#include<cmath>

using namespace std;

double f(const double x){
    return pow(x, 8) - 8 * pow(x, 7) + 28 * pow(x, 6) - 56 * pow(x, 5) + 70 * pow(x, 4) - 56 * pow(x, 3) + 28 * pow(x, 2) - 8 * x + 1;
}

double g(const double x){
    return (((((((x - 8) * x + 28) * x - 56) * x + 70) * x - 56) * x + 28) * x - 8) * x + 1;
}

double h(const double x){
    return pow(x - 1, 8);
}

void test(ofstream& osf, ofstream& osg, ofstream& osh){
    double x = 0.99;
    double interval = 0.02/101.0;
    for(int i = 0; i < 101; i++){
        osf << f(x) << " ";
        osg << g(x) << " ";
        osh << h(x) << " ";
        x += interval;
    }
}

int main(int argc, char *argv[])
{
    ofstream osf("output/ma1"),
        osg("output/ma2"),
        osh("output/ma3");
    test(osf, osg, osh);
    return 0;
}
