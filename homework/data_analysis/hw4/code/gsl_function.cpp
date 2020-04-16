#include"gsl_function.h"



template<>
double my_f</*Trigonometric*/Fractional>(double x, void* p){
    /**
     * calculate tanx - x
     *
     * x: position where f calculate
     * p: meanless in this function, but to meet requirement of gsl_function
     * return f value
     */
    /*
      assert(p == 0 &&
      "Trigonometric function shouldn't have params");
      double d = 0;
      d = tan(x) - x;
      return d;
    */
    assert(p != 0 &&
           "This is fractional function, should have params");
    my_f_params * params = static_cast<my_f_params*>(p);
    double d1 = 0, d2 = 0;
    std::vector<double> a = params->get()[0],b;
    x = x - params->getmid();
    if(params->get().size() > 1)
        b = params->get()[1];
    for_each(a.begin(), a.end(), [&d1, x] (int k) {d1 = d1 * x + k;});
    assert(d1 != 0 &&
           "Denominator shouldn't be 0.");
    d1 = 1/d1;
    if(b.empty())
        return d1;
    for_each(b.begin(), b.end(), [&d2, x] (int k) {d2 = d2 * x + k;});
    d1 = d2 * d1;
    return d1;
}


template<>
double my_f<Trigonometric>(double x, void* p){
    /**
     * calculate tanx - x
     *
     * x: position where f calculate
     * p: meanless in this function, but to meet requirement of gsl_function
     * return f value
     */
    //assert(p == 0 && "Trigonometric function shouldn't have params");
    my_f_params * params = static_cast<my_f_params *>(p);
    x = x - params->getmid();
    double d = 0;
    d = tan(x) - x;
    return d;
}