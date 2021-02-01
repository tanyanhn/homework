#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER


#include <cmath>
#include "real.h"
#include "coord.h"
#include <memory>

namespace solver_ty{
enum boundary_type{Diri =0 , Neum=1};
enum restrict_operator{full = 0, inject = 1};
enum interpolat_operator{linear = 0, quadratic = 1};
enum cycle_type{V_cycle = 0, full_multi_cycle = 1};
enum stop_criteria{maxstep = 0, relate_error = 1};
enum relax_type{Jacobi = 0, weightJacobi = 1, GaussSeidel = 2, RB_GaussSeidel = 3};

template<int dim>
class cycle;

template<int dim>
class equationsolver {
  protected:
    int n[dim];
    Real h[dim];
    std::unique_ptr<Real[]> fvalue;
    std::unique_ptr<Real[]> value;
    std::unique_ptr<Real[]> rvalue;
    std::unique_ptr<Real[]> dfvalue;
    std::unique_ptr<basecoord<dim>> pcoord;
    boundary_type btype = Diri;
    restrict_operator restrictop = full;
    interpolat_operator interpolatop = linear;
    cycle_type ctype = V_cycle;
    stop_criteria stopcondi = maxstep;
    relax_type rtype = GaussSeidel;
  public:
    equationsolver(int nn[dim], Real hh[dim],
                   std::unique_ptr<Real[]> f,
                   std::unique_ptr<Real[]> initialvalue = nullptr);
    equationsolver(const equationsolver&);
    equationsolver& operator=(const equationsolver&);
    bool is_coarsest() const;
    void plusvalue(const equationsolver& addv, bool boundary = true);
    void minusvalue(const equationsolver& addv, bool boundayr = true);
    void assignvalue(const equationsolver& addv, bool boundary = true);
    void setvaluezero(bool boundary = true);
    std::pair<std::unique_ptr<Real[]>, std::unique_ptr<Real[]>> getresult();
    bool is_boundary(int) const;
    bool is_homogeneous() const;
    int size() const{
        int size = 1;
        for(auto i = 0; i < dim; ++i){
            size *= (n[i] + 1);
        }
        return size;
    }
    virtual void copyflag(equationsolver<dim>& );
    virtual bool chooseoption(std::string, int);
    virtual std::unique_ptr<Real[]> test_getvalue(std::string) const;
    virtual void solver(int*, Real st = 15);
    virtual void Relax();
    virtual std::unique_ptr<equationsolver<dim>> Restrict();
    virtual std::unique_ptr<equationsolver<dim>> Interpolate();
    virtual void updatervalue();
    virtual ~equationsolver(){}
    friend class cycle<dim>;
};

}
#endif