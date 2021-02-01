#ifndef COORD
#define COORD


#include "real.h"
#include "coord.h"
#include <assert.h>
#include <cmath>
#include <memory>

namespace solver_ty{

template<int dim>
class basecoord  {
  protected:
    int n[dim];
  public:
    basecoord(int nn[dim]) {
        for(auto i = 0; i < dim; ++i){
            n[i] = nn[i];
        }
    }
    virtual std::unique_ptr<int[]> getnearcoord(int) const{}
    virtual std::unique_ptr<int[]> getdiagnolcoord(int) const{}
    virtual ~basecoord() = 0;
};

template<int dim>
basecoord<dim>::~basecoord(){};

template<int dim>
class coordinate : public basecoord<dim> {  
  public:
    explicit coordinate(int nn[dim]) : basecoord<dim>(nn) {}
    virtual std::unique_ptr<int[]> getnearcoord(int co) const override;
    virtual std::unique_ptr<int[]> getdiagnolcoord(int co) const override;
    ~coordinate(){}
};


}

#endif