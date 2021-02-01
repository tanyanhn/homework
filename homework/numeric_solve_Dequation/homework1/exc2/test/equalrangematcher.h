#ifndef EQUALRANGMATCHER_TY
#define EQUALRANGMATCHER_TY

#include "catch.hpp"
#include <algorithm>
#include <iterator>

namespace solver_ty{
// ...

template<typename Range>
struct EqualsRangeMatcher : Catch::MatcherBase<Range> {
    explicit EqualsRangeMatcher(Range const& range, double e = 1e-5):
            range{ range } , epsilon(e)
    {}
    
    bool match(Range const& other) const {
        using std::begin; using std::end;

        return std::equal(begin(range), end(range), begin(other), digitcompare(epsilon));
    }

    std::string describe() const override {
        return "Equals: " + Catch::rangeToString(range);
    }

    
    struct digitcompare {
        double epsilon;
        explicit digitcompare(double e) : epsilon(e){}
        template<typename T>
        bool operator()(T a, T b){
            return (a - b) < epsilon && (b - a) < epsilon;
        }
    };

  private:
    double epsilon;
    Range const& range;
};

template<typename Range>
auto EqualsRange(const Range& range, double e = 1e-5) -> EqualsRangeMatcher<Range> {
    return EqualsRangeMatcher<Range>{range, e};
}

}
#endif