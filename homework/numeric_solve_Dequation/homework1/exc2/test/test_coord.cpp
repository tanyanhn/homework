#include "../src/coord.h"
#include "catch.hpp"
#include "equalrangematcher.h"

using namespace solver_ty;

SCENARIO("getnearcoord and getdiagnolcoord test", "[coordinate]"){
    int i[1] = {8};
    coordinate<1> coord1(i);
    std::vector<int> coord1_near_3_aim = {2, 4};
    auto temp = coord1.getnearcoord(3);
    std::vector<int> coord1_near_3_out{temp[0], temp[1]};
    REQUIRE_THAT(coord1_near_3_out, EqualsRange(coord1_near_3_aim));
    temp.reset(coord1.getnearcoord(7).release());
    std::vector<int> coord1_near_7_aim = {6, 8},
            coord1_near_7_out{temp[0], temp[1]};
    REQUIRE_THAT(coord1_near_7_out, EqualsRange(coord1_near_7_aim));
    int d[2] = { 9, 7};
    coordinate<2> coord2(d);
    temp.reset(coord2.getnearcoord(35).release());
    auto coord2_near_35_aim = {27, 34, 36, 43},
            coord2_near_35_out = {temp[0], temp[1], temp[2], temp[3]};
    REQUIRE_THAT(coord2_near_35_out, EqualsRange(coord2_near_35_aim));
    temp.reset(coord2.getnearcoord(51).release());
    auto coord2_near_51_aim = {43, 50, 52, 59},
            coord2_near_51_out = {temp[0], temp[1], temp[2], temp[3]};
    REQUIRE_THAT(coord2_near_51_out, EqualsRange(coord2_near_51_aim));
    temp.reset(coord2.getdiagnolcoord(38).release());
    auto coord2_diag_38_aim = {29, 31, 45, 47},
            coord2_diag_38_out = {temp[0], temp[1], temp[2], temp[3]};
    REQUIRE_THAT(coord2_diag_38_out, EqualsRange(coord2_diag_38_aim));
    temp.reset(coord2.getdiagnolcoord(65).release());
    auto coord2_diag_65_aim = {56, 58, 72, 74},
            coord2_diag_65_out = {temp[0], temp[1], temp[2], temp[3]};
    REQUIRE_THAT(coord2_diag_65_out, EqualsRange(coord2_diag_65_aim));
    //coord2.getdiagnolcoord(73);
}
