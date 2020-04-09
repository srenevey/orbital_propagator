//
// Created by Sylvain  on 4/6/20.
//

#include "Area.h"


namespace Dimension {
    Area::Area(double a): BaseDimension(a) {}

    Area operator*(const Distance& d1, const Distance& d2) {
        return Area(d1.data() * d2.data());
    }
}

namespace Unit {
    // Base unit for areas used throughout the program is km^2.
    using namespace Dimension;
    Area operator "" _cm2(long double a) {
        return Area(a * 1E-10);
    }

    Area operator "" _m2(long double a) {
        return Area(a * 1E-6);
    }

    Area operator "" _km2(long double a) {
        return Area(a);
    }

    Area operator "" _in2(long double a) {
        return Area(a * 6.4516E-10);
    }

    Area operator "" _ft2(long double a) {
        return Area(a * 9.2903E-8);
    }

    Area operator "" _acre(long double a) {
        return Area(a * 0.00404686);
    }

    Area operator "" _ha(long double a) {
        return Area(a * 0.01);
    }

    Area operator "" _mi2(long double a) {
        return Area(a * 2.58999);
    }
}