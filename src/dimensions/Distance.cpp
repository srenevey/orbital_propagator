//
// Created by Sylvain  on 4/1/20.
//

#include "Distance.h"

namespace Dimension {
    Distance::Distance(double d): BaseDimension(d) {}
    Distance::Distance(const BaseDimension& a): BaseDimension(a) {}
}

namespace Unit {
    // Base unit for distances used throughout the program is km.
    using namespace Dimension;
    Distance operator "" _mm(long double d) {
        return Distance(1.0E-6 * d);
    }

    Distance operator "" _cm(long double d) {
        return Distance(1.0E-5 * d);
    }

    Distance operator "" _dm(long double d) {
        return Distance(1.0E-4 * d);
    }

    Distance operator "" _m(long double d) {
        return Distance(1.0E-3 * d);
    }

    Distance operator "" _km(long double d) {
        return Distance(d);
    }

    Distance operator "" _AU(long double d) {
        return Distance(constants::AU * d);
    }

    Distance operator "" _in(long double d) {
        return Distance(2.54E-5 * d);
    }

    Distance operator "" _ft(long double d) {
        return Distance(3.048E-4 * d);
    }

    Distance operator "" _mi(long double d) {
        return Distance(1.60934 * d);
    }
}