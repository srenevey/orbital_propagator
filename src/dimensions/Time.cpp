//
// Created by Sylvain  on 4/1/20.
//

#include "Time.h"

namespace Dimension {
    Time::Time(const double time): BaseDimension(time) {}
    Time::Time(const BaseDimension& a): BaseDimension(a) {}
}

namespace Unit {
    using namespace Dimension;
    Time operator "" _s(long double d) {
        return Time(d);
    }

    Time operator "" _min(long double d) {
        return Time(60.0 * d);
    }

    Time operator "" _h(long double d) {
        return Time(3600.0 * d);
    }

    Time operator "" _d(long double d) {
        return Time(24. * 3600. * d);
    }
}