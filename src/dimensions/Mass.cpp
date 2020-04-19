//
// Created by Sylvain  on 4/19/20.
//

#include "Mass.h"

namespace Dimension {
    Mass::Mass(double mass): BaseDimension(mass) {}
    Mass::Mass(const BaseDimension& a): BaseDimension(a) {}
}

namespace Unit {
    using namespace Dimension;

    Mass operator "" _g(long double m) {
        return Mass(1E-3 * m);
    }

    Mass operator "" _kg(long double m) {
        return Mass(m);
    }

    Mass operator "" _t(long double m) {
        return Mass(1E3 * m);
    }

    Mass operator "" _oz(long double m) {
        return Mass(28.349523125E-3 * m);
    }

    Mass operator "" _lb(long double m) {
        return Mass(0.45359237 * m);
    }

    Mass operator "" _ust(long double m) {
        return Mass(907.185 * m);
    }
}