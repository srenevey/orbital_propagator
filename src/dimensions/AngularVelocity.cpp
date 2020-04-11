//
// Created by Sylvain  on 4/1/20.
//

#include "AngularVelocity.h"

namespace Dimension {
    AngularVelocity::AngularVelocity(double d): BaseDimension(d) {}
    AngularVelocity::AngularVelocity(const BaseDimension &a): BaseDimension(a) {}
}

namespace Unit {
    using namespace Dimension;
    AngularVelocity operator "" _rads(long double d) {
        return AngularVelocity(d);
    }

    AngularVelocity operator "" _degs(long double d) {
        return AngularVelocity(constants::DEG_TO_RAD * d);
    }
}