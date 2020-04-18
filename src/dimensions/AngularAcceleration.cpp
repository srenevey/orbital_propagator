//
// Created by Sylvain  on 4/14/20.
//

#include "AngularAcceleration.h"

namespace Dimension {
    AngularAcceleration::AngularAcceleration(double acceleration): BaseDimension(acceleration) {}
    AngularAcceleration::AngularAcceleration(const BaseDimension& a): BaseDimension(a) {}
}

namespace Unit {
    // Base unit for angular accelerations used throughout the program is rad/s^2.
    using namespace Dimension;
    AngularAcceleration operator "" _rads2(long double a) {
        return AngularAcceleration(a);
    }
}