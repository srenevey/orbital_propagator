//
// Created by Sylvain  on 4/13/20.
//

#include "Acceleration.h"

namespace Dimension {
    Acceleration::Acceleration(double acceleration): BaseDimension(acceleration) {}
    Acceleration::Acceleration(const BaseDimension& a): BaseDimension(a) {}
}

namespace Unit {
    // Base unit for accelerations used throughout the program is km/s^2.
    using namespace Dimension;

    Acceleration operator "" _kms2(long double a) {
        return Acceleration(a);
    }
}