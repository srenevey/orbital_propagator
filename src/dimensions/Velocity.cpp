//
// Created by Sylvain  on 4/1/20.
//

#include "Velocity.h"

namespace Dimension {
    Velocity::Velocity(double velocity): BaseDimension(velocity) {}
    Velocity::Velocity(const BaseDimension& a): BaseDimension(a) {}
}

namespace Unit {
    // Base unit for velocities used throughout the program is km/s.
    using namespace Dimension;
    Velocity operator "" _mms(long double d) {
        return Velocity(1.0E-6 * d);
    }

    Velocity operator "" _cms(long double d) {
        return Velocity(1.0E-5 * d);
    }

    Velocity operator "" _ms(long double d) {
        return Velocity(1.0E-3 * d);
    }

    Velocity operator "" _kms(long double d) {
        return Velocity(d);
    }

    Velocity operator "" _kph(long double d) {
        return Velocity(d / 3600.0);
    }

    Velocity operator "" _mps(long double d) {
        return Velocity(1.60934 * d);
    }

    Velocity operator "" _mph(long double d) {
        return Velocity(4.4704E-4 * d);
    }

}