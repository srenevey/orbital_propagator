//
// Created by Sylvain  on 4/13/20.
//

#ifndef ORBITAL_PROPAGATOR_ACCELERATION_H
#define ORBITAL_PROPAGATOR_ACCELERATION_H

#include "BaseDimension.h"

namespace Dimension {
    /** Definition of acceleration quantities. The base unit used throughout the program is km/s<sup>2</sup>. */
    class Acceleration : public BaseDimension {
    public:
        /**
         * @param acceleration in kilometers per second
         */
        explicit Acceleration(double acceleration = 0.);
        Acceleration(const BaseDimension& a);
    };
}

namespace Unit {
    // Base unit for accelerations used throughout the program is km/s^2.
    using namespace Dimension;
    /** Creates an acceleration in kilometers per square second (km/s<sup>2</sup>). */
    Acceleration operator "" _kms2(long double a);
}

#endif //ORBITAL_PROPAGATOR_ACCELERATION_H
