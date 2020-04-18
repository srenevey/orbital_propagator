//
// Created by Sylvain  on 4/14/20.
//

#ifndef ORBITAL_PROPAGATOR_ANGULARACCELERATION_H
#define ORBITAL_PROPAGATOR_ANGULARACCELERATION_H

#include "BaseDimension.h"

namespace Dimension {
    /** Definition of angular acceleration quantities. The base unit used throughout the program is rad/s<sup>2</sup>. */
    class AngularAcceleration: public BaseDimension {
    public:
        /**
         * @param acceleration in radians per square second
         */
        explicit AngularAcceleration(double acceleration = 0.);
        explicit AngularAcceleration(const BaseDimension& a);
    };
}

namespace Unit {
    // Base unit for angular accelerations used throughout the program is rad/s^2.
    using namespace Dimension;
    /** Creates an angular acceleration in radians per square second (rad/s<sup>2</sup>). */
    AngularAcceleration operator "" _rads2(long double a);
}


#endif //ORBITAL_PROPAGATOR_ANGULARACCELERATION_H
