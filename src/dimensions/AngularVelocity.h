//
// Created by Sylvain  on 4/1/20.
//

#ifndef ORBITAL_PROPAGATOR_ANGULARVELOCITY_H
#define ORBITAL_PROPAGATOR_ANGULARVELOCITY_H

#include "BaseDimension.h"
#include "constants.h"

namespace Dimension {
    /** Definition of angular velocity quantities. The base unit used throughout the program is rad/s. */
    class AngularVelocity : public BaseDimension {
    public:
        /**
         * @param angular_velocity in radians per second
         */
        explicit AngularVelocity(double angular_velocity = 0.);
        AngularVelocity(const BaseDimension& a);
    };
}

namespace Unit {
    // Base unit for angular velocities used throughout the program is rad/s.
    using namespace Dimension;

    /** Creates an angular velocity in radians per second (rad/s). */
    AngularVelocity operator "" _rads(long double d);
    /** Creates an angular velocity in degrees per second (deg/s). */
    AngularVelocity operator "" _degs(long double d);
}


#endif //ORBITAL_PROPAGATOR_ANGULARVELOCITY_H
