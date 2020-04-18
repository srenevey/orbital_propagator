//
// Created by Sylvain  on 4/1/20.
//

#ifndef ORBITAL_PROPAGATOR_TIME_H
#define ORBITAL_PROPAGATOR_TIME_H

#include "BaseDimension.h"

namespace Dimension {
    /** Definition of time quantities. The base unit used throughout the program is second. */
    class Time: public BaseDimension {
    public:
        /**
         * @param time in seconds
         */
        explicit Time(double time = 0.);
        Time(const BaseDimension& a);
    };
}

namespace Unit {
    // Base unit for time used throughout the program is second.
    using namespace Dimension;
    /** Creates a time in seconds (s). */
    Time operator "" _s(long double d);
    /** Creates a time in minutes (min). */
    Time operator "" _min(long double d);
    /** Creates a time in hours (h). */
    Time operator "" _h(long double d);
    /** Creates a time in days (d). */
    Time operator "" _d(long double d);
}

#endif //ORBITAL_PROPAGATOR_TIME_H
