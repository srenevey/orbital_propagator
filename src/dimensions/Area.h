//
// Created by Sylvain  on 4/6/20.
//

#ifndef ORBITAL_PROPAGATOR_AREA_H
#define ORBITAL_PROPAGATOR_AREA_H

#include "BaseDimension.h"
#include "dimensions/Distance.h"

namespace Dimension {
    /** Definition of area quantities. The base unit used throughout the program is km<sup>2</sup>. */
    class Area : public BaseDimension {
    public:
        /**
         * @param area in square kilometers
         */
        explicit Area(double area);
    };
}

namespace Unit {
    // Base unit for areas used throughout the program is km^2.
    using namespace Dimension;
    /** Creates an area in square centimeters (cm<sup>2</sup>). */
    Area operator "" _cm2(long double a);
    /** Creates an area in square meters (m<sup>2</sup>). */
    Area operator "" _m2(long double a);
    /** Creates an area in square kilometers (km<sup>2</sup>). */
    Area operator "" _km2(long double a);
    /** Creates an area in square inches (in<sup>2</sup>). */
    Area operator "" _in2(long double a);
    /** Creates an area in square feet (ft<sup>2</sup>). */
    Area operator "" _ft2(long double a);
    /** Creates an area in acres. */
    Area operator "" _acre(long double a);
    /** Creates an area in hectares (ha). */
    Area operator "" _ha(long double a);
    /** Creates an area in square miles (mi<sup>2</sup>). */
    Area operator "" _mi2(long double a);
}


#endif //ORBITAL_PROPAGATOR_AREA_H
