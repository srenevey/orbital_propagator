//
// Created by Sylvain  on 3/28/20.
//

#ifndef ORBITAL_PROPAGATOR_REFERENCEFRAME_H
#define ORBITAL_PROPAGATOR_REFERENCEFRAME_H

/*!
 * @file
 * Contains definitions of reference frames and coordinate systems (currently not used).
 */

/** Definition of different reference frames used throught the program. */
enum class ReferenceFrame {
    J2000, /*!< Earth Centered Inertial */
    ITRF93, /*!< Earth Centered, Earth Fixed */
    ICRF /*!< International Celestial Reference Frame */
};


/** Definition of different coordinate systems used throught the program (currently not used). */
enum class CoordinateSystem {
    Cartesian, /*!< Cartesian coordinate system */
    Spherical, /*!< Spherical coordinate system */
};

#endif //ORBITAL_PROPAGATOR_REFERENCEFRAME_H
