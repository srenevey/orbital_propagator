//
// Created by Sylvain  on 4/9/20.
//

#ifndef ORBITAL_PROPAGATOR_ATMMODEL_H
#define ORBITAL_PROPAGATOR_ATMMODEL_H

/*!
 * @file
 * Atmosphere models definitions.
 */

/** Enumeration of the diffente atmospheric models implemented in the program. */
enum class AtmModel {
    EarthGRAM, /*!< EarthGRAM 2016 model. */
    Exponential, /*!< Exponentially decaying atmospheric model. */
    None /*!< No atmospheric model is used. */
};

#endif //ORBITAL_PROPAGATOR_ATMMODEL_H
