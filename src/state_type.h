//
// Created by Sylvain Renevey on 5/16/18.
//

#ifndef CPP_PROPAGATOR_STATE_TYPE_H
#define CPP_PROPAGATOR_STATE_TYPE_H

#include <array>

extern "C" {
    #include "SpiceUsr.h"
};

/*! \typedef state_type
 *  \brief Definition of the state_type type as a std array of doubles.
 */
typedef std::array<double, 6> state_type;

#endif //CPP_PROPAGATOR_STATE_TYPE_H
