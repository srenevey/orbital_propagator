//
// Created by Sylvain Renevey on 5/11/18.
//

#ifndef CPP_PROPAGATOR_EOM_H
#define CPP_PROPAGATOR_EOM_H

#include <string>
#include <array>
#include <Eigen/Dense>
#include "state_type.h"
#include "EnvironmentModel.h"
#include <iostream>
#include "Spacecraft.h"
extern "C" {
#include "SpiceUsr.h"
};

/*! \class EqMotion
 *  \brief Class defining the equations of motion.
 */

class EqMotion {
private:
    const EnvironmentModel& env_model_;
	SpiceDouble& initial_et_;
    const Spacecraft& spacecraft_;

public:
    /*!\fn EqMotion
    *  \brief Constructor of the class.
    *
    *  \param sc           Instance of the Spacecraft object containing the properties of the spacecraft.
    *  \param env_model    Instance of the EnvironmentModel class containing the parameters of the environment model.
    *  \param initial_et   Initial ephemeris time.
    */
    EqMotion(const Spacecraft &sc, const EnvironmentModel &env_model, SpiceDouble &initial_et);

    /*! \fn void operator()( const state_type &x , state_type &dxdt , double t )
     *  \brief Computes the time derivative of the state at a given time.
     *
     *  This function computes the time derivative of the six states of the spacecraft. It is used by the integrator to propagate the state.
     *
     *  \param x    Inertial state of the spacecraft in the J2000 frame for Earth bound orbits and Ecliptic J2000 otherwise.
     *  \param dxdt state_type instance where the time derivative of the state will be stored.
     *  \param t    Time at which the time derivative is taken
     */
    void operator()( const state_type &x , state_type &dxdt , double t );
};


#endif //CPP_PROPAGATOR_EOM_H
