//
// Created by Sylvain Renevey on 8/17/18.
//

/*! \namespace perturbations
 * \brief This namespace contains methods to compute the different perturbations.
 */


#ifndef CPP_PROPAGATOR_PERTURBATIONS_H
#define CPP_PROPAGATOR_PERTURBATIONS_H

#include <array>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <string>

#include "EnvironmentModel.h"
#include "constants.h"
#include "Spacecraft.h"

extern "C" {
    #include "SpiceUsr.h"
};


namespace perturbations {

    /*! \fn Eigen::Vector3d drag(const Spacecraft &sc, const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, double elapsed_time, const EnvironmentModel &env_model)
     *  \brief Compute the perturbations induced by the drag.
     *
     * This function computes the perturbations produced by the atmospheric drag. Based on the settings of the environment model, the atmospheric density
     * is retrieved either from an exponential model or from EarthGRAM2016.
     *
     *  \param sc				Instance of the Spacecraft object containing the properties of the spacecraft.
     *  \param state			State vector of the spacecraft expressed in the J2000 frame (ECI).
     *  \param et				Ephemeris time.
     *  \param elapsed_time     Time elapsed since the beginning of the integration (sec).
     *  \param env_model		Instance of the EnvironmentModel class containing the parameters of the environment model.
     */
    Eigen::Vector3d drag(const Spacecraft &sc, const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, double elapsed_time, const EnvironmentModel &env_model);



    /*! \fn Eigen::Vector3d geopotential(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const EnvironmentModel &env_model)
     *  \brief Compute the gravitational perturbation induced by the nonspherical geopotential.
     *
     * This function computes the perturbations produced by the nonspherical geopotential.
     * The algorithm described in Montenbruck, O., Gill, E.,
     * <em>Satellites Orbits</em>, Springer-Verlag Berlin Heidelberg, 2000, is used. This method uses the EGM2008 model.
     *
     *  \param state        State vector of the spacecraft expressed in the J2000 frame (ECI).
     *  \param et           Ephemeris time.
     *  \param env_model    Instance of the EnvironmentModel class containing the parameters of the environment model.
     */
    Eigen::Vector3d geopotential(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const EnvironmentModel &env_model);


    /*! \fn Eigen::Vector3d nbody(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const EnvironmentModel &env_model)
     *  \brief Compute the gravitational perturbations produced by celestial bodies.
     *
     * This function computes the perturbations produced by all other planets and/or moons taken into account.
     * The algorithm is taken from Vallado, D.A., <em>Fundamentals of Astrodynamics and Applications</em>, Fourth Edition, Space Technology Library, 2013.
     *
     *  \param state        State vector of the spacecraft expressed in the J2000 or Ecliptic J2000 frame.
     *  \param et			Ephemeris time
     *  \param env_model    Instance of the EnvironmentModel class containing the parameters of the environment model.
     *
     *  \return Acceleration resulting from the third-body perturbations expressed in the J2000
     *  frame if the central body is the Earth and in the Ecliptic J2000 otherwise.
     */
    Eigen::Vector3d nbody(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const EnvironmentModel &env_model);



    /*! \fn Eigen::Vector3d perturbation_from_body(const Eigen::Ref<const Eigen::VectorXd>& state, const std::string& third_body, const std::string& central_body, SpiceDouble et)
     *  \brief Compute the perturbation induced by a single body.
     *
     *  This function computes the perturbing acceleration produced by a single body on the spacecraft at a given date.
     *  JPL's ephemerides are used to retrieve the body's state vector.
     *
     *  \param state        State vector of the spacecraft expressed in the J2000 or Ecliptic J2000 frame.
     *  \param third_body   Name of the perturbing body.
     *  \param central_body Celestial body to which the spacecraft is bound.
     *  \param et			Ephemeris time.
     */
    Eigen::Vector3d perturbation_from_body(const Eigen::Ref<const Eigen::VectorXd>& state, const std::string& third_body, const std::string& central_body, SpiceDouble et);


    /*! \fn Eigen::Vector3d solar_radiation_pressure(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const EnvironmentModel &env_model, const Spacecraft &sc)
     *  \brief Compute the perturbations induced by the solar radiation pressure.
     *
     * This function computes the orbital perturbations induced by the solar radiation pressure. It is assumed that the
     * surface normal points in the direction of the Sun.
     *
     *  \param state			State vector of the spacecraft expressed in the J2000 or Ecliptic J2000 frame.
     *  \param et				Ephemeris time.
     *  \param env_model		Instance of the EnvironmentModel class containing the parameters of the environment model.
     *  \param sc				Instance of the Spacecraft object containing the properties of the spacecraft.
     */
    Eigen::Vector3d solar_radiation_pressure(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const EnvironmentModel &env_model, const Spacecraft &sc);
}

#endif //CPP_PROPAGATOR_PERTURBATIONS_H
