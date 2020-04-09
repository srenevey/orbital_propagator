//
// Created by Sylvain Renevey on 5/11/18.
//

#ifndef CPP_PROPAGATOR_EOM_H
#define CPP_PROPAGATOR_EOM_H

#include <string>
#include <array>
#include <Eigen/Dense>
#include <iostream>
#include "EnvironmentModel.h"
#include "Spacecraft.h"
#include "State.h"
#include "Quaternion.h"
extern "C" {
#include "SpiceUsr.h"
};

/** Class defining the equations of motion. */
class EqMotion {
private:
    const EnvironmentModel &m_env_model;
    SpiceDouble &m_initial_et;
    const Spacecraft &m_spacecraft;

public:
    /**
    *  @param sc           Instance of the Spacecraft object containing the properties of the spacecraft.
    *  @param env_model    Instance of the EnvironmentModel class containing the parameters of the environment model.
    *  @param initial_et   Initial ephemeris time.
    */
    EqMotion(const Spacecraft &sc, const EnvironmentModel &env_model, SpiceDouble &initial_et);

    /** Computes the time derivative of the state at a given time.
     *
     *  This function computes the time derivative of the thirteen states of the spacecraft. It is used by the integrator to propagate the states.
     *
     *  @param[in] state                Inertial state of the spacecraft.
     *  @param[in,out] state_derivative     State vector where the time derivatives of the states are stored.
     *  @param[in,out] t                    Time at which the time derivatives are taken
     */
    void operator()(const State &state, State &state_derivative, double t);

    /** Computes the gravitational perturbation induced by the nonspherical geopotential.
     *
     * This function computes the perturbations produced by the nonspherical geopotential.
     * The algorithm described in Montenbruck, O., Gill, E.,
     * <em>Satellites Orbits</em>, Springer-Verlag Berlin Heidelberg, 2000, is used. This method uses the EGM2008 model.
     *
     *  @param[out] acc              Reference to the acceleration vector.
     *  @param[out] torque           Reference to the torque vector.
     *  @param[in] state			State vector of the spacecraft.
     *  @param[in] et				Ephemeris time.
     */
    void geopotential(Eigen::Vector3d &acc, Eigen::Vector3d &torque, const State &state, SpiceDouble et) const;

    /** Computes the perturbations induced by the drag.
     *
     *  This function computes the perturbations produced by the atmospheric drag. Based on the settings of the environment model, the atmospheric density
     *  is retrieved either from an exponential model or from EarthGRAM2016.
     *
     *  @param[out] acc              Reference to the acceleration vector.
     *  @param[out] torque           Reference to the torque vector.
     *  @param[in] state			State vector of the spacecraft expressed in the J2000 frame (ECI).
     *  @param[in] et				Ephemeris time.
     *  @param[in] elapsed_time     Time elapsed since the beginning of the integration (sec).
     */
    void drag(Eigen::Vector3d &acc, Eigen::Vector3d &torque, const State &state, SpiceDouble et, double elapsed_time) const;

    /** Computes the perturbations induced by the solar radiation pressure.
     *
     *  @param[out] acc              Reference to the acceleration vector.
     *  @param[out] torque           Reference to the torque vector.
     *  @param[in] state		    State vector of the spacecraft.
     *  @param[in] et				Ephemeris time.
     */
    void solar_radiation_pressure(Eigen::Vector3d &acc, Eigen::Vector3d &torque, const State &state, SpiceDouble et) const;

    /** Computes the perturbations induced by the Earth magnetic field.
     *
     *  @param[out] torque           Reference to the torque vector.
     *  @param[in] state			State vector of the spacecraft expressed in the J2000 frame.
     *  @param[in] et				Ephemeris time.
     */
    void magnetic_perturbations(Eigen::Vector3d &torque, const State &state, SpiceDouble et) const;

    /** Computes the perturbations induced by third bodies.
     *
     *  @param[out] acc              Reference to the acceleration vector.
     *  @param[in] state			State vector of the spacecraft.
     *  @param[in] et				Ephemeris time.
     */
    void third_body_effect(Eigen::Vector3d &acc, const State &state, SpiceDouble et) const;
};

#endif //CPP_PROPAGATOR_EOM_H
