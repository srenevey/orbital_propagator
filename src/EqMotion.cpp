//
// Created by Sylvain Renevey on 5/11/18.
//

#include "EqMotion.h"

EqMotion::EqMotion(const Spacecraft &sc, const EnvironmentModel &env_model, SpiceDouble &initial_et): env_model_(env_model), initial_et_(initial_et), spacecraft_(sc) {};

void EqMotion::operator()( const state_type &x , state_type &dxdt , double t ) {
    Eigen::Matrix<double, 6, 1> state;
    state << x[0], x[1], x[2], x[3], x[4], x[5];

    // Compute the current ephemeris time
	SpiceDouble et = initial_et_ + t;

	// Two-body acceleration
    double r = state.head(3).norm();
    Eigen::Vector3d a_kepler = - env_model_.mu() / pow(r, 3) * state.head(3);


    // Perturbations
    Eigen::Vector3d a_geopot(0.0, 0.0, 0.0);
    if (env_model_.gp_degree() > 1) {
        a_geopot = env_model_.geopotential(state, et);
    }

    Eigen::Vector3d a_nbody(0.0, 0.0, 0.0);
    if (env_model_.is_third_body_flag()) {
        a_nbody = env_model_.nbody(state, et);
    }

    Eigen::Vector3d a_drag(0.0, 0.0, 0.0);
    if (env_model_.is_drag_flag() && spacecraft_.mass() != 0.0) {
        a_drag = env_model_.drag(spacecraft_, state, et, t);
    }

    Eigen::Vector3d a_srp(0.0, 0.0, 0.0);
    if (env_model_.is_srp_flag()) {
        a_srp = env_model_.solar_radiation_pressure(state, et, spacecraft_);
    }

    Eigen::Vector3d a_total = a_kepler + a_geopot + a_nbody + a_drag + a_srp;

    dxdt[0] = state(3);
    dxdt[1] = state(4);
    dxdt[2] = state(5);
    dxdt[3] = a_total(0);
    dxdt[4] = a_total(1);
    dxdt[5] = a_total(2);
}