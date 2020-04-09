//
// Created by Sylvain Renevey on 5/11/18.
//

#include "EqMotion.h"

EqMotion::EqMotion(const Spacecraft &sc, const EnvironmentModel &env_model, SpiceDouble &initial_et): m_env_model(env_model), m_initial_et(initial_et), m_spacecraft(sc) {};

void EqMotion::operator()(const State &state, State &state_derivative, double t ) {

    // Compute the current ephemeris time
	SpiceDouble et = m_initial_et + t;

	// Two-body acceleration
    Eigen::Vector3d a_kepler = - m_env_model.central_body()->mu() / pow(state.position().norm(), 3) * state.position();

    // Geopotential
    Eigen::Vector3d a_geopot(0., 0., 0.);
    Eigen::Vector3d t_gg(0., 0., 0.);
    if (m_env_model.gp_degree() > 1) {
        geopotential(a_geopot, t_gg, state, et);
    }

    // Third-body effect
    Eigen::Vector3d a_third_body(0., 0., 0.);
    if (!m_env_model.third_body().empty()) {
         third_body_effect(a_third_body, state, et);
    }

    // Aerodynamic perturbations
    Eigen::Vector3d a_aero(0., 0., 0.);
    Eigen::Vector3d t_aero(0., 0., 0.);
    if (m_env_model.is_drag() && m_spacecraft.mass() != 0.0) {
        drag(a_aero, t_aero, state, et, t);
    }

    // Solar radiation pressure
    Eigen::Vector3d a_srp(0., 0., 0.);
    Eigen::Vector3d t_srp(0., 0., 0.);
    if (m_env_model.is_srp_flag()) {
        solar_radiation_pressure(a_srp, t_srp, state, et);
    }

    // Magnetic perturbations
    Eigen::Vector3d t_mag(0., 0., 0.);
    magnetic_perturbations(t_mag, state, et);

    Eigen::Vector3d a_total = a_kepler + a_geopot + a_aero + a_third_body + a_srp;
    Eigen::Vector3d torques = t_aero + t_gg + t_srp + t_mag;

    state_derivative.set_position(state.velocity());
    state_derivative.set_velocity(a_total);

    // Attitude kinematics
    Quaternion q = state.orientation();
    q.normalize();
    Eigen::Vector3d ang_velocity = state.ang_velocity();
    Quaternion dq(0.5 * q.xi() * ang_velocity);
    state_derivative.set_orientation(dq);

    // Attitude dynamics
    Eigen::Vector3d dw = m_spacecraft.inertia_matrix().inverse() * (torques - ang_velocity.cross(m_spacecraft.inertia_matrix() * ang_velocity));
    state_derivative.set_ang_velocity(dw);

}

void EqMotion::geopotential(Eigen::Vector3d &acc, Eigen::Vector3d &torque, const State &state, SpiceDouble et) const {

    // Returns the index for degree n and order m.
    auto index = [](int n, int m) { return n * (n + 1) / 2 + m; };

    Eigen::MatrixXd cs_coeff = m_env_model.cs_coeffs();
    auto C = [&](int m, int n) { return cs_coeff(index(m, n), 0); };
    auto S = [&](int m, int n) { return cs_coeff(index(m, n), 1); };

    Eigen::MatrixXd vw_coeff = m_env_model.geopotential_harmonic_coeff(state, et);
    auto V = [&](int m, int n) { return vw_coeff(index(m, n), 0); };
    auto W = [&](int m, int n) { return vw_coeff(index(m, n), 1); };

    // Compute the acceleration
    SpiceDouble a_geopot[3] = {0.0, 0.0, 0.0};

    for (int n = 0; n < m_env_model.gp_degree() + 1; ++n) {
        for (int m = 0; m < n + 1; ++m) {
            double scale = 0;
            if (m == 0) {
                scale = sqrt((2 * n + 1));
                a_geopot[0] += scale * constants::MU_EARTH_EGM08 / (constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08) * (-C(n, 0) * V(n + 1, 1));
                a_geopot[1] += scale * constants::MU_EARTH_EGM08 / (constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08) * (-C(n, 0) * W(n + 1, 1));
            } else {
                scale = sqrt(2 * (2 * n + 1) * boost::math::factorial<double>(n - m) /
                             boost::math::factorial<double>(n + m));
                a_geopot[0] += scale * constants::MU_EARTH_EGM08 / (constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08) * 0.5 *
                               ((-C(n, m) * V(n + 1, m + 1) - S(n, m) * W(n + 1, m + 1)) +
                                boost::math::factorial<double>(n - m + 2) / boost::math::factorial<double>(n - m) *
                                (C(n, m) * V(n + 1, m - 1) + S(n, m) * W(n + 1, m - 1)));
                a_geopot[1] += scale * constants::MU_EARTH_EGM08 / (constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08) * 0.5 *
                               ((-C(n, m) * W(n + 1, m + 1) + S(n, m) * V(n + 1, m + 1)) +
                                boost::math::factorial<double>(n - m + 2) / boost::math::factorial<double>(n - m) *
                                (-C(n, m) * W(n + 1, m - 1) + S(n, m) * V(n + 1, m - 1)));
            }
            a_geopot[2] += scale * constants::MU_EARTH_EGM08 / (constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08) *
                           ((n - m + 1) * (-C(n, m) * V(n + 1, m) - S(n, m) * W(n + 1, m)));
        }
    }

    // Transform acceleration from ITRF93 (rotating) to J2000 (inertial)
    SpiceDouble rot[3][3];
    pxform_c("ITRF93", "J2000", et, rot);
    SpiceDouble a_geopot_eci_c[3];
    mxv_c(rot, a_geopot, a_geopot_eci_c);

    // Convert from SpiceDouble to Eigen::Vector3d
    Eigen::Vector3d a_geopot_eci(a_geopot_eci_c[0], a_geopot_eci_c[1], a_geopot_eci_c[2]);
    acc = a_geopot_eci;

    // Nadir-pointing vector in the body frame
    Eigen::Vector3d nadir = state.orientation().attitude_matrix().inverse() * state.position() / state.position().norm();
    torque = 3.0 * constants::MU_EARTH / pow(state.position().norm(), 3) * nadir.cross(m_spacecraft.inertia_matrix()*nadir);
}

void EqMotion::drag(Eigen::Vector3d& acc, Eigen::Vector3d& torque, const State& state, const SpiceDouble et, const double elapsed_time) const {

    double density = m_env_model.atm_density(state, elapsed_time, et);

    // Compute the relative velocity
    Eigen::Vector3d earth_ang_vel(0.0, 0.0, constants::EARTH_ANGULAR_VELOCITY);
    Eigen::Vector3d v_rel = state.velocity() - earth_ang_vel.cross(state.position());

    // Force and torque
    Eigen::Vector3d v_rel_bff = state.orientation().attitude_matrix() * v_rel;
    for (int i = 0; i < m_spacecraft.face_normals().size(); ++i) {
        double cos_theta = m_spacecraft.face_normals()[i].dot(v_rel_bff) / (v_rel.norm()*v_rel.norm());
        Eigen::Vector3d force = - 0.5 * density * m_spacecraft.drag_coefficient() * v_rel.norm() * v_rel_bff * m_spacecraft.face_areas()[i] * max(cos_theta, 0.);
        acc += force / m_spacecraft.mass();
        torque += m_spacecraft.face_cop_positions()[i].cross(force);
    }
    acc = state.orientation().attitude_matrix().inverse() * acc;
}

void EqMotion::solar_radiation_pressure(Eigen::Vector3d& acc, Eigen::Vector3d& torque, const State& state, const SpiceDouble et) const {
    Eigen::Vector3d r_sun2sc = m_env_model.sun_spacecraft_vector(state, et);
    double nu = m_env_model.in_shadow(state, et);
    double p_sr = (constants::SOLAR_PRESSURE * 1E-3) * constants::AU / r_sun2sc.norm(); // kg / (km s^2)

    Eigen::Vector3d r_sc2sun_bff = state.orientation().attitude_matrix() * (-r_sun2sc) / r_sun2sc.norm();
    for (int i = 0; i < m_spacecraft.face_normals().size(); ++i) {
        double cos_theta = m_spacecraft.face_normals()[i].dot(r_sc2sun_bff);
        Eigen::Vector3d force = - nu * p_sr * m_spacecraft.face_areas()[i] * (2.0 * (m_spacecraft.diffuse_reflection_coeff()[i] / 3.0 + m_spacecraft.specular_reflection_coeff()[i]*cos_theta)*m_spacecraft.face_normals()[i] + (1.0 - m_spacecraft.specular_reflection_coeff()[i])*r_sc2sun_bff) * max(cos_theta, 0.);
        acc += force / m_spacecraft.mass();
        torque += m_spacecraft.face_cop_positions()[i].cross(force);
    }
    acc = state.orientation().attitude_matrix().inverse() * acc;
}

void EqMotion::magnetic_perturbations(Eigen::Vector3d& torque, const State& state, SpiceDouble et) const {
    Eigen::Vector3d B_eci = m_env_model.magnetic_field(state, et);
    Eigen::Vector3d B_bff = state.orientation().attitude_matrix() * B_eci;

    // Test internal magnetic dipole
    double area = M_PI * 5E-5*5E-5; // km2
    double current = 0.1; // A
    double n_turns = 5;
    Eigen::Vector3d dipole = n_turns * current * area * Eigen::Vector3d(1., 0., 0.);
    torque = dipole.cross(B_bff);
}

void EqMotion::third_body_effect(Eigen::Vector3d& acc, const State& state, const SpiceDouble et) const {
    std::vector<BodyContainer*> bodies = m_env_model.third_body();
    for (auto body: bodies) {
        Eigen::Vector3d r_cb2tb = m_env_model.body_vector(body, et); // central body to third body
        Eigen::Vector3d r_cb2sc = state.position(); // central body to spacecraft
        Eigen::Vector3d r_sc2tb = r_cb2tb - r_cb2sc; // spacecraft to third body

        double numerator = r_cb2sc.norm() * r_cb2sc.norm() + r_cb2tb.norm() * r_cb2tb.norm() -
                           ((r_cb2tb - r_cb2sc).norm()) * ((r_cb2tb - r_cb2sc).norm());
        double denominator = 2 * r_cb2sc.norm() * r_cb2tb.norm();
        double cos_zeta = numerator / denominator;
        double h = r_cb2sc.norm() / r_cb2tb.norm();
        double B = 0.0;
        double err = 1.0;
        int i = 1;

        while (std::abs(err) > 1e-8 && i < 100) {
            double B_new = B + boost::math::legendre_p(i, cos_zeta) * pow(h, i);
            err = B_new - B;
            B = B_new;
            ++i;
        }

        double beta = 3.0 * B + 3.0 * B * B + B * B * B;
        acc += - body->mu() / pow(r_cb2tb.norm(), 3) * (r_cb2sc - beta * (r_cb2tb - r_cb2sc));
    }
}