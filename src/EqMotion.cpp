//
// Created by Sylvain Renevey on 5/11/18.
//

#include "EqMotion.h"

EqMotion::EqMotion(const Spacecraft &sc, const EnvironmentModel &env_model, double initial_et): m_env_model(env_model), m_initial_et(initial_et), m_spacecraft(sc) {};

void EqMotion::operator()(StateVector &state, StateVector &state_derivative, double t) {

    // Compute the current ephemeris time
    double et = m_initial_et + t;
    state.set_time(et);
    state.normalize_quaternion();

    // Two-body acceleration
    double factor = - m_env_model.central_body().mu() / pow(state.position().norm(), 3);
    Vector3d<Dimension::Acceleration> a_kepler(state.frame(), {factor*state[0], factor*state[1], factor*state[2]});

    // Geopotential
    Vector3d<Dimension::Acceleration> a_geopot(state.frame());
    Vector3d<double> t_gg(ReferenceFrame::BODY);
    if (m_env_model.gp_degree() > 1) {
        geopotential(a_geopot, t_gg, state, et);
    }

    // Third-body effect
    Vector3d<Dimension::Acceleration> a_third_body(state.frame());
    if (!m_env_model.third_bodies().empty()) {
        third_body_effect(a_third_body, state, et);
    }

    // Aerodynamic perturbations
    Vector3d<Dimension::Acceleration> a_aero(state.frame());
    Vector3d<double> t_aero(ReferenceFrame::BODY);
    if (m_env_model.is_drag() && m_spacecraft.mass() != 0.0) {
        drag(a_aero, t_aero, state, et, t);
    }

    // Solar radiation pressure
    Vector3d<Dimension::Acceleration> a_srp(state.frame());
    Vector3d<double> t_srp(ReferenceFrame::BODY);
    if (m_env_model.is_srp_flag()) {
        solar_radiation_pressure(a_srp, t_srp, state, et);
    }

    // Magnetic perturbations
    Vector3d<double> t_mag(ReferenceFrame::BODY);
    if (m_env_model.mag_flag()) {
        magnetic_perturbations(t_mag, state, et);
    }

    Vector3d<Dimension::Acceleration> a_total = a_kepler + a_geopot + a_third_body + a_aero + a_srp;
    Vector3d<double> torques = t_gg + t_aero + t_srp + t_mag;

    state_derivative.set_position_derivative(state.velocity());
    state_derivative.set_velocity_derivative(a_total);


    // Attitude kinematics
    Quaternion q = state.orientation();
    q.normalize();
    Quaternion dq(0.5 * q.xi() * state.ang_velocity());
    state_derivative.set_orientation(dq);

    // Attitude dynamics
    Vector3d<Dimension::AngularAcceleration> ang_acceleration(m_spacecraft.inertia_matrix().inverse() * (torques - state.ang_velocity().cross(m_spacecraft.inertia_matrix() * state.ang_velocity())));
    state_derivative.set_ang_velocity_derivative(ang_acceleration);
}

void EqMotion::geopotential(Vector3d<Dimension::Acceleration> &acc, Vector3d<double> &torque, const StateVector &state, double et) const {
    // Computes the index for degree n and order m.
    auto index = [](int n, int m) { return n * (n + 1) / 2 + m; };

    Eigen::MatrixXd cs_coeff = m_env_model.cs_coeffs(); // TODO: no need to use Eigen here, return native array instead
    auto C = [&](int m, int n) { return cs_coeff(index(m, n), 0); };
    auto S = [&](int m, int n) { return cs_coeff(index(m, n), 1); };

    Eigen::MatrixXd vw_coeff = m_env_model.geopotential_harmonic_coeff(state, et); // TODO: no need to use Eigen here, return native array instead
    auto V = [&](int m, int n) { return vw_coeff(index(m, n), 0); };
    auto W = [&](int m, int n) { return vw_coeff(index(m, n), 1); };

    // Computes the acceleration
    Vector3d<Dimension::Acceleration> a_geopot(ReferenceFrame::ITRF93);

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

    // Rotates back to the original frame. Note: a_geopot already accounts for Coriolis, centripetal, etc. accelerations.
    acc = transformations::rotate_to_frame(a_geopot, state.frame(), et);

    // Nadir-pointing vector in the body frame
    Vector3d<Dimension::Distance> nadir = state.position() / state.position().norm();
    Vector3d<Dimension::Distance> nadir_bff = transformations::rotate_to_frame(nadir, ReferenceFrame::BODY, et, state.orientation());
    torque = 3.0 * constants::MU_EARTH / pow(state.position().norm(), 3) * nadir_bff.cross(m_spacecraft.inertia_matrix()*nadir_bff);
}

void EqMotion::drag(Vector3d<Dimension::Acceleration> &acc, Vector3d<double> &torque, const StateVector &state, double et, double elapsed_time) const {

    // Retrieves atmospheric density
    double density = m_env_model.atm_density(state, elapsed_time, et);

    // Computes the velocity relative to Earth's atmosphere
    Vector3d<Dimension::AngularVelocity> earth_ang_vel(ReferenceFrame::J2000, {0., 0., constants::EARTH_ANGULAR_VELOCITY});
    Vector3d<Dimension::Velocity> earth_velocity(earth_ang_vel.cross(state.position()));
    Vector3d<Dimension::Velocity> v_rel = state.velocity() - earth_velocity;
    Vector3d<Dimension::Velocity> v_rel_bff = transformations::rotate_to_frame(v_rel, ReferenceFrame::BODY, et, state.orientation());

    // Force and torque
    Vector3d<Dimension::Acceleration> acc_bff(ReferenceFrame::BODY);
    for (int i = 0; i < m_spacecraft.face_normals().size(); ++i) {
        double cos_theta = m_spacecraft.face_normals()[i].dot(v_rel_bff) / (v_rel.norm()*v_rel.norm());
        Vector3d<double> force_bff(-0.5 * density * m_spacecraft.drag_coefficient() * v_rel.norm() * v_rel_bff * m_spacecraft.face_areas()[i] * max(cos_theta, 0.));
        acc_bff += force_bff/m_spacecraft.mass();
        torque += m_spacecraft.face_cop_positions()[i].cross(force_bff);
    }
    acc = transformations::rotate_to_frame(acc_bff, state.frame(), et, state.orientation());
}

void EqMotion::solar_radiation_pressure(Vector3d<Dimension::Acceleration>& acc, Vector3d<double>& torque, const StateVector& state, const double et) const {
    Vector3d<Dimension::Distance> r_sun2sc = m_env_model.sun_spacecraft_vector(state, et);
    double nu = m_env_model.in_shadow(state, et);
    double p_sr = (constants::SOLAR_PRESSURE * 1E-3) * constants::AU / r_sun2sc.norm(); // kg / (km s^2)

    Vector3d<Dimension::Distance> r_sc2sun_bff = transformations::rotate_to_frame((-r_sun2sc) / r_sun2sc.norm(), ReferenceFrame::BODY, et, state.orientation());
    Vector3d<Dimension::Acceleration> acc_bff(ReferenceFrame::BODY);
    for (int i = 0; i < m_spacecraft.face_normals().size(); ++i) {
        double cos_theta = m_spacecraft.face_normals()[i].dot(r_sc2sun_bff);
        Vector3d<double> force_bff(-nu * p_sr * m_spacecraft.face_areas()[i] * (2.0 * (m_spacecraft.diffuse_reflection_coeff()[i] / 3.0 + m_spacecraft.specular_reflection_coeff()[i]*cos_theta)*m_spacecraft.face_normals()[i] + (1.0 - m_spacecraft.specular_reflection_coeff()[i])*r_sc2sun_bff) * max(cos_theta, 0.));
        acc_bff += force_bff/m_spacecraft.mass();
        torque += m_spacecraft.face_cop_positions()[i].cross(force_bff);
    }
    acc = transformations::rotate_to_frame(acc_bff, state.frame(), et);
}

void EqMotion::magnetic_perturbations(Vector3d<double>& torque, const StateVector& state, double et) const {

    Vector3d<double> B = m_env_model.magnetic_field(state, et);
    Vector3d<double> B_bff = transformations::rotate_to_frame(B, ReferenceFrame::BODY, et, state.orientation());

    // Test internal magnetic dipole
    double area = M_PI * 5E-5*5E-5; // km2
    double current = 0.1; // A
    double n_turns = 5;
    Vector3d<double> dipole = n_turns * current * area * Vector3d<double>(ReferenceFrame::BODY, {1., 0., 0.});
    torque = dipole.cross(B_bff);
}

void EqMotion::third_body_effect(Vector3d<Dimension::Acceleration>& acc, const StateVector& state, const double et) const {
    auto bodies = m_env_model.third_bodies();

    for (auto body: bodies) {
        Vector3d<Dimension::Distance> r_cb2tb = m_env_model.body_vector(body, et, state.frame()); // central body to third body
        Vector3d<Dimension::Distance> r_cb2sc = state.position(); // central body to spacecraft
        Vector3d<Dimension::Distance> r_sc2tb = r_cb2tb - r_cb2sc; // spacecraft to third body

        double numerator = r_cb2sc.norm() * r_cb2sc.norm() + r_cb2tb.norm() * r_cb2tb.norm() -
                           ((r_cb2tb - r_cb2sc).norm()) * ((r_cb2tb - r_cb2sc).norm());
        double denominator = 2 * r_cb2sc.norm() * r_cb2tb.norm();
        double cos_zeta = numerator / denominator;
        double h = r_cb2sc.norm() / r_cb2tb.norm();
        double B = 0.0;
        double err = 1.0;
        int i = 1;

        while (std::abs(err) > 1E-8 && i < 100) {
            double B_new = B + boost::math::legendre_p(i, cos_zeta) * pow(h, i);
            err = B_new - B;
            B = B_new;
            ++i;
        }

        double beta = 3.0 * B + 3.0 * B * B + B * B * B;
        acc += -body.mu()/pow(r_cb2tb.norm(), 3) * (r_cb2sc - beta * (r_cb2tb - r_cb2sc));
    }
}