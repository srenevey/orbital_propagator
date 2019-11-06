//
// Created by Sylvain Renevey on 8/17/18.
//

#include "perturbations.h"

namespace perturbations {

    Eigen::Vector3d drag(const Spacecraft &sc, const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, double elapsed_time, const EnvironmentModel &env_model) {

        // Get atmospheric parameters
        double density = 0.0;
        SpiceDouble inertial_position[3];
        //SpiceDouble inertial_velocity[3];
        for (int i = 0; i < 3; ++i) {
            inertial_position[i] = state[i];
            //inertial_velocity[i + 3] = state[i + 3];
        }
        env_model.getAtm_parameters(inertial_position, density, elapsed_time, et);
        density = density * 1E9; // transform from kg/m^3 to kg/km^3

        // Compute the relative velocity
        Eigen::Vector3d earth_ang_vel(0.0, 0.0, constants::kEarthAngularVelocity);
        Eigen::Vector3d r = state.head(3);
        Eigen::Vector3d v = state.tail(3);        
        Eigen::Vector3d v_rel = v - earth_ang_vel.cross(r);

        double cd = sc.getDragCoefficient();
        double reference_area = sc.getReferenceArea() * 1E-6; // transform from m^2 to km^2
        double mass = sc.getMass();
        Eigen::Vector3d a_drag = - 0.5 * cd * reference_area / mass * density * v_rel.norm() * v_rel;

        return a_drag;
    }



    Eigen::Vector3d geopotential(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const EnvironmentModel &env_model) {

        // Transform the position from J2000 (inertial) to ITRF93 (rotating).
		SpiceDouble rot[3][3];
		pxform_c("J2000", "ITRF93", et, rot);
		
		SpiceDouble ecef_position[3];
		SpiceDouble s_eci_position[3];
		for (int i = 0; i < 3; ++i)
			s_eci_position[i] = state(i);
		mxv_c(rot, s_eci_position, ecef_position);
		

        // Returns the index for degree n and order m.
        auto index = [](int n, int m) { return n * (n + 1) / 2 + m; };


        // Lambda expressions used as helpers to get C and S coefficients
        Eigen::MatrixXd CS_coeffs = env_model.getCS_coeffs();
        auto C = [&](int m, int n) { return CS_coeffs(index(m, n), 0); };
        auto S = [&](int m, int n) { return CS_coeffs(index(m, n), 1); };


        // Compute V and W coefficients
        int degree = env_model.getGp_degree();
        int num_rows = (degree + 2) * (degree + 3) / 2;
        double r = vnorm_c(ecef_position);
        Eigen::MatrixXd VW_coeffs = Eigen::MatrixXd::Zero(num_rows, 2);
        VW_coeffs(0, 0) = constants::kEarthRadiusEgm08 / r;
        VW_coeffs(0, 1) = 0;

        // Zonal terms
        VW_coeffs(index(1, 0), 0) = ecef_position[2] * constants::kEarthRadiusEgm08 / (r * r) * VW_coeffs(index(0, 0), 0);
        VW_coeffs(index(1, 0), 1) = 0;
        for (int n = 2; n < degree + 2; ++n) {
            VW_coeffs(index(n, 0), 0) =
                    (2 * n - 1) * ecef_position[2] * constants::kEarthRadiusEgm08 / (n * r * r) * VW_coeffs(index(n - 1, 0), 0) -
                    (n - 1) * constants::kEarthRadiusEgm08 * constants::kEarthRadiusEgm08 / (n * r * r) * VW_coeffs(index(n - 2, 0), 0);
            VW_coeffs(index(n, 0), 1) = 0;
        }

        // Tesseral terms
        for (int n = 1; n < degree + 2; ++n) {
            VW_coeffs(index(n, n), 0) = (2 * n - 1) *
                                        (ecef_position[0] * constants::kEarthRadiusEgm08 / (r * r) * VW_coeffs(index(n - 1, n - 1), 0) -
                                            ecef_position[1] * constants::kEarthRadiusEgm08 / (r * r) * VW_coeffs(index(n - 1, n - 1), 1));
            VW_coeffs(index(n, n), 1) = (2 * n - 1) *
                                        (ecef_position[0] * constants::kEarthRadiusEgm08 / (r * r) * VW_coeffs(index(n - 1, n - 1), 1) +
                                            ecef_position[1] * constants::kEarthRadiusEgm08 / (r * r) * VW_coeffs(index(n - 1, n - 1), 0));
        }

        // Remaining terms
        for (int n = 1; n < degree + 2; ++n) {
            for (int m = 1; m < n; ++m) {
                if (n == m + 1) {
                    VW_coeffs(index(n, m), 0) = (2 * n - 1) * ecef_position[2] * constants::kEarthRadiusEgm08 / ((n - m) * r * r) *
                                                VW_coeffs(index(n - 1, m), 0);
                    VW_coeffs(index(n, m), 1) = (2 * n - 1) * ecef_position[2] * constants::kEarthRadiusEgm08 / ((n - m) * r * r) *
                                                VW_coeffs(index(n - 1, m), 1);
                } else {
                    VW_coeffs(index(n, m), 0) = (2 * n - 1) * ecef_position[2] * constants::kEarthRadiusEgm08 / ((n - m) * r * r) *
                                                VW_coeffs(index(n - 1, m), 0) -
                                                (n + m - 1) * constants::kEarthRadiusEgm08 * constants::kEarthRadiusEgm08 / ((n - m) * r * r) *
                                                VW_coeffs(index(n - 2, m), 0);
                    VW_coeffs(index(n, m), 1) = (2 * n - 1) * ecef_position[2] * constants::kEarthRadiusEgm08 / ((n - m) * r * r) *
                                                VW_coeffs(index(n - 1, m), 1) -
                                                (n + m - 1) * constants::kEarthRadiusEgm08 * constants::kEarthRadiusEgm08 / ((n - m) * r * r) *
                                                VW_coeffs(index(n - 2, m), 1);
                }
            }
        }

        // Lambda expressions used as helpers to get V and W coefficients
        auto V = [&](int m, int n) { return VW_coeffs(index(m, n), 0); };
        auto W = [&](int m, int n) { return VW_coeffs(index(m, n), 1); };


        // Compute the acceleration
        SpiceDouble a_geopot[3];
        for (int i = 0; i < 3; ++i)
            a_geopot[i] = 0.0;

        for (int n = 0; n < degree + 1; ++n) {
            for (int m = 0; m < n + 1; ++m) {
                double scale = 0;
                if (m == 0) {
                    scale = sqrt((2 * n + 1));
                    a_geopot[0] += scale * constants::kMuEarthEgm08 / (constants::kEarthRadiusEgm08 * constants::kEarthRadiusEgm08) * (-C(n, 0) * V(n + 1, 1));
                    a_geopot[1] += scale * constants::kMuEarthEgm08 / (constants::kEarthRadiusEgm08 * constants::kEarthRadiusEgm08) * (-C(n, 0) * W(n + 1, 1));
                } else {
                    scale = sqrt(2 * (2 * n + 1) * boost::math::factorial<double>(n - m) /
                                    boost::math::factorial<double>(n + m));
                    a_geopot[0] += scale * constants::kMuEarthEgm08 / (constants::kEarthRadiusEgm08 * constants::kEarthRadiusEgm08) * 0.5 *
                                    ((-C(n, m) * V(n + 1, m + 1) - S(n, m) * W(n + 1, m + 1)) +
                                    boost::math::factorial<double>(n - m + 2) / boost::math::factorial<double>(n - m) *
                                    (C(n, m) * V(n + 1, m - 1) + S(n, m) * W(n + 1, m - 1)));
                    a_geopot[1] += scale * constants::kMuEarthEgm08 / (constants::kEarthRadiusEgm08 * constants::kEarthRadiusEgm08) * 0.5 *
                                    ((-C(n, m) * W(n + 1, m + 1) + S(n, m) * V(n + 1, m + 1)) +
                                    boost::math::factorial<double>(n - m + 2) / boost::math::factorial<double>(n - m) *
                                    (-C(n, m) * W(n + 1, m - 1) + S(n, m) * V(n + 1, m - 1)));
                }
                a_geopot[2] += scale * constants::kMuEarthEgm08 / (constants::kEarthRadiusEgm08 * constants::kEarthRadiusEgm08) *
                                ((n - m + 1) * (-C(n, m) * V(n + 1, m) - S(n, m) * W(n + 1, m)));
            }
        }

        // Transform acceleration from ITRF93 (rotating) to J2000 (inertial)
		pxform_c("ITRF93", "J2000", et, rot);
        SpiceDouble a_geopot_eci_c[3];
        mxv_c(rot, a_geopot, a_geopot_eci_c);

        // Convert from SpiceDouble to Eigen::Vector3d
        Eigen::Vector3d a_geopot_eci(a_geopot_eci_c[0], a_geopot_eci_c[1], a_geopot_eci_c[2]);

        return a_geopot_eci;
    }






    Eigen::Vector3d nbody(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const EnvironmentModel &env_model) {

        std::string central_body = env_model.getCentral_body();
        std::array<int,10> third_body_flags = env_model.getThird_body_flags();

        // Compute perturbations due to each body
        Eigen::Vector3d a_sun(0.0, 0.0, 0.0);
        if ( third_body_flags[0] && central_body != "Sun") {
            a_sun = perturbation_from_body(state, "Sun", central_body, et);
        }

        Eigen::Vector3d a_mercury(0.0, 0.0, 0.0);
        if (third_body_flags[1] && central_body != "Mercury") {
            a_mercury = perturbation_from_body(state, "Mercury", central_body, et);
        }

        Eigen::Vector3d a_venus(0.0, 0.0, 0.0);
        if (third_body_flags[2] && central_body != "Venus") {
            a_venus = perturbation_from_body(state, "Venus", central_body, et);
        }

        Eigen::Vector3d a_earth(0.0, 0.0, 0.0);
        if (third_body_flags[3] && central_body != "Earth") {
            a_earth = perturbation_from_body(state, "Earth", central_body, et);
        }

        Eigen::Vector3d a_moon(0.0, 0.0, 0.0);
        if (third_body_flags[4] && central_body != "Moon") {
            a_moon = perturbation_from_body(state, "Moon", central_body, et);
        }

        Eigen::Vector3d a_mars(0.0, 0.0, 0.0);
        if (third_body_flags[5] && central_body != "Mars") {
            a_mars = perturbation_from_body(state, "Mars", central_body, et);
        }

        Eigen::Vector3d a_jupiter(0.0, 0.0, 0.0);
        if (third_body_flags[6] && central_body != "Jupiter") {
            a_jupiter = perturbation_from_body(state, "Jupiter", central_body, et);
        }

        Eigen::Vector3d a_saturn(0.0, 0.0, 0.0);
        if (third_body_flags[7] && central_body != "Saturn") {
            a_saturn = perturbation_from_body(state, "Saturn", central_body, et);
        }

        Eigen::Vector3d a_uranus(0.0, 0.0, 0.0);
        if (third_body_flags[8] && central_body != "Uranus") {
            a_uranus = perturbation_from_body(state, "Uranus", central_body, et);
        }

        Eigen::Vector3d a_neptune(0.0, 0.0, 0.0);
        if (third_body_flags[9] && central_body != "Neptune") {
            a_neptune = perturbation_from_body(state, "Neptune", central_body, et);
        }

        Eigen::Vector3d a_nbody(0.0, 0.0, 0.0);
        a_nbody = a_sun + a_mercury + a_venus + a_earth + a_moon + a_mars + a_jupiter + a_saturn + a_uranus + a_neptune;
        return a_nbody;
    }



    Eigen::Vector3d perturbation_from_body(const Eigen::Ref<const Eigen::VectorXd>& state, const std::string &third_body, const std::string &central_body, SpiceDouble et) {

		// Change names to map SPICE's nomenclature
		std::string cb(central_body);
		std::string tb(third_body);
		if (central_body != std::string("Earth") && central_body != std::string("Sun") && central_body != std::string("Moon"))
			cb += " barycenter";
		if (third_body != std::string("Earth") && third_body != std::string("Sun") && third_body != std::string("Moon"))
			tb += " barycenter";

		// Get the state of the third body relative to the central body
		SpiceDouble third_body_state[6];
		SpiceDouble lt = 0;
		if (central_body == std::string("Earth"))
			spkezr_c(tb.c_str(), et, "J2000", "NONE", cb.c_str(), third_body_state, &lt);
		else
			spkezr_c(tb.c_str(), et, "ECLIPJ2000", "NONE", cb.c_str(), third_body_state, &lt);


		Eigen::Vector3d r_cb2tb = {third_body_state[0], third_body_state[1], third_body_state[2]}; // central body to third body
		Eigen::Vector3d r_cb2sc = state.head(3);            // central body to spacecraft
        Eigen::Vector3d r_sc2tb = r_cb2tb - r_cb2sc;        // spacecraft to third body


		// Compute perturbation
		double numerator = r_cb2sc.norm() * r_cb2sc.norm() + r_cb2tb.norm() * r_cb2tb.norm() - ((r_cb2tb - r_cb2sc).norm()) * ((r_cb2tb - r_cb2sc).norm());
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
            i++;
        }

        double beta = 3 * B + 3 * B * B + B * B * B;

        // Retrieve the gravitational parameter of the third body
        SpiceDouble mu_third_body = 0;
        SpiceInt n = 0;
        bodvrd_c(tb.c_str(), "GM", 1, &n, &mu_third_body);

        // Compute the acceleration due to the third body
        Eigen::Vector3d a_body(0.0, 0.0, 0.0);
        a_body = - mu_third_body / pow(r_cb2tb.norm(), 3) * (r_cb2sc - beta * (r_cb2tb - r_cb2sc));

        return a_body;
    }

    Eigen::Vector3d solar_radiation_pressure(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const EnvironmentModel &env_model, const Spacecraft &sc) {
        Eigen::Vector3d a_srp(0.0, 0.0, 0.0);
        std::string central_body = env_model.getCentral_body();

        // Get vector from spacecraft to the Sun
        Eigen::Vector3d r_sun2sc(0.0, 0.0, 0.0);
        if (central_body == std::string("Sun")) {
            r_sun2sc = state.head(3);
        }
        else {
            std::string cb(central_body);
            if (central_body != std::string("Earth") && central_body != std::string("Sun") && central_body != std::string("Moon"))
                cb += " barycenter";
            SpiceDouble central_body_state[6];
            SpiceDouble lt = 0.0;
            if (central_body == std::string("Earth"))
                spkezr_c(cb.c_str(), et, "J2000", "NONE", "Sun", central_body_state, &lt);
            else
                spkezr_c(cb.c_str(), et, "ECLIPJ2000", "NONE", "Sun", central_body_state, &lt);

            r_sun2sc = {central_body_state[0] + state(0), central_body_state[1] + state(1), central_body_state[2] + state(2)};
        }

        // Shadow function
        double a = asin(constants::radius("Sun") / r_sun2sc.norm());
        double b = asin(constants::radius(central_body) / state.head(3).norm());
        double c = acos( state.head(3).dot(r_sun2sc) / (state.head(3).norm() * r_sun2sc.norm()));

        double nu = 0.0;
        if (abs(a-b) < c && c < a + b) {
            double x = (c * c + a * a - b * b) / (2.0 * c);
            double y = sqrt(a * a - x * x);
            double A = a * a * acos(x / a) + b * b * acos((c - x) / b) - c * y;
            nu = 1.0 - A / (constants::kPi * a * a);
        } else if (c >= a + b) {
            nu = 0.0;
        } else if (c < b - a) {
            nu = 1.0;
        }



        // Compute acceleration
        double p_sr = (constants::kSolarPressure * 1E-3) * constants::kAstronomicalUnit / r_sun2sc.norm(); // kg / (km s^2)
        double reflectivity = sc.getReflectivity();
        double area = sc.getSunExposedArea() * 1E-6; // transform from m^2 to km^2
        a_srp = - nu * p_sr * (1 + reflectivity) * area / sc.getMass() * r_sun2sc / r_sun2sc.norm();

        return a_srp;
    }

}