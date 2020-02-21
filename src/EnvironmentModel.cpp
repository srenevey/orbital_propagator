//
// Created by Sylvain Renevey on 1/26/18.
//



#include "EnvironmentModel.h"

EnvironmentModel::EnvironmentModel():
        mu_(0.0),
        gp_degree_(0),
        third_body_flags_({0,0,0,0,0,0,0,0,0,0}),
        third_body_flag_(false),
        drag_flag_(false),
        earth_gram_atm_model_(nullptr),
        srp_flag_(false)
{}


EnvironmentModel::EnvironmentModel(std::string earthgram_path, std::string central_body, int gp_degree, std::string geopot_model_path, std::array<int, 10> third_body_flags, bool drag_flag, std::string atm_model, bool srp_flag, std::string epoch):
        central_body_(central_body),
        mu_(constants::mu(central_body)),
        gp_degree_(gp_degree),
        third_body_flags_(third_body_flags),
        drag_flag_(drag_flag),
        atm_model_(atm_model),
        srp_flag_(srp_flag),
        earth_gram_atm_model_(nullptr)
{
    // Load geopotential coefficients from model file
    if (gp_degree > 1 && !geopot_model_path.empty())
        CS_coeffs_ = load_coefficients(gp_degree, geopot_model_path);

    // Check if some third body flags are on
    third_body_flag_ = (std::accumulate(third_body_flags.begin(), third_body_flags.end(), 0) != 0);

    // If the drag flag is on, initialize EarthGRAM2016 or exponential atmospheric model
    if (drag_flag && central_body == "Earth") {
        if (atm_model == "EarthGRAM") {
            if (!epoch.empty()) {
                earth_gram_atm_model_ = new Atm1;
                earth_gram_atm_model_->initdata(earthgram_path, "NameRef.txt", epoch);
            } else {
                std::cout << "The epoch must be passed in as argument." << std::endl;
            }
        } else if (atm_model == "exp") {
            init_exp_model();
        } else {
            std::cout << "Wrong atmosphere model name. Drag effects were turned off." << std::endl;
            drag_flag_ = false;
        }
    } else if (drag_flag && central_body != "Earth") {
        std::cout << "The drag effects are currently only supported for the Earth." << std::endl;
        drag_flag_ = false;
    }
}


EnvironmentModel::~EnvironmentModel() {
    delete earth_gram_atm_model_;
}


double EnvironmentModel::mu() const {
    return mu_;
}


std::array<int, 10> EnvironmentModel::third_body_flags() const {
    return third_body_flags_;
}


bool EnvironmentModel::is_third_body_flag() const {
    return third_body_flag_;
}


int EnvironmentModel::gp_degree() const {
    return gp_degree_;
}


const Eigen::MatrixXd &EnvironmentModel::cs_coeffs() const {
    return CS_coeffs_;
}


std::string EnvironmentModel::central_body() const {
    return central_body_;
}


bool EnvironmentModel::is_drag_flag() const {
    return drag_flag_;
}

bool EnvironmentModel::is_srp_flag() const {
    return srp_flag_;
}

Eigen::MatrixXd EnvironmentModel::load_coefficients(int degree, std::string model_file_name) {
    // Open file containing the geopotential coefficients
    ifstream data_file;
    data_file.open(model_file_name);

    // Initialize matrix containing the coefficients
    int num_entries = (degree + 1) * (degree + 2) / 2;

    Eigen::MatrixXd coefficients = Eigen::MatrixXd::Zero(num_entries, 2);
    coefficients.row(0) << 0, 0;
    coefficients.row(1) << 0, 0;
    coefficients.row(2) << 0, 0;

    if (data_file.is_open()) {

        // Read required lines
        for (int i = 3; i < num_entries; ++i) {
            string line;
            getline(data_file, line);

            // Parse the line, convert to double and store into matrix
            unsigned long n = line.find('D');
            while (n != string::npos) {
                line.replace(n, 1, "e");
                n = line.find('D');
            }
            std::istringstream iss(line);
            std::vector<std::string> results(std::istream_iterator<std::string>{iss},std::istream_iterator<std::string>());
            for (int j = 0; j < 2; ++j) {
                coefficients(i, j) = stod(results[j + 2]);
            }
        }

        data_file.close();
    } else {
    	std::cerr << "Error opening the file containing the geopotential coefficients." << std::endl;
    }

    return coefficients;
}

// void EnvironmentModel::getAtm_parameters(const Eigen::Ref<const Eigen::VectorXd>& state, double& density, double elapsed_time, SpiceDouble et) const {
void EnvironmentModel::get_atm_parameters(const SpiceDouble *inertial_position, double& density, double elapsed_time, SpiceDouble et) const {

    if (atm_model_ == "EarthGRAM") {

        // Transform the state from J2000 (inertial) to ITRF93 (rotating)
        SpiceDouble rot[3][3];
        pxform_c("J2000", "ITRF93", et, rot);
		
		SpiceDouble ecef_position[3];
		mxv_c(rot, inertial_position, ecef_position);


		// Compute Earth's flattening
		SpiceInt dim = 0;
		SpiceDouble radii[3];
		bodvrd_c("EARTH", "RADII", 3, &dim, radii);
		SpiceDouble f = (radii[0] - radii[2]) / radii[1];

		// Transform from rectangular to geodetic coordinates
		SpiceDouble lon = 0.0, lat = 0.0, alt = 0.0;
		recgeo_c(ecef_position, radii[0], f, &lon, &lat, &alt);

        double pm1, dm1, tm1, um1, vm1, wm1, pp1, dp1, tp1, up1, vp1, wp1, ps1, ds1, ts1,
                us1, vs1, ws1, psmall, dsmall, tsmall, usmall, vsmall, wsmall, sos, sosp;
        
        // Convert from radians to degrees
        lat = lat * dpr_c();
        lon = lon * dpr_c();

        int iupdate = 1;
        int initonce = 0;

        if (elapsed_time == 0.0)
            initonce = 1;

        earth_gram_atm_model_->traj(alt, lat, lon, elapsed_time, iupdate, initonce, &dm1, &pm1, &tm1, &um1, &vm1, &wm1,
                                     &dp1, &pp1, &tp1, &up1, &vp1, &wp1, &ds1, &ps1, &ts1, &us1, &vs1, &ws1, &dsmall,
                                     &psmall, &tsmall, &usmall, &vsmall, &wsmall, &sos, &sosp);
        density = dm1;

    } else if (atm_model_ == "exp") {

        // Retrieve Earth's radius
        SpiceInt dim = 0;
		SpiceDouble radius[3];
		bodvrd_c("EARTH", "RADII", 3, &dim, radius);

		// Compute the altitude above the ellipsoid. There is no atmosphere above 1100 km.
        double h_ell = vnorm_c(inertial_position) - radius[0];
        if (h_ell <= 1100.0) {
            for (int i = 0; i < 28; ++i) {
                if (exp_atm_model_[i][0] <= h_ell && exp_atm_model_[i][1] > h_ell) {
                    double h0 = exp_atm_model_[i][0];
                    double rho0 = exp_atm_model_[i][2];
                    double H = exp_atm_model_[i][3];
                    density = rho0 * exp(- (h_ell - h0) / H);
                }
            }
        } else {
            density = 0.0;
        }
    }

}


void EnvironmentModel::print_parameters() const {
    // Print model parameters
    std::cout << "------------------" << std::endl;
    std::cout << "Environment model:" << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << "Central body: " << central_body_ << std::endl;
    std::cout << "Geopotential degree: " << gp_degree_ << std::endl;
    std::cout << "Third body flag: " << third_body_flag_ << std::endl;
    std::cout << "Third body:" << std::endl;
    std::cout << "\t Sun: " << third_body_flags_[0] << std::endl;
    std::cout << "\t Mercury: " << third_body_flags_[1] << std::endl;
    std::cout << "\t Venus: " << third_body_flags_[2] << std::endl;
    std::cout << "\t Earth: " << third_body_flags_[3] << std::endl;
    std::cout << "\t Moon: " << third_body_flags_[4] << std::endl;
    std::cout << "\t Mars: " << third_body_flags_[5] << std::endl;
    std::cout << "\t Jupiter: " << third_body_flags_[6] << std::endl;
    std::cout << "\t Saturn: " << third_body_flags_[7] << std::endl;
    std::cout << "\t Uranus: " << third_body_flags_[8] << std::endl;
    std::cout << "\t Neptune: " << third_body_flags_[9] << std::endl;
    std::cout << "Drag flag: " << drag_flag_ << std::endl;
    std::cout << "Atm model: " << atm_model_ << std::endl;
}

Eigen::Vector3d EnvironmentModel::drag(const Spacecraft &sc, const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, double elapsed_time) const {

    // Get atmospheric parameters
    double density = 0.0;
    SpiceDouble inertial_position[3];
    //SpiceDouble inertial_velocity[3];
    for (int i = 0; i < 3; ++i) {
        inertial_position[i] = state[i];
        //inertial_velocity[i + 3] = state[i + 3];
    }
    get_atm_parameters(inertial_position, density, elapsed_time, et);
    density = density * 1E9; // transform from kg/m^3 to kg/km^3

    // Compute the relative velocity
    Eigen::Vector3d earth_ang_vel(0.0, 0.0, constants::kEarthAngularVelocity);
    Eigen::Vector3d r = state.head(3);
    Eigen::Vector3d v = state.tail(3);
    Eigen::Vector3d v_rel = v - earth_ang_vel.cross(r);

    double cd = sc.drag_coefficient();
    double reference_area = sc.reference_area() * 1E-6; // transform from m^2 to km^2
    double mass = sc.mass();
    Eigen::Vector3d a_drag = - 0.5 * cd * reference_area / mass * density * v_rel.norm() * v_rel;

    return a_drag;
}



Eigen::Vector3d EnvironmentModel::geopotential(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et) const {

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
    Eigen::MatrixXd CS_coeffs = cs_coeffs();
    auto C = [&](int m, int n) { return CS_coeffs(index(m, n), 0); };
    auto S = [&](int m, int n) { return CS_coeffs(index(m, n), 1); };


    // Compute V and W coefficients
    int num_rows = (gp_degree_ + 2) * (gp_degree_ + 3) / 2;
    double r = vnorm_c(ecef_position);
    Eigen::MatrixXd VW_coeffs = Eigen::MatrixXd::Zero(num_rows, 2);
    VW_coeffs(0, 0) = constants::kEarthRadiusEgm08 / r;
    VW_coeffs(0, 1) = 0;

    // Zonal terms
    VW_coeffs(index(1, 0), 0) = ecef_position[2] * constants::kEarthRadiusEgm08 / (r * r) * VW_coeffs(index(0, 0), 0);
    VW_coeffs(index(1, 0), 1) = 0;
    for (int n = 2; n < gp_degree_ + 2; ++n) {
        VW_coeffs(index(n, 0), 0) =
                (2 * n - 1) * ecef_position[2] * constants::kEarthRadiusEgm08 / (n * r * r) * VW_coeffs(index(n - 1, 0), 0) -
                (n - 1) * constants::kEarthRadiusEgm08 * constants::kEarthRadiusEgm08 / (n * r * r) * VW_coeffs(index(n - 2, 0), 0);
        VW_coeffs(index(n, 0), 1) = 0;
    }

    // Tesseral terms
    for (int n = 1; n < gp_degree_ + 2; ++n) {
        VW_coeffs(index(n, n), 0) = (2 * n - 1) *
                                    (ecef_position[0] * constants::kEarthRadiusEgm08 / (r * r) * VW_coeffs(index(n - 1, n - 1), 0) -
                                     ecef_position[1] * constants::kEarthRadiusEgm08 / (r * r) * VW_coeffs(index(n - 1, n - 1), 1));
        VW_coeffs(index(n, n), 1) = (2 * n - 1) *
                                    (ecef_position[0] * constants::kEarthRadiusEgm08 / (r * r) * VW_coeffs(index(n - 1, n - 1), 1) +
                                     ecef_position[1] * constants::kEarthRadiusEgm08 / (r * r) * VW_coeffs(index(n - 1, n - 1), 0));
    }

    // Remaining terms
    for (int n = 1; n < gp_degree_ + 2; ++n) {
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

    for (int n = 0; n < gp_degree_ + 1; ++n) {
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


Eigen::Vector3d EnvironmentModel::nbody(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et) const {


    // Compute perturbations due to each body
    Eigen::Vector3d a_sun(0.0, 0.0, 0.0);
    if (third_body_flags_[0] && central_body_ != "Sun") {
        a_sun = perturbation_from_body(state, "Sun", central_body_, et);
    }

    Eigen::Vector3d a_mercury(0.0, 0.0, 0.0);
    if (third_body_flags_[1] && central_body_ != "Mercury") {
        a_mercury = perturbation_from_body(state, "Mercury", central_body_, et);
    }

    Eigen::Vector3d a_venus(0.0, 0.0, 0.0);
    if (third_body_flags_[2] && central_body_ != "Venus") {
        a_venus = perturbation_from_body(state, "Venus", central_body_, et);
    }

    Eigen::Vector3d a_earth(0.0, 0.0, 0.0);
    if (third_body_flags_[3] && central_body_ != "Earth") {
        a_earth = perturbation_from_body(state, "Earth", central_body_, et);
    }

    Eigen::Vector3d a_moon(0.0, 0.0, 0.0);
    if (third_body_flags_[4] && central_body_ != "Moon") {
        a_moon = perturbation_from_body(state, "Moon", central_body_, et);
    }

    Eigen::Vector3d a_mars(0.0, 0.0, 0.0);
    if (third_body_flags_[5] && central_body_ != "Mars") {
        a_mars = perturbation_from_body(state, "Mars", central_body_, et);
    }

    Eigen::Vector3d a_jupiter(0.0, 0.0, 0.0);
    if (third_body_flags_[6] && central_body_ != "Jupiter") {
        a_jupiter = perturbation_from_body(state, "Jupiter", central_body_, et);
    }

    Eigen::Vector3d a_saturn(0.0, 0.0, 0.0);
    if (third_body_flags_[7] && central_body_ != "Saturn") {
        a_saturn = perturbation_from_body(state, "Saturn", central_body_, et);
    }

    Eigen::Vector3d a_uranus(0.0, 0.0, 0.0);
    if (third_body_flags_[8] && central_body_ != "Uranus") {
        a_uranus = perturbation_from_body(state, "Uranus", central_body_, et);
    }

    Eigen::Vector3d a_neptune(0.0, 0.0, 0.0);
    if (third_body_flags_[9] && central_body_ != "Neptune") {
        a_neptune = perturbation_from_body(state, "Neptune", central_body_, et);
    }

    Eigen::Vector3d a_nbody(0.0, 0.0, 0.0);
    a_nbody = a_sun + a_mercury + a_venus + a_earth + a_moon + a_mars + a_jupiter + a_saturn + a_uranus + a_neptune;
    return a_nbody;
}



Eigen::Vector3d EnvironmentModel::perturbation_from_body(const Eigen::Ref<const Eigen::VectorXd>& state, const std::string &third_body, const std::string &central_body, SpiceDouble et) const {

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

Eigen::Vector3d EnvironmentModel::solar_radiation_pressure(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const Spacecraft &sc) const {
    Eigen::Vector3d a_srp(0.0, 0.0, 0.0);

    // Get vector from spacecraft to the Sun
    Eigen::Vector3d r_sun2sc(0.0, 0.0, 0.0);
    if (central_body_ == std::string("Sun")) {
        r_sun2sc = state.head(3);
    }
    else {
        std::string cb(central_body_);
        if (central_body_ != std::string("Earth") && central_body_ != std::string("Sun") && central_body_ != std::string("Moon"))
            cb += " barycenter";
        SpiceDouble central_body_state[6];
        SpiceDouble lt = 0.0;
        if (central_body_ == std::string("Earth"))
            spkezr_c(cb.c_str(), et, "J2000", "NONE", "Sun", central_body_state, &lt);
        else
            spkezr_c(cb.c_str(), et, "ECLIPJ2000", "NONE", "Sun", central_body_state, &lt);

        r_sun2sc = {central_body_state[0] + state(0), central_body_state[1] + state(1), central_body_state[2] + state(2)};
    }

    // Shadow function
    double a = asin(constants::radius("Sun") / r_sun2sc.norm());
    double b = asin(constants::radius(central_body_) / state.head(3).norm());
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
    double reflectivity = sc.reflectivity();
    double area = sc.sun_exposed_area() * 1E-6; // transform from m^2 to km^2
    a_srp = - nu * p_sr * (1 + reflectivity) * area / sc.mass() * r_sun2sc / r_sun2sc.norm() * 1E-3;

    return a_srp;
}


void EnvironmentModel::init_exp_model() {
    // Second dimension is:
    // - h_ellp (km) lower bound (which corresponds to base altitude h0)
    // - h_ellp (km) upper bound
    // - nominal density (kg/m^3)
    // - scale height H (km)

    exp_atm_model_[0][0] = 0.0;
    exp_atm_model_[0][1] = 25.0;
    exp_atm_model_[0][2] = 1.225;
    exp_atm_model_[0][3] = 7.249;

    exp_atm_model_[1][0] = 25.0;
    exp_atm_model_[1][1] = 30.0;
    exp_atm_model_[1][2] = 3.899E-2;
    exp_atm_model_[1][3] = 6.349;

    exp_atm_model_[2][0] = 30.0;
    exp_atm_model_[2][1] = 40.0;
    exp_atm_model_[2][2] = 1.774E-2;
    exp_atm_model_[2][3] = 6.682;

    exp_atm_model_[3][0] = 40.0;
    exp_atm_model_[3][1] = 50.0;
    exp_atm_model_[3][2] = 3.972E-3;
    exp_atm_model_[3][3] = 7.554;

    exp_atm_model_[4][0] = 50.0;
    exp_atm_model_[4][1] = 60.0;
    exp_atm_model_[4][2] = 1.057E-3;
    exp_atm_model_[4][3] = 8.382;

    exp_atm_model_[5][0] = 60.0;
    exp_atm_model_[5][1] = 70.0;
    exp_atm_model_[5][2] = 3.206E-4;
    exp_atm_model_[5][3] = 7.714;

    exp_atm_model_[6][0] = 70.0;
    exp_atm_model_[6][1] = 80.0;
    exp_atm_model_[6][2] = 8.770E-5;
    exp_atm_model_[6][3] = 6.549;

    exp_atm_model_[7][0] = 80.0;
    exp_atm_model_[7][1] = 90.0;
    exp_atm_model_[7][2] = 1.905E-5;
    exp_atm_model_[7][3] = 5.799;

    exp_atm_model_[8][0] = 90.0;
    exp_atm_model_[8][1] = 100.0;
    exp_atm_model_[8][2] = 3.396E-6;
    exp_atm_model_[8][3] = 5.382;

    exp_atm_model_[9][0] = 100.0;
    exp_atm_model_[9][1] = 110.0;
    exp_atm_model_[9][2] = 5.297E-7;
    exp_atm_model_[9][3] = 5.877;

    exp_atm_model_[10][0] = 110.0;
    exp_atm_model_[10][1] = 120.0;
    exp_atm_model_[10][2] = 9.661E-8;
    exp_atm_model_[10][3] = 7.263;

    exp_atm_model_[11][0] = 120.0;
    exp_atm_model_[11][1] = 130.0;
    exp_atm_model_[11][2] = 2.438E-8;
    exp_atm_model_[11][3] = 9.473;

    exp_atm_model_[12][0] = 130.0;
    exp_atm_model_[12][1] = 140.0;
    exp_atm_model_[12][2] = 8.484E-9;
    exp_atm_model_[12][3] = 12.636;

    exp_atm_model_[13][0] = 140.0;
    exp_atm_model_[13][1] = 150.0;
    exp_atm_model_[13][2] = 3.845E-9;
    exp_atm_model_[13][3] = 16.149;

    exp_atm_model_[14][0] = 150.0;
    exp_atm_model_[14][1] = 180.0;
    exp_atm_model_[14][2] = 2.070E-9;
    exp_atm_model_[14][3] = 22.523;

    exp_atm_model_[15][0] = 180.0;
    exp_atm_model_[15][1] = 200.0;
    exp_atm_model_[15][2] = 5.464E-10;
    exp_atm_model_[15][3] = 29.740;

    exp_atm_model_[16][0] = 200.0;
    exp_atm_model_[16][1] = 250.0;
    exp_atm_model_[16][2] = 2.789E-10;
    exp_atm_model_[16][3] = 37.105;

    exp_atm_model_[17][0] = 250.0;
    exp_atm_model_[17][1] = 300.0;
    exp_atm_model_[17][2] = 7.248E-11;
    exp_atm_model_[17][3] = 45.546;

    exp_atm_model_[18][0] = 300.0;
    exp_atm_model_[18][1] = 350.0;
    exp_atm_model_[18][2] = 2.418E-11;
    exp_atm_model_[18][3] = 53.628;

    exp_atm_model_[19][0] = 350.0;
    exp_atm_model_[19][1] = 400.0;
    exp_atm_model_[19][2] = 9.518E-12;
    exp_atm_model_[19][3] = 53.298;

    exp_atm_model_[20][0] = 400.0;
    exp_atm_model_[20][1] = 450.0;
    exp_atm_model_[20][2] = 3.725E-12;
    exp_atm_model_[20][3] = 58.515;

    exp_atm_model_[21][0] = 450.0;
    exp_atm_model_[21][1] = 500.0;
    exp_atm_model_[21][2] = 1.585E-12;
    exp_atm_model_[21][3] = 60.828;

    exp_atm_model_[22][0] = 500.0;
    exp_atm_model_[22][1] = 600.0;
    exp_atm_model_[22][2] = 6.967E-13;
    exp_atm_model_[22][3] = 63.822;

    exp_atm_model_[23][0] = 600.0;
    exp_atm_model_[23][1] = 700.0;
    exp_atm_model_[23][2] = 1.454E-13;
    exp_atm_model_[23][3] = 71.835;

    exp_atm_model_[24][0] = 700.0;
    exp_atm_model_[24][1] = 800.0;
    exp_atm_model_[24][2] = 3.614E-14;
    exp_atm_model_[24][3] = 88.667;

    exp_atm_model_[25][0] = 800.0;
    exp_atm_model_[25][1] = 900.0;
    exp_atm_model_[25][2] = 1.170E-14;
    exp_atm_model_[25][3] = 124.64;

    exp_atm_model_[26][0] = 900.0;
    exp_atm_model_[26][1] = 1000.0;
    exp_atm_model_[26][2] = 5.245E-15;
    exp_atm_model_[26][3] = 181.05;

    exp_atm_model_[27][0] = 1000.0;
    exp_atm_model_[27][1] = 1100.0;
    exp_atm_model_[27][2] = 3.019E-15;
    exp_atm_model_[27][3] = 268.0;
}