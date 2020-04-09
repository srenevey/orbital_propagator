//
// Created by Sylvain Renevey on 1/26/18.
//



#include "EnvironmentModel.h"

EnvironmentModel::EnvironmentModel():
        m_gp_degree(0),
        m_atm_model(AtmModel::None),
        m_third_body({}),
        m_earth_gram_atm_model(nullptr),
        m_srp_flag(false)
{
    m_central_body = new Earth();
}


EnvironmentModel::EnvironmentModel(Body central_body, int gp_degree, const std::string& geopot_model_path, const std::vector<Body>& third_body, AtmModel atm_model, const std::string& earthgram_path, bool srp_flag, const std::string& epoch):
        m_gp_degree(gp_degree),
        m_atm_model(atm_model),
        m_srp_flag(srp_flag),
        m_earth_gram_atm_model(nullptr)
{
    m_central_body = BodyContainer::create_body(central_body);
    for (auto body: third_body) {
        m_third_body.push_back(BodyContainer::create_body(body));
    }

    // Load geopotential coefficients from model file
    if (gp_degree > 1 && !geopot_model_path.empty())
        m_cs_coeffs = load_coefficients(gp_degree, geopot_model_path);

    // If an atmospheric model is selected, initialize coefficients
    if (m_atm_model == AtmModel::EarthGRAM && m_central_body->name() == Body::Earth) {
        if (!epoch.empty()) {
            m_earth_gram_atm_model = new Atm1;
            m_earth_gram_atm_model->initdata(earthgram_path, "NameRef.txt", epoch);
        } else {
            std::cerr << "The epoch must be passed in as argument." << std::endl;
        }
    } else if (m_atm_model == AtmModel::Exponential && m_central_body->name() == Body::Earth) {
        init_exp_model();
    } else if (m_atm_model != AtmModel::None && m_central_body->name() != Body::Earth) {
        m_atm_model = AtmModel::None;
        std::cerr << "The drag effects are currently only supported for the Earth." << std::endl;
    }
}


EnvironmentModel::~EnvironmentModel() {
    delete m_central_body;
    for (auto body: m_third_body)
        delete body;
    delete m_earth_gram_atm_model;
}


std::vector<BodyContainer*> EnvironmentModel::third_body() const {
    return m_third_body;
}


int EnvironmentModel::gp_degree() const {
    return m_gp_degree;
}


const Eigen::MatrixXd &EnvironmentModel::cs_coeffs() const {
    return m_cs_coeffs;
}


BodyContainer* EnvironmentModel::central_body() const {
    return m_central_body;
}

bool EnvironmentModel::is_drag() const {
    return m_atm_model != AtmModel::None;
}


bool EnvironmentModel::is_srp_flag() const {
    return m_srp_flag;
}

Eigen::MatrixXd EnvironmentModel::load_coefficients(const int degree, const std::string model_file_name) {
    // Initialize matrix containing the coefficients
    int num_entries = (degree + 1) * (degree + 2) / 2;

    Eigen::MatrixXd coefficients = Eigen::MatrixXd::Zero(num_entries, 2);
    coefficients.row(0) << 0, 0;
    coefficients.row(1) << 0, 0;
    coefficients.row(2) << 0, 0;

    // Open file containing the geopotential coefficients
    ifstream data_file;
    data_file.open(model_file_name);
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

double EnvironmentModel::atm_density(const State& state, const double elapsed_time, const SpiceDouble et) const {

    double density = 0.0;
    SpiceDouble inertial_position[3];
    for (int i = 0; i < 3; ++i)
        inertial_position[i] = state.position()[i];

    if (m_atm_model == AtmModel::EarthGRAM) {

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

        m_earth_gram_atm_model->traj(alt, lat, lon, elapsed_time, iupdate, initonce, &dm1, &pm1, &tm1, &um1, &vm1, &wm1,
                                     &dp1, &pp1, &tp1, &up1, &vp1, &wp1, &ds1, &ps1, &ts1, &us1, &vs1, &ws1, &dsmall,
                                     &psmall, &tsmall, &usmall, &vsmall, &wsmall, &sos, &sosp);
        density = dm1;

    } else if (m_atm_model == AtmModel::Exponential) {

        // Retrieve Earth's radius
        SpiceInt dim = 0;
		SpiceDouble radius[3];
		bodvrd_c("EARTH", "RADII", 3, &dim, radius);

		// Compute the altitude above the ellipsoid. There is no atmosphere above 1100 km.
        double h_ell = vnorm_c(inertial_position) - radius[0];
        if (h_ell <= 1100.0) {
            for (int i = 0; i < 28; ++i) {
                if (m_exp_atm_model[i][0] <= h_ell && m_exp_atm_model[i][1] > h_ell) {
                    double h0 = m_exp_atm_model[i][0];
                    double rho0 = m_exp_atm_model[i][2];
                    double H = m_exp_atm_model[i][3];
                    density = rho0 * exp(- (h_ell - h0) / H);
                }
            }
        } else {
            density = 0.0;
        }
    }
    density = density * 1E9; // convert from kg/m^3 to kg/km^3
    return density;
}

Eigen::Vector3d EnvironmentModel::body_vector(const BodyContainer* body, const SpiceDouble et) const {

    std::string cb = std::to_string(m_central_body->naif_id());
    std::string tb = std::to_string(body->naif_id());
    /*
    if (m_central_body != std::string("Earth") && m_central_body != std::string("Sun") && central_body != std::string("Moon"))
        cb += " barycenter";
    if (third_body != std::string("Earth") && third_body != std::string("Sun") && third_body != std::string("Moon"))
        tb += " barycenter";
    */

    // Get the state of the third body relative to the central body
    SpiceDouble third_body_state[6];
    SpiceDouble lt = 0;
    if (m_central_body->name() == Body::Earth)
        spkezr_c(tb.c_str(), et, "J2000", "NONE", cb.c_str(), third_body_state, &lt);
    else
        spkezr_c(tb.c_str(), et, "ECLIPJ2000", "NONE", cb.c_str(), third_body_state, &lt);

   return Eigen::Vector3d(third_body_state[0], third_body_state[1], third_body_state[2]); // central body to third body
}

Eigen::MatrixXd EnvironmentModel::geopotential_harmonic_coeff(const State& state, SpiceDouble et) const {
    // Transform the position from J2000 (inertial) to ITRF93 (rotating).
    SpiceDouble rot[3][3];
    pxform_c("J2000", "ITRF93", et, rot);

    SpiceDouble ecef_position[3];
    SpiceDouble s_eci_position[3];
    for (int i = 0; i < 3; ++i)
        s_eci_position[i] = state.position()[i];
    mxv_c(rot, s_eci_position, ecef_position);


    // Returns the index for degree n and order m.
    auto index = [](int n, int m) { return n * (n + 1) / 2 + m; };

    // Lambda expressions used as helpers to get C and S coefficients
    Eigen::MatrixXd CS_coeffs = cs_coeffs();
    auto C = [&](int m, int n) { return CS_coeffs(index(m, n), 0); };
    auto S = [&](int m, int n) { return CS_coeffs(index(m, n), 1); };


    // Compute V and W coefficients
    int num_rows = (m_gp_degree + 2) * (m_gp_degree + 3) / 2;
    double r = vnorm_c(ecef_position);
    Eigen::MatrixXd VW_coeffs = Eigen::MatrixXd::Zero(num_rows, 2);
    VW_coeffs(0, 0) = constants::R_EARTH_EGM08 / r;
    VW_coeffs(0, 1) = 0;

    // Zonal terms
    VW_coeffs(index(1, 0), 0) = ecef_position[2] * constants::R_EARTH_EGM08 / (r * r) * VW_coeffs(index(0, 0), 0);
    VW_coeffs(index(1, 0), 1) = 0;
    for (int n = 2; n < m_gp_degree + 2; ++n) {
        VW_coeffs(index(n, 0), 0) =
                (2 * n - 1) * ecef_position[2] * constants::R_EARTH_EGM08 / (n * r * r) * VW_coeffs(index(n - 1, 0), 0) -
                (n - 1) * constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08 / (n * r * r) * VW_coeffs(index(n - 2, 0), 0);
        VW_coeffs(index(n, 0), 1) = 0;
    }

    // Tesseral terms
    for (int n = 1; n < m_gp_degree + 2; ++n) {
        VW_coeffs(index(n, n), 0) = (2 * n - 1) *
                                    (ecef_position[0] * constants::R_EARTH_EGM08 / (r * r) * VW_coeffs(index(n - 1, n - 1), 0) -
                                     ecef_position[1] * constants::R_EARTH_EGM08 / (r * r) * VW_coeffs(index(n - 1, n - 1), 1));
        VW_coeffs(index(n, n), 1) = (2 * n - 1) *
                                    (ecef_position[0] * constants::R_EARTH_EGM08 / (r * r) * VW_coeffs(index(n - 1, n - 1), 1) +
                                     ecef_position[1] * constants::R_EARTH_EGM08 / (r * r) * VW_coeffs(index(n - 1, n - 1), 0));
    }

    // Remaining terms
    for (int n = 1; n < m_gp_degree + 2; ++n) {
        for (int m = 1; m < n; ++m) {
            if (n == m + 1) {
                VW_coeffs(index(n, m), 0) = (2 * n - 1) * ecef_position[2] * constants::R_EARTH_EGM08 / ((n - m) * r * r) *
                                            VW_coeffs(index(n - 1, m), 0);
                VW_coeffs(index(n, m), 1) = (2 * n - 1) * ecef_position[2] * constants::R_EARTH_EGM08 / ((n - m) * r * r) *
                                            VW_coeffs(index(n - 1, m), 1);
            } else {
                VW_coeffs(index(n, m), 0) = (2 * n - 1) * ecef_position[2] * constants::R_EARTH_EGM08 / ((n - m) * r * r) *
                                            VW_coeffs(index(n - 1, m), 0) -
                                            (n + m - 1) * constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08 / ((n - m) * r * r) *
                                            VW_coeffs(index(n - 2, m), 0);
                VW_coeffs(index(n, m), 1) = (2 * n - 1) * ecef_position[2] * constants::R_EARTH_EGM08 / ((n - m) * r * r) *
                                            VW_coeffs(index(n - 1, m), 1) -
                                            (n + m - 1) * constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08 / ((n - m) * r * r) *
                                            VW_coeffs(index(n - 2, m), 1);
            }
        }
    }
    return VW_coeffs;
}


Eigen::Vector3d EnvironmentModel::sun_spacecraft_vector(const State& state, SpiceDouble et) const {
    Eigen::Vector3d r_sun2sc(0.0, 0.0, 0.0);
    if (m_central_body->name() == Body::Sun) {
        r_sun2sc = state.position();
    }
    else {
        std::string cb = std::to_string(m_central_body->naif_id());

        //if (m_central_body->name() != Body::Earth && m_central_body->name() != Body::Sun && m_central_body->name() != Body::Moon)
        //    cb += " barycenter";

        SpiceDouble central_body_state[6];
        SpiceDouble lt = 0.0;
        if (m_central_body->name() == Body::Earth) // TODO: update to return vector in spacecraft ref frame
            spkezr_c(cb.c_str(), et, "J2000", "NONE", "Sun", central_body_state, &lt);
        else
            spkezr_c(cb.c_str(), et, "ECLIPJ2000", "NONE", "Sun", central_body_state, &lt);

        r_sun2sc = {central_body_state[0] + state.position()[0], central_body_state[1] + state.position()[1], central_body_state[2] + state.position()[2]};
    }
    return r_sun2sc;
}

double EnvironmentModel::in_shadow(const State& state, SpiceDouble et) const {
    Eigen::Vector3d r_sun2sc = sun_spacecraft_vector(state, et);

    double a = asin(Sun::RADIUS / r_sun2sc.norm());
    double b = asin(m_central_body->radius() / state.position().norm());
    double c = acos( state.position().dot(r_sun2sc) / (state.position().norm() * r_sun2sc.norm()));

    double nu = 0.0;
    if (abs(a-b) < c && c < a + b) {
        double x = (c * c + a * a - b * b) / (2.0 * c);
        double y = sqrt(a * a - x * x);
        double A = a * a * acos(x / a) + b * b * acos((c - x) / b) - c * y;
        nu = 1.0 - A / (constants::PI * a * a);
    } else if (c >= a + b) {
        nu = 0.0;
    } else if (c < b - a) {
        nu = 1.0;
    }
    return nu;
}

void EnvironmentModel::init_exp_model() {
    // Second dimension is:
    // - h_ellp (km) lower bound (which corresponds to base altitude h0)
    // - h_ellp (km) upper bound
    // - nominal density (kg/m^3)
    // - scale height H (km)

    m_exp_atm_model[0][0] = 0.0;
    m_exp_atm_model[0][1] = 25.0;
    m_exp_atm_model[0][2] = 1.225;
    m_exp_atm_model[0][3] = 7.249;

    m_exp_atm_model[1][0] = 25.0;
    m_exp_atm_model[1][1] = 30.0;
    m_exp_atm_model[1][2] = 3.899E-2;
    m_exp_atm_model[1][3] = 6.349;

    m_exp_atm_model[2][0] = 30.0;
    m_exp_atm_model[2][1] = 40.0;
    m_exp_atm_model[2][2] = 1.774E-2;
    m_exp_atm_model[2][3] = 6.682;

    m_exp_atm_model[3][0] = 40.0;
    m_exp_atm_model[3][1] = 50.0;
    m_exp_atm_model[3][2] = 3.972E-3;
    m_exp_atm_model[3][3] = 7.554;

    m_exp_atm_model[4][0] = 50.0;
    m_exp_atm_model[4][1] = 60.0;
    m_exp_atm_model[4][2] = 1.057E-3;
    m_exp_atm_model[4][3] = 8.382;

    m_exp_atm_model[5][0] = 60.0;
    m_exp_atm_model[5][1] = 70.0;
    m_exp_atm_model[5][2] = 3.206E-4;
    m_exp_atm_model[5][3] = 7.714;

    m_exp_atm_model[6][0] = 70.0;
    m_exp_atm_model[6][1] = 80.0;
    m_exp_atm_model[6][2] = 8.770E-5;
    m_exp_atm_model[6][3] = 6.549;

    m_exp_atm_model[7][0] = 80.0;
    m_exp_atm_model[7][1] = 90.0;
    m_exp_atm_model[7][2] = 1.905E-5;
    m_exp_atm_model[7][3] = 5.799;

    m_exp_atm_model[8][0] = 90.0;
    m_exp_atm_model[8][1] = 100.0;
    m_exp_atm_model[8][2] = 3.396E-6;
    m_exp_atm_model[8][3] = 5.382;

    m_exp_atm_model[9][0] = 100.0;
    m_exp_atm_model[9][1] = 110.0;
    m_exp_atm_model[9][2] = 5.297E-7;
    m_exp_atm_model[9][3] = 5.877;

    m_exp_atm_model[10][0] = 110.0;
    m_exp_atm_model[10][1] = 120.0;
    m_exp_atm_model[10][2] = 9.661E-8;
    m_exp_atm_model[10][3] = 7.263;

    m_exp_atm_model[11][0] = 120.0;
    m_exp_atm_model[11][1] = 130.0;
    m_exp_atm_model[11][2] = 2.438E-8;
    m_exp_atm_model[11][3] = 9.473;

    m_exp_atm_model[12][0] = 130.0;
    m_exp_atm_model[12][1] = 140.0;
    m_exp_atm_model[12][2] = 8.484E-9;
    m_exp_atm_model[12][3] = 12.636;

    m_exp_atm_model[13][0] = 140.0;
    m_exp_atm_model[13][1] = 150.0;
    m_exp_atm_model[13][2] = 3.845E-9;
    m_exp_atm_model[13][3] = 16.149;

    m_exp_atm_model[14][0] = 150.0;
    m_exp_atm_model[14][1] = 180.0;
    m_exp_atm_model[14][2] = 2.070E-9;
    m_exp_atm_model[14][3] = 22.523;

    m_exp_atm_model[15][0] = 180.0;
    m_exp_atm_model[15][1] = 200.0;
    m_exp_atm_model[15][2] = 5.464E-10;
    m_exp_atm_model[15][3] = 29.740;

    m_exp_atm_model[16][0] = 200.0;
    m_exp_atm_model[16][1] = 250.0;
    m_exp_atm_model[16][2] = 2.789E-10;
    m_exp_atm_model[16][3] = 37.105;

    m_exp_atm_model[17][0] = 250.0;
    m_exp_atm_model[17][1] = 300.0;
    m_exp_atm_model[17][2] = 7.248E-11;
    m_exp_atm_model[17][3] = 45.546;

    m_exp_atm_model[18][0] = 300.0;
    m_exp_atm_model[18][1] = 350.0;
    m_exp_atm_model[18][2] = 2.418E-11;
    m_exp_atm_model[18][3] = 53.628;

    m_exp_atm_model[19][0] = 350.0;
    m_exp_atm_model[19][1] = 400.0;
    m_exp_atm_model[19][2] = 9.518E-12;
    m_exp_atm_model[19][3] = 53.298;

    m_exp_atm_model[20][0] = 400.0;
    m_exp_atm_model[20][1] = 450.0;
    m_exp_atm_model[20][2] = 3.725E-12;
    m_exp_atm_model[20][3] = 58.515;

    m_exp_atm_model[21][0] = 450.0;
    m_exp_atm_model[21][1] = 500.0;
    m_exp_atm_model[21][2] = 1.585E-12;
    m_exp_atm_model[21][3] = 60.828;

    m_exp_atm_model[22][0] = 500.0;
    m_exp_atm_model[22][1] = 600.0;
    m_exp_atm_model[22][2] = 6.967E-13;
    m_exp_atm_model[22][3] = 63.822;

    m_exp_atm_model[23][0] = 600.0;
    m_exp_atm_model[23][1] = 700.0;
    m_exp_atm_model[23][2] = 1.454E-13;
    m_exp_atm_model[23][3] = 71.835;

    m_exp_atm_model[24][0] = 700.0;
    m_exp_atm_model[24][1] = 800.0;
    m_exp_atm_model[24][2] = 3.614E-14;
    m_exp_atm_model[24][3] = 88.667;

    m_exp_atm_model[25][0] = 800.0;
    m_exp_atm_model[25][1] = 900.0;
    m_exp_atm_model[25][2] = 1.170E-14;
    m_exp_atm_model[25][3] = 124.64;

    m_exp_atm_model[26][0] = 900.0;
    m_exp_atm_model[26][1] = 1000.0;
    m_exp_atm_model[26][2] = 5.245E-15;
    m_exp_atm_model[26][3] = 181.05;

    m_exp_atm_model[27][0] = 1000.0;
    m_exp_atm_model[27][1] = 1100.0;
    m_exp_atm_model[27][2] = 3.019E-15;
    m_exp_atm_model[27][3] = 268.0;
}

Eigen::Vector3d EnvironmentModel::magnetic_field(const State& state, SpiceDouble et) const {
    Eigen::Vector3d m = pow(constants::IGRF_2020_RADIUS, 3) * Eigen::Vector3d(constants::IGRF_2020_G11, constants::IGRF_2020_H11, constants::IGRF_2020_G01);

    // Transform the position from J2000 (inertial) to ITRF93 (rotating).
    SpiceDouble rot[3][3];
    pxform_c("J2000", "ITRF93", et, rot);
    SpiceDouble ecef_position[3];
    SpiceDouble s_eci_position[3];
    for (int i = 0; i < 3; ++i)
        s_eci_position[i] = state.position()[i];
    mxv_c(rot, s_eci_position, ecef_position);
    Eigen::Vector3d r_ecef(ecef_position[0], ecef_position[1], ecef_position[2]);

    Eigen::Vector3d B = 3.0 * (m.dot(r_ecef) * r_ecef - r_ecef.norm()*r_ecef.norm()*m) / pow(r_ecef.norm(), 5);
    SpiceDouble B_ecef[3]{B[0], B[1], B[2]};

    // Transform from ITRF93 (rotating) to J2000 (inertial)
    pxform_c("ITRF93", "J2000", et, rot);
    SpiceDouble B_eci[3];
    mxv_c(rot, B_ecef, B_eci);

    return Eigen::Vector3d(B_eci[0], B_eci[1], B_eci[2]);
}