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