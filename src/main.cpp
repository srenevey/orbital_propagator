#include <iostream>
#include <array>
#include <vector>
#include "dimensions/Dimensions.h"
#include "algebra/Vector.h"
#include "algebra/Matrix.h"
#include "algebra/Quaternion.h"
#include "StateVector.h"
#include "EnvironmentModel.h"
#include "Body.h"
#include "Sim.h"
#include "Magnetometer.h"

using namespace Unit;

int main() {

    // Define paths to the meta-kernel, geopotential model and EarthGRAM lib
    std::string meta_kernel(R"(/Users/sylvain/Developer/orbital_propagator/assets/kernels.tm)");
    std::string geopot_model_path(R"(/Users/sylvain/Developer/orbital_propagator/assets/EGM2008_to2190_TideFree)");
    std::string earthgram_path(R"(/Users/sylvain/Developer/orbital_propagator/libs/earthGRAM2016/)");

    // Create a Sim object to load the SPICE kernels
    Sim sim(meta_kernel);

    // Create the initial state of the spacecraft
    Vector3d<Dimension::Distance> position(ReferenceFrame::J2000, {-7009.22977_km, 0.0_km, 0.0_km});
    Vector3d<Dimension::Velocity> velocity(ReferenceFrame::J2000, {0.0_kms, -7.0033_kms, 3.2657_kms});
    Quaternion orientation(ReferenceFrame::J2000, 0., 0., 0., 1.);
    Vector3d<Dimension::AngularVelocity> angular_velocity(ReferenceFrame::BODY, {0.2_degs, 0.0_degs, 0.1_degs});
    StateVector initial_state(121000.0, position, velocity, orientation, angular_velocity);


    // Model of a 1U CubeSat
    // Spacecraft are represented as a collection of N flat plates.
    double mass = 1.33;
    Dimension::Distance side = 0.1_m;
    double I_ii = mass * side * side / 6.;
    Matrix3d inertia_matrix({{{I_ii, 0., 0.}, {0., I_ii, 0.}, {0., 0., I_ii}}});
    std::vector<Vector3d<double>> face_normals;
    face_normals.push_back(Vector3d<double>(ReferenceFrame::BODY, {1., 0., 0.}));
    face_normals.push_back(Vector3d<double>(ReferenceFrame::BODY, {0., 1., 0.}));
    face_normals.push_back(Vector3d<double>(ReferenceFrame::BODY, {0., 0., 1.}));
    face_normals.push_back(Vector3d<double>(ReferenceFrame::BODY, {-1., 0., 0.}));
    face_normals.push_back(Vector3d<double>(ReferenceFrame::BODY, {0., -1., 0.}));
    face_normals.push_back(Vector3d<double>(ReferenceFrame::BODY, {0., 0., -1.}));

    Dimension::Area face_area = side * side;
    std::vector<Dimension::Area> face_areas{face_area, face_area, face_area, face_area, face_area, face_area};

    std::vector<Vector3d<Dimension::Distance>> cop_positions;
    cop_positions.push_back(Vector3d<Dimension::Distance>(ReferenceFrame::BODY, {side/2., 0._m, 0._m}));
    cop_positions.push_back(Vector3d<Dimension::Distance>(ReferenceFrame::BODY, {0._m, side/2., 0._m}));
    cop_positions.push_back(Vector3d<Dimension::Distance>(ReferenceFrame::BODY, {0._m, 0._m, side/2.}));
    cop_positions.push_back(Vector3d<Dimension::Distance>(ReferenceFrame::BODY, {-side/2., 0._m, 0._m}));
    cop_positions.push_back(Vector3d<Dimension::Distance>(ReferenceFrame::BODY, {0._m, -side/2., 0._m}));
    cop_positions.push_back(Vector3d<Dimension::Distance>(ReferenceFrame::BODY, {0._m, 0._m, -side/2.}));

    std::vector<double> specular_reflection_coeff{0.042, 0.042, 0.184, 0.042, 0.042, 0.184}; // Solar panel: 0.042. Gold foil: 0.184.
    std::vector<double> diffuse_reflection_coeff{0.168, 0.168, 0.736, 0.168, 0.168, 0.736}; // Solar panel: 0.168. Gold foil: 0.736.

    // Sensors (work in progress)
    Vector3d<double> bias({1.0E-7, 1.0E-7, 1.0E-7}); // Tesla
    Vector3d<double> std_dev({5.0E-6, 5.0E-6, 5.0E-6});
    Matrix3d rot_bff2sff({{{0., 0.2588, 0.9659}, {-1., 0., 0.}, {0., -0.9659, 0.2588}}}); // rotation matrix from body-fixed frame to sensor-fixed frame
    double max_value = 2.0E-4; // Tesla
    double min_value = 0.0; // Tesla
    std::shared_ptr<Magnetometer> mag1 = std::make_shared<Magnetometer>(bias, std_dev, rot_bff2sff, max_value, min_value);
    std::vector<std::shared_ptr<Sensor>> sensors;
    sensors.push_back(std::move(mag1));

    double drag_coeff = 2.1;

    // Create a spacecraft object. The name is used to save the output of the simulation in a file "name.dat".
    Spacecraft sc("sc", initial_state, mass, drag_coeff, specular_reflection_coeff, diffuse_reflection_coeff, inertia_matrix, face_normals, face_areas, cop_positions, sensors);

    // Define the initial epoch, final time, and step size.
    std::string epoch = "2020-04-05 16:00:00";
    Dimension::Time t0 = 0._s;
    Dimension::Time t1 = 24._h;
    Dimension::Time dt = 5._s;

    // Create the environment model
    EnvironmentModel env_model(Body::Earth, 6, geopot_model_path, {Body::Sun, Body::Moon}, AtmModel::EarthGRAM, earthgram_path, epoch, true, true);

    // Integrate the equations of motion.
    sim.integrate(sc, env_model, epoch, t0, t1, dt);
    std::cout << "Simulation results saved in \"" << sc.name() << ".dat\"." << std::endl;


    return 0;
}