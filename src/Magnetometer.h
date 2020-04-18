//
// Created by Sylvain  on 4/7/20.
//

#ifndef ORBITAL_PROPAGATOR_MAGNETOMETER_H
#define ORBITAL_PROPAGATOR_MAGNETOMETER_H

#include <array>
#include <random>
#include "Sensor.h"
#include "StateVector.h"
#include "algebra/Vector.h"
#include "algebra/Matrix.h"

/** Model of a 3-axis magnetometer (work in progress). */
class Magnetometer: public Sensor {
public:
    /**
     * @param bias                  Biases on the three axes
     * @param stddev                Standard deviation on the three axes
     * @param rotation_matrix       Rotation matrix to go from body-fixed to sensor-fixed frame
     * @param max_value             Maximum value that the magnetometer can read
     * @param min_value             Minimum value that the magnetometer can read
     */
    Magnetometer(Vector3d<double> bias, Vector3d<double> stddev, Matrix3d rotation_matrix, double max_value, double min_value);
    ~Magnetometer();

    /** Measures the local magnetic field.
     *
     * @param state         State vector of the spacecraft
     * @param et            Ephemeris time
     * @param env_model     Reference to the environment model
     * @return              Magnetic field in the sensor-fixed frame (in Tesla).
     */
    [[nodiscard]] Vector3d<double> measure(const StateVector& state, double et, const EnvironmentModel& env_model);
    [[nodiscard]] Vector3d<double> bias() const;
    [[nodiscard]] Vector3d<double> stddev() const;

private:
    Matrix3d m_rot_bff2sff; // rotation matrix from body-fixed frame to sensor-fixed frame
    double m_max_value;
    double m_min_value;
    std::default_random_engine m_generator;
    std::array<std::normal_distribution<double>, 3> m_distribution;
};


#endif //ORBITAL_PROPAGATOR_MAGNETOMETER_H
