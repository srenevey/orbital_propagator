//
// Created by Sylvain  on 4/7/20.
//

#include "Magnetometer.h"

Magnetometer::Magnetometer(Vector3d<double> bias, Vector3d<double> stddev, Matrix3d rotation_matrix, double max_value, double min_value): m_rot_bff2sff(rotation_matrix), m_max_value(max_value), m_min_value(min_value) {
    for (int i = 0; i < 3; ++i) {
        m_distribution[i] = std::normal_distribution<double>(bias[i], stddev[i]);
    }
}

Magnetometer::~Magnetometer() {}

Vector3d<double> Magnetometer::measure(const StateVector& state, double et, const EnvironmentModel& env_model) {

    Vector3d<double> B_eci = env_model.magnetic_field(state, et);
    Vector3d<double> noise(state.frame(), {m_distribution[0](m_generator), m_distribution[1](m_generator), m_distribution[2](m_generator)});
    Vector3d<double> B_meas_eci = B_eci + noise;
    Vector3d<double> B_meas_bff = transformations::rotate_to_frame(B_meas_eci, ReferenceFrame::BODY, et, state.orientation());
    Vector3d<double> B_meas_sff = m_rot_bff2sff * B_meas_bff;
    return B_meas_sff;
}

Vector3d<double> Magnetometer::bias() const {
    return Vector3d<double>({m_distribution[0].mean(), m_distribution[1].mean(), m_distribution[2].mean()});
}

Vector3d<double> Magnetometer::stddev() const {
    return Vector3d<double>({m_distribution[0].stddev(), m_distribution[1].stddev(), m_distribution[2].stddev()});
}