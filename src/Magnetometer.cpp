//
// Created by Sylvain  on 4/7/20.
//

#include "Magnetometer.h"

Magnetometer::Magnetometer(std::array<double, 3> bias, std::array<double, 3> stddev, std::array<double, 9> rotation_matrix, double max_value, double min_value): m_max_value(max_value), m_min_value(min_value) {
    for (int i = 0; i < 3; ++i) {
        m_distribution[i] = std::normal_distribution<double>(bias[i], stddev[i]);
        for (int j = 0; j < 3; ++j)
            m_rot_bff2sff(i,j) = rotation_matrix[3*i+j];
    }
}

Magnetometer::~Magnetometer() {}

Eigen::Vector3d Magnetometer::measure(const State& state, double et, const EnvironmentModel& env_model) {
    Eigen::Vector3d B_eci = env_model.magnetic_field(state, et);
    Eigen::Vector3d noise(m_distribution[0](m_generator), m_distribution[1](m_generator), m_distribution[2](m_generator));
    Eigen::Vector3d B_meas_eci = B_eci + noise;
    Eigen::Vector3d B_meas_bff = state.orientation().attitude_matrix() * B_meas_eci;
    Eigen::Vector3d B_meas_sff = m_rot_bff2sff * B_meas_bff;
    return B_meas_sff;
}

Eigen::Vector3d Magnetometer::bias() const {
    return Eigen::Vector3d(m_distribution[0].mean(), m_distribution[1].mean(), m_distribution[2].mean());
}

Eigen::Vector3d Magnetometer::stddev() const {
    return Eigen::Vector3d(m_distribution[0].stddev(), m_distribution[1].stddev(), m_distribution[2].stddev());
}