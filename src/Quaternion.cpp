//
// Created by Sylvain  on 4/3/20.
//

#include "Quaternion.h"

#include <utility>

Quaternion::Quaternion(double x, double y, double z, double w) {
    m_q = Eigen::Vector4d(x, y, z, w);
}

Quaternion::Quaternion(Eigen::Vector4d  q): m_q(std::move(q)) {}

Eigen::Vector4d Quaternion::data() const {
    return m_q;
}

Eigen::Matrix3d Quaternion::cross_product_matrix() const {
    Eigen::Matrix3d skew;
    skew << 0., -m_q(2), m_q(1), m_q(2), 0., -m_q(0), -m_q(1), m_q(0), 0.;
    return skew;
}

Eigen::Matrix<double, 4, 3> Quaternion::xi() const {
    Eigen::Matrix<double, 4, 3> xi;
    xi.topRows(3) = m_q(3) * Eigen::Matrix3d::Identity() + this->cross_product_matrix();
    xi.bottomRows(1) = -m_q.head(3).transpose();
    return xi;
}

Eigen::Matrix<double, 4, 3> Quaternion::psi() const {
    Eigen::Matrix<double, 4, 3> psi;
    psi.topRows(3) = m_q(3) * Eigen::Matrix3d::Identity() - this->cross_product_matrix();
    psi.bottomRows(1) = -m_q.head(3).transpose();
    return psi;
}


Eigen::Matrix3d Quaternion::attitude_matrix() const {
    return xi().transpose() * psi();
}

double Quaternion::norm() const {
    return m_q.norm();
}

double Quaternion::norm_inf() const {
    return m_q.lpNorm<Eigen::Infinity>();
}

Quaternion Quaternion::abs() const {
    return Quaternion(m_q.array().abs());
}

Quaternion Quaternion::conjugate() const {
    return Quaternion(-m_q[0], -m_q[1], -m_q[2], m_q[3]);
}

Quaternion Quaternion::inverse() const {
    return conjugate() / (norm()*norm());
}

void Quaternion::normalize() {
    m_q /= m_q.norm();
}

double& Quaternion::operator[](int i) {
    return m_q[i];
}

const double& Quaternion::operator[](int i) const {
    return m_q[i];
}

Quaternion& Quaternion::operator+=(const Quaternion& q) {
    for (int i = 0; i < 4; ++i)
        m_q[i] += q.data()[i];
    return *this;
}

Quaternion& Quaternion::operator*=(double a) {
    for (int i = 0; i < 4; ++i)
        m_q[i] *= a;
    return *this;
}

Quaternion operator+(const Quaternion& q1, const Quaternion& q2) {
    return Quaternion(q1.m_q + q2.m_q);
}

Quaternion operator*(double a, const Quaternion& q) {
    return Quaternion(a * q.m_q);
}

Quaternion operator/(const Quaternion& q, double a) {
    return Quaternion(q.m_q / a);
}

Quaternion identity() {
    return Quaternion{0., 0., 0., 1.};
}