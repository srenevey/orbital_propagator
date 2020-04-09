//
// Created by Sylvain  on 4/3/20.
//

#ifndef ORBITAL_PROPAGATOR_QUATERNION_H
#define ORBITAL_PROPAGATOR_QUATERNION_H

#include "Eigen/Dense"

/** Representation of a quaternion.

    The first three elements of the quaternion represent the vector component and the fourth element is the scalar component.
*/
class Quaternion {
    Eigen::Vector4d m_q;
public:
    Quaternion(double x, double y, double z, double w);
    explicit Quaternion(Eigen::Vector4d q);

    [[nodiscard]] Eigen::Vector4d data() const;
    [[nodiscard]] Eigen::Matrix3d cross_product_matrix() const;
    [[nodiscard]] Eigen::Matrix<double, 4, 3> xi() const;
    [[nodiscard]] Eigen::Matrix<double, 4, 3> psi() const;

    /** Computes the rotation matrix to go from inertial to body-fixed frame.
     *
     *  If A = q.attitude_matrix(), then x_body = A * x_inertial and x_inertial = A.inverse() * x_body.
     */
    [[nodiscard]] Eigen::Matrix3d attitude_matrix() const;

    [[nodiscard]] double norm() const;
    [[nodiscard]] double norm_inf() const;
    [[nodiscard]] Quaternion abs() const;
    [[nodiscard]] Quaternion conjugate() const;
    [[nodiscard]] Quaternion inverse() const;

    void normalize();

    double& operator[](int i);
    const double& operator[](int i) const;
    Quaternion& operator+=(const Quaternion& q);
    Quaternion& operator*=(double a);
    friend Quaternion operator+(const Quaternion& q1, const Quaternion& q2);
    friend Quaternion operator*(double a, const Quaternion& q);
    friend Quaternion operator/(const Quaternion& q, double a);

    /** Returns the identity quaternion. */
    static Quaternion identity();
};


#endif //ORBITAL_PROPAGATOR_QUATERNION_H
