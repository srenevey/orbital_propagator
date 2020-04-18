//
// Created by Sylvain  on 4/17/20.
//

#ifndef ORBITAL_PROPAGATOR_QUATERNION_H
#define ORBITAL_PROPAGATOR_QUATERNION_H

#include <array>
#include "Vector.h"
#include "Matrix.h"
#include "ReferenceFrame.h"

template <class T, std::size_t N>
class Vector;

/** Representation of a quaternion.

    The first three elements of the quaternion represent the vector component and the fourth element is the scalar component.
*/
class Quaternion : public Vector<double, 4> {
public:
    Quaternion();
    Quaternion(ReferenceFrame base_frame, double x, double y, double z, double w);
    Quaternion(ReferenceFrame base_frame, std::array<double, 4> q);
    Quaternion(const Vector<double, 4>& q);

    [[nodiscard]] ReferenceFrame base_frame() const;
    [[nodiscard]] Matrix3d cross_product_matrix() const;
    [[nodiscard]] Matrix<double, 4, 3> xi() const;
    [[nodiscard]] Matrix<double, 4, 3> psi() const;

    /** Computes the rotation matrix to go from the base frame to the body-fixed frame.
     *
     *  If A = q.attitude_matrix(), then x_body = A * x_inertial and x_inertial = A.inverse() * x_body.
     */
    [[nodiscard]] Matrix3d attitude_matrix() const;
    [[nodiscard]] Quaternion conjugate() const;
    [[nodiscard]] Quaternion inverse() const;

    /** Returns the identity quaternion. */
    static Quaternion identity();
};

#endif //ORBITAL_PROPAGATOR_QUATERNION_H
