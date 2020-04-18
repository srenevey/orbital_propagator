//
// Created by Sylvain  on 4/12/20.
//

#ifndef ORBITAL_PROPAGATOR_STATEVECTOR_H
#define ORBITAL_PROPAGATOR_STATEVECTOR_H

#include <exception>
#include <iostream>
#include <fstream>
#include "ReferenceFrame.h"
#include "dimensions/Dimensions.h"
#include "algebra/Vector.h"
#include "algebra/Quaternion.h"
#include <boost/numeric/odeint.hpp>
extern "C" {
#include "SpiceUsr.h"
};

/** State vector of a spacecraft containing its position, velocity, orientation, and angular velocity.
 *
 * The state vector is composed of 13 elements: 3 for the position, 3 for the velocity, 4 for the orientation (quaternion), and 3 for the angular velocity.
 */
class StateVector {
public:
    StateVector();

    /**
     * @param ephemeris_time    Ephemeris time.
     * @param position          Vector containing the position in km.
     * @param velocity          Vector containing the velocity in km/s.
     * @param orientation       Quaternion describing the rotation from a base frame to the body-fixed frame.
     * @param ang_velocity      Vector containing the angular velocity in rad/s.
     */
    StateVector(double ephemeris_time, Vector3d<Dimension::Distance> position, Vector3d<Dimension::Velocity> velocity, Quaternion orientation, Vector3d<Dimension::AngularVelocity> ang_velocity);

    [[nodiscard]] ReferenceFrame frame() const;
    [[nodiscard]] double time() const;
    [[nodiscard]] Vector3d<Dimension::Distance> position() const;
    [[nodiscard]] Vector3d<Dimension::Velocity> velocity() const;
    [[nodiscard]] Quaternion orientation() const;
    [[nodiscard]] Vector3d<Dimension::AngularVelocity> ang_velocity() const;

    void set_frame(ReferenceFrame frame);
    void set_time(double time);
    void set_position(Vector3d<Dimension::Distance> position);
    void set_velocity(Vector3d<Dimension::Velocity> velocity);
    void set_orientation(Quaternion orientation);
    void set_ang_velocity(Vector3d<Dimension::AngularVelocity> ang_velocity);

    void set_position_derivative(Vector3d<Dimension::Velocity> velocity);
    void set_velocity_derivative(Vector3d<Dimension::Acceleration> acceleration);
    void set_ang_velocity_derivative(Vector3d<Dimension::AngularAcceleration> ang_acceleration);

    /** Transforms the position and velocity vector into a given reference frame.
     *
     * @param frame ReferenceFrame to transform into.
     * @return StateVector where the position and velocity are expressed in the given frame.
     */
    [[nodiscard]] StateVector to_frame(ReferenceFrame frame) const;
    [[nodiscard]] double norm_inf() const;

    /** Normalizes the quaternion describing the orientation of the spacecraft. */
    void normalize_quaternion();

    double operator[](int i) const;
    StateVector& operator+=(const StateVector& state);
    StateVector& operator*=(double a);
    friend StateVector operator+(const StateVector& state1, const StateVector& state2);
    friend StateVector operator*(double a, const StateVector& state);
    friend StateVector operator/(const StateVector& state1, const StateVector& state2);
    friend StateVector abs(const StateVector &state);
    friend std::ofstream& operator<<(std::ofstream &out, const StateVector& state);

private:
    ReferenceFrame m_frame;
    double m_ephemeris_time;
    Vector3d<Dimension::Distance> m_position;
    Vector3d<Dimension::Velocity> m_velocity;
    Quaternion m_orientation;
    Vector3d<Dimension::AngularVelocity> m_ang_velocity;
};


// Template specialization to be compatible with boost odeint
namespace boost::numeric::odeint {
    template<>
    struct vector_space_norm_inf< StateVector >
    {
        typedef double result_type;
        double operator()( const StateVector &state ) const
        {
            return std::max(std::max(state.position().norm_inf(), state.velocity().norm_inf()), std::max( state.orientation().norm_inf(), state.ang_velocity().norm_inf()));
        }
    };
}

#endif //ORBITAL_PROPAGATOR_STATEVECTOR_H
