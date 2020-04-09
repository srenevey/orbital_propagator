//
// Created by Sylvain  on 3/28/20.
//

#ifndef ORBITAL_PROPAGATOR_STATE_H
#define ORBITAL_PROPAGATOR_STATE_H

#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include "ReferenceFrame.h"
#include "Quaternion.h"
#include "dimensions/Dimensions.h"
#include <boost/numeric/odeint.hpp>
#include <string>

/** Definition of the state vector.
 *
 *  The state vector contains 3 elements for the position, 3 elements for the velocity, 4 elements for the orientation, and 3 elements for the angular velocity.
 *  The ephemeris time and a reference frame are associated with the state.
 */
class State {
    double m_et;
    ReferenceFrame m_frame;
    std::array<double, 13> m_state;

public:
    State();
    State(const State& s);
    State(const double val);
    State(double et, ReferenceFrame frame, std::array<double, 13> state);
    State(double et, ReferenceFrame frame, Eigen::Vector3d position, Eigen::Vector3d velocity, Quaternion orientation, Eigen::Vector3d ang_velocity);
    State(double et, ReferenceFrame frame, std::array<double, 3> position, std::array<double, 3> velocity, std::array<double, 4> orientation, std::array<double, 3> ang_velocity);
    State(double et, ReferenceFrame frame, std::array<Dimension::Distance, 3> position, std::array<Dimension::Velocity, 3> velocity, Quaternion orientation, std::array<Dimension::AngularVelocity, 3> ang_velocity);

    [[nodiscard]] double et() const;
    [[nodiscard]] ReferenceFrame frame() const;
    [[nodiscard]] std::array<double, 13> state() const;
    [[nodiscard]] Eigen::Vector3d position() const;
    [[nodiscard]] Eigen::Vector3d velocity() const;
    [[nodiscard]] Quaternion orientation() const;
    [[nodiscard]] Eigen::Vector3d ang_velocity() const;

    void set_et(double et);
    void set_position(const Eigen::Vector3d& position);
    void set_velocity(const Eigen::Vector3d& velocity);
    void set_orientation(const Quaternion& orientation);
    void set_ang_velocity(const Eigen::Vector3d& ang_velocity);
    void set_position(const std::array<double, 3>& position);
    void set_velocity(const std::array<double, 3>& velocity);

    void normalize_quaternion();

    State& operator+=(const State& state);
    State& operator*=(double a);
    friend State operator+(const State& state1, const State& state2);
    friend State operator*(double a, const State& state);
    friend State operator/(const State& state1, const State& state2);
    friend State abs(const State &state);
    friend std::ofstream& operator<<(std::ofstream &out, const State &state);
};

std::ostream& operator<<( std::ostream &out , const State &state);

// Template specialization to be compatible with boost odeint
namespace boost::numeric::odeint {
    template<>
    struct vector_space_norm_inf< State >
    {
        typedef double result_type;
        double operator()( const State &state ) const
        {
            return std::max(std::max(state.position().lpNorm<Eigen::Infinity>(), state.velocity().lpNorm<Eigen::Infinity>()), std::max( state.orientation().norm_inf(), state.ang_velocity().lpNorm<Eigen::Infinity>()));
        }
    };
}

#endif //ORBITAL_PROPAGATOR_STATE_H
