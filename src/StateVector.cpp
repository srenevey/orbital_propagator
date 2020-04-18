//
// Created by Sylvain  on 4/12/20.
//

#include "StateVector.h"

StateVector::StateVector(): m_frame(ReferenceFrame::NONE), m_ephemeris_time(0.0), m_position(Vector3d<Dimension::Distance>()), m_velocity(Vector3d<Dimension::Velocity>()), m_orientation(Quaternion()), m_ang_velocity(Vector3d<Dimension::AngularVelocity>()) {}

StateVector::StateVector(double ephemeris_time, Vector3d<Dimension::Distance> position, Vector3d<Dimension::Velocity> velocity, Quaternion orientation, Vector3d<Dimension::AngularVelocity> ang_velocity) :
    m_ephemeris_time(ephemeris_time),
    m_orientation(orientation)
{
    try {
        if (position.frame() != velocity.frame())
            throw std::invalid_argument("StateVector: The position and velocity must be in the same reference frame.");
        if (ang_velocity.frame() != ReferenceFrame::BODY)
            throw std::invalid_argument("StateVector: The angular velocity must be specified in the body-fixed frame.");

        m_frame = position.frame();
        m_position = position;
        m_velocity = velocity;
        m_ang_velocity = ang_velocity;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
}

ReferenceFrame StateVector::frame() const {
    return m_frame;
}

double StateVector::time() const {
    return m_ephemeris_time;
}

Vector3d<Dimension::Distance> StateVector::position() const {
    return m_position;
}

Vector3d<Dimension::Velocity> StateVector::velocity() const {
    return m_velocity;
}

Quaternion StateVector::orientation() const {
    return m_orientation;
}

Vector3d<Dimension::AngularVelocity> StateVector::ang_velocity() const {
    return m_ang_velocity;
}

void StateVector::set_frame(ReferenceFrame frame) {
    m_frame = frame;
}

void StateVector::set_time(double time) {
    m_ephemeris_time = time;
}

void StateVector::set_position(Vector3d<Dimension::Distance> position) {
    try {
        if (m_frame != position.frame())
            throw std::invalid_argument("StateVector: Cannot set position in the state vector: the reference frames are different.");

        m_position = position;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
}

void StateVector::set_velocity(Vector3d<Dimension::Velocity> velocity) {
    try {
        if (m_frame != velocity.frame())
            throw std::invalid_argument("StateVector: Cannot set velocity in the state vector: the reference frames are different.");

        m_velocity = velocity;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
}

void StateVector::set_orientation(Quaternion orientation) {
    m_orientation = orientation;
}

void StateVector::set_ang_velocity(Vector3d<Dimension::AngularVelocity> ang_velocity) {
    try {
        if (ang_velocity.frame() != ReferenceFrame::BODY)
            throw std::invalid_argument("StateVector: The angular velocity must be specified in the body-fixed frame.");

        m_ang_velocity = ang_velocity;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    m_ang_velocity = ang_velocity;
}

void StateVector::set_position_derivative(Vector3d<Dimension::Velocity> velocity) {
    m_position = Vector3d<Dimension::Distance>(velocity.frame(), {velocity[0], velocity[1], velocity[2]});
}

void StateVector::set_velocity_derivative(Vector3d<Dimension::Acceleration> acceleration) {
    m_velocity = Vector3d<Dimension::Velocity>(acceleration.frame(), {acceleration[0], acceleration[1], acceleration[2]});
}

void StateVector::set_ang_velocity_derivative(Vector3d<Dimension::AngularAcceleration> ang_acceleration) {
    m_ang_velocity = Vector3d<Dimension::AngularVelocity>(ang_acceleration.frame(), {ang_acceleration[0], ang_acceleration[1], ang_acceleration[2]});
}

StateVector StateVector::to_frame(ReferenceFrame frame) const {
    if (m_frame == frame)
        return *this;
    if (m_frame == ReferenceFrame::NONE || frame == ReferenceFrame::NONE)
        return *this;

    // Transform the position and velocity vectors
    std::string in_frame(m_frame);
    std::string out_frame(frame);
    double rot[6][6];
    sxform_c(in_frame.c_str(), out_frame.c_str(), m_ephemeris_time, rot);

    double in[6] = {m_position[0], m_position[1], m_position[2], m_velocity[0], m_velocity[1], m_velocity[2]};
    double out[6];
    mxvg_c(rot, in, 6, 6, out);

    Vector3d<Dimension::Distance> pos_out(frame, {out[0], out[1], out[2]});
    Vector3d<Dimension::Velocity> vel_out(frame, {Dimension::Velocity(out[3]), Dimension::Velocity(out[4]), Dimension::Velocity(out[5])});
    return StateVector(m_ephemeris_time, pos_out, vel_out, m_orientation, m_ang_velocity);
}

double StateVector::norm_inf() const {
    return std::max(std::max(m_position.norm_inf(), m_velocity.norm_inf()), std::max(m_orientation.norm_inf(), m_ang_velocity.norm_inf()));
}

void StateVector::normalize_quaternion() {
    m_orientation.normalize();
}


double StateVector::operator[](int i) const {
    if (i < 3)
        return m_position[i];
    else if (i < 6)
        return m_velocity[i-3];
    else if (i < 10)
        return m_orientation[i-6];
    else
        return m_ang_velocity[i-10];
}

StateVector& StateVector::operator+=(const StateVector& state) {
    m_position += state.position();
    m_velocity += state.velocity();
    m_orientation += state.orientation();
    m_ang_velocity += state.ang_velocity();
    return *this;
}

StateVector& StateVector::operator*=(double a) {
    m_position *= a;
    m_velocity *= a;
    m_orientation *= a;
    m_ang_velocity *= a;
    return *this;
}


StateVector operator+(const StateVector& state1, const StateVector& state2) {
    return StateVector(state1.time(), state1.position()+state2.position(), state1.velocity()+state2.velocity(), state1.orientation()+state2.orientation(), state1.ang_velocity()+state2.ang_velocity());
}

StateVector operator*(double a, const StateVector& state) {
    return StateVector(state.time(), a*state.position(), a*state.velocity(), a*state.orientation(), a*state.ang_velocity());
}

StateVector operator/(const StateVector& state1, const StateVector& state2) {
    return StateVector(state1.time(), state1.position()/state2.position(), state1.velocity()/state2.velocity(), state1.orientation()/state2.orientation(), state1.ang_velocity()/state2.ang_velocity());
}

StateVector abs(const StateVector &state) {
    return StateVector(state.time(), state.position().abs(), state.velocity().abs(), state.orientation().abs(), state.ang_velocity().abs());
}

std::ofstream& operator<<(std::ofstream &out, const StateVector& state) {
    std::stringstream ss;
    ss << state.m_ephemeris_time;
    ss << std::fixed << std::setprecision(6);
    for (int i = 0; i < 13; ++i) {
        ss << ", ";
        ss << state[i];
    }
    out << ss.rdbuf();
    return out;
}