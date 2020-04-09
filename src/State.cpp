//
// Created by Sylvain  on 3/28/20.
//

#include "State.h"

State::State(): m_et(0.), m_frame(ReferenceFrame::J2000), m_state({}) {}

State::State(const State& s): m_et(s.et()), m_frame(s.frame()), m_state(s.state()) {}

State::State(const double val): m_et(0.), m_frame(ReferenceFrame::J2000), m_state({val, val, val, val, val, val, val, val, val, val, val, val, val}) {}

State::State(const double et, const ReferenceFrame frame, const std::array<double, 13> state):
    m_et(et),
    m_frame(frame),
    m_state(state)
{}

State::State(const double et, const ReferenceFrame frame, const Eigen::Vector3d position, const Eigen::Vector3d velocity, const Quaternion orientation, const Eigen::Vector3d ang_velocity):
        m_et(et),
        m_frame(frame)
{
    for (int i = 0; i < 3; ++i) {
        m_state[i] = position[i];
        m_state[i+3] = velocity[i];
        m_state[i+6] = orientation[i];
        m_state[i+10] = ang_velocity[i];
    }
    m_state[9] = orientation[3];
}

State::State(const double et, const ReferenceFrame frame, const std::array<double, 3> position, const std::array<double, 3> velocity, const std::array<double, 4> orientation, const std::array<double, 3> ang_velocity):
    m_et(et),
    m_frame(frame)
{
    for (int i = 0; i < 3; ++i) {
        m_state[i] = position[i];
        m_state[i+3] = velocity[i];
        m_state[i+6] = orientation[i];
        m_state[i+10] = ang_velocity[i];
    }
    m_state[9] = orientation[3];
}

State::State(double et, ReferenceFrame frame, const std::array<Dimension::Distance, 3> position, const std::array<Dimension::Velocity, 3> velocity, const Quaternion orientation, const std::array<Dimension::AngularVelocity, 3> ang_velocity):
    m_et(et),
    m_frame(frame)
{
    for (int i = 0; i < 3; ++i) {
        m_state[i] = position[i].data();
        m_state[i+3] = velocity[i].data();
        m_state[i+6] = orientation[i];
        m_state[i+10] = ang_velocity[i].data();
    }
    m_state[9] = orientation[3];
}

double State::et() const {
    return m_et;
}

std::array<double, 13> State::state() const {
    return m_state;
}

ReferenceFrame State::frame() const {
    return m_frame;
}

void State::set_et(const double et) {
    m_et = et;
}

void State::set_position(const Eigen::Vector3d& position) {
    for (int i = 0; i < 3; ++i)
        m_state[i] = position[i];
}

void State::set_velocity(const Eigen::Vector3d& velocity) {
    for (int i = 0; i < 3; ++i)
        m_state[i+3] = velocity[i];
}

void State::set_orientation(const Quaternion& orientation) {
    for (int i = 0; i < 4; ++i)
        m_state[i+6] = orientation[i];
}

void State::set_ang_velocity(const Eigen::Vector3d& ang_velocity) {
    for (int i = 0; i < 3; ++i)
        m_state[i+10] = ang_velocity[i];
}


void State::set_position(const std::array<double, 3>& position) {
    for (int i = 0; i < 3; ++i)
        m_state[i] = position[i];
}

void State::set_velocity(const std::array<double, 3>& velocity) {
    for (int i = 0; i < 3; ++i)
        m_state[i+3] = velocity[i];
}

void State::normalize_quaternion() {
    double norm = sqrt(m_state[6]*m_state[6] + m_state[7]*m_state[7] + m_state[8]*m_state[8] + m_state[9]*m_state[9]);
    for (int i = 6; i < 10; ++i)
        m_state[i] /= norm;
}

Eigen::Vector3d State::position() const {
    return Eigen::Vector3d(m_state[0], m_state[1], m_state[2]);
}

Eigen::Vector3d State::velocity() const {
    return Eigen::Vector3d(m_state[3], m_state[4], m_state[5]);
}

Quaternion State::orientation() const {
    return Quaternion(m_state[6], m_state[7], m_state[8], m_state[9]);
}

Eigen::Vector3d State::ang_velocity() const {
    return Eigen::Vector3d(m_state[10], m_state[11], m_state[12]);
}

State& State::operator+=(const State &state) {
    for (int i = 0; i < 13; ++i)
        m_state[i] += state.state()[i];
    return *this;
}

State& State::operator*=(const double a) {
    for (int i = 0; i < 13; ++i)
        m_state[i] *= a;
    return *this;
}

State operator+(const State& state1, const State& state2) {
    std::array<double, 13> state = {};
    for (int i = 0; i < 13; ++i)
        state[i] = state1.m_state[i] + state2.m_state[i];
    return State(state1.m_et, state1.m_frame, state);
}

State operator/(const State& state1, const State& state2) {
    std::array<double, 13> state = {};
    for (int i = 0; i < 13; ++i)
        state[i] = state1.m_state[i] / state2.m_state[i];
    return State(state1.m_et, state1.m_frame, state);
}

State operator*(const double a, const State& state) {
    return State(state.m_et, state.m_frame, a * state.position(), a * state.velocity(), a * state.orientation(), a * state.ang_velocity());
}


State abs(const State& state) {
    return State(state.m_et, state.m_frame, state.position().array().abs(), state.velocity().array().abs(), state.orientation().abs(), state.ang_velocity().array().abs());
}

std::ostream& operator<<( std::ostream &out , const State& state) {
    std::stringstream ss;
    ss << state.et();
    ss << std::fixed << std::setprecision(6);
    for (int i = 0; i < 13; ++i) {
        ss << ", ";
        ss << state.state()[i];
    }
    out << ss.rdbuf();
    return out;

    //out << state.et() << '\n' << state.position() << '\n' << state.velocity() << '\n' << state.orientation() << '\n' << state.ang_velocity();
    //return out;
}

std::ofstream& operator<<(std::ofstream &out, const State& state) {
    std::stringstream ss;
    ss << state.m_et;
    ss << std::fixed << std::setprecision(6);
    for (int i = 0; i < 13; ++i) {
        ss << ", ";
        ss << state.m_state[i];
    }
    out << ss.rdbuf();
    return out;
}