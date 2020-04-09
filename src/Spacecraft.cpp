//
// Created by Sylvain Renevey on 8/21/18.
//

#include "Spacecraft.h"

#include <utility>

Spacecraft::Spacecraft():
        m_name("sc"),
        m_state(State()),
        m_wet_mass(0.0),
        m_drag_coefficient(0.0),
        m_specular_reflection({}),
        m_diffuse_reflection({}),
        m_sensors({})
{}

Spacecraft::Spacecraft(
        const std::string& name,
        const State& initial_state,
        double mass,
        double drag_coefficient,
        std::vector<double> specular_reflection,
        std::vector<double> diffuse_reflection,
        const std::array<double, 9>& inertia_matrix,
        const std::vector<std::array<double, 3>>& face_normals,
        std::vector<double> face_areas,
        const std::vector<std::array<double, 3>>& face_cop_positions,
        std::vector<std::shared_ptr<Sensor>> sensors
):
        m_name(name),
        m_state(initial_state),
        m_wet_mass(mass),
        m_drag_coefficient(drag_coefficient),
        m_specular_reflection(std::move(specular_reflection)),
        m_diffuse_reflection(std::move(diffuse_reflection)),
        m_face_areas(std::move(face_areas)),
        m_sensors(sensors)
{
    Eigen::Matrix3d inertia;
    inertia << inertia_matrix[0], inertia_matrix[1], inertia_matrix[2],
        inertia_matrix[3], inertia_matrix[4], inertia_matrix[5],
        inertia_matrix[6], inertia_matrix[7], inertia_matrix[8];
    m_inertia_matrix = inertia;

    std::vector<Eigen::Vector3d> normals;
    normals.reserve(face_normals.size());
    for (auto normal: face_normals) {
        normals.emplace_back(Eigen::Vector3d(normal[0], normal[1], normal[2]));
    }
    m_face_normals = normals;

    std::vector<Eigen::Vector3d> cop;
    cop.reserve(face_cop_positions.size());
    for (auto position: face_cop_positions)
        cop.emplace_back(Eigen::Vector3d(position[0], position[1], position[2]));
    m_face_cop_positions = cop;
}

std::string Spacecraft::name() const {
    return m_name;
}

State Spacecraft::state() const {
    return m_state;
}


void Spacecraft::set_state(const State& state) {
    m_state = state;
}

double Spacecraft::mass() const {
    return m_wet_mass;
}

double Spacecraft::drag_coefficient() const {
    return m_drag_coefficient;
}

std::vector<double> Spacecraft::specular_reflection_coeff() const {
    return m_specular_reflection;
}

std::vector<double> Spacecraft::diffuse_reflection_coeff() const {
    return m_diffuse_reflection;
}

Eigen::Matrix3d Spacecraft::inertia_matrix() const {
    return m_inertia_matrix;
}

std::vector<Eigen::Vector3d> Spacecraft::face_normals() const {
    return m_face_normals;
}

std::vector<double> Spacecraft::face_areas() const {
    return m_face_areas;
}

std::vector<Eigen::Vector3d> Spacecraft::face_cop_positions() const {
    return m_face_cop_positions;
}

std::vector<std::shared_ptr<Sensor>> Spacecraft::sensors() const {
    return m_sensors;
}

void Spacecraft::operator() (const State& state, const double et) {
    m_state = state;
    m_state.normalize_quaternion();
    m_state.set_et(et);

    std::string filename = m_name + ".dat";
    {
        std::ofstream f(filename, std::ios::app);
        f << m_state << '\n';
    }
}