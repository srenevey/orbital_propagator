//
// Created by Sylvain Renevey on 8/21/18.
//

#include "Spacecraft.h"

Spacecraft::Spacecraft():
        m_name("sc"),
        m_state(StateVector()),
        m_wet_mass(0.0),
        m_drag_coefficient(0.0),
        m_specular_reflection({}),
        m_diffuse_reflection({}),
        m_sensors({})
{}

Spacecraft::Spacecraft(
        const std::string& name,
        const StateVector& initial_state,
        double mass,
        double drag_coefficient,
        std::vector<double> specular_reflection,
        std::vector<double> diffuse_reflection,
        const Matrix3d& inertia_matrix,
        const std::vector<Vector3d<double>>& face_normals,
        std::vector<Dimension::Area> face_areas,
        const std::vector<Vector3d<Dimension::Distance>>& face_cop_positions,
        std::vector<std::shared_ptr<Sensor>> sensors
):
        m_name(name),
        m_state(initial_state),
        m_wet_mass(mass),
        m_drag_coefficient(drag_coefficient),
        m_specular_reflection(std::move(specular_reflection)),
        m_diffuse_reflection(std::move(diffuse_reflection)),
        m_inertia_matrix(inertia_matrix),
        m_face_normals(face_normals),
        m_face_areas(std::move(face_areas)),
        m_face_cop_positions(face_cop_positions),
        m_sensors(sensors)
{}

std::string Spacecraft::name() const {
    return m_name;
}

StateVector Spacecraft::state() const {
    return m_state;
}


void Spacecraft::set_state(const StateVector& state) {
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

Matrix3d Spacecraft::inertia_matrix() const {
    return m_inertia_matrix;
}

std::vector<Vector3d<double>> Spacecraft::face_normals() const {
    return m_face_normals;
}

std::vector<Dimension::Area> Spacecraft::face_areas() const {
    return m_face_areas;
}

std::vector<Vector3d<Dimension::Distance>> Spacecraft::face_cop_positions() const {
    return m_face_cop_positions;
}

std::vector<std::shared_ptr<Sensor>> Spacecraft::sensors() const {
    return m_sensors;
}

void Spacecraft::operator() (const StateVector& state, const double et) {
    m_state = state;
    m_state.orientation().normalize();
    m_state.set_time(et);

    std::string filename = m_name + ".dat";
    {
        std::ofstream f(filename, std::ios::app);
        f << m_state << '\n';
    }
}