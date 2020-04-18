//
// Created by Sylvain Renevey on 8/21/18.
//

#ifndef CPP_PROPAGATOR_SPACECRAFT_H
#define CPP_PROPAGATOR_SPACECRAFT_H

#include <string>
#include <fstream>
#include "StateVector.h"
#include "Sensor.h"
#include "algebra/Vector.h"
#include "algebra/Matrix.h"

class Integrator;

/** Contains the state vector and the physical properties of a spacecraft.
 *
 *  A spacecraft is modeled as a collection of N flat plates with given area, surface normal unit vector, and position of the center of pressure w.r.t. the center of mass.
*/
class Spacecraft {
    friend class Integrator;

public:
    Spacecraft();

    /**
     *  @param name                 Name of the spacecraft, used as filename to save simulation output
     *  @param initial_state        Initial state of the spacecraft
     *  @param mass                 Mass of the spacecraft (kg)
     *  @param drag_coefficient     Drag coefficient of the spacecraft
     *  @param specular_reflection  Specular reflection coefficient
     *  @param diffuse_reflection   Diffuse reflection coefficient
     *  @param inertia_matrix       Inertia matrix w.r.t. the center of mass expressed in the body-frame (kg*km<sup>2</sup>)
     *  @param face_normals         Outward normal unit vector of each flat plate
     *  @param face_areas           Area of each flat plate (km<sup>2</sup>)
     *  @param face_cop_positions   Vectors from the center of mass to the center of pressure of each flat plate (km)
     *  @param sensors              Vector containing the sensors of the spacecraft (work in progress)
     */
    Spacecraft(
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
            std::vector<std::shared_ptr<Sensor>> sensors);

    [[nodiscard]] std::string name() const;
    [[nodiscard]] StateVector state() const;
    [[nodiscard]] double mass() const;
    [[nodiscard]] double drag_coefficient() const;
    [[nodiscard]] std::vector<double> specular_reflection_coeff() const;
    [[nodiscard]] std::vector<double> diffuse_reflection_coeff() const;
    [[nodiscard]] Matrix3d inertia_matrix() const;
    [[nodiscard]] std::vector<Vector3d<double>> face_normals() const;
    [[nodiscard]] std::vector<Dimension::Area> face_areas() const;
    [[nodiscard]] std::vector<Vector3d<Dimension::Distance>> face_cop_positions() const;
    [[nodiscard]] std::vector<std::shared_ptr<Sensor>> sensors() const;

    void set_state(const StateVector& state);

    /** This method is called after each integration step to save the state of the spacecraft.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et			Ephemeris time.
     */
    void operator() (const StateVector& state, double et);

private:
    std::string m_name; /*!< Name of the spacecraft, used to save data. */
    double m_wet_mass; /*!< Wet mass (kg). */
    StateVector m_state; /*!< State composed of 13 elements. */
    double m_drag_coefficient; /*!< Drag coefficient. */
    std::vector<double> m_specular_reflection; /*!< Specular reflection coefficient. */
    std::vector<double> m_diffuse_reflection; /*!< Diffuse reflection coefficient. */
    Matrix3d m_inertia_matrix; /*!< Inertia tensor of the spacecraft (kg*km<sup>2</sup>). */
    std::vector<Vector3d<double>> m_face_normals; /*!< Outward normal unit vector of each flat plate. */
    std::vector<Dimension::Area> m_face_areas; /*!< Area of each flat plate (km<sup>2</sup>). */
    std::vector<Vector3d<Dimension::Distance>> m_face_cop_positions; /*!< Position of the center of pressure of each flat plate w.r.t. the center of mass (km). */
    std::vector<std::shared_ptr<Sensor>> m_sensors; /*!< Sensors of the spacecraft. */
};


#endif //CPP_PROPAGATOR_SPACECRAFT_H
