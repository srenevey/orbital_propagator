//
// Created by Sylvain Renevey on 8/21/18.
//

#ifndef CPP_PROPAGATOR_SPACECRAFT_H
#define CPP_PROPAGATOR_SPACECRAFT_H

#include "State.h"
#include "Sensor.h"
#include <string>
#include <fstream>
#include <Eigen/Dense>

class Integrator;
class Sensor;

/** Contains the physical properties of the spacecraft.
 *
 *  A spacecraft is modeled as a collection of N flat plates with given area, surface normal unit vector, and position of the center of pressure w.r.t. the center of mass.
*/
class Spacecraft {
    friend class Integrator;

    std::string m_name; /*!< Name of the spacecraft, used to save data. */
    State m_state; /*!< State composed of 13 elements. */
    double m_wet_mass; /*!< Wet mass (kg). */
    double m_drag_coefficient; /*!< Drag coefficient. */
    std::vector<double> m_specular_reflection; /*!< Specular reflection coefficient. */
    std::vector<double> m_diffuse_reflection; /*!< Diffuse reflection coefficient. */
    Eigen::Matrix3d m_inertia_matrix; /*!< Inertia tensor of the spacecraft (kg*km<sup>2</sup>). */
    std::vector<Eigen::Vector3d> m_face_normals; /*!< Outward normal unit vector of each flat plate. */
    std::vector<double> m_face_areas; /*!< Area of each flat plate (km<sup>2</sup>). */
    std::vector<Eigen::Vector3d> m_face_cop_positions; /*!< Position of the center of pressure of each flat plate w.r.t. the center of mass (km). */
    std::vector<std::shared_ptr<Sensor>> m_sensors; /*!< Sensors of the spacecraft. */

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
            const State& initial_state,
            double mass,
            double drag_coefficient,
            std::vector<double> specular_reflection,
            std::vector<double> diffuse_reflection,
            const std::array<double, 9>& inertia_matrix,
            const std::vector<std::array<double, 3>>& face_normals,
            std::vector<double> face_areas,
            const std::vector<std::array<double, 3>>& face_cop_positions,
            std::vector<std::shared_ptr<Sensor>> sensors);

    [[nodiscard]] std::string name() const;
    [[nodiscard]] State state() const;
    [[nodiscard]] double mass() const;
    [[nodiscard]] double drag_coefficient() const;
    [[nodiscard]] std::vector<double> specular_reflection_coeff() const;
    [[nodiscard]] std::vector<double> diffuse_reflection_coeff() const;
    [[nodiscard]] Eigen::Matrix3d inertia_matrix() const;
    [[nodiscard]] std::vector<Eigen::Vector3d> face_normals() const;
    [[nodiscard]] std::vector<double> face_areas() const;
    [[nodiscard]] std::vector<Eigen::Vector3d> face_cop_positions() const;
    [[nodiscard]] std::vector<std::shared_ptr<Sensor>> sensors() const;

    void set_state(const State& state);

    /** This method is called after each integration step to save the state of the spacecraft.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et			Ephemeris time.
     */
    void operator() (const State& state, double et);
};


#endif //CPP_PROPAGATOR_SPACECRAFT_H
