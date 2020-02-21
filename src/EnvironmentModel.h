//
// Created by Sylvain Renevey on 1/26/18.
//

#ifndef PROPAGATION_ENVIRONMENTMODEL_H
#define PROPAGATION_ENVIRONMENTMODEL_H


#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <array>
#include <numeric>
#include <iterator>
#include <iostream>
#include "constants.h"
#include "Atmod1.h"
#include "Spacecraft.h"

extern "C" {
    #include "SpiceUsr.h"
};

/*! \class EnvironmentModel
    \brief This class represents the environment model used for the simulation.

    Acceleration perturbations due to several disturbing effects can be turned on and off.
*/
class EnvironmentModel {

protected:
    double              mu_;                       /*!< Gravitational parameter of the central body. */
    std::string         central_body_;             /*!< String containing the name of the central body. */
    std::array<int,10>  third_body_flags_;         /*!< Array of flags indicating if the perturbation due to a third body is active. */
    bool                third_body_flag_;          /*!< Flag indicating if the perturbations due to a third body are active. */
    int                 gp_degree_;                /*!< Degree of expansion of the Earth geopotential. 0 if not set */
    Eigen::MatrixXd     CS_coeffs_;                /*!< Matrix containing the C and S normalized coefficients used to compute the geopotential effect. */
    bool                drag_flag_;                /*!< Flag indicating if the perturbations due to drag are active. */
    std::string         atm_model_;                /*!< String containing the name of the atmospheric model to be used. Current options are exp (default) and EarthGRAM. */
    Atm1*               earth_gram_atm_model_;     /*!< Pointer to an instance of EarthGRAM 2016 model's Atm1 object. */
    double              exp_atm_model_[28][4];     /*!< Array containing the exponential atmospheric model */
    bool                srp_flag_;                 /*!< Flag indicating if the perturbations due to solar pressure radiation are active. */

public:
    EnvironmentModel();

    /*! \fn explicit EnvironmentModel(std::string earthgram_path, std::string central_body, int gp_degree = 0, std::string geopot_model_path = "", std::array<int, 10> third_body_flags = std::array<int, 10> {0,0,0,0,0,0,0,0,0,0}, bool drag_flag = false, std::string atm_model = "exp", bool srp_flag = false, std::string epoch = "")
     *
     * @param earthgram_path        String containing the path to EarthGRAM model
     * @param central_body          Name of the central body
     * @param gp_degree             Degree of expansion of the geopotential
     * @param geopot_model_path     Path to the geopotential model
     * @param third_body_flags      Flags indicating which third body are accounted for
     * @param drag_flag             Flag indicating if the drag is accounted for
     * @param atm_model             Atmospheric model used for drag computation. Can be "exp" or "EarthGRAM"
     * @param srp_flag              Flag indicating if the solar radiation pressure is accounted for
     * @param epoch                 Initial epoch with format "YYYY-MM-DD hh:mm:ss"
     */
    explicit EnvironmentModel(
            std::string earthgram_path,
            std::string central_body,
            int gp_degree = 0,
            std::string geopot_model_path = "",
            std::array<int, 10> third_body_flags = std::array<int, 10> {0,0,0,0,0,0,0,0,0,0},
            bool drag_flag = false,
            std::string atm_model = "exp",
            bool srp_flag = false,
            std::string epoch = "");
    virtual ~EnvironmentModel();

    double mu() const;
    std::string central_body() const;
    std::array<int, 10> third_body_flags() const;
    bool is_third_body_flag() const;
    int gp_degree() const;
    const Eigen::MatrixXd &cs_coeffs() const;
    bool is_drag_flag() const;
    bool is_srp_flag() const;


    /*! \fn void getAtm_parameters(const SpiceDouble& inertial_position, double& density, double elapsed_time, SpiceDouble et) const
     *  \brief Retrieve the atmosphere density from the model (kg/m<sup>3</sup>).
     *
     *  \param inertial_position    Inertial position of the spacecraft.
     *  \param density              Variable where the density will be stored. Units are kg/m<sup>3</sup>.
     *  \param elapsed_time         Time elapsed since the beginning of the integration (sec)
     *  \param et                   Ephemeris time.
     */
    void get_atm_parameters(const SpiceDouble* inertial_position, double& density, double elapsed_time, SpiceDouble et) const;


    /*! \fn void print_parameters() const;
     *  \brief Display the parameters of the environment model.
     */
    void print_parameters() const;


    /*! \fn Eigen::Vector3d drag(const Spacecraft &sc, const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, double elapsed_time, const EnvironmentModel &env_model)
     *  \brief Compute the perturbations induced by the drag.
     *
     * This function computes the perturbations produced by the atmospheric drag. Based on the settings of the environment model, the atmospheric density
     * is retrieved either from an exponential model or from EarthGRAM2016.
     *
     *  \param sc				Instance of the Spacecraft object containing the properties of the spacecraft.
     *  \param state			State vector of the spacecraft expressed in the J2000 frame (ECI).
     *  \param et				Ephemeris time.
     *  \param elapsed_time     Time elapsed since the beginning of the integration (sec).
     */
    Eigen::Vector3d drag(const Spacecraft &sc, const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, double elapsed_time) const;



    /*! \fn Eigen::Vector3d geopotential(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const EnvironmentModel &env_model)
     *  \brief Compute the gravitational perturbation induced by the nonspherical geopotential.
     *
     * This function computes the perturbations produced by the nonspherical geopotential.
     * The algorithm described in Montenbruck, O., Gill, E.,
     * <em>Satellites Orbits</em>, Springer-Verlag Berlin Heidelberg, 2000, is used. This method uses the EGM2008 model.
     *
     *  \param state        State vector of the spacecraft expressed in the J2000 frame (ECI).
     *  \param et           Ephemeris time.
     */
    Eigen::Vector3d geopotential(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et) const;


    /*! \fn Eigen::Vector3d nbody(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const EnvironmentModel &env_model)
     *  \brief Compute the gravitational perturbations produced by celestial bodies.
     *
     * This function computes the perturbations produced by all other planets and/or moons taken into account.
     * The algorithm is taken from Vallado, D.A., <em>Fundamentals of Astrodynamics and Applications</em>, Fourth Edition, Space Technology Library, 2013.
     *
     *  \param state        State vector of the spacecraft expressed in the J2000 or Ecliptic J2000 frame.
     *  \param et			Ephemeris time
     *
     *  \return Acceleration resulting from the third-body perturbations expressed in the J2000
     *  frame if the central body is the Earth and in the Ecliptic J2000 otherwise.
     */
    Eigen::Vector3d nbody(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et) const;


    /*! \fn Eigen::Vector3d solar_radiation_pressure(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const EnvironmentModel &env_model, const Spacecraft &sc)
     *  \brief Compute the perturbations induced by the solar radiation pressure.
     *
     * This function computes the orbital perturbations induced by the solar radiation pressure. It is assumed that the
     * surface normal points in the direction of the Sun.
     *
     *  \param state			State vector of the spacecraft expressed in the J2000 or Ecliptic J2000 frame.
     *  \param et				Ephemeris time.
     *  \param sc				Instance of the Spacecraft object containing the properties of the spacecraft.
     */
    Eigen::Vector3d solar_radiation_pressure(const Eigen::Ref<const Eigen::VectorXd>& state, SpiceDouble et, const Spacecraft &sc) const;

private:
    Eigen::MatrixXd load_coefficients(int degree, std::string model_file_name);
    void init_exp_model();

    /*! \fn Eigen::Vector3d perturbation_from_body(const Eigen::Ref<const Eigen::VectorXd>& state, const std::string& third_body, const std::string& central_body, SpiceDouble et)
     *  \brief Compute the perturbation induced by a single body.
     *
     *  This function computes the perturbing acceleration produced by a single body on the spacecraft at a given date.
     *  JPL's ephemerides are used to retrieve the body's state vector.
     *
     *  \param state        State vector of the spacecraft expressed in the J2000 or Ecliptic J2000 frame.
     *  \param third_body   Name of the perturbing body.
     *  \param central_body Celestial body to which the spacecraft is bound.
     *  \param et			Ephemeris time.
     */
    Eigen::Vector3d perturbation_from_body(const Eigen::Ref<const Eigen::VectorXd>& state, const std::string& third_body, const std::string& central_body, SpiceDouble et) const;
};


#endif //PROPAGATION_ENVIRONMENTMODEL_H
