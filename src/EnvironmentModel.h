//
// Created by Sylvain Renevey on 1/26/18.
//

#ifndef PROPAGATION_ENVIRONMENTMODEL_H
#define PROPAGATION_ENVIRONMENTMODEL_H


#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <fstream>
#include <string>
#include <array>
#include <memory>
#include <numeric>
#include <iterator>
#include <iostream>
#include <vector>
#include "constants.h"
#include "AtmModel.h"
#include "Atmod1.h"
#include "Body.h"
#include "StateVector.h"
#include "ReferenceFrame.h"
#include "algebra/transformations.h"

extern "C" {
    #include "SpiceUsr.h"
};


/** Environment model used for the simulation. */
class EnvironmentModel {
public:
    EnvironmentModel();

    /**
     * @param central_body          Central body
     * @param gp_degree             Degree of expansion of the geopotential
     * @param geopot_model_path     Path to the geopotential model
     * @param third_bodies          List of celestial bodies that are accounted for
     * @param atm_model             Atmospheric model used for drag computation. Current options are None, Exponential, and EarthGRAM.
     * @param earthgram_path        Path to the EarthGRAM 2016 root directory
     * @param epoch                 Initial epoch with format "YYYY-MM-DD hh:mm:ss"
     * @param srp_flag              Flag indicating if the solar radiation pressure is accounted for
     * @param mag_flag              Flag indicating if the magnetic perturbations are accounted for
     */
    EnvironmentModel(
            Body central_body,
            int gp_degree = 0,
            const std::string& geopot_model_path = "",
            const std::vector<Body>& third_bodies = std::vector<Body> {},
            AtmModel atm_model = AtmModel::None,
            const std::string& earthgram_path = "",
            const std::string& epoch = "",
            bool srp_flag = false,
            bool mag_flag = false);

    [[nodiscard]] Body central_body() const;
    [[nodiscard]] std::vector<Body> third_bodies() const;
    [[nodiscard]] int gp_degree() const;
    [[nodiscard]] bool is_drag() const;
    [[nodiscard]] bool is_srp_flag() const;
    [[nodiscard]] bool mag_flag() const;


    /** Returns the normalized C coefficient for the given degree and order. */
    [[nodiscard]] double c_coeff(int degree, int order) const;
    /** Returns the normalized S coefficient for the given degree and order. */
    [[nodiscard]] double s_coeff(int degree, int order) const;


    /** Retrieves the density from the atmosphere model specified in the EnvironmentModel.
     *
     *  @param state                State vector of the spacecraft.
     *  @param elapsed_time         Time elapsed since the beginning of the integration (sec)
     *  @param et                   Ephemeris time.
     *  @return                     Atmosphere density in kg/km<sup>3</sup>.
     */
    [[nodiscard]] double atm_density(const StateVector& state, double elapsed_time, double et) const;

    /** Computes the V and W coefficients required to compute the acceleration.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et           Ephemeris time.
     */
    [[nodiscard]] std::vector<std::array<double, 2>> geopotential_harmonic_coeff(const StateVector& state, double et) const;


    /** Computes the position vector from the central body to the given body.
     *
     *  @param body         Body to retrieve the position vector for.
     *  @param et			Ephemeris time.
     *  @param frame        Reference frame in which the vector is expressed.
     */
    [[nodiscard]] Vector3d<Dimension::Distance> body_vector(const Body& body, double et, ReferenceFrame frame) const;

    /** Computes the position vector from the sun to the spacecraft.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et			Ephemeris time.
     *  @return             Position vector expressed in the reference frame of the state vector.
     */
    [[nodiscard]] Vector3d<Dimension::Distance> sun_spacecraft_vector(const StateVector& state, double et) const;

    /** Computes the shadow condition on the spacecraft.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et			Ephemeris time.
     *  @return             Decimal number with value 0 if no shadow, 1 if in umbra, and between 0 and 1 if in penumbra.
     */
    [[nodiscard]] double in_shadow(const StateVector& state, double et) const;


    /** Computes the Earth magnetic field at the given position. The computation is based on a dipole model.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et			Ephemeris time.
     *  @return             Local magnetic field vector in Tesla, expressed in the reference frame of the state vector.
     */
    [[nodiscard]] Vector3d<double> magnetic_field(const StateVector& state, double et) const;

private:
    Body m_central_body; /*!< Central body. */
    int m_gp_degree; /*!< Degree of expansion of the Earth geopotential. */
    std::unique_ptr<double[]> m_c_coeffs; /*!< Array containing the C normalized coefficients used to compute the geopotential effect. */
    std::unique_ptr<double[]> m_s_coeffs; /*!< Array containing the S normalized coefficients used to compute the geopotential effect. */
    std::vector<Body> m_third_bodies; /*!< Vector containing the celestial bodies accounted for in the computation of third body perturbations. */
    AtmModel m_atm_model; /*!< Variant of the atmospheric model to be used. Current options are None, Exponential, and EarthGRAM. */
    std::unique_ptr<Atm1> m_earth_gram_atm_model; /*!< Pointer to an instance of EarthGRAM 2016 model's Atm1 object. */
    double m_exp_atm_model[28][4]; /*!< Array containing the exponential atmospheric model */
    bool m_srp_flag; /*!< Flag indicating if the perturbations due to solar pressure radiation are active. */
    bool m_mag_flag; /*!< Flag indicating if the magnetic perturbations are active. */

    // Methods
    void load_coefficients(int degree, std::string model_file_name);
    void init_exp_model();
    int index(int degree, int order) const;
};



#endif //PROPAGATION_ENVIRONMENTMODEL_H
