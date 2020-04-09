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
#include <vector>
#include "constants.h"
#include "AtmModel.h"
#include "State.h"
#include "Atmod1.h"
#include "Body.h"

extern "C" {
    #include "SpiceUsr.h"
};


/** Environment model used for the simulation. */
class EnvironmentModel {
    BodyContainer* m_central_body; /*!< String containing the name of the central body. */
    int m_gp_degree; /*!< Degree of expansion of the Earth geopotential. 0 if not set */
    Eigen::MatrixXd m_cs_coeffs; /*!< Matrix containing the C and S normalized coefficients used to compute the geopotential effect. */
    std::vector<BodyContainer*> m_third_body; /*!< Vector containing the celestial bodies accounted for in the computation of third body perturbations. */
    AtmModel m_atm_model; /*!< String containing the name of the atmospheric model to be used. Current options are exp (default) and EarthGRAM. */
    Atm1* m_earth_gram_atm_model; /*!< Pointer to an instance of EarthGRAM 2016 model's Atm1 object. */
    double m_exp_atm_model[28][4]; /*!< Array containing the exponential atmospheric model */
    bool m_srp_flag; /*!< Flag indicating if the perturbations due to solar pressure radiation are active. */

public:
    EnvironmentModel();

    /**
     * @param central_body          Body the spacecraft is orbiting
     * @param gp_degree             Degree of expansion of the geopotential
     * @param geopot_model_path     Path to the geopotential model
     * @param third_body            List of bodies that are accounted for
     * @param atm_model             Atmospheric model used for drag computation
     * @param earthgram_path        Path to EarthGRAM 2016
     * @param srp_flag              Flag indicating if the solar radiation pressure is accounted for
     * @param epoch                 Initial epoch with format "YYYY-MM-DD hh:mm:ss"
     */
    explicit EnvironmentModel(
            Body central_body,
            int gp_degree = 0,
            const std::string& geopot_model_path = "",
            const std::vector<Body>& third_body = std::vector<Body> {},
            AtmModel atm_model = AtmModel::None,
            const std::string& earthgram_path = "",
            bool srp_flag = false,
            const std::string& epoch = "");
    virtual ~EnvironmentModel();

    [[nodiscard]] BodyContainer* central_body() const;
    [[nodiscard]] std::vector<BodyContainer*> third_body() const;
    [[nodiscard]] int gp_degree() const;
    [[nodiscard]] const Eigen::MatrixXd &cs_coeffs() const;
    [[nodiscard]] bool is_drag() const;
    [[nodiscard]] bool is_srp_flag() const;


    /** Retrieves the density from the atmosphere model specified in the EnvironmentModel.
     *
     *  @param state                State vector of the spacecraft.
     *  @param elapsed_time         Time elapsed since the beginning of the integration (sec)
     *  @param et                   Ephemeris time.
     *  @return                     Atmosphere density in kg/km<sup>3</sup>.
     */
    [[nodiscard]] double atm_density(const State& state, double elapsed_time, SpiceDouble et) const;

    /** Computes the V and W coefficients required to compute the acceleration.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et           Ephemeris time.
     */
    [[nodiscard]] Eigen::MatrixXd geopotential_harmonic_coeff(const State& state, SpiceDouble et) const;

    /** Computes the position from the central body to the third body.
     *
     *  @param body         Third body.
     *  @param et			Ephemeris time.
     */
    Eigen::Vector3d body_vector(const BodyContainer* body, SpiceDouble et) const;

    /** Computes the vector from the sun to the spacecraft in the J2000 or ECLIPJ2000 reference frame.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et			Ephemeris time.
     */
    [[nodiscard]] Eigen::Vector3d sun_spacecraft_vector(const State& state, SpiceDouble et) const;


    /** Computes the shadow condition on the spacecraft.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et			Ephemeris time.
     *  @return             Decimal number with value 0 if no shadow, 1 if in umbra, and between 0 and 1 if in penumbra.
     */
    [[nodiscard]] double in_shadow(const State& state, SpiceDouble et) const;


    /** Computes the Earth magnetic field at the given position. The computation is based on a dipole model.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et			Ephemeris time.
     *  @return             Local magnetic field vector expressed in the J2000 reference frame (in Tesla).
     */
    [[nodiscard]] Eigen::Vector3d magnetic_field(const State& state, SpiceDouble et) const;

private:
    Eigen::MatrixXd load_coefficients(int degree, std::string model_file_name);
    void init_exp_model();
};



#endif //PROPAGATION_ENVIRONMENTMODEL_H
