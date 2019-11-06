//
// Created by Sylvain Renevey on 8/21/18.
//

#ifndef CPP_PROPAGATOR_SPACECRAFT_H
#define CPP_PROPAGATOR_SPACECRAFT_H

/*! \class Spacecraft
    \brief This class contains some properties of the spacecraft.
*/

class Spacecraft {
protected:
    double mass_;					/*!< Mass (kg). */
    double drag_coefficient_;       /*!< Drag coefficient. */
    double reference_area_;			/*!< Reference area (m<sup>2</sup>). */
    double reflectivity_;           /*!< Coefficient of reflectivity. */
    double sun_exposed_area_;       /*!< Area exposed to the Sun (m<sup>2</sup>). */

public:
    Spacecraft();

    /*! \fn Spacecraft(double mass, double drag_coefficient, double reference_area, double reflectivity, double sun_exposed_area)
     *  \brief Constructor.
     *
     *  \param mass                 Mass of the spacecraft (kg)
     *  \param drag_coefficient     Drag coefficient of the spacecraft
     *  \param reference_area       Reference area of the spacecraft (m<sup>2</sup>)
     *  \param reflectivity         Coefficient of reflectivity
     *  \param sun_exposed_area     Area exposed to the Sun (m<sup>2</sup>)
     */
    Spacecraft(double mass, double drag_coefficient, double reference_area, double reflectivity, double sun_exposed_area);

    /*! \fn double getMass() const
     *  \brief Return the mass of the spacecraft.
     */
    double getMass() const;

    /*! \fn double getDragCoefficient() const
     *  \brief Return the drag coefficient of the spacecraft.
     */
    double getDragCoefficient() const;

    /*! \fn double getReferenceArea() const
     *  \brief Return the reference area of the spacecraft.
     */
    double getReferenceArea() const;

    /*! \fn double getReflectivity() const
     * \brief Return the reflectivity of the spacecraft.
     */
    double getReflectivity() const;

    /*! \fn double getSunExposedArea() const
     *  \brief Return the area of the spacecraft which is exposed to the Sun.
     */
    double getSunExposedArea() const;
};


#endif //CPP_PROPAGATOR_SPACECRAFT_H
