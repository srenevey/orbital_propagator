//
// Created by Sylvain Renevey on 8/21/18.
//

#include "Spacecraft.h"

Spacecraft::Spacecraft():
    mass_(0.0),
	drag_coefficient_(0.0),
    reference_area_(0.0),
    reflectivity_(0.0),
    sun_exposed_area_(0.0)
{}

Spacecraft::Spacecraft(double mass, double drag_coefficient, double reference_area, double reflectivity, double sun_exposed_area):
    mass_(mass),
	drag_coefficient_(drag_coefficient),
    reference_area_(reference_area),
    reflectivity_(reflectivity),
    sun_exposed_area_(sun_exposed_area)
{}

double Spacecraft::getMass() const {
    return mass_;
}

double Spacecraft::getDragCoefficient() const {
    return drag_coefficient_;
}

double Spacecraft::getReferenceArea() const {
    return reference_area_;
}

double Spacecraft::getReflectivity() const {
    return reflectivity_;
}

double Spacecraft::getSunExposedArea() const {
    return sun_exposed_area_;
}