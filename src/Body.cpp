//
// Created by Sylvain  on 3/19/20.
//

#include "Body.h"

BodyContainer::BodyContainer(const double mu, const double radius, const Body name, const int naif_id): m_mu(mu), m_radius(radius), m_name(name), m_naif_id(naif_id) {}

Body BodyContainer::name() const {
    return m_name;
}

double BodyContainer::mu() const {
    return m_mu;
}

double BodyContainer::radius() const {
    return m_radius;
}

int BodyContainer::naif_id() const {
    return m_naif_id;
}

BodyContainer* BodyContainer::create_body(const Body body) {
    BodyContainer* body_container;
    switch (body) {
        case Body::Sun:
            body_container = new Sun();
            break;
        case Body::Mercury:
            body_container = new Mercury();
            break;
        case Body::Venus:
            body_container = new Venus();
            break;
        case Body::Earth:
            body_container = new Earth();
            break;
        case Body::Moon:
            body_container = new Moon();
            break;
        case Body::Mars:
            body_container = new Mars();
            break;
        case Body::Jupiter:
            body_container = new Jupiter();
            break;
        case Body::Saturn:
            body_container = new Saturn();
            break;
        case Body::Uranus:
            body_container = new Uranus();
            break;
        case Body::Neptune:
            body_container = new Neptune();
            break;
    }
    return body_container;
}

Sun::Sun(): BodyContainer(constants::MU_SUN, constants::R_SUN, Body::Sun, 10) {}
Mercury::Mercury(): BodyContainer(constants::MU_MERCURY, constants::R_MERCURY, Body::Mercury, 199) {}
Venus::Venus(): BodyContainer(constants::MU_VENUS, constants::R_VENUS, Body::Venus, 299) {}
Earth::Earth(): BodyContainer(constants::MU_EARTH, constants::R_EARTH, Body::Earth, 399), m_obliquity(constants::EARTH_OBLIQUITY), m_polar_radius(constants::EARTH_POLAR_RADIUS), m_angular_velocity(constants::EARTH_ANGULAR_VELOCITY) {}
Moon::Moon(): BodyContainer(constants::MU_MOON, constants::R_MOON, Body::Moon, 301) {}
Mars::Mars(): BodyContainer(constants::MU_MARS, constants::R_MARS, Body::Mars, 499) {}
Jupiter::Jupiter(): BodyContainer(constants::MU_JUPITER, constants::R_JUPITER, Body::Jupiter, 599) {}
Saturn::Saturn(): BodyContainer(constants::MU_SATURN, constants::R_SATURN, Body::Saturn, 699) {}
Uranus::Uranus(): BodyContainer(constants::MU_URANUS, constants::R_URANUS, Body::Uranus, 799) {}
Neptune::Neptune(): BodyContainer(constants::MU_NEPTUNE, constants::R_NEPTUNE, Body::Neptune, 899) {}