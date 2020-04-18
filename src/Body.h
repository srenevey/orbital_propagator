//
// Created by Sylvain  on 4/18/20.
//

#ifndef ORBITAL_PROPAGATOR_BODY_H
#define ORBITAL_PROPAGATOR_BODY_H

#include "constants.h"

/** Wrapper to easily retrieve the physical constants of the different celestial bodies. */
class Body {
public:
    enum CelestialBody {
        Sun,
        Mercury,
        Venus,
        Earth,
        Moon,
        Mars,
        Jupiter,
        Saturn,
        Uranus,
        Neptune
    };

    Body()= default;
    constexpr Body(CelestialBody body) : m_body(body) { }
    constexpr bool operator==(Body a) const { return m_body == a.m_body; }
    constexpr bool operator!=(Body a) const { return m_body != a.m_body; }
    explicit operator bool() = delete;

    /** Returns the gravitational parameter of the body in km<sup>3</sup>/s<sup>2</sup>. */
    [[nodiscard]] constexpr double mu() const {
        switch (m_body) {
            case Sun:
                return constants::MU_SUN;
            case Mercury:
                return constants::MU_MERCURY;
            case Venus:
                return constants::MU_VENUS;
            case Earth:
                return constants::MU_EARTH;
            case Moon:
                return constants::MU_MOON;
            case Mars:
                return constants::MU_MARS;
            case Jupiter:
                return constants::MU_JUPITER;
            case Saturn:
                return constants::MU_SATURN;
            case Uranus:
                return constants::MU_URANUS;
            case Neptune:
                return constants::MU_NEPTUNE;
        }
    }

    /** Returns the equatorial radius of the body in km. */
    [[nodiscard]] constexpr double radius() const {
        switch (m_body) {
            case Sun:
                return constants::R_SUN;
            case Mercury:
                return constants::R_MERCURY;
            case Venus:
                return constants::R_VENUS;
            case Earth:
                return constants::R_EARTH;
            case Moon:
                return constants::R_MOON;
            case Mars:
                return constants::R_MARS;
            case Jupiter:
                return constants::R_JUPITER;
            case Saturn:
                return constants::R_SATURN;
            case Uranus:
                return constants::R_URANUS;
            case Neptune:
                return constants::R_NEPTUNE;
        }
    }

    /** Returns the NAIF identification number. */
    [[nodiscard]] constexpr int naif_id() const {
        switch (m_body) {
            case Sun:
                return 10;
            case Mercury:
                return 199;
            case Venus:
                return 299;
            case Earth:
                return 399;
            case Moon:
                return 301;
            case Mars:
                return 499;
            case Jupiter:
                return 599;
            case Saturn:
                return 699;
            case Uranus:
                return 799;
            case Neptune:
                return 899;
        }
    }

    /** Returns the name of the body as a string. */
    explicit operator std::string() const {
        switch (m_body) {
            case Sun:
                return "Sun";
            case Mercury:
                return "Mercury";
            case Venus:
                return "Venus";
            case Earth:
                return "Earth";
            case Moon:
                return "Moon";
            case Mars:
                return "Mars";
            case Jupiter:
                return "Jupiter";
            case Saturn:
                return "Saturn";
            case Uranus:
                return "Uranus";
            case Neptune:
                return "Neptune";
        }
    }

private:
    CelestialBody m_body;

};


#endif //ORBITAL_PROPAGATOR_BODY_H
