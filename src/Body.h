//
// Created by Sylvain  on 3/19/20.
//

#ifndef ORBITAL_PROPAGATOR_BODY_H
#define ORBITAL_PROPAGATOR_BODY_H

#include "constants.h"

/** @file
 */

/** Names of the bodies implemented in the program. */
enum class Body {
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


/** Base class to define celestial bodies. */
class BodyContainer {
protected:
    Body m_name;
    double m_mu;
    double m_radius;
    int m_naif_id;

public:
    BodyContainer(const double mu, const double radius, const Body name, const int naif_id);
    virtual ~BodyContainer() = default;

    Body name() const;
    double mu() const;
    double radius() const;
    static BodyContainer* create_body(const Body body);

    /** Returns the NAIF identification number. */
    int naif_id() const;
};

class Sun: public BodyContainer {
public:
    static constexpr double MU { 132712440041.939400 } ;
    static constexpr double RADIUS { 696000.0 };
    Sun();
};

class Mercury: public BodyContainer {
public:
    Mercury();
};

class Venus: public BodyContainer {
public:
    Venus();
};

class Earth: public BodyContainer {
    double m_obliquity;
    double m_polar_radius;
    double m_angular_velocity;
public:
    Earth();
};

class Moon: public BodyContainer {
public:
    Moon();
};

class Mars: public BodyContainer {
public:
    Mars();
};

class Jupiter: public BodyContainer {
public:
    Jupiter();
};

class Saturn: public BodyContainer {
public:
    Saturn();
};

class Uranus: public BodyContainer {
public:
    Uranus();
};

class Neptune: public BodyContainer {
public:
    Neptune();
};

#endif //ORBITAL_PROPAGATOR_BODY_H
