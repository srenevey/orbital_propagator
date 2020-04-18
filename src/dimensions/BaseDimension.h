//
// Created by Sylvain  on 4/1/20.
//

#ifndef ORBITAL_PROPAGATOR_BASEDIMENSION_H
#define ORBITAL_PROPAGATOR_BASEDIMENSION_H

/** Contains definitions of physical dimensions. */
namespace Dimension {
    /** Base class used to define physical dimensions. */
    class BaseDimension {
    public:
        BaseDimension(double d);
        virtual ~BaseDimension();

        /** Returns the underlying data of the dimension. */
        [[nodiscard]] double data() const;

        /** Converts the dimension to double. */
        operator double() const;

        //const BaseDimension& operator-();
        BaseDimension operator-();
        BaseDimension& operator+=(const BaseDimension& a);
        friend BaseDimension operator+(const BaseDimension& a, const BaseDimension& b);
        friend BaseDimension operator+(const BaseDimension& a, const double& b);
        friend BaseDimension operator-(const BaseDimension& a, const BaseDimension& b);
        BaseDimension& operator*=(const double& a);
        friend BaseDimension operator*(const double& a, const BaseDimension& b);
        friend BaseDimension operator*(const int& a, const BaseDimension& b);
        BaseDimension& operator/=(const double& a);
        friend BaseDimension operator/(const BaseDimension& a, const double& b);

    protected:
        double m_data;
    };
}

/** Namespace containing the definitions of literals to specify the unit of each dimension. */
namespace Unit {}

#endif //ORBITAL_PROPAGATOR_BASEDIMENSION_H
