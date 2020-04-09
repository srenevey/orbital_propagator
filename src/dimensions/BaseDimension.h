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
        explicit BaseDimension(double d);
        virtual ~BaseDimension();

        /** Returns the underlying data of the dimension. */
        [[nodiscard]] double data() const;

        /** Converts the dimension to double. */
        operator double() const;

        friend BaseDimension operator+(const BaseDimension& a, const BaseDimension& b);
        friend BaseDimension operator-(const BaseDimension& a, const BaseDimension& b);
        friend BaseDimension operator*(const double& a, const BaseDimension& b);
        friend BaseDimension operator/(const BaseDimension& a, const double& b);

    protected:
        double m_data;
    };
}

/*! @brief Namespace containing the definitions of literals to specify the unit of each dimension.
 */
namespace Unit {}

#endif //ORBITAL_PROPAGATOR_BASEDIMENSION_H
