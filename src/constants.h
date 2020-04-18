//
// Created by Sylvain Renevey on 2/27/19.
//

#ifndef CPP_PROPAGATOR_CONSTANTS_H
#define CPP_PROPAGATOR_CONSTANTS_H

namespace constants {

	// Gravitational parameters

	/** \defgroup GRAV_PARAM Gravitational parameters
	 *
	 * Units: km<sup>3</sup> / s<sup>2</sup>.<br>
	 * The values are taken from Folkner, W.M. et al., <em>The Planetary and Lunar Ephemerides DE430 and DE431</em>, IPN Progress Report 42-196, 2014.
	 *
	 * @{
	*/
	constexpr static double MU_SUN = 132712440041.939400;
	constexpr static double MU_MERCURY = 22031.780000;
	constexpr static double MU_VENUS = 324858.592000;
	constexpr static double MU_EARTH = 398600.435436;
	constexpr static double MU_MOON = 4902.800066;
	constexpr static double MU_MARS = 42828.375214;
	constexpr static double MU_JUPITER = 126712764.800000;
	constexpr static double MU_SATURN = 37940585.200000;
	constexpr static double MU_URANUS = 5794548.600000;
	constexpr static double MU_NEPTUNE = 6836527.100580;
	/** @} */


	/** \defgroup EQ_RADIUS Equatorial radii
	 *
	 * Units: km.<br>
	 * The values are taken from Archinal, B.A. et al., <em>Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009</em>, Celestial Mechanics and Dynamical Astronomy, Vol.109-2, Springer, 2011.
	 *
	 * @{
	*/
	constexpr static double R_SUN = 696000.0;
	constexpr static double R_MERCURY = 2439.7;
	constexpr static double R_VENUS = 6051.8;
	constexpr static double R_EARTH = 6378.1366;
	constexpr static double R_MOON = 1737.4;
	constexpr static double R_MARS = 3396.19;
	constexpr static double R_JUPITER = 71492;
	constexpr static double R_SATURN = 60268;
	constexpr static double R_URANUS = 25559;
	constexpr static double R_NEPTUNE = 24764;
	/** @} */


	/** \defgroup MISC Misc
	 * @{
	 */
	/** Obliquity of the ecliptic (rad) */
	constexpr static double EARTH_OBLIQUITY = 0.409106;
	/** Earth polar radius (km) */
	constexpr static double EARTH_POLAR_RADIUS = 6356.7519;
	/** Earth angular velocity (rad) */
	constexpr static double EARTH_ANGULAR_VELOCITY = 7.272205216643040E-05;
	/** Astronomical Unit (km) */
	constexpr static double AU = 149597870.700;
	/** Pi value **/
	constexpr static double PI = 3.14159265359;
	/** Factor to convert from degrees to radians **/
	constexpr static double DEG_TO_RAD = PI / 180.0;
	/** Factor to convert from radians to degrees **/
	constexpr static double RAD_TO_DEG = 180.0 / PI;
	/** Solar irradiance (W/km<sup>2</sup>) (Kopp, G.; Lean, J. L., <em>A new, lower value of total solar irradiance: Evidence and climate significance</em>. Geophysical Research Letters, 38, 2011, doi:10.1029/2010GL045777.) */
	constexpr static double SOLAR_IRRADIANCE = 1360.8E6;
	/** Speed of light (km/s) **/
	constexpr static double SPEED_OF_LIGHT = 299792.458;
	/** Solar pressure (N/km<sup>2</sup>) **/
	constexpr static double SOLAR_PRESSURE = SOLAR_IRRADIANCE / SPEED_OF_LIGHT;
	/** @} */


	/** \defgroup ZONAL_COEFF Zonal coefficients
	 *
	 * Zonal coefficients used for the expansion of the geopotential field up to degree 6. <br>
	 * The values are taken from Markley, F.L, Crassidis, J.L., <em>Fundamentals of Spacecraft Attitude Determination and Control</em>, Space Technology Library, Springer, New York, 2014.
	 *
	 * @{
	 */
	constexpr static double EARTH_J2 = 1.08262668355E-03;
	constexpr static double EARTH_J3 = -2.53265648533E-06;
	constexpr static double EARTH_J4 = -1.61962159137E-06;
	constexpr static double EARTH_J5 = -2.27296082869E-07;
	constexpr static double EARTH_J6 = 5.40681239107E-07;
	/** @} */


	/** \defgroup MU_EARTH_EGM08 Earth Gravitational Parameter EGM08
	 *  \brief Gravitational parameter in km<sup>3</sup> / s<sup>2</sup>.
	 *
	 *  This value is specific to the EGM2008 geopotential model.
	 */
	constexpr static double MU_EARTH_EGM08 = 398600.4415;


	/** \defgroup R_EARTH_EGM08 Earth Radius EGM08
	 *  \brief Value of the Earth's radius in km.
	 *
	 *  This value is specific to the EGM2008 geopotential model.
	 */
	constexpr static double R_EARTH_EGM08 = 6378.1363;


    /** \defgroup IGRF_COEFF Magnetic Field Coefficients
     *
     * Coefficients used for the Earth's magnetic field modeled as a magnetic dipole.
     *
     * @{
     */
    constexpr static double IGRF_2020_RADIUS = 6371.2; // km
	constexpr static double IGRF_2020_G01 = -29404.8E-9; // T
	constexpr static double IGRF_2020_G11 = -1450.9E-9; // T
	constexpr static double IGRF_2020_H11 = 4652.5E-9; // T
	/** @} */
}

#endif //CPP_PROPAGATOR_CONSTANTS_H