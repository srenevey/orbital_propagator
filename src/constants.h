//
// Created by Sylvain Renevey on 2/27/19.
//

#ifndef CPP_PROPAGATOR_CONSTANTS_H
#define CPP_PROPAGATOR_CONSTANTS_H

#include <string>
#include <map>

namespace constants {

	// Gravitational parameters

	/*! \defgroup GRAV_PARAM Gravitational parameters
	 *
	 * Units: km<sup>3</sup> / s<sup>2</sup>.<br>
	 * The values are taken from Folkner, W.M. et al., <em>The Planetary and Lunar Ephemerides DE430 and DE431</em>, IPN Progress Report 42-196, 2014.
	 *
	 * @{
	*/
	constexpr static double kMuSun = 132712440041.939400;
	constexpr static double kMuMercury = 22031.780000;
	constexpr static double kMuVenus = 324858.592000;
	constexpr static double kMuEarth = 398600.435436;
	constexpr static double kMuMoon = 4902.800066;
	constexpr static double kMuMars = 42828.375214;
	constexpr static double kMuJupiter = 126712764.800000;
	constexpr static double kMuSaturn = 37940585.200000;
	constexpr static double kMuUranus = 5794548.600000;
	constexpr static double kMuNeptune = 6836527.100580;
	/** @} */


	/*! \defgroup EQ_RADIUS Equatorial radii
	 *
	 * Units: km.<br>
	 * The values are taken from Archinal, B.A. et al., <em>Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009</em>, Celestial Mechanics and Dynamical Astronomy, Vol.109-2, Springer, 2011.
	 *
	 * @{
	*/
	constexpr static double kSunRadius = 696000.0;
	constexpr static double kMercuryRadius = 2439.7;
	constexpr static double kVenusRadius = 6051.8;
	constexpr static double kEarthRadius = 6378.1366;
	constexpr static double kMoonRadius = 1737.4;
	constexpr static double kMarsRadius = 3396.19;
	constexpr static double kJupiterRadius = 71492;
	constexpr static double kSaturnRadius = 60268;
	constexpr static double kUranusRadius = 25559;
	constexpr static double kNeptuneRadius = 24764;
	/** @} */


	/*! \defgroup MISC Misc
	 * @{
	 */
	/** Obliquity of the ecliptic (rad) */
	constexpr static double kObliquity = 0.409106;
	/** Earth polar radius (km) */
	constexpr static double kEarthPolarRadius = 6356.7519;
	/** Earth angular velocity (rad) */
	constexpr static double kEarthAngularVelocity = 7.272205216643040E-05;
	/** Astronomical Unit (km) */
	constexpr static double kAstronomicalUnit = 149597870.700;
	/** Pi value **/
	constexpr static double kPi = 3.14159265359;
	/** Factor to convert from degrees to radians **/
	constexpr static double kDegToRad = kPi / 180.0;
	/** Factor to convert from radians to degrees **/
	constexpr static double kRadToDeg = 180.0 / kPi;
	/** Solar irradiance (W/km<sup>2</sup>) (Kopp, G.; Lean, J. L., <em>A new, lower value of total solar irradiance: Evidence and climate significance</em>. Geophysical Research Letters, 38, 2011, doi:10.1029/2010GL045777.) */
	constexpr static double kSolarIrradiance = 1360.8E6;
	/** Speed of light (km/s) **/
	constexpr static double kSpeedOfLight = 299792.458;
	/** Solar pressure (N/km<sup>2</sup>) **/
	constexpr static double kSolarPressure = kSolarIrradiance / kSpeedOfLight;
	/** @} */


	/*! \defgroup ZONAL_COEFF Zonal coefficients
	 *
	 * Zonal coefficients used for the expansion of the geopotential field up to degree 6. <br>
	 * The values are taken from Markley, F.L, Crassidis, J.L., <em>Fundamentals of Spacecraft Attitude Determination and Control</em>, Space Technology Library, Springer, New York, 2014.
	 *
	 * @{
	 */
	constexpr static double kEarthJ2 = 1.08262668355E-03;
	constexpr static double kEarthJ3 = -2.53265648533E-06;
	constexpr static double kEarthJ4 = -1.61962159137E-06;
	constexpr static double kEarthJ5 = -2.27296082869E-07;
	constexpr static double kEarthJ6 = 5.40681239107E-07;
	/** @} */


	/*! \defgroup MU_EARTH_EGM08 Earth Gravitational Parameter EGM08
	 *  \brief Gravitational parameter in km<sup>3</sup> / s<sup>2</sup>.
	 *
	 *  This value is specific to the EGM2008 geopotential model.
	 */
	constexpr static double kMuEarthEgm08 = 398600.4415;


	/*! \defgroup R_EARTH_EGM08 Earth Radius EGM08
	 *  \brief Value of the Earth's radius in km.
	 *
	 *  This value is specific to the EGM2008 geopotential model.
	 */
	constexpr static double kEarthRadiusEgm08 = 6378.1363;


	/*!\fn inline double mu(const std::string &body_name)
	 *  \brief Returns the gravitational parameter of the body.
	 *
	 *  \param body_name    Name of the body for which the gravitational parameter is returned.
	 */
	inline double mu(const std::string &body_name) {
		std::map<std::string, double> muLookupTable;
		muLookupTable["Sun"] = kMuSun;
		muLookupTable["Mercury"] = kMuMercury;
		muLookupTable["Venus"] = kMuVenus;
		muLookupTable["Earth"] = kMuEarth;
		muLookupTable["Moon"] = kMuMoon;
		muLookupTable["Mars"] = kMuMars;
		muLookupTable["Jupiter"] = kMuJupiter;
		muLookupTable["Saturn"] = kMuSaturn;
		muLookupTable["Uranus"] = kMuUranus;
		muLookupTable["Neptune"] = kMuNeptune;

		return muLookupTable[body_name];
	}


	/*!\fn inline double radius(const std::string &body_name)
	 *  \brief Returns the equatorial radius of the body.
	 *
	 *  \param body_name    Name of the body for which the equatorial radius is returned.
	 */
	inline double radius(const std::string &body_name) {
		std::map<std::string, double> radiusLookupTable;
		radiusLookupTable["Sun"] = kSunRadius;
		radiusLookupTable["Mercury"] = kMercuryRadius;
		radiusLookupTable["Venus"] = kVenusRadius;
		radiusLookupTable["Earth"] = kEarthRadius;
		radiusLookupTable["Moon"] = kMoonRadius;
		radiusLookupTable["Mars"] = kMarsRadius;
		radiusLookupTable["Jupiter"] = kJupiterRadius;
		radiusLookupTable["Saturn"] = kSaturnRadius;
		radiusLookupTable["Uranus"] = kUranusRadius;
		radiusLookupTable["Neptune"] = kNeptuneRadius;

		return radiusLookupTable[body_name];
	}
}

#endif //CPP_PROPAGATOR_CONSTANTS_H