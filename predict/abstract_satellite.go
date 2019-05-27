/**
  predict4java: An SDP4 / SGP4 library for satellite orbit predictions

  Copyright (C)  2004-2010  David A. B. Johnson, G4DPZ.

  This class is a Java port of one of the core elements of
  the Predict program, Copyright John A. Magliacane,
  KD2BD 1991-2003: http://www.qsl.net/kd2bd/predict.html

  Dr. T.S. Kelso is the author of the SGP4/SDP4 orbital models,
  originally written in Fortran and Pascal, and released into the
  public domain through his website (http://www.celestrak.com/).
  Neoklis Kyriazis, 5B4AZ, later re-wrote Dr. Kelso's code in C,
  and released it under the GNU GPL in 2002.
  PREDICT's core is based on 5B4AZ's code translation efforts.

  Author: David A. B. Johnson, G4DPZ <dave@g4dpz.me.uk>

  Comments, questions and bugreports should be submitted via
  http://sourceforge.net/projects/websat/
  More details can be found at the project home page:

  http://websat.sourceforge.net

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, visit http://www.fsf.org/
*/
package predict

import (
	"math"
	"time"
)

const (
	MINS_PER_DAY                             = 1.44E3
	PI_OVER_TWO                              = math.Pi / 2.0
	SECS_PER_DAY                             = 8.6400E4
	FLATTENING_FACTOR                        = 3.35281066474748E-3
	CK4                                      = 6.209887E-7
	EARTH_GRAVITATIONAL_CONSTANT             = 3.986008E5
	S                                        = 1.012229
	QOMS2T                                   = 1.880279E-09
	EARTH_ROTATIONS_PER_SIDERIAL_DAY         = 1.00273790934
	EARTH_ROTATIONS_RADIANS_PER_SIDERIAL_DAY = EARTH_ROTATIONS_PER_SIDERIAL_DAY * TWO_PI
	RHO                                      = 1.5696615E-1
	MFACTOR                                  = 7.292115E-5
	SOLAR_RADIUS_KM                          = 6.96000E5
	ASTRONOMICAL_UNIT                        = 1.49597870691E8
	SPEED_OF_LIGHT                           = 2.99792458E8
	PERIGEE_156_KM                           = 156.0
	/* WGS 84 Earth radius km */
	EARTH_RADIUS = 6.378137E3
	/* Solar radius - km (IAU 76) */
	SOLAR_RADIUS = 6.96000E5
)

var _ Satellite = &AbstractSatellite{}

type AbstractSatellite struct {
	S4           float64
	Qoms24       float64
	perigee      float64
	tle          TLE
	eclipseDepth float64

	position Vector4
	velocity Vector4
	julUTC   float64
	satPos   SatPos
	julEpoch float64

	calculateSDP4orSGP4 func(float64)
}

func NewAbstractSatellite(tle TLE) *AbstractSatellite {
	return &AbstractSatellite{
		tle:      tle,
		julEpoch: julianDateOfEpoch(tle.Epoch),
	}
}

/**
 * Calculate the geodetic position of an object given its ECI position pos
 * and time. It is intended to be used to determine the ground track of a
 * satellite. The calculations assume the earth to be an oblate spheroid as
 * defined in WGS '72.
 *
 * Reference: The 1992 Astronomical Almanac, page K12.
 *
 * @param time
 *            the time
 */
func (a *AbstractSatellite) calculateLatLonAltFromTime(t float64) {
	calculateLatLonAlt(t, &a.satPos, &a.position)
}

func (a *AbstractSatellite) GetPosition(gsPos *GroundStationPosition, t time.Time) SatPos {

	/* This is the stuff we need to do repetitively while tracking. */
	a.satPos = SatPos{}

	a.julUTC = calcCurrentDaynum(t) + 2444238.5

	/* Convert satellite'S epoch time to Julian */
	/* and calculate time since epoch in minutes */
	tsince := (a.julUTC - a.julEpoch) * MINS_PER_DAY

	a.calculateSDP4orSGP4(tsince)

	/* Scale position and velocity vectors to km and km/sec */
	convertSatState(&a.position, &a.velocity)

	/* Calculate velocity of satellite */
	magnitude(&a.velocity)

	squintVector := &Vector4{}

	/** All angles in rads. Distance in km. Velocity in km/S **/
	/* Calculate satellite Azi, Ele, Range and Range-rate */
	a.calculateObs(a.julUTC, &a.position, &a.velocity, gsPos, squintVector)

	/* Calculate satellite Lat North, Lon East and Alt. */
	a.calculateLatLonAltFromTime(a.julUTC)

	a.satPos.Time = t

	a.satPos.Eclipsed = a.isEclipsed()
	a.satPos.EclipseDepth = a.eclipseDepth

	return a.satPos
}

/**
 * Calculate_User_PosVel() passes the user'S observer position and the time
 * of interest and returns the ECI position and velocity of the observer.
 * The velocity calculation assumes the observer position is stationary
 * relative to the earth'S surface.
 *
 * Reference: The 1992 Astronomical Almanac, page K11.
 *
 * @param time
 *            the time
 * @param gsPos
 *            the ground station position
 * @param obsPos
 *            the position of the observer
 * @param obsVel
 *            the velocity of the observer
 */
func calculateUserPosVel(t float64, gsPos *GroundStationPosition, gsPosTheta *float64, obsPos *Vector4, obsVel *Vector4) {
	*gsPosTheta = mod2PI(thetaGJD(t) + DEG2RAD*gsPos.Longitude)
	c := invert(math.Sqrt(1.0 + FLATTENING_FACTOR*(FLATTENING_FACTOR-2)*sqr(math.Sin(DEG2RAD*gsPos.Latitude))))
	sq := sqr(1.0-FLATTENING_FACTOR) * c
	achcp := (EARTH_RADIUS_KM*c + (gsPos.HeightAMSL / 1000.0)) * math.Cos(DEG2RAD*gsPos.Latitude)
	obsPos.setXYZ(
		achcp*math.Cos(*gsPosTheta),
		achcp*math.Sin(*gsPosTheta),
		(EARTH_RADIUS_KM*sq+(gsPos.HeightAMSL/1000.0))*math.Sin(DEG2RAD*gsPos.Latitude))
	obsVel.setXYZ(-MFACTOR*obsPos.Y, MFACTOR*obsPos.X, 0)
	magnitude(obsPos)
	magnitude(obsVel)
}

/**
 * The procedures Calculate_Obs and Calculate_RADec calculate thetopocentric
 * coordinates of the object with ECI position, {pos}, and velocity, {vel},
 * from location {geodetic} at {time}. The {obs_set} returned for
 * Calculate_Obs consists of azimuth, elevation, range, and range rate (in
 * that order) with units of radians, radians, kilometers, and
 * kilometers/second, respectively. The WGS '72 geoid is used and the effect
 * of atmospheric refraction (under standard temperature and pressure) is
 * incorporated into the elevation calculation; the effect of atmospheric
 * refraction on range and range rate has not yet been quantified.
 *
 * The {obs_set} for Calculate_RADec consists of right ascension and
 * declination (in that order) in radians. Again, calculations are based
 * ontopocentric position using the WGS '72 geoid and incorporating
 * atmospheric refraction.
 *
 * @param julianUTC
 *            Julian date of UTC
 * @param positionVector
 *            the position vector
 * @param velocityVector
 *            the velocity vector
 * @param gsPos
 *            the ground tstation position
 * @param squintVector
 *            the squint vector
 * @param satellitePosition
 *            the satellite position
 *
 */
func (a *AbstractSatellite) calculateObs(julianUTC float64, positionVector *Vector4, velocityVector *Vector4, gsPos *GroundStationPosition, squintVector *Vector4) {
	obsPos := &Vector4{}
	obsVel := &Vector4{}
	Range := &Vector4{}
	rgvel := &Vector4{}

	var gsPosTheta float64
	calculateUserPosVel(julianUTC, gsPos, &gsPosTheta, obsPos, obsVel)

	Range.setXYZ(positionVector.X-obsPos.X, positionVector.Y-obsPos.Y, positionVector.Z-obsPos.Z)

	/* Save these values globally for calculating squint angles later... */
	squintVector.setXYZ(Range.X, Range.Y, Range.Z)

	rgvel.setXYZ(velocityVector.X-obsVel.X, velocityVector.Y-obsVel.Y, velocityVector.Z-obsVel.Z)

	magnitude(Range)

	sinLat := math.Sin(DEG2RAD * gsPos.Latitude)
	cosLat := math.Cos(DEG2RAD * gsPos.Latitude)
	sinTheta := math.Sin(gsPosTheta)
	cosTheta := math.Cos(gsPosTheta)
	topS := sinLat*cosTheta*Range.X + sinLat*sinTheta*Range.Y - cosLat*Range.Z
	topE := -sinTheta*Range.X + cosTheta*Range.Y
	topZ := cosLat*cosTheta*Range.X + cosLat*sinTheta*Range.Y + sinLat*Range.Z
	azim := math.Atan(-topE / topS)

	if topS > 0.0 {
		azim = azim + math.Pi
	}

	if azim < 0.0 {
		azim = azim + TWO_PI
	}

	a.satPos.Azimuth = azim
	a.satPos.Elevation = math.Asin(topZ / Range.W)
	a.satPos.Range = Range.W
	a.satPos.RangeRate = dot(Range, rgvel) / Range.W

	sector := int(a.satPos.Azimuth / TWO_PI * 360.0 / 10.0)

	elevation := (a.satPos.Elevation / TWO_PI) * 360.0

	if elevation > 90 {
		elevation = 180 - elevation
	}

	a.satPos.AboveHorizon = (elevation - float64(gsPos.HorizonElevations[sector])) > EPSILON
}

func (a *AbstractSatellite) WillBeSeen(qth *GroundStationPosition) bool {
	if a.tle.Meanmo < 1e-8 {
		return false
	} else {
		lin := a.tle.Incl

		if lin >= 90.0 {
			lin = 180.0 - lin
		}

		sma := 331.25 * math.Exp(math.Log(1440.0/a.tle.Meanmo)*(2.0/3.0))
		apogee := sma*(1.0+a.tle.Eccn) - EARTH_RADIUS_KM

		return (math.Acos(EARTH_RADIUS_KM/(apogee+EARTH_RADIUS_KM)) + (lin * DEG2RAD)) > math.Abs(qth.Latitude*DEG2RAD)
	}

}

/**
 * Checks and adjusts the calculation if the perigee is less tan 156KM.
 */
func (a *AbstractSatellite) checkPerigee() {
	a.S4 = S
	a.Qoms24 = QOMS2T

	if a.perigee < PERIGEE_156_KM {
		if a.perigee <= 98.0 {
			a.S4 = 20.0
		} else {
			a.S4 = a.perigee - 78.0
		}

		a.Qoms24 = math.Pow((120-a.S4)/EARTH_RADIUS_KM, 4)
		a.S4 = a.S4/EARTH_RADIUS_KM + 1.0
	}
}

/**
 * Sets perigee and checks and adjusts the calculation if the perigee is
 * less tan 156KM.
 *
 * @param perigee
 *            the perigee to set
 */
func (a *AbstractSatellite) setPerigee(perigee float64) {
	a.perigee = perigee
	a.checkPerigee()
}

func (a *AbstractSatellite) CalculateSatelliteVectors(t time.Time) {
	// Re-initialize, object can contain data from previous calculations
	a.satPos = SatPos{}

	// Date/time for which the satellite position and velocity are
	// calculated
	a.julUTC = calcCurrentDaynum(t) + 2444238.5

	// Calculate time since epoch in minutes
	tsince := (a.julUTC - a.julEpoch) * MINS_PER_DAY

	// Calculations of satellite position, no ground stations involved here yet
	a.calculateSDP4orSGP4(tsince)

	// Scale position and velocity vectors to km and km/s
	convertSatState(&a.position, &a.velocity)

	// Calculate the magnitude of the velocity of satellite
	magnitude(&a.velocity)

	a.satPos.Eclipsed = a.isEclipsed()
	a.satPos.EclipseDepth = a.eclipseDepth
	a.satPos.Time = t
}

func (a *AbstractSatellite) CalculateSatelliteVectorsAndReturnPosition(t time.Time) SatPos {
	a.CalculateSatelliteVectors(t)
	return a.satPos
}

func (a *AbstractSatellite) CalculateSatelliteGroundTrack() SatPos {
	a.calculateLatLonAltFromTime(a.julUTC)

	return a.satPos
}

func (a *AbstractSatellite) CalculateSatPosForGroundStation(gsPos *GroundStationPosition) SatPos {
	squintVector := &Vector4{}

	// All angles in rads. Distance in km. Velocity in km/s
	// Calculate satellite Azi, Ele, Range and Range-rate
	a.calculateObs(a.julUTC, &a.position, &a.velocity, gsPos, squintVector)

	return a.satPos
}

/**
 * Determines if the satellite is in sunlight.
 */
func (a *AbstractSatellite) isEclipsed() bool {
	sunVector := a.calculateSunVector()

	/* Calculates stellite's eclipse status and depth */
	/* Determine partial eclipse */

	sdEarth := math.Asin(EARTH_RADIUS / a.position.W)
	rho := subtract(sunVector, &a.position)
	sdSun := math.Asin(SOLAR_RADIUS / rho.W)
	earth := a.position.scalarMultiply(-1)
	delta := angle(sunVector, earth)
	a.eclipseDepth = sdEarth - sdSun - delta

	if sdEarth < sdSun {
		return false
	} else {
		return a.eclipseDepth >= 0
	}
}

func (a *AbstractSatellite) calculateSunVector() *Vector4 {
	mjd := a.julUTC - 2415020.0
	year := 1900 + mjd/365.25
	solTime := (mjd + deltaEt(year)/SECS_PER_DAY) / 36525.0

	m := radians(modulus(358.47583+modulus(35999.04975*solTime, 360.0)-(0.000150+0.0000033*solTime)*sqr(solTime), 360.0))
	l := radians(modulus(279.69668+modulus(36000.76892*solTime, 360.0)+0.0003025*sqr(solTime), 360.0))
	e := 0.01675104 - (0.0000418+0.000000126*solTime)*solTime
	c := radians((1.919460-(0.004789+0.000014*solTime)*solTime)*math.Sin(m) + (0.020094-0.000100*solTime)*math.Sin(2*m) + 0.000293*math.Sin(3*m))
	o := radians(modulus(259.18-1934.142*solTime, 360.0))
	lsa := modulus(l+c-radians(0.00569-0.00479*math.Sin(o)), TWO_PI)
	nu := modulus(m+c, TWO_PI)
	r := 1.0000002 * (1.0 - sqr(e)) / (1.0 + e*math.Cos(nu))
	eps := radians(23.452294 - (0.0130125+(0.00000164-0.000000503*solTime)*solTime)*solTime + 0.00256*math.Cos(o))
	r = ASTRONOMICAL_UNIT * r

	return &Vector4{
		W: r,
		X: r * math.Cos(lsa),
		Y: r * math.Sin(lsa) * math.Cos(eps),
		Z: r * math.Sin(lsa) * math.Sin(eps),
	}
}

func (a *AbstractSatellite) calculatePhase(xlt float64, xnode float64, omegadf float64) {
	/* Phase in radians */
	phaseValue := xlt - xnode - omegadf + TWO_PI

	if phaseValue < 0.0 {
		phaseValue += TWO_PI
	}

	a.satPos.Phase = mod2PI(phaseValue)
}

func (a *AbstractSatellite) calculatePositionAndVelocity(rk float64, uk float64, xnodek float64, xinck float64, rdotk float64, rfdotk float64) {
	/* Orientation vectors */
	sinuk := math.Sin(uk)
	cosuk := math.Cos(uk)
	sinik := math.Sin(xinck)
	cosik := math.Cos(xinck)
	sinnok := math.Sin(xnodek)
	cosnok := math.Cos(xnodek)
	xmx := -sinnok * cosik
	xmy := cosnok * cosik
	ux := xmx*sinuk + cosnok*cosuk
	uy := xmy*sinuk + sinnok*cosuk
	uz := sinik * sinuk
	vx := xmx*cosuk - cosnok*sinuk
	vy := xmy*cosuk - sinnok*sinuk
	vz := sinik * cosuk

	/* Position and velocity */
	a.position.setXYZ(ux, uy, uz)
	a.position.multiply(rk)
	a.velocity.setXYZ(
		rdotk*ux+rfdotk*vx,
		rdotk*uy+rfdotk*vy,
		rdotk*uz+rfdotk*vz)
}

func invert(value float64) float64 {
	return 1.0 / value
}

/**
 * The function Julian_Date_of_Epoch returns the Julian Date of an epoch
 * specified in the format used in the NORAD two-line element sets. It has
 * been modified to support dates beyond the year 1999 assuming that
 * two-digit years in the range 00-56 correspond to 2000-2056. Until the
 * two-line element set format is changed, it is only valid for dates
 * through 2056 December 31.
 *
 * @param epoch
 *            the Epoch
 * @return The Julian date of the Epoch
 */
func julianDateOfEpoch(epoch float64) float64 {
	/* Modification to support Y2K */
	/* Valid 1957 through 2056 */
	year := math.Floor(epoch * 1e-3)
	day := (epoch*1e-3 - year) * 1000.0

	if year < 57 {
		year += 2000
	} else {
		year += 1900
	}

	return julianDateOfYear(year) + day
}

/**
 * Calculates the Julian Day of the Year.
 *
 * The function Julian_Date_of_Year calculates the Julian Date of Day 0.0 of
 * {year}. This function is used to calculate the Julian Date of any date by
 * using Julian_Date_of_Year, DOY, and Fraction_of_Day.
 *
 * Astronomical Formulae for Calculators, Jean Meeus, pages 23-25. Calculate
 * Julian Date of 0.0 Jan aYear
 *
 * @param theYear
 *            the year
 * @return the Julian day number
 */
func julianDateOfYear(theYear float64) float64 {
	aYear := theYear - 1
	i := int64(math.Floor(aYear / 100))
	a := i
	i = a / 4
	b := 2 - a + i
	i = int64(math.Floor(365.25 * aYear))
	i += int64(428) //int64(float64(30.6001 * 14))

	return float64(i) + 1720994.5 + float64(b)
}

var sgp4Epoch = time.Date(1979, 12, 31, 0, 0, 0, 0, time.UTC)

/**
 * Read the system clock and return the number of days since 31Dec79
 * 00:00:00 UTC (daynum 0).
 *
 * @param date
 *            the date we wan to get the offset for
 * @return the number of days offset
 */
func calcCurrentDaynum(date time.Time) float64 {
	now := float64(date.UnixNano()) / 1000000
	then := float64(sgp4Epoch.UnixNano()) / 1000000
	millis := now - then
	return millis / 1000.0 / 60.0 / 60.0 / 24.0
}

/**
 * Returns the square of a double.
 *
 * @param arg
 *            the value for which to get the double
 * @return the arg squared
 */
func sqr(arg float64) float64 {
	return arg * arg
}

/**
 * Calculates scalar magnitude of a vector4 argument.
 *
 * @param v
 *            the vector were measuring
 *
 */
func magnitude(v *Vector4) {
	v.W = math.Sqrt(sqr(v.X) + sqr(v.Y) + sqr(v.Z))
}

/**
 * Multiplies the vector v1 by the scalar k.
 *
 * @param k
 *            the multiplier
 * @param v
 *            the vector
 */
func scaleVector(k float64, v *Vector4) {
	v.multiply(k)
	magnitude(v)
}

/**
 * Calculates the dot product of two vectors.
 *
 * @param v1
 *            vector 1
 * @param v2
 *            vector 2
 * @return the dot product
 */

func dot(v1 *Vector4, v2 *Vector4) float64 {
	return v1.X*v2.X + v1.Y*v2.Y + v1.Z*v2.Z
}

type Vector4 struct {
	W float64
	X float64
	Y float64
	Z float64
}

func (v *Vector4) multiply(multiplier float64) {
	v.X *= multiplier
	v.Y *= multiplier
	v.Z *= multiplier
}

func (v *Vector4) setXYZ(xValue float64, yValue float64, zValue float64) {
	v.X = xValue
	v.Y = yValue
	v.Z = zValue
}

func (v *Vector4) subtract(vector *Vector4) *Vector4 {
	return &Vector4{
		W: v.W - vector.W,
		X: v.X - vector.X,
		Y: v.Y - vector.Y,
		Z: v.Z - vector.Z,
	}
}

func (v *Vector4) scalarMultiply(multiplier float64) *Vector4 {
	return &Vector4{
		W: v.W * math.Abs(multiplier),
		X: v.X * multiplier,
		Y: v.Y * multiplier,
		Z: v.Z * multiplier,
	}
}

func angle(v1 *Vector4, v2 *Vector4) float64 {
	magnitude(v1)
	magnitude(v2)
	return math.Acos(dot(v1, v2) / (v1.W * v2.W))
}

func subtract(v1 *Vector4, v2 *Vector4) *Vector4 {
	v3 := &Vector4{
		X: v1.X - v2.X,
		Y: v1.Y - v2.Y,
		Z: v1.Z - v2.Z,
	}
	magnitude(v3)
	return v3
}

/**
 * Gets the modulus of a double value.
 *
 * @param arg1
 *            the value to be tested
 * @param arg2
 *            the divisor
 * @return the remainder
 */
func modulus(arg1 float64, arg2 float64) float64 {
	/* Returns arg1 mod arg2 */

	returnValue := arg1

	i := int(math.Floor(returnValue / arg2))
	returnValue -= float64(i) * arg2

	if returnValue < 0.0 {
		returnValue += arg2
	}

	return returnValue
}

func frac(arg float64) float64 {
	/* Returns fractional part of double argument */
	return arg - math.Floor(arg)
}

func thetaGJD(theJD float64) float64 {
	/* Reference: The 1992 Astronomical Almanac, page B6. */

	ut := frac(theJD + 0.5)
	aJD := theJD - ut
	tu := (aJD - 2451545.0) / 36525.0
	gmst := 24110.54841 + tu*(8640184.812866+tu*(0.093104-tu*6.2E-6))
	gmst = modulus(gmst+SECS_PER_DAY*EARTH_ROTATIONS_PER_SIDERIAL_DAY*ut, SECS_PER_DAY)

	return TWO_PI * gmst / SECS_PER_DAY
}

/**
 * Calculates the modulus of 2 * PI.
 *
 * @param testValue
 *            the value under test
 * @return the modulus
 */
func mod2PI(testValue float64) float64 {
	/* Returns mod 2PI of argument */
	retVal := testValue
	i := int(retVal / TWO_PI)
	retVal -= float64(i) * TWO_PI

	if retVal < 0.0 {
		retVal += TWO_PI
	}

	return retVal
}

func calculateLatLonAlt(t float64, satPos *SatPos, position *Vector4) {
	satPos.Theta = math.Atan2(position.Y, position.X)
	satPos.Longitude = mod2PI(satPos.Theta - thetaGJD(t))
	r := math.Sqrt(sqr(position.X) + sqr(position.Y))
	e2 := FLATTENING_FACTOR * (2.0 - FLATTENING_FACTOR)
	satPos.Latitude = math.Atan2(position.Z, r)

	var phi float64
	var c float64
	i := 0
	converged := false

	for {
		phi = satPos.Latitude
		c = invert(math.Sqrt(1.0 - e2*sqr(math.Sin(phi))))
		satPos.Latitude = math.Atan2(position.Z+EARTH_RADIUS_KM*c*e2*math.Sin(phi), r)

		converged = math.Abs(satPos.Latitude-phi) < EPSILON
		i += 1
		if i >= 10 || converged {
			break
		}
	}

	satPos.Altitude = r/math.Cos(satPos.Latitude) - EARTH_RADIUS_KM*c

	temp := satPos.Latitude

	if temp > PI_OVER_TWO {
		temp -= TWO_PI
		satPos.Latitude = temp
	}
}

/**
 * Converts the satellite'S position and velocity vectors from normalized
 * values to km and km/sec.
 *
 * @param pos
 *            the position
 * @param vel
 *            the velocity
 */
func convertSatState(pos *Vector4, vel *Vector4) {
	/* Converts the satellite'S position and velocity */
	/* vectors from normalized values to km and km/sec */
	scaleVector(EARTH_RADIUS_KM, pos)
	scaleVector(EARTH_RADIUS_KM*MINS_PER_DAY/SECS_PER_DAY, vel)
}

/**
 * The function Delta_ET has been added to allow calculations on the
 * position of the sun. It provides the difference between UT (approximately
 * the same as UTC) and ET (now referred to as TDT) This function is based
 * on a least squares fit of data from 1950 to 1991 and will need to be
 * updated periodically.
 *
 * Values determined using data from 1950-1991 in the 1990 Astronomical
 * Almanac. See DELTA_ET.WQ1 for details.
 */
func deltaEt(year float64) float64 {
	return 26.465 + 0.747622*(year-1950) + 1.886913*math.Sin(TWO_PI*(year-1975)/33)
}

/**
 * Returns angle in radians from argument in degrees.
 */
func radians(degrees float64) float64 {
	return degrees * DEG2RAD
}

/**
 * Solves Keplers' Equation.
 *
 * @param temp
 *            an array of temporary values we pass around as part of the
 *            orbit calculation.
 * @param axn
 * @param ayn
 * @param capu
 */
func converge(temp []float64, axn float64, ayn float64, capu float64) {
	converged := false
	i := 0

	for {
		temp[7] = math.Sin(temp[2])
		temp[8] = math.Cos(temp[2])
		temp[3] = axn * temp[7]
		temp[4] = ayn * temp[8]
		temp[5] = axn * temp[8]
		temp[6] = ayn * temp[7]
		epw := (capu-temp[4]+temp[3]-temp[2])/(1.0-temp[5]-temp[6]) + temp[2]

		if math.Abs(epw-temp[2]) <= EPSILON {
			converged = true
		} else {
			temp[2] = epw
		}

		i += 1
		if i >= 10 || converged {
			break
		}
	}
}
