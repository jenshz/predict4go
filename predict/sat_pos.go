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
	"fmt"
	"math"
	"time"
)

const (
	DEG_CR = " deg.\n"
	R0     = 6378.16
)

/**
 *
 * @author g4dpz
 *
 */

type SatPos struct {
	Azimuth   float64
	Elevation float64
	Latitude  float64
	Longitude float64
	Time      time.Time
	Range     float64
	RangeRate float64
	Phase     float64
	Altitude  float64
	Theta     float64

	EclipseDepth float64
	Eclipsed     bool

	AboveHorizon bool
}

// 	/**
// 	 * Constructs a Satellite Position.
// 	 *
// 	 * @param azimuth
// 	 *            the Azimuth
// 	 * @param elevation
// 	 *            the Elevation
// 	 * @param theTime
// 	 *            the Time
// 	 */
func NewSatPos(azimuth float64, elevation float64, theTime time.Time) *SatPos {
	return &SatPos{
		Azimuth:   azimuth,
		Elevation: elevation,
		Time:      theTime,
	}
}

func ftoa(f float64) string {
	return fmt.Sprintf("%f", f)
}

/**
 * @return a pretty printed version of the Satellite Position
 */
func (s *SatPos) String() string {
	return "Azimuth:      " + ftoa(s.Azimuth/(math.Pi*2.0)*360) + DEG_CR +
		"Elevation:    " + ftoa(s.Elevation/(math.Pi*2.0)*360) + DEG_CR +
		"Latitude:     " + ftoa(s.Latitude/(math.Pi*2.0)*360) + DEG_CR +
		"Longitude:    " + ftoa(s.Longitude/(math.Pi*2.0)*360) + DEG_CR +
		"Date:         " + s.Time.String() + "\n" +
		"Range:        " + ftoa(s.Range) + " km.\n" +
		"Range rate:   " + ftoa(s.RangeRate) + " m/S.\n" +
		"Phase:        " + ftoa(s.Phase) + " /(256)\n" +
		"Altitude:     " + ftoa(s.Altitude) + " km\n" +
		"Theta:        " + ftoa(s.Theta) + " rad/sec\n" +
		"Eclipsed:     " + fmt.Sprintf("%v", s.Eclipsed) + "\n" +
		"Eclipse depth:" + ftoa(s.EclipseDepth) + " radians\n"
}

func (s *SatPos) ToShortString() string {
	return "Elevation: " + fmt.Sprintf("%.0f", s.Elevation/(math.Pi*2.0)*360) + DEG_CR +
		"Azimuth: " + fmt.Sprintf("%.0f", s.Azimuth/(math.Pi*2.0)*360) + DEG_CR +
		"Latitude: " + fmt.Sprintf("%.2f", s.Latitude/(math.Pi*2.0)*360) + DEG_CR +
		"Longitude: " + fmt.Sprintf("%.2f", s.Longitude/(math.Pi*2.0)*360) + DEG_CR +
		"Range: " + fmt.Sprintf("%.0f", s.Range) + " Km"
}

/**
 * Calculates the footprint range circle using the given increment. TODO
 * where is first point, give heading.
 *
 * @param incrementDegrees
 * @return
 */
func (s *SatPos) getRangeCircle(incrementDegrees float64) []Position {
	return calculateRangeCirclePoints(s, incrementDegrees)
}

/**
 * Calculates the footprint range circle using an increment of 1.0 degrees.
 *
 * @param pos
 * @return a list of {@link Position}
 */
func (s *SatPos) getRangeCircleWithIncrement() []Position {
	return s.getRangeCircle(1.0)
}

func (s *SatPos) getRangeCircleRadiusKm() float64 {
	return 0.5 * (12756.33 * math.Acos(EARTH_RADIUS_KM/(EARTH_RADIUS_KM+s.Altitude)))
}

/**
* Calculates the footprint range circle using the given increment.
*
* @param pos
* @return a list of {@link Position}
 */
func calculateRangeCirclePoints(pos *SatPos, incrementDegrees float64) []Position {
	radiusKm := pos.getRangeCircleRadiusKm()
	latitude := pos.Latitude
	longitude := pos.Longitude
	beta := radiusKm / R0
	var result []Position
	for azi := 0.0; azi < 360.0; azi += incrementDegrees {
		azimuth := (azi / 360.0) * 2.0 * math.Pi
		rangeLat := math.Asin(math.Sin(latitude)*math.Cos(beta) + math.Cos(azimuth)*math.Sin(beta)*math.Cos(latitude))
		num := math.Cos(beta) - (math.Sin(latitude) * math.Sin(rangeLat))
		den := math.Cos(latitude) * math.Cos(rangeLat)
		var rangeLon float64

		if azi == 0 && (beta > ((math.Pi / 2.0) - latitude)) {
			rangeLon = longitude + math.Pi
		} else if azi == 180 && (beta > ((math.Pi / 2.0) - latitude)) {
			rangeLon = longitude + math.Pi
		} else if math.Abs(num/den) > 1.0 {
			rangeLon = longitude
		} else {
			if (180 - azi) >= 0 {
				rangeLon = longitude - math.Acos(num/den)
			} else {
				rangeLon = longitude + math.Acos(num/den)
			}
		}

		for rangeLon < 0.0 {
			rangeLon += math.Pi * 2.0
		}

		for rangeLon > math.Pi*2.0 {
			rangeLon -= math.Pi * 2.0
		}

		rangeLat = (rangeLat / (2.0 * math.Pi)) * 360.0
		rangeLon = (rangeLon / (2.0 * math.Pi)) * 360.0

		result = append(result, Position{Lat: rangeLat, Lon: rangeLon})
	}

	return result
}
