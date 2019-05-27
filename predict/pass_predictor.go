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
	"errors"
	"math"
	"time"
)

const (
	UTC           = "UTC"
	SOUTH         = "south"
	NORTH         = "north"
	TWOPI         = math.Pi * 2.0
	DEADSPOT_NONE = "none"
)

/**
 * Class which provides Pass Prediction.
 *
 * @author David A. B. Johnson, g4dpz
 *
 */
type PassPredictor struct {
	tle            *TLE
	qth            *GroundStationPosition
	sat            Satellite
	iterationCount int
}

/**
 * Constructor.
 *
 * @param tle
 *            the Three Line Elements
 * @param qth
 *            the ground station position
 */
func NewPassPredictor(tle *TLE, qth *GroundStationPosition) (*PassPredictor, error) {
	result := &PassPredictor{
		tle: tle,
		qth: qth,
		sat: NewSatellite(tle),
	}

	if !result.sat.WillBeSeen(result.qth) {
		return nil, errors.New("Satellite will never appear above the horizon.")
	}

	return result, nil
}

/**
 * Gets the downlink frequency corrected for doppler.
 *
 * @param freq
 *            the original frequency in Hz
 * @return the doppler corrected frequency in Hz
 * @throws InvalidTleException
 *             bad TLE passed in
 * @throws SatNotFoundException
 */
func (p *PassPredictor) GetDownlinkFreq(freq float64, t time.Time) int64 {
	// get the current position
	satPos := p.getSatPos(t)
	rangeRate := satPos.RangeRate
	return int64(float64(freq * (SPEED_OF_LIGHT - rangeRate*1000.0) / SPEED_OF_LIGHT))
}

func (p *PassPredictor) GetUplinkFreq(freq float64, t time.Time) int64 {
	satPos := p.getSatPos(t)
	rangeRate := satPos.RangeRate
	return int64(freq * (SPEED_OF_LIGHT + rangeRate*1000.0) / SPEED_OF_LIGHT)
}

func (p *PassPredictor) getSatPos(t time.Time) SatPos {
	p.iterationCount += 1
	return p.sat.GetPosition(p.qth, t)
}

func (p *PassPredictor) nextSatPass(t time.Time, targetElevation float64) *SatPassTime {
	return p.nextSatPassInternal(t, false, targetElevation)
}

func (p *PassPredictor) nextSatPassInternal(t time.Time, windBack bool, targetElevation float64) *SatPassTime {
	aosAzimuth := 0
	losAzimuth := 0
	maxElevation := 0.0

	polePassed := DEADSPOT_NONE

	cal := t
	stopSearch := cal.Add(time.Duration(24) * time.Hour)

	// wind back time 1/4 of an orbit
	if windBack {
		cal = cal.Add(time.Minute * time.Duration(int(-24.0*60.0/p.tle.Meanmo/4.0)))
	}

	satPos := p.getSatPos(cal)
	prevPos := satPos

	// test for the elevation being above the horizon
	if satPos.Elevation > targetElevation {
		// move time forward in 30 second intervals until the sat goes below
		// the horizon
		for {
			cal, satPos = p.getPosition(cal, 60)
			if satPos.Elevation >= targetElevation {
				break
			}
			if !cal.Before(stopSearch) {
				return nil
			}
		}

		// move time forward 3/4 orbit
		cal = cal.Add(time.Minute * time.Duration(p.threeQuarterOrbitMinutes()))
	}

	var tca time.Time

	// now find the next time it comes above the horizon
	findCrossing := func(jump int, aos bool, checkPolePassed bool) bool {
		for {
			cal, satPos = p.getPosition(cal, jump)

			if checkPolePassed {
				currPolePassed := p.getPolePassed(prevPos, satPos)
				if currPolePassed != DEADSPOT_NONE {
					polePassed = currPolePassed
				}
			}

			if satPos.Elevation > maxElevation {
				maxElevation = satPos.Elevation
				tca = cal
			}
			prevPos = satPos

			if satPos.Elevation >= targetElevation && aos {
				break // aos
			} else if satPos.Elevation <= targetElevation && !aos {
				break // los
			}
			if !cal.Before(stopSearch) {
				return false
			}
		}
		return true
	}

	if !findCrossing(60, true, false) {
		return nil
	}

	// refine it to 1 seconds
	cal = cal.Add(time.Second * time.Duration(-60))
	if !findCrossing(1, true, false) {
		return nil
	}

	startDate := satPos.Time

	aosAzimuth = int((satPos.Azimuth / (2.0 * math.Pi)) * 360.0)

	// now find when it goes below
	if !findCrossing(30, false, true) {
		return nil
	}

	// refine it to 1 seconds
	cal = cal.Add(time.Second * time.Duration(-30))
	if !findCrossing(1, false, false) {
		return nil
	}

	endDate := satPos.Time
	losAzimuth = int((satPos.Azimuth / (2.0 * math.Pi)) * 360.0)

	return NewSatPassTime(startDate, endDate, tca, polePassed, aosAzimuth, losAzimuth, (maxElevation/(2.0*math.Pi))*360.0)

}

func (p *PassPredictor) getPosition(t time.Time, offset int) (time.Time, SatPos) {
	t = t.Add(time.Second * time.Duration(offset))
	return t, p.getSatPos(t)
}

/**
 * Gets a list of SatPassTime
 */
func (p *PassPredictor) GetPasses(start time.Time, hoursAhead int, windBack bool, targetElevation float64) []*SatPassTime {
	p.iterationCount = 0
	windBackTime := windBack

	var passes []*SatPassTime

	trackStartTime := start.Truncate(time.Second)
	trackEndTime := trackStartTime.Add(time.Duration(hoursAhead) * time.Hour)

	var lastAOS time.Time

	count := 0

	for lastAOS.Before(trackEndTime) {
		if count > 0 {
			windBackTime = false
		}

		pass := p.nextSatPassInternal(trackStartTime, windBackTime, DEG2RAD*targetElevation)
		if pass == nil {
			break
		}

		lastAOS = pass.StartTime

		passes = append(passes, pass)

		trackStartTime = pass.EndTime.Add(time.Duration(p.threeQuarterOrbitMinutes()) * time.Minute)

		count += 1
	}

	return passes
}

// 	/**
// 	 * Returns the iterationCount. @VisibleForTesting
// 	 *
// 	 * @return the iterationCount
// 	 */
// 	final int getIterationCount() {
// 		return iterationCount;
// 	}

/**
 * @return time in mS for 3/4 of an orbit
 */
func (s *PassPredictor) threeQuarterOrbitMinutes() int {
	return (int)(24.0 * 60.0 / s.tle.Meanmo * 0.75)
}

func (s *PassPredictor) getPolePassed(prevPos SatPos, satPos SatPos) string {
	polePassed := DEADSPOT_NONE

	az1 := prevPos.Azimuth / TWOPI * 360.0
	az2 := satPos.Azimuth / TWOPI * 360.0

	if az1 > az2 {
		// we may be moving from 350 or greateer thru north
		if az1 > 350 && az2 < 10 {
			polePassed = NORTH
		} else {
			// we may be moving from 190 or greateer thru south
			if az1 > 180 && az2 < 180 {
				polePassed = SOUTH
			}
		}
	} else {
		// we may be moving from 10 or less through north
		if az1 < 10 && az2 > 350 {
			polePassed = NORTH
		} else {
			// we may be moving from 170 or more through south
			if az1 < 180 && az2 > 180 {
				polePassed = SOUTH
			}
		}
	}

	return polePassed
}

/**
 * Calculates positions of satellite for a given point in time, time range
 * and step increment.
 */
func (p *PassPredictor) GetPositions(t time.Time, incrementSeconds int, minutesBefore int, minutesAfter int) []SatPos {
	trackTime := t.Add(-(time.Duration(minutesBefore) * time.Minute))
	endTime := t.Add(time.Duration(minutesAfter) * time.Minute)

	var positions []SatPos

	for trackTime.Before(endTime) {
		positions = append(positions, p.getSatPos(trackTime))
		trackTime = trackTime.Add(time.Duration(incrementSeconds) * time.Second)
	}

	return positions
}
