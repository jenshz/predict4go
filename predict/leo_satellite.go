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
)

type LEOSatellite struct {
	super *AbstractSatellite
	aodp float64
	aycof float64
	c1 float64
	c4 float64
	c5 float64
	cosio float64
	d2 float64
	d3 float64
	d4 float64
	delmo float64
	omgcof float64
	eta float64
	omgdot float64
	sinio float64
	xnodp float64
	sinmo float64
	t2cof float64
	t3cof float64
	t4cof float64
	t5cof float64
	x1mth2 float64
	x3thm1 float64
	x7thm1 float64
	xmcof float64
	xmdot float64
	xnodcf float64
	xnodot float64
	xlcof float64

	sgp4Simple bool
}

/**
 * 
 * @author g4dpz
 * 
 */
func NewLEOSatellite(tle TLE) (Satellite) {
	sat := &LEOSatellite{
		super: NewAbstractSatellite(tle),
	}

	/* Recover original mean motion (xnodp) and */
	/* semimajor axis (aodp) from input elements. */

	a1 := math.Pow(XKE / tle.Xno, TWO_THIRDS)
	sat.cosio = math.Cos(tle.Xincl)
	theta2 := sqr(sat.cosio)
	sat.x3thm1 = 3.0 * theta2 - 1.0
	eo := tle.Eo
	eosq := sqr(eo)
	betao2 := 1.0 - eosq
	betao := math.Sqrt(betao2)
	del1 := 1.5 * CK2 * sat.x3thm1 / (sqr(a1) * betao * betao2)
	ao := a1 * (1.0 - del1 * (0.5 * TWO_THIRDS + del1 * (1.0 + 134.0 / 81.0 * del1)))
	delo := 1.5 * CK2 * sat.x3thm1 / (sqr(ao) * betao * betao2)
	sat.xnodp = tle.Xno / (1.0 + delo)
	sat.aodp = ao / (1.0 - delo)

	/* For perigee less than 220 kilometers, the "simple" */
	/* flag is set and the equations are truncated to linear */
	/* variation in sqrt a and quadratic variation in mean */
	/* anomaly. Also, the c3 term, the delta omega term, and */
	/* the delta m term are dropped. */

	sat.sgp4Simple = (sat.aodp * (1.0 - eo)) < (220 / EARTH_RADIUS_KM + 1.0)

	/* For perigees below 156 km, the */
	/* values of S and QOMS2T are altered. */
	sat.super.setPerigee((sat.aodp * (1.0 - eo) - 1.0) * EARTH_RADIUS_KM)

	pinvsq := invert(sqr(sat.aodp) * sqr(betao2))
	tsi := invert(sat.aodp - sat.super.S4)
	sat.eta = sat.aodp * eo * tsi
	etasq := sat.eta * sat.eta
	eeta := eo * sat.eta
	psisq := math.Abs(1.0 - etasq)
	coef := sat.super.Qoms24 * math.Pow(tsi, 4)
	coef1 := coef / math.Pow(psisq, 3.5)
	bstar := tle.Bstar
	c2 := coef1 * sat.xnodp * (sat.aodp * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) + 0.75 * CK2 * tsi / psisq * sat.x3thm1 * (8.0 + 3.0 * etasq * (8.0 + etasq)))
	sat.c1 = bstar * c2
	sat.sinio = math.Sin(tle.Xincl)
	a3ovk2 := -J3_HARMONIC / CK2
	c3 := coef * tsi * a3ovk2 * sat.xnodp * sat.sinio / eo
	sat.x1mth2 = 1.0 - theta2

	omegao := tle.Omegao

	sat.c4 = 2 * sat.xnodp * coef1 * sat.aodp * betao2 * (sat.eta * (2.0 + 0.5 * etasq) + eo * (0.5 + 2 * etasq) - 2 * CK2 * tsi / (sat.aodp * psisq) * (-3 * sat.x3thm1 * (1.0 - 2 * eeta + etasq * (1.5 - 0.5 * eeta)) + 0.75 * sat.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * math.Cos(2.0 * omegao)))
	sat.c5 = 2.0 * coef1 * sat.aodp * betao2 * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq)

	theta4 := sqr(theta2)
	temp1 := 3.0 * CK2 * pinvsq * sat.xnodp
	temp2 := temp1 * CK2 * pinvsq
	temp3 := 1.25 * CK4 * pinvsq * pinvsq * sat.xnodp
	sat.xmdot = sat.xnodp + 0.5 * temp1 * betao * sat.x3thm1 + 0.0625 * temp2 * betao * (13.0 - 78.0 * theta2 + 137.0 * theta4)
	x1m5th := 1.0 - 5.0 * theta2
	sat.omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4) + temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4)
	xhdot1 := -temp1 * sat.cosio
	sat.xnodot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 * (3.0 - 7.0 * theta2)) * sat.cosio
	sat.omgcof = bstar * c3 * math.Cos(omegao)
	sat.xmcof = -TWO_THIRDS * coef * bstar / eeta
	sat.xnodcf = 3.5 * betao2 * xhdot1 * sat.c1
	sat.t2cof = 1.5 * sat.c1
	sat.xlcof = 0.125 * a3ovk2 * sat.sinio * (3.0 + 5 * sat.cosio) / (1.0 + sat.cosio)
	sat.aycof = 0.25 * a3ovk2 * sat.sinio
	xmo := tle.Xmo
	sat.delmo = math.Pow(1.0 + sat.eta * math.Cos(xmo), 3)
	sat.sinmo = math.Sin(xmo)
	sat.x7thm1 = 7.0 * theta2 - 1

	if (!sat.sgp4Simple) {
		c1sq := sqr(sat.c1)
		sat.d2 = 4.0 * sat.aodp * tsi * c1sq
		temp := sat.d2 * tsi * sat.c1 / 3.0
		sat.d3 = (17 * sat.aodp + sat.super.S4) * temp
		sat.d4 = 0.5 * temp * sat.aodp * tsi * (221 * sat.aodp + 31 * sat.super.S4) * sat.c1
		sat.t3cof = sat.d2 + 2 * c1sq
		sat.t4cof = 0.25 * (3.0 * sat.d3 + sat.c1 * (12 * sat.d2 + 10 * c1sq))
		sat.t5cof = 0.2 * (3.0 * sat.d4 + 12 * sat.c1 * sat.d3 + 6 * sat.d2 * sat.d2 + 15 * c1sq * (2.0 * sat.d2 + c1sq))
	} else {
		sat.d2 = 0
		sat.d3 = 0
		sat.d4 = 0
		sat.t3cof = 0
		sat.t4cof = 0
		sat.t5cof = 0
	}

	sat.super.calculateSDP4orSGP4 = sat.calculateSGP4

	return sat.super
}

func (sat *LEOSatellite) calculateSGP4(tsince float64) {
	temp := make([]float64, 9)

	/* Update for secular gravity and atmospheric drag. */
	xmdf := sat.super.tle.Xmo + sat.xmdot * tsince
	omgadf := sat.super.tle.Omegao + sat.omgdot * tsince
	xnoddf := sat.super.tle.Xnodeo + sat.xnodot * tsince
	omega := omgadf
	xmp := xmdf
	tsq := sqr(tsince)
	xnode := xnoddf + sat.xnodcf * tsq
	bstar := sat.super.tle.Bstar
	tempa := 1.0 - sat.c1 * tsince
	tempe := bstar * sat.c4 * tsince
	templ := sat.t2cof * tsq

	if (!sat.sgp4Simple) {
		delomg := sat.omgcof * tsince
		delm := sat.xmcof * (math.Pow(1.0 + sat.eta * math.Cos(xmdf), 3) - sat.delmo)
		temp[0] = delomg + delm
		xmp = xmdf + temp[0]
		omega = omgadf - temp[0]
		tcube := tsq * tsince
		tfour := tsince * tcube
		tempa = tempa - sat.d2 * tsq - sat.d3 * tcube - sat.d4 * tfour
		tempe = tempe + bstar * sat.c5 * (math.Sin(xmp) - sat.sinmo)
		templ = templ + sat.t3cof * tcube + tfour * (sat.t4cof + tsince * sat.t5cof)
	}

	a := sat.aodp * math.Pow(tempa, 2)
	eo := sat.super.tle.Eo
	e := eo - tempe
	xl := xmp + omega + xnode + sat.xnodp * templ
	beta := math.Sqrt(1.0 - e * e)
	xn := XKE / math.Pow(a, 1.5)

	/* Long period periodics */
	axn := e * math.Cos(omega)
	temp[0] = invert(a * sqr(beta))
	xll := temp[0] * sat.xlcof * axn
	aynl := temp[0] * sat.aycof
	xlt := xl + xll
	ayn := e * math.Sin(omega) + aynl

	/* Solve Kepler'S Equation */
	capu := mod2PI(xlt - xnode)
	temp[2] = capu

	converge(temp, axn, ayn, capu)

	sat.calculatePositionAndVelocity(temp, xnode, a, xn, axn, ayn)

	sat.super.calculatePhase(xlt, xnode, omgadf)
}


func (sat *LEOSatellite) calculatePositionAndVelocity(temp []float64, xnode float64, a float64, xn float64, axn float64, ayn float64) {
	ecose := temp[5] + temp[6]
	esine := temp[3] - temp[4]
	elsq := sqr(axn) + sqr(ayn)
	temp[0] = 1.0 - elsq
	pl := a * temp[0]
	r := a * (1.0 - ecose)
	temp[1] = invert(r)
	rdot := XKE * math.Sqrt(a) * esine * temp[1]
	rfdot := XKE * math.Sqrt(pl) * temp[1]
	temp[2] = a * temp[1]
	betal := math.Sqrt(temp[0])
	temp[3] = invert(1.0 + betal)
	cosu := temp[2] * (temp[8] - axn + ayn * esine * temp[3])
	sinu := temp[2] * (temp[7] - ayn - axn * esine * temp[3])
	u := math.Atan2(sinu, cosu)
	sin2u := 2.0 * sinu * cosu
	cos2u := 2.0 * cosu * cosu - 1
	temp[0] = invert(pl)
	temp[1] = CK2 * temp[0]
	temp[2] = temp[1] * temp[0]

	/* Update for short periodics */
	rk := r * (1.0 - 1.5 * temp[2] * betal * sat.x3thm1) + 0.5 * temp[1] * sat.x1mth2 * cos2u
	uk := u - 0.25 * temp[2] * sat.x7thm1 * sin2u
	xnodek := xnode + 1.5 * temp[2] * sat.cosio * sin2u
	xinck := sat.super.tle.Xincl + 1.5 * temp[2] * sat.cosio * sat.sinio * cos2u
	rdotk := rdot - xn * temp[1] * sat.x1mth2 * sin2u
	rfdotk := rfdot + xn * temp[1] * (sat.x1mth2 * cos2u + 1.5 * sat.x3thm1)

	sat.super.calculatePositionAndVelocity(rk, uk, xnodek, xinck, rdotk, rfdotk)
}
