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

const (
	ZSINIS = 3.9785416E-1
	ZSINGS = -9.8088458E-1
	ZNS    = 1.19459E-5
	C1SS   = 2.9864797E-6
	ZES    = 1.675E-2
	ZNL    = 1.5835218E-4
	C1L    = 4.7968065E-7
	ZEL    = 5.490E-2
	ROOT22 = 1.7891679E-6
	ROOT32 = 3.7393792E-7
	ROOT44 = 7.3636953E-9
	ROOT52 = 1.1428639E-7
	ROOT54 = 2.1765803E-9
	THDT   = 4.3752691E-3
	Q22    = 1.7891679E-6
	Q31    = 2.1460748E-6
	Q33    = 2.2123015E-7
	G22    = 5.7686396
	G32    = 9.5240898E-1
	G44    = 1.8014998
	G52    = 1.0508330
	G54    = 4.4108898
)

type DeepSpaceSatellite struct {
	super AbstractSatellite

	c1     float64
	c4     float64
	x1mth2 float64
	x3thm1 float64
	xlcof  float64
	xnodcf float64
	t2cof  float64
	aycof  float64
	x7thm1 float64

	dsv  DeepSpaceValueObject
	deep *DeepSpaceCalculator
}

/**
 * DeepSpaceSatellite Constructor.
 *
 * @param tle
 *            the three line elements
 */
func NewDeepSpaceSatellite(tle TLE) Satellite {
	this := &DeepSpaceSatellite{
		super: *NewAbstractSatellite(tle),
	}

	// ////////////////////////////
	// initSDP4
	// ////////////////////////////

	/* Recover original mean motion (xnodp) and */
	/* semimajor axis (aodp) from input elements. */

	// recoverMeanMotionAndSemiMajorAxis

	a1 := math.Pow(XKE/tle.Xno, TWO_THIRDS)
	this.dsv.cosio = math.Cos(tle.Xincl)
	this.dsv.theta2 = this.dsv.cosio * this.dsv.cosio
	this.x3thm1 = 3.0*this.dsv.theta2 - 1
	this.dsv.eosq = tle.Eo * tle.Eo
	this.dsv.betao2 = 1.0 - this.dsv.eosq
	this.dsv.betao = math.Sqrt(this.dsv.betao2)
	del1 := 1.5 * CK2 * this.x3thm1 / (a1 * a1 * this.dsv.betao * this.dsv.betao2)
	ao := a1 * (1.0 - del1*(0.5*TWO_THIRDS+del1*(1.0+134/81*del1)))
	delo := 1.5 * CK2 * this.x3thm1 / (ao * ao * this.dsv.betao * this.dsv.betao2)
	this.dsv.xnodp = tle.Xno / (1.0 + delo)
	this.dsv.aodp = ao / (1.0 - delo)

	/* For perigee below 156 km, the values */
	/* of S and QOMS2T are altered. */
	this.super.setPerigee((this.dsv.aodp*(1.0-tle.Eo) - 1.0) * EARTH_RADIUS_KM)

	pinvsq := invert(this.dsv.aodp * this.dsv.aodp * this.dsv.betao2 * this.dsv.betao2)
	this.dsv.sing = math.Sin(tle.Omegao)
	this.dsv.cosg = math.Cos(tle.Omegao)
	tsi := invert(this.dsv.aodp - this.super.S4)
	eta := this.dsv.aodp * tle.Eo * tsi
	etasq := eta * eta
	eeta := tle.Eo * eta
	psisq := math.Abs(1.0 - etasq)
	coef := this.super.Qoms24 * math.Pow(tsi, 4)
	coef1 := coef / math.Pow(psisq, 3.5)
	c2 := coef1 * this.dsv.xnodp * (this.dsv.aodp*(1.0+1.5*etasq+eeta*(4.0+etasq)) + 0.75*CK2*tsi/psisq*this.x3thm1*(8.0+3.0*etasq*(8.0+etasq)))
	this.c1 = tle.Bstar * c2
	this.dsv.sinio = math.Sin(tle.Xincl)
	a3ovk2 := -J3_HARMONIC / CK2
	this.x1mth2 = 1.0 - this.dsv.theta2
	this.c4 = 2 * this.dsv.xnodp * coef1 * this.dsv.aodp * this.dsv.betao2 * (eta*(2.0+0.5*etasq) + tle.Eo*(0.5+2*etasq) - 2*CK2*tsi/(this.dsv.aodp*psisq)*(-3*this.x3thm1*(1.0-2*eeta+etasq*(1.5-0.5*eeta))+0.75*this.x1mth2*(2.0*etasq-eeta*(1.0+etasq))*math.Cos(2.0*tle.Omegao)))

	theta4 := this.dsv.theta2 * this.dsv.theta2
	temp1 := 3.0 * CK2 * pinvsq * this.dsv.xnodp
	temp2 := temp1 * CK2 * pinvsq
	temp3 := 1.25 * CK4 * pinvsq * pinvsq * this.dsv.xnodp
	this.dsv.xmdot = this.dsv.xnodp + 0.5*temp1*this.dsv.betao*this.x3thm1 + 0.0625*temp2*this.dsv.betao*(13-78*this.dsv.theta2+137*theta4)
	x1m5th := 1.0 - 5*this.dsv.theta2
	this.dsv.omgdot = -0.5*temp1*x1m5th + 0.0625*temp2*(7.0-114*this.dsv.theta2+395*theta4) + temp3*(3.0-36*this.dsv.theta2+49*theta4)
	xhdot1 := -temp1 * this.dsv.cosio
	this.dsv.xnodot = xhdot1 + (0.5*temp2*(4.0-19*this.dsv.theta2)+2*temp3*(3.0-7*this.dsv.theta2))*this.dsv.cosio
	this.xnodcf = 3.5 * this.dsv.betao2 * xhdot1 * this.c1
	this.t2cof = 1.5 * this.c1
	this.xlcof = 0.125 * a3ovk2 * this.dsv.sinio * (3.0 + 5*this.dsv.cosio) / (1.0 + this.dsv.cosio)
	this.aycof = 0.25 * a3ovk2 * this.dsv.sinio
	this.x7thm1 = 7.0*this.dsv.theta2 - 1

	/* initialize Deep() */
	this.deep = NewDeepSpaceCalculator(tle, &this.dsv)

	this.super.calculateSDP4orSGP4 = this.calculateSDP4

	return &this.super
}

// 	/**
// 	 * This function is used to calculate the position and velocity of
// 	 * deep-space (period > 225 minutes) satellites. tsince is time since epoch
// 	 * in minutes, tle is a pointer to a tle_t structure with Keplerian orbital
// 	 * elements and pos and vel are vector_t structures returning ECI satellite
// 	 * position and velocity. Use Convert_Sat_State() to convert to km and km/S.
// 	 *
// 	 * @param tsince
// 	 *            time since the epoch
// 	 * @param position
// 	 *            the position
// 	 * @param velocity
// 	 *            the velocity
// 	 * @param satPos
// 	 *            the position of the satellite
// 	 */
func (this *DeepSpaceSatellite) calculateSDP4(tsince float64) {
	temp := make([]float64, 12)
	xmdf := this.super.tle.Xmo + this.dsv.xmdot*tsince
	tsq := tsince * tsince
	templ := this.t2cof * tsq
	this.dsv.xll = xmdf + this.dsv.xnodp*templ

	this.dsv.omgadf = this.super.tle.Omegao + this.dsv.omgdot*tsince
	xnoddf := this.super.tle.Xnodeo + this.dsv.xnodot*tsince
	this.dsv.xnode = xnoddf + this.xnodcf*tsq
	tempa := 1.0 - this.c1*tsince
	tempe := this.super.tle.Bstar * this.c4 * tsince
	this.dsv.xn = this.dsv.xnodp

	this.dsv.t = tsince

	this.deep.dpsec(this.super.tle)

	a := math.Pow(XKE/this.dsv.xn, TWO_THIRDS) * tempa * tempa
	this.dsv.em = this.dsv.em - tempe
	this.deep.dpper()

	xl := this.dsv.xll + this.dsv.omgadf + this.dsv.xnode
	beta := math.Sqrt(1.0 - this.dsv.em*this.dsv.em)
	this.dsv.xn = XKE / math.Pow(a, 1.5)

	/* Long period periodics */
	axn := this.dsv.em * math.Cos(this.dsv.omgadf)
	temp[0] = invert(a * beta * beta)
	xll := temp[0] * this.xlcof * axn
	aynl := temp[0] * this.aycof
	xlt := xl + xll
	ayn := this.dsv.em*math.Sin(this.dsv.omgadf) + aynl

	/* Solve Kepler'S Equation */
	capu := mod2PI(xlt - this.dsv.xnode)
	temp[2] = capu

	converge(temp, axn, ayn, capu)
	this.calculatePositionAndVelocity(temp, a, axn, ayn)
	this.super.calculatePhase(xlt, this.dsv.xnode, this.dsv.omgadf)
}

func (this *DeepSpaceSatellite) calculatePositionAndVelocity(temp []float64, a float64, axn float64, ayn float64) {
	ecose := temp[5] + temp[6]
	esine := temp[3] - temp[4]
	elsq := axn*axn + ayn*ayn
	temp[0] = 1.0 - elsq
	pl := a * temp[0]
	temp[9] = a * (1.0 - ecose)
	temp[1] = invert(temp[9])
	temp[10] = XKE * math.Sqrt(a) * esine * temp[1]
	temp[11] = XKE * math.Sqrt(pl) * temp[1]
	temp[2] = a * temp[1]
	betal := math.Sqrt(temp[0])
	temp[3] = invert(1.0 + betal)
	cosu := temp[2] * (temp[8] - axn + ayn*esine*temp[3])
	sinu := temp[2] * (temp[7] - ayn - axn*esine*temp[3])
	u := math.Atan2(sinu, cosu)
	sin2u := 2.0 * sinu * cosu
	cos2u := 2.0*cosu*cosu - 1
	temp[0] = invert(pl)
	temp[1] = CK2 * temp[0]
	temp[2] = temp[1] * temp[0]

	/* Update for short periodics */
	rk := temp[9]*(1.0-1.5*temp[2]*betal*this.x3thm1) + 0.5*temp[1]*this.x1mth2*cos2u
	uk := u - 0.25*temp[2]*this.x7thm1*sin2u
	xnodek := this.dsv.xnode + 1.5*temp[2]*this.dsv.cosio*sin2u
	xinck := this.dsv.xinc + 1.5*temp[2]*this.dsv.cosio*this.dsv.sinio*cos2u
	rdotk := temp[10] - this.dsv.xn*temp[1]*this.x1mth2*sin2u
	rfdotk := temp[11] + this.dsv.xn*temp[1]*(this.x1mth2*cos2u+1.5*this.x3thm1)

	this.super.calculatePositionAndVelocity(rk, uk, xnodek, xinck, rdotk, rfdotk)
}

type DeepSpaceCalculator struct {
	/* This function is used by SDP4 to add lunar and solar */
	/* perturbation effects to deep-space orbit objects. */

	thgr   float64
	xnq    float64
	xqncl  float64
	omegaq float64
	zmol   float64
	zmos   float64

	savtsn float64
	ee2    float64
	e3     float64
	xi2    float64
	xl2    float64
	xl3    float64
	xl4    float64
	xgh2   float64
	xgh3   float64
	xgh4   float64
	xh2    float64
	xh3    float64
	sse    float64
	ssi    float64
	ssg    float64
	xi3    float64
	se2    float64
	si2    float64
	sl2    float64
	sgh2   float64
	sh2    float64
	se3    float64
	si3    float64
	sl3    float64
	sgh3   float64
	sh3    float64
	sl4    float64
	sgh4   float64
	ssl    float64
	ssh    float64
	d3210  float64
	d3222  float64
	d4410  float64
	d4422  float64
	d5220  float64
	d5232  float64
	d5421  float64
	d5433  float64
	del1   float64
	del2   float64
	del3   float64
	fasx2  float64
	fasx4  float64
	fasx6  float64
	xlamo  float64
	xfact  float64
	xni    float64
	atime  float64
	stepp  float64
	stepn  float64
	step2  float64
	preep  float64
	pl     float64
	sghs   float64
	xli    float64
	d2201  float64
	d2211  float64
	sghl   float64
	sh1    float64
	pinc   float64
	pe     float64
	shs    float64
	zsingl float64
	zcosgl float64
	zsinhl float64
	zcoshl float64
	zsinil float64
	zcosil float64

	a1     float64
	a2     float64
	a3     float64
	a4     float64
	a5     float64
	a6     float64
	a7     float64
	a8     float64
	a9     float64
	a10    float64
	ainv2  float64
	alfdp  float64
	aqnv   float64
	sgh    float64
	sini2  float64
	sinis  float64
	sinok  float64
	sh     float64
	si     float64
	sil    float64
	day    float64
	betdp  float64
	dalf   float64
	bfact  float64
	c      float64
	cc     float64
	cosis  float64
	cosok  float64
	cosq   float64
	ctem   float64
	f322   float64
	zx     float64
	zy     float64
	dbet   float64
	dls    float64
	eoc    float64
	eq     float64
	f2     float64
	f220   float64
	f221   float64
	f3     float64
	f311   float64
	f321   float64
	xnoh   float64
	f330   float64
	f441   float64
	f442   float64
	f522   float64
	f523   float64
	f542   float64
	f543   float64
	g200   float64
	g201   float64
	g211   float64
	pgh    float64
	ph     float64
	s1     float64
	s2     float64
	s3     float64
	s4     float64
	s5     float64
	s6     float64
	s7     float64
	se     float64
	sel    float64
	ses    float64
	xls    float64
	g300   float64
	g310   float64
	g322   float64
	g410   float64
	g422   float64
	g520   float64
	g521   float64
	g532   float64
	g533   float64
	gam    float64
	sinq   float64
	sinzf  float64
	sis    float64
	sl     float64
	sll    float64
	sls    float64
	stem   float64
	temp   float64
	temp1  float64
	x1     float64
	x2     float64
	x2li   float64
	x2omi  float64
	x3     float64
	x4     float64
	x5     float64
	x6     float64
	x7     float64
	x8     float64
	xl     float64
	xldot  float64
	xmao   float64
	xnddt  float64
	xndot  float64
	xno2   float64
	xnodce float64
	xnoi   float64
	xomi   float64
	xpidot float64
	z1     float64
	z11    float64
	z12    float64
	z13    float64
	z2     float64
	z21    float64
	z22    float64
	z23    float64
	z3     float64
	z31    float64
	z32    float64
	z33    float64
	ze     float64
	zf     float64
	zm     float64
	zn     float64
	zsing  float64
	zsinh  float64
	zsini  float64
	zcosg  float64
	zcosh  float64
	zcosi  float64
	delt   float64
	ft     float64

	resonance    bool
	synchronous  bool
	doLoop       bool
	epochRestart bool
	dsv          *DeepSpaceValueObject
}

func NewDeepSpaceCalculator(tle TLE, dsv *DeepSpaceValueObject) (this *DeepSpaceCalculator) {
	this = &DeepSpaceCalculator{
		dsv:    dsv,
		thgr:   dsv.thetaG(tle.Refepoch),
		eq:     tle.Eo,
		xnq:    dsv.xnodp,
		aqnv:   invert(dsv.aodp),
		xqncl:  tle.Xincl,
		xmao:   tle.Xmo,
		xpidot: dsv.omgdot + dsv.xnodot,
		sinq:   math.Sin(tle.Xnodeo),
		cosq:   math.Cos(tle.Xnodeo),
		omegaq: tle.Omegao,
		/* Initialize lunar solar terms */
		/* Days since 1900 Jan 0.5 */
		day: dsv.ds50 + 18261.5,
	}

	if math.Abs(this.day-this.preep) > 1.0E-6 {
		this.preep = this.day
		this.xnodce = 4.5236020 - 9.2422029E-4*this.day
		this.stem = math.Sin(this.xnodce)
		this.ctem = math.Cos(this.xnodce)
		this.zcosil = 0.91375164 - 0.03568096*this.ctem
		this.zsinil = math.Sqrt(1.0 - this.zcosil*this.zcosil)
		this.zsinhl = 0.089683511 * this.stem / this.zsinil
		this.zcoshl = math.Sqrt(1.0 - this.zsinhl*this.zsinhl)
		this.c = 4.7199672 + 0.22997150*this.day
		this.gam = 5.8351514 + 0.0019443680*this.day
		this.zmol = mod2PI(this.c - this.gam)
		this.zx = 0.39785416 * this.stem / this.zsinil
		this.zy = this.zcoshl*this.ctem + 0.91744867*this.zsinhl*this.stem
		this.zx = math.Atan2(this.zx, this.zy)
		this.zx = this.gam + this.zx - this.xnodce
		this.zcosgl = math.Cos(this.zx)
		this.zsingl = math.Sin(this.zx)
		this.zmos = mod2PI(6.2565837 + 0.017201977*this.day)
	} else {
		this.zmol = 0
		this.zmos = 0
	}

	/* Do solar terms */

	this.doSolarTerms()

	/* Geopotential resonance initialization for 12 hour orbits */
	this.resonance = false
	this.synchronous = false

	if !((this.xnq < 0.0052359877) && (this.xnq > 0.0034906585)) {
		if (this.xnq < 0.00826) || (this.xnq > 0.00924) {
			return
		}

		if this.eq < 0.5 {
			return
		}

		// calculateResonance

		this.resonance = true
		this.eoc = this.eq * this.dsv.eosq
		this.g201 = -0.306 - (this.eq-0.64)*0.440

		if this.eq <= 0.65 {
			this.g211 = 3.616 - 13.247*this.eq + 16.290*this.dsv.eosq
			this.g310 = -19.302 + 117.390*this.eq - 228.419*this.dsv.eosq + 156.591*this.eoc
			this.g322 = -18.9068 + 109.7927*this.eq - 214.6334*this.dsv.eosq + 146.5816*this.eoc
			this.g410 = -41.122 + 242.694*this.eq - 471.094*this.dsv.eosq + 313.953*this.eoc
			this.g422 = -146.407 + 841.880*this.eq - 1629.014*this.dsv.eosq + 1083.435*this.eoc
			this.g520 = -532.114 + 3017.977*this.eq - 5740*this.dsv.eosq + 3708.276*this.eoc
		} else {
			this.g211 = -72.099 + 331.819*this.eq - 508.738*this.dsv.eosq + 266.724*this.eoc
			this.g310 = -346.844 + 1582.851*this.eq - 2415.925*this.dsv.eosq + 1246.113*this.eoc
			this.g322 = -342.585 + 1554.908*this.eq - 2366.899*this.dsv.eosq + 1215.972*this.eoc
			this.g410 = -1052.797 + 4758.686*this.eq - 7193.992*this.dsv.eosq + 3651.957*this.eoc
			this.g422 = -3581.69 + 16178.11*this.eq - 24462.77*this.dsv.eosq + 12422.52*this.eoc

			if this.eq <= 0.715 {
				this.g520 = 1464.74 - 4664.75*this.eq + 3763.64*this.dsv.eosq
			} else {
				this.g520 = -5149.66 + 29936.92*this.eq - 54087.36*this.dsv.eosq + 31324.56*this.eoc
			}
		}

		if this.eq < 0.7 {
			this.g533 = -919.2277 + 4988.61*this.eq - 9064.77*this.dsv.eosq + 5542.21*this.eoc
			this.g521 = -822.71072 + 4568.6173*this.eq - 8491.4146*this.dsv.eosq + 5337.524*this.eoc
			this.g532 = -853.666 + 4690.25*this.eq - 8624.77*this.dsv.eosq + 5341.4*this.eoc
		} else {
			this.g533 = -37995.78 + 161616.52*this.eq - 229838.2*this.dsv.eosq + 109377.94*this.eoc
			this.g521 = -51752.104 + 218913.95*this.eq - 309468.16*this.dsv.eosq + 146349.42*this.eoc
			this.g532 = -40023.88 + 170470.89*this.eq - 242699.48*this.dsv.eosq + 115605.82*this.eoc
		}

		this.sini2 = this.dsv.sinio * this.dsv.sinio
		this.f220 = 0.75 * (1.0 + 2*dsv.cosio + dsv.theta2)
		this.f221 = 1.5 * this.sini2
		this.f321 = 1.875 * this.dsv.sinio * (1.0 - 2*this.dsv.cosio - 3.0*this.dsv.theta2)
		this.f322 = -1.875 * this.dsv.sinio * (1.0 + 2*this.dsv.cosio - 3.0*this.dsv.theta2)
		this.f441 = 35 * this.sini2 * this.f220
		this.f442 = 39.3750 * this.sini2 * this.sini2
		this.f522 = 9.84375 * this.dsv.sinio * (this.sini2*(1.0-2*this.dsv.cosio-5*this.dsv.theta2) + 0.33333333*(-2+4*this.dsv.cosio+6*this.dsv.theta2))
		this.f523 = this.dsv.sinio * (4.92187512*this.sini2*(-2-4*this.dsv.cosio+10*this.dsv.theta2) + 6.56250012*(1.0+2*this.dsv.cosio-3.0*this.dsv.theta2))
		this.f542 = 29.53125 * this.dsv.sinio * (2.0 - 8*this.dsv.cosio + this.dsv.theta2*(-12+8*this.dsv.cosio+10*this.dsv.theta2))
		this.f543 = 29.53125 * this.dsv.sinio * (-2 - 8*this.dsv.cosio + this.dsv.theta2*(12+8*this.dsv.cosio-10*this.dsv.theta2))
		this.xno2 = this.xnq * this.xnq
		this.ainv2 = this.aqnv * this.aqnv
		this.temp1 = 3.0 * this.xno2 * this.ainv2
		this.temp = this.temp1 * ROOT22
		this.d2201 = this.temp * this.f220 * this.g201
		this.d2211 = this.temp * this.f221 * this.g211
		this.temp1 = this.temp1 * this.aqnv
		this.temp = this.temp1 * ROOT32
		this.d3210 = this.temp * this.f321 * this.g310
		this.d3222 = this.temp * this.f322 * this.g322
		this.temp1 = this.temp1 * this.aqnv
		this.temp = 2.0 * this.temp1 * ROOT44
		this.d4410 = this.temp * this.f441 * this.g410
		this.d4422 = this.temp * this.f442 * this.g422
		this.temp1 = this.temp1 * this.aqnv
		this.temp = this.temp1 * ROOT52
		this.d5220 = this.temp * this.f522 * this.g520
		this.d5232 = this.temp * this.f523 * this.g532
		this.temp = 2.0 * this.temp1 * ROOT54
		this.d5421 = this.temp * this.f542 * this.g521
		this.d5433 = this.temp * this.f543 * this.g533
		this.xlamo = this.xmao + tle.Xnodeo + tle.Xnodeo - this.thgr - this.thgr
		this.bfact = this.dsv.xmdot + this.dsv.xnodot + this.dsv.xnodot - THDT - THDT
		this.bfact = this.bfact + this.ssl + this.ssh + this.ssh
	} else {

		// initSynchronousResonanceTerms

		this.resonance = true
		this.synchronous = true

		this.g200 = 1.0 + this.dsv.eosq*(-2.5+0.8125*this.dsv.eosq)
		this.g310 = 1.0 + 2*this.dsv.eosq
		this.g300 = 1.0 + this.dsv.eosq*(-6+6.60937*this.dsv.eosq)
		this.f220 = 0.75 * (1.0 + this.dsv.cosio) * (1.0 + this.dsv.cosio)
		this.f311 = 0.9375*this.dsv.sinio*this.dsv.sinio*(1.0+3.0*this.dsv.cosio) - 0.75*(1.0+this.dsv.cosio)
		this.f330 = 1.0 + this.dsv.cosio
		this.f330 = 1.875 * this.f330 * this.f330 * this.f330
		this.del1 = 3.0 * this.xnq * this.xnq * this.aqnv * this.aqnv
		this.del2 = 2.0 * this.del1 * this.f220 * this.g200 * Q22
		this.del3 = 3.0 * this.del1 * this.f330 * this.g300 * Q33 * this.aqnv
		this.del1 = this.del1 * this.f311 * this.g310 * Q31 * this.aqnv
		this.fasx2 = 0.13130908
		this.fasx4 = 2.8843198
		this.fasx6 = 0.37448087
		this.xlamo = this.xmao + tle.Xnodeo + tle.Omegao - this.thgr
		this.bfact = this.dsv.xmdot + this.xpidot - THDT
		this.bfact = this.bfact + this.ssl + this.ssg + this.ssh
	}

	this.xfact = this.bfact - this.xnq

	/* Initialize integrator */
	this.xli = this.xlamo
	this.xni = this.xnq
	this.atime = 0
	this.stepp = 720
	this.stepn = -720
	this.step2 = 259200

	return // implicit return of "this"
}

func (this *DeepSpaceCalculator) doSolarTerms() {
	this.savtsn = 1E20
	this.zcosg = 1.945905E-1
	this.zsing = ZSINGS
	this.zcosi = 9.1744867E-1
	this.zsini = ZSINIS
	this.zcosh = this.cosq
	this.zsinh = this.sinq
	this.cc = C1SS
	this.zn = ZNS
	this.ze = ZES
	this.xnoi = invert(this.xnq)

	/* Solar terms done again after Lunar terms are done */
	this.calculateSolarTerms()

	/* Do lunar terms */
	this.calculateLunarTerms()

	this.calculateSolarTerms()

	this.sse = this.sse + this.se
	this.ssi = this.ssi + this.si
	this.ssl = this.ssl + this.sl
	this.ssg = this.ssg + this.sgh - this.dsv.cosio/this.dsv.sinio*this.sh
	this.ssh = this.ssh + this.sh/this.dsv.sinio
}

/**
 *
 */
func (this *DeepSpaceCalculator) calculateLunarTerms() {
	this.sse = this.se
	this.ssi = this.si
	this.ssl = this.sl
	this.ssh = this.sh / this.dsv.sinio
	this.ssg = this.sgh - this.dsv.cosio*this.ssh
	this.se2 = this.ee2
	this.si2 = this.xi2
	this.sl2 = this.xl2
	this.sgh2 = this.xgh2
	this.sh2 = this.xh2
	this.se3 = this.e3
	this.si3 = this.xi3
	this.sl3 = this.xl3
	this.sgh3 = this.xgh3
	this.sh3 = this.xh3
	this.sl4 = this.xl4
	this.sgh4 = this.xgh4
	this.zcosg = this.zcosgl
	this.zsing = this.zsingl
	this.zcosi = this.zcosil
	this.zsini = this.zsinil
	this.zcosh = this.zcoshl*this.cosq + this.zsinhl*this.sinq
	this.zsinh = this.sinq*this.zcoshl - this.cosq*this.zsinhl
	this.zn = ZNL
	this.cc = C1L
	this.ze = ZEL
}

/**
 *
 */
func (this *DeepSpaceCalculator) calculateSolarTerms() {
	this.a1 = this.zcosg*this.zcosh + this.zsing*this.zcosi*this.zsinh
	this.a3 = -this.zsing*this.zcosh + this.zcosg*this.zcosi*this.zsinh
	this.a7 = -this.zcosg*this.zsinh + this.zsing*this.zcosi*this.zcosh
	this.a8 = this.zsing * this.zsini
	this.a9 = this.zsing*this.zsinh + this.zcosg*this.zcosi*this.zcosh
	this.a10 = this.zcosg * this.zsini
	this.a2 = this.dsv.cosio*this.a7 + this.dsv.sinio*this.a8
	this.a4 = this.dsv.cosio*this.a9 + this.dsv.sinio*this.a10
	this.a5 = -this.dsv.sinio*this.a7 + this.dsv.cosio*this.a8
	this.a6 = -this.dsv.sinio*this.a9 + this.dsv.cosio*this.a10
	this.x1 = this.a1*this.dsv.cosg + this.a2*this.dsv.sing
	this.x2 = this.a3*this.dsv.cosg + this.a4*this.dsv.sing
	this.x3 = -this.a1*this.dsv.sing + this.a2*this.dsv.cosg
	this.x4 = -this.a3*this.dsv.sing + this.a4*this.dsv.cosg
	this.x5 = this.a5 * this.dsv.sing
	this.x6 = this.a6 * this.dsv.sing
	this.x7 = this.a5 * this.dsv.cosg
	this.x8 = this.a6 * this.dsv.cosg
	this.z31 = 12*this.x1*this.x1 - 3.0*this.x3*this.x3
	this.z32 = 24*this.x1*this.x2 - 6*this.x3*this.x4
	this.z33 = 12*this.x2*this.x2 - 3.0*this.x4*this.x4
	this.z1 = 3.0*(this.a1*this.a1+this.a2*this.a2) + this.z31*this.dsv.eosq
	this.z2 = 6.0*(this.a1*this.a3+this.a2*this.a4) + this.z32*this.dsv.eosq
	this.z3 = 3.0*(this.a3*this.a3+this.a4*this.a4) + this.z33*this.dsv.eosq
	this.z11 = -6*this.a1*this.a5 + this.dsv.eosq*(-24*this.x1*this.x7-6*this.x3*this.x5)
	this.z12 = -6*(this.a1*this.a6+this.a3*this.a5) + this.dsv.eosq*(-24*(this.x2*this.x7+this.x1*this.x8)-6*(this.x3*this.x6+this.x4*this.x5))
	this.z13 = -6*this.a3*this.a6 + this.dsv.eosq*(-24*this.x2*this.x8-6*this.x4*this.x6)
	this.z21 = 6.0*this.a2*this.a5 + this.dsv.eosq*(24*this.x1*this.x5-6*this.x3*this.x7)
	this.z22 = 6.0*(this.a4*this.a5+this.a2*this.a6) + this.dsv.eosq*(24*(this.x2*this.x5+this.x1*this.x6)-6*(this.x4*this.x7+this.x3*this.x8))
	this.z23 = 6.0*this.a4*this.a6 + this.dsv.eosq*(24*this.x2*this.x6-6*this.x4*this.x8)
	this.z1 = this.z1 + this.z1 + this.dsv.betao2*this.z31
	this.z2 = this.z2 + this.z2 + this.dsv.betao2*this.z32
	this.z3 = this.z3 + this.z3 + this.dsv.betao2*this.z33
	this.s3 = this.cc * this.xnoi
	this.s2 = -0.5 * this.s3 / this.dsv.betao
	this.s4 = this.s3 * this.dsv.betao
	this.s1 = -15 * this.eq * this.s4
	this.s5 = this.x1*this.x3 + this.x2*this.x4
	this.s6 = this.x2*this.x3 + this.x1*this.x4
	this.s7 = this.x2*this.x4 - this.x1*this.x3
	this.se = this.s1 * this.zn * this.s5
	this.si = this.s2 * this.zn * (this.z11 + this.z13)
	this.sl = -this.zn * this.s3 * (this.z1 + this.z3 - 14 - 6*this.dsv.eosq)
	this.sgh = this.s4 * this.zn * (this.z31 + this.z33 - 6)
	this.sh = -this.zn * this.s2 * (this.z21 + this.z23)

	if this.xqncl < 5.2359877E-2 {
		this.sh = 0
	}

	this.ee2 = 2.0 * this.s1 * this.s6
	this.e3 = 2.0 * this.s1 * this.s7
	this.xi2 = 2.0 * this.s2 * this.z12
	this.xi3 = 2.0 * this.s2 * (this.z13 - this.z11)
	this.xl2 = -2 * this.s3 * this.z2
	this.xl3 = -2 * this.s3 * (this.z3 - this.z1)
	this.xl4 = -2 * this.s3 * (-21 - 9*this.dsv.eosq) * this.ze
	this.xgh2 = 2.0 * this.s4 * this.z32
	this.xgh3 = 2.0 * this.s4 * (this.z33 - this.z31)
	this.xgh4 = -18 * this.s4 * this.ze
	this.xh2 = -2 * this.s2 * this.z22
	this.xh3 = -2 * this.s2 * (this.z23 - this.z21)
}

/**
 * Entrance for deep space secular effects.
 *
 * @param tle
 *            the three line elements
 * @param dsv
 *            the deep space values
 */
func (this *DeepSpaceCalculator) dpsec(tle TLE) {
	this.dsv.xll = this.dsv.xll + this.ssl*this.dsv.t
	this.dsv.omgadf = this.dsv.omgadf + this.ssg*this.dsv.t
	this.dsv.xnode = this.dsv.xnode + this.ssh*this.dsv.t
	this.dsv.em = tle.Eo + this.sse*this.dsv.t
	this.dsv.xinc = tle.Xincl + this.ssi*this.dsv.t

	if this.dsv.xinc < 0 {
		this.dsv.xinc = -this.dsv.xinc
		this.dsv.xnode = this.dsv.xnode + math.Pi
		this.dsv.omgadf = this.dsv.omgadf - math.Pi
	}

	if !this.resonance {
		return
	}

	for {
		this.processEpochRestartLoop()
		if !this.doLoop || !this.epochRestart {
			break
		}
	}

	this.dsv.xn = this.xni + this.xndot*this.ft + this.xnddt*this.ft*this.ft*0.5
	this.xl = this.xli + this.xldot*this.ft + this.xndot*this.ft*this.ft*0.5
	this.temp = -this.dsv.xnode + this.thgr + this.dsv.t*THDT

	if this.synchronous {
		this.dsv.xll = this.xl - this.dsv.omgadf + this.temp
	} else {
		this.dsv.xll = this.xl + this.temp + this.temp
	}
}

/**
 *
 */
func (this *DeepSpaceCalculator) processEpochRestartLoop() {
	if (this.atime == 0) || ((this.dsv.t >= 0) && (this.atime < 0)) || ((this.dsv.t < 0) && (this.atime >= 0)) {
		/* Epoch restart */

		this.calclateDelt()

		this.atime = 0
		this.xni = this.xnq
		this.xli = this.xlamo
	} else if math.Abs(this.dsv.t) >= math.Abs(this.atime) {
		this.calclateDelt()
	}

	this.processNotEpochRestartLoop()
}

func (this *DeepSpaceCalculator) calclateDelt() {
	if this.dsv.t < 0 {
		this.delt = this.stepn
	} else {
		this.delt = this.stepp
	}
}

/**
*
 */
func (this *DeepSpaceCalculator) processNotEpochRestartLoop() {
	for {
		if math.Abs(this.dsv.t-this.atime) >= this.stepp {
			this.doLoop = true
			this.epochRestart = false
		} else {
			this.ft = this.dsv.t - this.atime
			this.doLoop = false
		}

		if math.Abs(this.dsv.t) < math.Abs(this.atime) {
			if this.dsv.t >= 0 {
				this.delt = this.stepn
			} else {
				this.delt = this.stepp
			}

			this.doLoop = this.doLoop || this.epochRestart
		}

		/* Dot terms calculated */
		if this.synchronous {
			this.xndot = this.del1*math.Sin(this.xli-this.fasx2) + this.del2*math.Sin(2.0*(this.xli-this.fasx4)) + this.del3*math.Sin(3.0*(this.xli-this.fasx6))
			this.xnddt = this.del1*math.Cos(this.xli-this.fasx2) + 2*this.del2*math.Cos(2.0*(this.xli-this.fasx4)) + 3.0*this.del3*math.Cos(3.0*(this.xli-this.fasx6))
		} else {
			this.xomi = this.omegaq + this.dsv.omgdot*this.atime
			this.x2omi = this.xomi + this.xomi
			this.x2li = this.xli + this.xli
			this.xndot = this.d2201*math.Sin(this.x2omi+this.xli-G22) + this.d2211*math.Sin(this.xli-G22) + this.d3210*math.Sin(this.xomi+this.xli-G32) + this.d3222*math.Sin(-this.xomi+this.xli-G32) + this.d4410*math.Sin(this.x2omi+this.x2li-G44) + this.d4422*math.Sin(this.x2li-G44) + this.d5220*math.Sin(this.xomi+this.xli-G52) + this.d5232*math.Sin(-this.xomi+this.xli-G52) + this.d5421*math.Sin(this.xomi+this.x2li-G54) + this.d5433*math.Sin(-this.xomi+this.x2li-G54)
			this.xnddt = this.d2201*math.Cos(this.x2omi+this.xli-G22) + this.d2211*math.Cos(this.xli-G22) + this.d3210*math.Cos(this.xomi+this.xli-G32) + this.d3222*math.Cos(-this.xomi+this.xli-G32) + this.d5220*math.Cos(this.xomi+this.xli-G52) + this.d5232*math.Cos(-this.xomi+this.xli-G52) + 2*(this.d4410*math.Cos(this.x2omi+this.x2li-G44)+this.d4422*math.Cos(this.x2li-G44)+this.d5421*math.Cos(this.xomi+this.x2li-G54)+this.d5433*math.Cos(-this.xomi+this.x2li-G54))
		}

		this.xldot = this.xni + this.xfact
		this.xnddt = this.xnddt * this.xldot

		if this.doLoop {
			this.xli = this.xli + this.xldot*this.delt + this.xndot*this.step2
			this.xni = this.xni + this.xndot*this.delt + this.xnddt*this.step2
			this.atime = this.atime + this.delt
		}

		if !this.doLoop || this.epochRestart {
			break
		}
	}
}

/**
 * Entrance for lunar-solar periodics.
 *
 * @param tle
 *            the three line elements
 * @param dsv
 *            the deep space values
 */
func (this *DeepSpaceCalculator) dpper() {
	this.sinis = math.Sin(this.dsv.xinc)
	this.cosis = math.Cos(this.dsv.xinc)

	if math.Abs(this.savtsn-this.dsv.t) >= 30 {
		this.savtsn = this.dsv.t
		this.zm = this.zmos + ZNS*this.dsv.t
		this.zf = this.zm + 2*ZES*math.Sin(this.zm)
		this.sinzf = math.Sin(this.zf)
		this.f2 = 0.5*this.sinzf*this.sinzf - 0.25
		this.f3 = -0.5 * this.sinzf * math.Cos(this.zf)
		this.ses = this.se2*this.f2 + this.se3*this.f3
		this.sis = this.si2*this.f2 + this.si3*this.f3
		this.sls = this.sl2*this.f2 + this.sl3*this.f3 + this.sl4*this.sinzf
		this.sghs = this.sgh2*this.f2 + this.sgh3*this.f3 + this.sgh4*this.sinzf
		this.shs = this.sh2*this.f2 + this.sh3*this.f3
		this.zm = this.zmol + ZNL*this.dsv.t
		this.zf = this.zm + 2*ZEL*math.Sin(this.zm)
		this.sinzf = math.Sin(this.zf)
		this.f2 = 0.5*this.sinzf*this.sinzf - 0.25
		this.f3 = -0.5 * this.sinzf * math.Cos(this.zf)
		this.sel = this.ee2*this.f2 + this.e3*this.f3
		this.sil = this.xi2*this.f2 + this.xi3*this.f3
		this.sll = this.xl2*this.f2 + this.xl3*this.f3 + this.xl4*this.sinzf
		this.sghl = this.xgh2*this.f2 + this.xgh3*this.f3 + this.xgh4*this.sinzf
		this.sh1 = this.xh2*this.f2 + this.xh3*this.f3
		this.pe = this.ses + this.sel
		this.pinc = this.sis + this.sil
		this.pl = this.sls + this.sll
	}

	this.pgh = this.sghs + this.sghl
	this.ph = this.shs + this.sh1
	this.dsv.xinc = this.dsv.xinc + this.pinc
	this.dsv.em = this.dsv.em + this.pe

	if this.xqncl >= 0.2 {
		/* Apply periodics directly */
		this.ph = this.ph / this.dsv.sinio
		this.pgh = this.pgh - this.dsv.cosio*this.ph
		this.dsv.omgadf = this.dsv.omgadf + this.pgh
		this.dsv.xnode = this.dsv.xnode + this.ph
		this.dsv.xll = this.dsv.xll + this.pl
	} else {
		this.applyPeriodics()

		/* This is a patch to Lyddane modification */
		/* suggested by Rob Matson. */
		if math.Abs(this.xnoh-this.dsv.xnode) > math.Pi {
			if this.dsv.xnode < this.xnoh {
				this.dsv.xnode += TWO_PI
			} else {
				this.dsv.xnode -= TWO_PI
			}
		}

		this.dsv.xll = this.dsv.xll + this.pl
		this.dsv.omgadf = this.xls - this.dsv.xll - math.Cos(this.dsv.xinc)*this.dsv.xnode
	}
}

/**
 * Apply periodics with Lyddane modification.
 *
 * @param dsv
 *            the space values
 */
func (this *DeepSpaceCalculator) applyPeriodics() {
	this.sinok = math.Sin(this.dsv.xnode)
	this.cosok = math.Cos(this.dsv.xnode)
	this.alfdp = this.sinis * this.sinok
	this.betdp = this.sinis * this.cosok
	this.dalf = this.ph*this.cosok + this.pinc*this.cosis*this.sinok
	this.dbet = -this.ph*this.sinok + this.pinc*this.cosis*this.cosok
	this.alfdp = this.alfdp + this.dalf
	this.betdp = this.betdp + this.dbet
	this.dsv.xnode = mod2PI(this.dsv.xnode)
	this.xls = this.dsv.xll + this.dsv.omgadf + this.cosis*this.dsv.xnode
	this.dls = this.pl + this.pgh - this.pinc*this.dsv.xnode*this.sinis
	this.xls = this.xls + this.dls
	this.xnoh = this.dsv.xnode
	this.dsv.xnode = math.Atan2(this.alfdp, this.betdp)
}

type DeepSpaceValueObject struct {
	eosq   float64
	sinio  float64
	cosio  float64
	betao  float64
	aodp   float64
	theta2 float64
	sing   float64
	cosg   float64
	betao2 float64
	xmdot  float64
	omgdot float64
	xnodot float64
	xnodp  float64

	/* Used by dpsec and dpper parts of Deep() */
	xll    float64
	omgadf float64
	xnode  float64
	em     float64
	xinc   float64
	xn     float64
	t      float64

	/* Used by thetg and Deep() */
	ds50 float64
}

/**
* The function ThetaG calculates the Greenwich Mean Sidereal Time for
* an epoch specified in the format used in the NORAD two-line element
* sets. It has now been adapted for dates beyond the year 1999, as
* described above. The function ThetaG_JD provides the same calculation
* except that it is based on an input in the form of a Julian Date.
*
* Reference: The 1992 Astronomical Almanac, page B6.
*
* @param epoch
*            the epach
* @param dsv
*            the deep space values
* @return the Greenwich Mean Sidereal Time
 */
func (dsv *DeepSpaceValueObject) thetaG(epoch float64) float64 {
	/* Modification to support Y2K */
	/* Valid 1957 through 2056 */
	year := math.Floor(epoch * 1E-3)
	dayOfYear := (epoch*1E-3 - year) * 1000.0

	if year < 57 {
		year = year + 2000
	} else {
		year = year + 1900
	}

	dayFloor := math.Floor(dayOfYear)
	dayFraction := dayOfYear - dayFloor
	dayOfYear = dayFloor

	jd := julianDateOfYear(year) + dayOfYear
	dsv.ds50 = jd - 2433281.5 + dayFraction

	return mod2PI(6.3003880987*dsv.ds50 + 1.72944494)
}
