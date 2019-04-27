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
    "strings"
    "strconv"
    "errors"
)

const (
    THREELINES = 3
    MINS_PERDAY = 1.44E3
)

/**
 * TLE representation to aid SGP4 calculations. Instances of this class are
 * immutable and thus thread safe.
 */
type TLE struct {
    Catnum int
    Name string
    Setnum int
    Year int
    Refepoch float64
    Incl float64
    Raan float64
    Eccn float64
    Argper float64
    Meanan float64
    Meanmo float64
    Drag float64
    Nddot6 float64
    Bstar float64
    Orbitnum int
    Epoch float64
    Xndt2o float64
    Xincl float64
    Xnodeo float64
    Eo float64
    Omegao float64
    Xmo float64
    Xno float64
    Deepspace bool
}

func NewTLE(tle []string) (*TLE, error) {
    if len(tle) != THREELINES {
        return nil, errors.New("TLE not 3 lines long.")
    }

    parseInt := func(str string) int {
        i, err := strconv.ParseInt(strings.TrimSpace(str), 10, 32)
        if err != nil {
            panic(err)
        }
        return int(i)
    }

    parseDouble := func(str string) float64 {
        f, err := strconv.ParseFloat(strings.TrimSpace(str), 64)
        if err != nil {
            panic(err)
        }
        return f
    }

    x := &TLE{}

    x.Catnum = parseInt(tle[1][2:7])
    x.Name = strings.TrimSpace(tle[0])
    x.Setnum = parseInt(tle[1][64:68])
    x.Year = parseInt(tle[1][18:20])
    x.Refepoch = parseDouble(tle[1][20:32])
    x.Incl = parseDouble(tle[2][8:16])
    x.Raan = parseDouble(tle[2][17:25])
    x.Eccn = 1.0e-07 * parseDouble(tle[2][26:33])
    x.Argper = parseDouble(tle[2][34:42])
    x.Meanan = parseDouble(tle[2][43:51])
    x.Meanmo = parseDouble(tle[2][52:63])
    x.Drag = parseDouble(tle[1][33:43])
    tempnum := 1.0e-5 * parseDouble(tle[1][44:50]);
    x.Nddot6 = tempnum / math.Pow(10.0, parseDouble(tle[1][51:52]))

    tempnum = 1.0e-5 * parseDouble(tle[1][53:59])
    x.Bstar = tempnum / math.Pow(10.0, parseDouble(tle[1][60:61]))
    x.Orbitnum = parseInt(strings.TrimSpace(tle[2][63:68]))

    /* reassign the values to these which get used in calculations */
    x.Epoch = (1000.0 * float64(x.Year)) + x.Refepoch

    temp := x.Incl
    temp *= DEG2RAD
    x.Xincl = temp

    temp = x.Raan
    temp *= DEG2RAD
    x.Xnodeo = temp

    x.Eo = x.Eccn

    temp = x.Argper
    temp *= DEG2RAD
    x.Omegao = temp

    temp = x.Meanan
    temp *= DEG2RAD
    x.Xmo = temp

    /* Preprocess tle set */
    {
        var temp float64
        temp = TWO_PI / MINS_PERDAY / MINS_PERDAY
        x.Xno = x.Meanmo * temp * MINS_PERDAY
        x.Xndt2o = x.Drag * temp

        dd1 := XKE / x.Xno
        a1 := math.Pow(dd1, TWO_THIRDS)
        r1 := math.Cos(x.Xincl)
        dd1 = 1.0 - x.Eo * x.Eo
        temp = CK2 * 1.5 * (r1 * r1 * 3.0 - 1.0) / math.Pow(dd1, 1.5)
        del1 := temp / (a1 * a1);
        ao := a1 * (1.0 - del1 * (TWO_THIRDS * 0.5 + del1 * (del1 * 1.654320987654321 + 1.0)))
        delo := temp / (ao * ao)
        xnodp := x.Xno / (delo + 1.0)

        /* Select a deep-space/near-earth ephemeris */
        x.Deepspace = TWO_PI / xnodp / MINS_PERDAY >= 0.15625
    }

    return x, nil
}
