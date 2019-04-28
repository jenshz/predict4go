package main

import (
	"./predict"
	"fmt"
	"time"
	"math"
)

func main() {
	hiber1_tle, _ := predict.NewTLE([]string{
		"HIBER-1                 ",
		"1 43744U 18096AB  19115.19815699  .00002003  00000-0  78676-4 0  9994",
		"2 43744  97.4641 185.2907 0018688 163.4737 196.7173 15.26755683 22421",
	})

	fmt.Println(hiber1_tle)

	hiber1 := predict.NewSatellite(hiber1_tle)

	tstart := time.Now()

	fmt.Println(tstart, tstart.Unix())

	hiber1.CalculateSatelliteVectors(tstart)
	pos := hiber1.CalculateSatelliteGroundTrack()
	fmt.Println(pos.String())
	fmt.Printf("%f,%f\n", pos.Latitude / (math.Pi * 2.0) * 360, math.Mod(pos.Longitude / (math.Pi * 2.0) * 360 + 180, 360.0) - 180)
}
