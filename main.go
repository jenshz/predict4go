package main

import (
	"./predict"
	"fmt"
	"time"
	_ "math"
)

func main() {
	hiber1_tle, _ := predict.NewTLE([]string{
		"HIBER-1                 ",
		"1 43744U 18096AB  19115.19815699  .00002003  00000-0  78676-4 0  9994",
		"2 43744  97.4641 185.2907 0018688 163.4737 196.7173 15.26755683 22421",
	})

	hiber1 := predict.NewSatellite(hiber1_tle)

	tstart := time.Now()

	for i := 0; i < 60; i++ {
		hiber1.CalculateSatelliteVectors(tstart.Add(time.Duration(i) * time.Minute))
		pos := hiber1.CalculateSatelliteGroundTrack()
		fmt.Println(pos.String())
		//fmt.Printf("%f,%f,%d,#0000%02x\n", pos.Latitude / (math.Pi * 2.0) * 360, pos.Longitude / (math.Pi * 2.0) * 360, i, i)
	}
}
