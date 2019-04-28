package main

import (
	"github.com/jenshz/predict4go/predict"
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

	hiber1 := predict.NewSatellite(hiber1_tle)

	tstart := time.Now()

	// Calculate current position
	hiber1.CalculateSatelliteVectors(tstart)
	pos := hiber1.CalculateSatelliteGroundTrack()
	fmt.Println(pos.String())
	fmt.Printf("%f,%f\n", pos.Latitude / (math.Pi * 2.0) * 360, math.Mod(pos.Longitude / (math.Pi * 2.0) * 360 + 180, 360.0) - 180)


	qth := &predict.GroundStationPosition{
		Latitude: 52.3547321,
		Longitude: 4.8284119,
		HeightAMSL: 10.0,
		HorizonElevations: [36]int{},
		Name: "Amsterdam",
	}

	// Calculate next passes over Amsterdam
	pp, err := predict.NewPassPredictor(hiber1_tle, qth)
	if err != nil {
		panic(err)
	}

	passes := pp.GetPasses(tstart, 48, false, 0)
	for _, pass := range passes {
		fmt.Println("================")
		fmt.Println(pass.String())
	}
}
