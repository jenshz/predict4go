package predict

import (
	"fmt"
)

type Position struct {
	Lat float64
	Lon float64
}

func (p *Position) String() string {
	return fmt.Sprintf("Position [lat=%f, lon=%f]", p.Lat, p.Lon)
}
