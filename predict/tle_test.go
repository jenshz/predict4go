package predict

import (
	"testing"
	"fmt"
)

func TestNewTLE(t *testing.T) {
	hiber1, _ := NewTLE([]string{
		"HIBER-1                 ",
		"1 43744U 18096AB  19115.19815699  .00002003  00000-0  78676-4 0  9994",
		"2 43744  97.4641 185.2907 0018688 163.4737 196.7173 15.26755683 22421",
	})
	fmt.Println(hiber1)
	hiber2, _ := NewTLE([]string{
		"HIBER-2                 ",
		"1 43774U 18099S   19115.19835931  .00000287  00000-0  31551-4 0  9994",
		"2 43774  97.7511 187.8310 0013991 141.9695 218.2515 14.94862146 21274",
	})
	fmt.Println(hiber2)
}
