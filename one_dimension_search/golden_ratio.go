package one_dimension_search

import (
	"fmt"
	"math"
)

var GLDNRT = (3 - math.Sqrt(5)) / float64(2)

type GoldenRatioSearch struct {
	aStart     float64
	bStart     float64
	precision  float64
	targetFunc func(x float64) float64
	k          int
	r          float64
}

func (gr *GoldenRatioSearch) Init(aStart float64, bStart float64, precision float64, targetFunc func(x float64) float64) {
	gr.aStart = aStart
	gr.bStart = bStart
	gr.precision = precision
	gr.targetFunc = targetFunc
}

func (gr *GoldenRatioSearch) Solve() (float64, float64) {
	var k int
	var delta float64
	var y, z float64
	var fY, fZ float64
	a := gr.aStart
	b := gr.bStart
	y = a + GLDNRT*(b-a)
	z = a + b - y
	for {
		fY = gr.targetFunc(y)
		fZ = gr.targetFunc(z)
		if checkDecreasing(fZ, fY) {
			a = y
			y = z
			z = a + b - z
		} else {
			b = z
			z = y
			y = a + b - y
		}
		delta = math.Abs(a - b)
		if delta <= gr.precision {
			//fmt.Printf("k value: %d\n", k)
			gr.k = k
			xMid := (a + b) / 2
			return xMid, gr.targetFunc(xMid)
		} else {
			k++
		}
	}
}

func (gr *GoldenRatioSearch) CountConvergence() (float64, error) {
	if gr.k != 0 {
		return math.Pow((1+math.Sqrt(5))/float64(2)-1, float64(gr.k)-1), nil
	}
	return 0, fmt.Errorf("the algorithm hasn't been run")
}
