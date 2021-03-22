package one_dimension_search

import (
	"fmt"
	"math"
)

type BreakInTwoSearch struct {
	aStart     float64
	bStart     float64
	precision  float64
	targetFunc func(x float64) float64
	k          int
	r          float64
}

func (bit *BreakInTwoSearch) Init(aStart float64, bStart float64, precision float64, targetFunc func(x float64) float64) {
	bit.aStart = aStart
	bit.bStart = bStart
	bit.precision = precision
	bit.targetFunc = targetFunc
}

func (bit *BreakInTwoSearch) Solve() (float64, float64) {
	var k int
	var y, z float64
	var fMid, fY, fZ float64
	a := bit.aStart
	b := bit.bStart
	xMid := (a + b) / 2
	length := b - a
	for {
		fMid = bit.targetFunc(xMid)
		y = a + length/4
		z = b - length/4
		fY = bit.targetFunc(y)
		if checkDecreasing(fY, fMid) {
			b = xMid
			xMid = y
		} else {
			fZ = bit.targetFunc(z)
			if checkDecreasing(fZ, fMid) {
				a = xMid
				xMid = z
			} else {
				a = y
				b = z
			}
		}
		length = b - a
		if length <= bit.precision {
			//fmt.Printf("k value: %d\n", k)
			bit.k = k
			return xMid, bit.targetFunc(xMid)
		} else {
			k++
		}
	}
}

func (bit *BreakInTwoSearch) CountConvergence() (float64, error) {
	if bit.k != 0 {
		return float64(1) / math.Pow(2, float64(bit.k)/2), nil
	}
	return 0, fmt.Errorf("the algorithm hasn't been run")
}
