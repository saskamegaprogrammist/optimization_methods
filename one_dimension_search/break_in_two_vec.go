package one_dimension_search

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"math"
)

type BreakInTwoSearchVector struct {
	aStart     []float64
	bStart     []float64
	precision  float64
	targetFunc func(x []float64) float64
	k          int
	r          float64
}

func (bit *BreakInTwoSearchVector) Init(aStart []float64, bStart []float64, precision float64, targetFunc func(xs []float64) float64) {
	bit.aStart = aStart
	bit.bStart = bStart
	bit.precision = precision
	bit.targetFunc = targetFunc
}

func (bit *BreakInTwoSearchVector) Solve() ([]float64, float64, error) {
	var err error
	var k int
	var fMid, fY, fZ float64
	var aVec, bVec la_methods.Vector
	var yVec, zVec la_methods.Vector
	err = aVec.InitWithPoints(len(bit.aStart), bit.aStart)
	if err != nil {
		return []float64{}, 0, fmt.Errorf("error initializing vector: %v", err)
	}
	err = bVec.InitWithPoints(len(bit.bStart), bit.bStart)
	if err != nil {
		return []float64{}, 0, fmt.Errorf("error initializing vector: %v", err)
	}
	sumAB, err := aVec.Add(bVec)
	if err != nil {
		return []float64{}, 0, fmt.Errorf("error adding vector: %v", err)
	}
	sumAB.MulOnValue(float64(1) / float64(2))
	xMid := sumAB
	AB, err := bVec.Sub(aVec)
	if err != nil {
		return []float64{}, 0, fmt.Errorf("error substracting vector: %v", err)
	}
	length := AB.Len()
	for {
		fMid = bit.targetFunc(xMid.Points)
		yVec = aVec.AddK(length / 4)
		zVec = bVec.SubK(length / 4)
		fY = bit.targetFunc(yVec.Points)
		if checkDecreasing(fY, fMid) {
			bVec = xMid
			xMid = yVec
		} else {
			fZ = bit.targetFunc(zVec.Points)
			if checkDecreasing(fZ, fMid) {
				aVec = xMid
				xMid = zVec
			} else {
				aVec = yVec
				bVec = zVec
			}
		}
		AB, err := bVec.Sub(aVec)
		if err != nil {
			return []float64{}, 0, fmt.Errorf("error substracting vector: %v", err)
		}
		length := AB.Len()
		if length <= bit.precision {
			//fmt.Printf("k value: %d\n", k)
			bit.k = k
			return xMid.Points, bit.targetFunc(xMid.Points), nil
		} else {
			k++
		}
	}
}

func (bit *BreakInTwoSearchVector) CountConvergence() (float64, error) {
	if bit.k != 0 {
		return float64(1) / math.Pow(2, float64(bit.k)/2), nil
	}
	return 0, fmt.Errorf("the algorithm hasn't been run")
}
