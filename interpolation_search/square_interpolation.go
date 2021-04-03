package interpolation_search

import (
	"math"
)

type SquareInterpolation struct {
	startPoint     float64
	step           float64
	precisionDelta float64
	precisionEps   float64
	targetFunc     func(x float64) float64
}

func (sqInt *SquareInterpolation) Init(startPoint float64, step float64, precisionDelta float64, precisionEps float64, targetFunc func(x float64) float64) {
	sqInt.startPoint = startPoint
	sqInt.step = step
	sqInt.precisionDelta = precisionDelta
	sqInt.precisionEps = precisionEps
	sqInt.targetFunc = targetFunc
}

func (sqInt *SquareInterpolation) Solve() (float64, float64) {
	var alph1, alph2, alph3 float64
	var f1, f2, f3 float64
	var min, minAlph, interpol, interpolAlph float64
	var isNotZero bool
	alph1 = sqInt.startPoint
OUTER:
	for {
		alph2 = alph1 + sqInt.step
		f1 = sqInt.targetFunc(alph1)
		f2 = sqInt.targetFunc(alph2)
		if checkLess(f2, f1) {
			alph3 = alph1 + 2*sqInt.step
		} else {
			alph3 = alph1 - 2*sqInt.step
		}
		f3 = sqInt.targetFunc(alph3)
		for {
			min, minAlph = sqInt.findMin(f1, f2, f3, alph1, alph2, alph3)
			isNotZero, interpol, interpolAlph = sqInt.interpolate(f1, f2, f3, alph1, alph2, alph3)
			if !isNotZero {
				alph1 = minAlph
				continue OUTER
			}
			if sqInt.checkFirstCond(min, interpol) && sqInt.checkSecondCond(minAlph, interpolAlph) {
				return interpolAlph, interpol
			} else if checkInInterval(interpolAlph, alph1, alph3) {
				alph1, f1, alph2, f2, alph3, f3 = sqInt.findBest(alph1, alph2, alph3,
					min, interpol,
					minAlph, interpolAlph)
			} else {
				alph1 = interpolAlph
				continue OUTER
			}
		}
	}
}

func checkLess(f1 float64, f2 float64) bool {
	return f1 < f2
}

func checkLessOrEq(f1 float64, f2 float64) bool {
	return f1 <= f2
}

func checkInInterval(alph float64, alph1 float64, alph2 float64) bool {
	return alph >= alph1 && alph <= alph2
}

func (sqInt *SquareInterpolation) checkFirstCond(min float64, interpol float64) bool {
	return math.Abs((min-interpol)/interpol) < sqInt.precisionEps
}

func (sqInt *SquareInterpolation) checkSecondCond(minAlph float64, interpolAlph float64) bool {
	return math.Abs((minAlph-interpolAlph)/interpolAlph) < sqInt.precisionDelta
}

func (sqInt *SquareInterpolation) interpolate(f1 float64, f2 float64, f3 float64, alph1 float64, alph2 float64, alph3 float64) (bool, float64, float64) {
	numerator := (math.Pow(alph2, 2)-math.Pow(alph3, 2))*f1 + (math.Pow(alph3, 2)-math.Pow(alph1, 2))*f2 + (math.Pow(alph1, 2)-math.Pow(alph2, 2))*f3
	denominator := 2 * ((alph2-alph3)*f1 + (alph3-alph1)*f2 + (alph1-alph2)*f3)
	if denominator == 0 {
		return false, 0, 0
	}
	alphInterpolate := numerator / denominator
	return true, sqInt.targetFunc(alphInterpolate), alphInterpolate
}

func (sqInt *SquareInterpolation) findMin(f1 float64, f2 float64, f3 float64, alph1 float64, alph2 float64, alph3 float64) (float64, float64) {
	var minInterm, minAlphInterm float64
	if f1 < f2 {
		minInterm = f1
		minAlphInterm = alph1
	} else {
		minInterm = f2
		minAlphInterm = alph2
	}
	if f3 < minInterm {
		minInterm = f3
		minAlphInterm = alph3
	}
	return minInterm, minAlphInterm
}

func (sqInt *SquareInterpolation) findBest(
	alph1 float64, alph2 float64, alph3 float64,
	min float64, interpol float64,
	minAlph float64, interpolAlph float64) (float64, float64, float64, float64, float64, float64) {
	var opt, optAlph float64
	var left, leftAlph float64
	var right, rightAlph float64
	var dists []float64
	var xs []float64
	if min < interpol {
		opt = min
		optAlph = minAlph
		dists = []float64{math.Abs(minAlph - alph1), math.Abs(minAlph - alph2), math.Abs(minAlph - alph3), math.Abs(minAlph - interpolAlph)}
		xs = []float64{alph1, alph2, alph3, interpolAlph}

	} else {
		opt = interpol
		optAlph = interpolAlph
		dists = []float64{math.Abs(minAlph - alph1), math.Abs(minAlph - alph2), math.Abs(minAlph - alph3)}
		xs = []float64{alph1, alph2, alph3}
	}
	minDist0 := dists[0]
	minDist1 := dists[0]
	indexes := make([]int, 2)
	for i, dist := range dists {
		if dist != 0 && dist <= minDist0 {
			if xs[i] < optAlph {
				minDist0 = dist
				indexes[0] = i
			}
		}
		if dist != 0 && dist <= minDist1 {
			if xs[i] > optAlph {
				minDist1 = dist
				indexes[1] = i
			}
		}
	}
	left, leftAlph = xs[indexes[0]], sqInt.targetFunc(xs[indexes[0]])
	right, rightAlph = xs[indexes[1]], sqInt.targetFunc(xs[indexes[1]])
	return left, leftAlph, opt, optAlph, right, rightAlph
}
