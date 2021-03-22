package one_dimension_search

import (
	"fmt"
	"math"
)

type Svenn struct {
	tStep      float64
	startPoint float64
	targetFunc func(x float64) float64
}

func (sv *Svenn) Init(tStep float64, startPoint float64, targetFunc func(x float64) float64) {
	sv.startPoint = startPoint
	sv.tStep = tStep
	sv.targetFunc = targetFunc
}

func (sv *Svenn) SetStartPoint(startPoint float64) {
	sv.startPoint = startPoint
}

func (sv *Svenn) Solve() (float64, float64, error) {
	var k float64
	x := sv.startPoint
	a := sv.startPoint - sv.tStep
	b := sv.startPoint + sv.tStep
	f1 := sv.targetFunc(a)
	f2 := sv.targetFunc(x)
	f3 := sv.targetFunc(b)
	if checkStopCond(f1, f2, f3) {
		return a, b, nil
	}
	if checkInvalidCond(f1, f2, f3) {
		return a, b, fmt.Errorf("function is not unimodal, select another start point")
	}
	var delta float64
	k = 1
	if checkFirstCond(f1, f2, f3) {
		delta = sv.tStep
		a = x
		x += sv.tStep
	}
	if checkSecondCond(f1, f2, f3) {
		delta = -sv.tStep
		b = x
		x -= sv.tStep
	}
	for {
		xNext := x + math.Pow(2, k)*delta
		fNext := sv.targetFunc(xNext)
		f := sv.targetFunc(x)
		if checkDecreasing(fNext, f) {
			if delta > 0 {
				a = x
			} else {
				b = x
			}
			x = xNext
		} else {
			if delta > 0 {
				b = xNext
			} else {
				a = xNext
			}
			break
		}
		k++
	}
	//fmt.Printf("k value: %v\n", k)
	return a, b, nil
}

func checkDecreasing(f1 float64, f2 float64) bool {
	return f1 < f2
}

func checkStopCond(f1 float64, f2 float64, f3 float64) bool {
	return f2 <= f1 && f2 <= f3
}

func checkInvalidCond(f1 float64, f2 float64, f3 float64) bool {
	return f2 >= f1 && f2 >= f3
}

func checkFirstCond(f1 float64, f2 float64, f3 float64) bool {
	return f2 <= f1 && f2 >= f3
}

func checkSecondCond(f1 float64, f2 float64, f3 float64) bool {
	return f2 >= f1 && f2 <= f3
}
