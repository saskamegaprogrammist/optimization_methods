package main

import "math"

func firstForIdeal(xs []float64) float64 {
	return 3*math.Pow(xs[0], 2) + 2*xs[0]*xs[1] + math.Pow(xs[1], 2) - xs[0] - xs[1]
}

func firstForIdealGradFirst(xs []float64) float64 {
	return 3*2*xs[0] + 2*xs[1] + math.Pow(xs[1], 2) - 1
}

func firstForIdealGradSecond(xs []float64) float64 {
	return 2*xs[0] + 2*xs[1] - 1
}

func firstForIdealGradThird(xs []float64) float64 {
	return 0
}

func secondForIdeal(xs []float64) float64 {
	return 4*math.Pow(xs[0]-3, 2) + math.Pow(xs[1]-1, 2) + 6*math.Pow(xs[2], 2)
}

func secondForIdealGradFirst(xs []float64) float64 {
	return 2 * 4 * (xs[0] - 3)
}

func secondForIdealGradSecond(xs []float64) float64 {
	return 2 * (xs[1] - 1)
}

func secondForIdealGradThird(xs []float64) float64 {
	return 6 * 2 * xs[2]
}

func funcForIdeal(xs []float64, ideal []float64, ws []float64) float64 {
	return ws[0]*math.Pow(firstForIdeal(xs)-ideal[0], 2) + ws[1]*math.Pow(secondForIdeal(xs)-ideal[1], 2)
}

func funcForIdealGradFirst(xs []float64, ideal []float64, ws []float64) float64 {
	return ws[0]*2*(firstForIdeal(xs)-ideal[0])*firstForIdealGradFirst(xs) + ws[1]*2*(secondForIdeal(xs)-ideal[1])*secondForIdealGradFirst(xs)
}

func funcForIdealGradSecond(xs []float64, ideal []float64, ws []float64) float64 {
	return ws[0]*2*(firstForIdeal(xs)-ideal[0])*firstForIdealGradSecond(xs) + ws[1]*2*(secondForIdeal(xs)-ideal[1])*secondForIdealGradSecond(xs)
}

func funcForIdealGradThird(xs []float64, ideal []float64, ws []float64) float64 {
	return ws[0]*2*(firstForIdeal(xs)-ideal[0])*firstForIdealGradThird(xs) + ws[1]*2*(secondForIdeal(xs)-ideal[1])*secondForIdealGradThird(xs)
}
