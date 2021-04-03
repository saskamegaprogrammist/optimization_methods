package main

import "math"

func firstConstraint(xs []float64) float64 {
	return math.Pow(xs[0], 2) + math.Pow(xs[1], 2) - 10
}

func secondConstraint(xs []float64) float64 {
	return -xs[0]
}

func thirdConstraint(xs []float64) float64 {
	return -xs[1]
}

func firstConstraintGradFirst(xs []float64) float64 {
	return 2 * xs[0]
}

func firstConstraintGradSecond(xs []float64) float64 {
	return 2 * xs[1]
}

func secondConstraintGradFirst(xs []float64) float64 {
	return -1
}

func secondConstraintGradSecond(xs []float64) float64 {
	return 0
}

func thirdConstraintGradFirst(xs []float64) float64 {
	return 0
}

func thirdConstraintGradSecond(xs []float64) float64 {
	return -1
}
