package main

import "github.com/saskamegaprogrammist/optimization_methods/la_methods"

func A(rows int, columns int) func(xs []float64) la_methods.Matrix {
	return func(xs []float64) la_methods.Matrix {
		var A la_methods.Matrix
		A.Init(rows, columns)
		A.Points[0][0] = firstConstraintGradFirst(xs)
		A.Points[0][1] = firstConstraintGradSecond(xs)
		A.Points[1][0] = secondConstraintGradFirst(xs)
		A.Points[1][1] = secondConstraintGradSecond(xs)
		A.Points[2][0] = thirdConstraintGradFirst(xs)
		A.Points[2][1] = thirdConstraintGradSecond(xs)
		return A
	}
}
