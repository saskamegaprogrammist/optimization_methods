package main

import (
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"math"
)

func constraintInt1Func(constraint1 func(xs []float64) float64, constraint2 func(xs []float64) float64, constraint3 func(xs []float64) float64) func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r * float64(1) /
			(constraint1(xs) + constraint2(xs) + constraint3(xs))
	}
}

func constraintInt2Func(constraint1 func(xs []float64) float64, constraint2 func(xs []float64) float64, constraint3 func(xs []float64) float64) func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r *
			(math.Log(-constraint1(xs)) + math.Log(-constraint2(xs)) + math.Log(-constraint3(xs)))
	}
}

func firstInt1ConstraintFirstGrad(xs []float64) float64 {
	return -math.Pow(firstConstraint(xs), -2) * 2 * xs[0]
}

func firstInt1ConstraintFirstGradFirst(xs []float64) float64 {
	return -2 * (math.Pow(firstConstraint(xs), -2) - 4*math.Pow(firstConstraint(xs), -3)*math.Pow(xs[0], 2))
}

func firstInt1ConstraintFirstGradSecond(xs []float64) float64 {
	return 8 * xs[0] * xs[1] * math.Pow(firstConstraint(xs), -3)
}

func firstInt1ConstraintSecondGrad(xs []float64) float64 {
	return -math.Pow(firstConstraint(xs), -2) * 2 * xs[1]
}

func firstInt1ConstraintSecondGradFirst(xs []float64) float64 {
	return 8 * xs[0] * xs[1] * math.Pow(firstConstraint(xs), -3)
}

func firstInt1ConstraintSecondGradSecond(xs []float64) float64 {
	return -2 * (math.Pow(firstConstraint(xs), -2) - 4*math.Pow(firstConstraint(xs), -3)*math.Pow(xs[1], 2))
}

func secondInt1ConstraintFirstGrad(xs []float64) float64 {
	return math.Pow(xs[0], -2)
}

func secondInt1ConstraintFirstGradFirst(xs []float64) float64 {
	return -2 * math.Pow(xs[0], -3)
}

func thirdInt1ConstraintSecondGrad(xs []float64) float64 {
	return math.Pow(xs[1], -2)
}

func thirdInt1ConstraintSecondGradSecond(xs []float64) float64 {
	return -2 * math.Pow(xs[1], -3)
}

func constraintInt1FuncFirstGrad() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r * (firstInt1ConstraintFirstGrad(xs) + secondInt1ConstraintFirstGrad(xs))
	}
}

func constraintInt1FuncSecondGrad() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r * (firstInt1ConstraintSecondGrad(xs) + thirdInt1ConstraintSecondGrad(xs))
	}
}

func constraintInt1FuncFirstGradFirst() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r * (firstInt1ConstraintFirstGradFirst(xs) + secondInt1ConstraintFirstGradFirst(xs))
	}
}

func constraintInt1FuncFirstGradSecond() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r * (firstInt1ConstraintFirstGradSecond(xs))
	}
}

func constraintInt1FuncSecondGradFirst() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r * (firstInt1ConstraintSecondGradFirst(xs))
	}
}

func constraintInt1FuncSecondGradSecond() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r * (firstInt1ConstraintSecondGradSecond(xs) + thirdInt1ConstraintSecondGradSecond(xs))
	}
}

func hessianConstraintInt1(xs []float64, r float64) la_methods.Matrix {
	var hess la_methods.Matrix
	dim := len(xs)
	hess.Init(dim, dim)
	hess.Points[0][0] = constraintInt1FuncFirstGradFirst()(xs, r)
	hess.Points[0][1] = constraintInt1FuncFirstGradSecond()(xs, r)
	hess.Points[1][0] = constraintInt1FuncSecondGradFirst()(xs, r)
	hess.Points[1][1] = constraintInt1FuncSecondGradSecond()(xs, r)
	return hess
}

func hessConstraintInt1() func(xs []float64, r float64) la_methods.Matrix {
	return func(xs []float64, r float64) la_methods.Matrix {
		return hessianConstraintInt1(xs, r)
	}
}

func firstInt2ConstraintFirstGrad(xs []float64) float64 {
	return 2 * xs[0] / firstConstraint(xs)
}

func firstInt2ConstraintFirstGradFirst(xs []float64) float64 {
	return (2*firstConstraint(xs) - 4*math.Pow(xs[0], 2)) / math.Pow(firstConstraint(xs), 2)
}

func firstInt2ConstraintFirstGradSecond(xs []float64) float64 {
	return (4 * xs[0] * xs[1]) / math.Pow(firstConstraint(xs), 2)
}

func firstInt2ConstraintSecondGrad(xs []float64) float64 {
	return 2 * xs[1] / firstConstraint(xs)
}

func firstInt2ConstraintSecondGradFirst(xs []float64) float64 {
	return (4 * xs[0] * xs[1]) / math.Pow(firstConstraint(xs), 2)
}

func firstInt2ConstraintSecondGradSecond(xs []float64) float64 {
	return (2*firstConstraint(xs) - 4*math.Pow(xs[1], 2)) / math.Pow(firstConstraint(xs), 2)
}

func secondInt2ConstraintFirstGrad(xs []float64) float64 {
	return float64(1) / xs[0]
}

func secondInt2ConstraintFirstGradFirst(xs []float64) float64 {
	return -1 * math.Pow(xs[0], -2)
}

func thirdInt2ConstraintSecondGrad(xs []float64) float64 {
	return float64(1) / xs[1]
}

func thirdInt2ConstraintSecondGradSecond(xs []float64) float64 {
	return -1 * math.Pow(xs[1], -2)
}

func constraintInt2FuncFirstGrad() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r * (firstInt2ConstraintFirstGrad(xs) + secondInt2ConstraintFirstGrad(xs))
	}
}

func constraintInt2FuncSecondGrad() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r * (firstInt2ConstraintSecondGrad(xs) + thirdInt2ConstraintSecondGrad(xs))
	}
}

func constraintInt2FuncFirstGradFirst() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r * (firstInt2ConstraintFirstGradFirst(xs) + secondInt2ConstraintFirstGradFirst(xs))
	}
}

func constraintInt2FuncFirstGradSecond() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r * (firstInt2ConstraintFirstGradSecond(xs))
	}
}

func constraintInt2FuncSecondGradFirst() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r * (firstInt2ConstraintSecondGradFirst(xs))
	}
}

func constraintInt2FuncSecondGradSecond() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return -r * (firstInt2ConstraintSecondGradSecond(xs) + thirdInt2ConstraintSecondGradSecond(xs))
	}
}

func hessianConstraintInt2(xs []float64, r float64) la_methods.Matrix {
	var hess la_methods.Matrix
	dim := len(xs)
	hess.Init(dim, dim)
	hess.Points[0][0] = constraintInt2FuncFirstGradFirst()(xs, r)
	hess.Points[0][1] = constraintInt2FuncFirstGradSecond()(xs, r)
	hess.Points[1][0] = constraintInt2FuncSecondGradFirst()(xs, r)
	hess.Points[1][1] = constraintInt2FuncSecondGradSecond()(xs, r)
	return hess
}

func hessConstraintInt2() func(xs []float64, r float64) la_methods.Matrix {
	return func(xs []float64, r float64) la_methods.Matrix {
		return hessianConstraintInt2(xs, r)
	}
}
