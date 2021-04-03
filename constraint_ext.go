package main

import (
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"math"
)

func firstExtConstraintFirstGrad(xs []float64) float64 {
	return 2 * firstConstraint(xs) * 2 * xs[0]
}

func firstExtConstraintFirstGradMod(xs []float64) float64 {
	if firstConstraint(xs) > 0 {
		return firstExtConstraintFirstGrad(xs)
	}
	return 0
}

func firstExtConstraintSecondGrad(xs []float64) float64 {
	return 2 * firstConstraint(xs) * 2 * xs[1]
}

func firstExtConstraintSecondGradMod(xs []float64) float64 {
	if firstConstraint(xs) > 0 {
		return firstExtConstraintSecondGrad(xs)
	}
	return 0
}

func firstExtConstraintFirstGradFirst(xs []float64) float64 {
	return 2 * 2 * (firstConstraint(xs) + 2*math.Pow(xs[0], 2))
}

func firstExtConstraintFirstGradFirstMod(xs []float64) float64 {
	if firstConstraint(xs) > 0 {
		return firstExtConstraintFirstGradFirst(xs)
	}
	return 0
}

func firstExtConstraintFirstGradSecond(xs []float64) float64 {
	return 2 * 2 * 2 * xs[0] * xs[1]
}

func firstExtConstraintFirstGradSecondMod(xs []float64) float64 {
	if firstConstraint(xs) > 0 {
		return firstExtConstraintFirstGradSecond(xs)
	}
	return 0
}

func firstExtConstraintSecondGradFirst(xs []float64) float64 {
	return 2 * 2 * 2 * xs[0] * xs[1]
}

func firstExtConstraintSecondGradFirstMod(xs []float64) float64 {
	if firstConstraint(xs) > 0 {
		return firstExtConstraintSecondGradFirst(xs)
	}
	return 0
}

func firstExtConstraintSecondGradSecond(xs []float64) float64 {
	return 2 * 2 * (firstConstraint(xs) + 2*math.Pow(xs[1], 2))
}

func firstExtConstraintSecondGradSecondMod(xs []float64) float64 {
	if firstConstraint(xs) > 0 {
		return firstExtConstraintSecondGradSecond(xs)
	}
	return 0
}

func secondExtConstraintFirstGrad(xs []float64) float64 {
	return 2 * xs[0]
}

func secondExtConstraintFirstGradFirst(xs []float64) float64 {
	return 2
}

func secondExtConstraintFirstGradMod(xs []float64) float64 {
	if secondConstraint(xs) > 0 {
		return secondExtConstraintFirstGrad(xs)
	}
	return 0
}

func secondExtConstraintFirstGradFirstMod(xs []float64) float64 {
	if secondConstraint(xs) > 0 {
		return secondExtConstraintFirstGradFirst(xs)
	}
	return 0
}

func thirdExtConstraintSecondGrad(xs []float64) float64 {
	return 2 * xs[1]
}

func thirdExtConstraintSecondGradSecond(xs []float64) float64 {
	return 2
}

func thirdExtConstraintSecondGradMod(xs []float64) float64 {
	val := thirdConstraint(xs)
	if val > 0 {
		return thirdExtConstraintSecondGrad(xs)
	}
	return 0
}

func thirdExtConstraintSecondGradSecondMod(xs []float64) float64 {
	if thirdConstraint(xs) > 0 {
		return thirdExtConstraintSecondGradSecond(xs)
	}
	return 0
}

func firstExtConstraintMod(xs []float64) float64 {
	val := firstConstraint(xs)
	if val > 0 {
		return val
	}
	return 0
}

func secondExtConstraintMod(xs []float64) float64 {
	val := secondConstraint(xs)
	if val > 0 {
		return val
	}
	return 0
}

func thirdExtConstraintMod(xs []float64) float64 {
	val := thirdConstraint(xs)
	if val > 0 {
		return val
	}
	return 0
}

func constraintExtFunc(constraint1 func(xs []float64) float64, constraint2 func(xs []float64) float64, constraint3 func(xs []float64) float64) func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return r / float64(2) * (math.Pow(constraint1(xs), 2) + math.Pow(constraint2(xs), 2) + math.Pow(constraint3(xs), 2))
	}
}

func constraintExtFuncFirstGrad() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return r / float64(2) * (firstExtConstraintFirstGradMod(xs) + secondExtConstraintFirstGradMod(xs))
	}
}

func constraintExtFuncSecondGrad() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return r / float64(2) * (firstExtConstraintSecondGradMod(xs) + thirdExtConstraintSecondGradMod(xs))
	}
}

func constraintExtFuncFirstGradFirst() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return r / float64(2) * (firstExtConstraintFirstGradFirstMod(xs) + secondExtConstraintFirstGradFirstMod(xs))
	}
}

func constraintExtFuncFirstGradSecond() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return r / float64(2) * (firstExtConstraintFirstGradSecondMod(xs))
	}
}

func constraintExtFuncSecondGradFirst() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return r / float64(2) * (firstExtConstraintSecondGradFirstMod(xs))
	}
}

func constraintExtFuncSecondGradSecond() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return r / float64(2) * (firstExtConstraintSecondGradSecondMod(xs) + thirdExtConstraintSecondGradSecondMod(xs))
	}
}

func hessianConstraintExt(xs []float64, r float64) la_methods.Matrix {
	var hess la_methods.Matrix
	dim := len(xs)
	hess.Init(dim, dim)
	hess.Points[0][0] = constraintExtFuncFirstGradFirst()(xs, r)
	hess.Points[0][1] = constraintExtFuncFirstGradSecond()(xs, r)
	hess.Points[1][0] = constraintExtFuncSecondGradFirst()(xs, r)
	hess.Points[1][1] = constraintExtFuncSecondGradSecond()(xs, r)
	return hess
}

func hessConstraintExt() func(xs []float64, r float64) la_methods.Matrix {
	return func(xs []float64, r float64) la_methods.Matrix {
		return hessianConstraintExt(xs, r)
	}
}
