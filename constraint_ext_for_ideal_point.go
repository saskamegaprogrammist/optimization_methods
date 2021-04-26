package main

import "math"

func firstConstraintForIdeal(xs []float64) float64 {
	return 2*xs[0] + xs[1] - 12
}

func firstConstraintForIdealMod(xs []float64) float64 {
	if firstConstraintForIdeal(xs) > 0 {
		return firstConstraintForIdeal(xs)
	}
	return 0
}

func secondConstraintForIdeal(xs []float64) float64 {
	return -2*xs[0] + xs[1] - 10
}

func secondConstraintForIdealMod(xs []float64) float64 {
	if secondConstraintForIdeal(xs) > 0 {
		return secondConstraintForIdeal(xs)
	}
	return 0
}

func firstConstraintForIdealGradFirst(xs []float64) float64 {
	return 2 * 2 * (2*xs[0] + xs[1])
}

func firstConstraintForIdealGradFirstMod(xs []float64) float64 {
	if firstConstraintForIdeal(xs) > 0 {
		return firstConstraintForIdealGradFirst(xs)
	}
	return 0
}

func firstConstraintForIdealGradSecond(xs []float64) float64 {
	return 2 * (2*xs[0] + xs[1])
}

func firstConstraintForIdealGradSecondMod(xs []float64) float64 {
	if firstConstraintForIdeal(xs) > 0 {
		return firstConstraintForIdealGradSecond(xs)
	}
	return 0
}

func secondConstraintForIdealGradFirst(xs []float64) float64 {
	return -2 * 2 * (-2*xs[0] + xs[1])
}

func secondConstraintForIdealGradFirstMod(xs []float64) float64 {
	if secondConstraintForIdeal(xs) > 0 {
		return secondConstraintForIdealGradFirst(xs)
	}
	return 0
}

func secondConstraintForIdealGradSecond(xs []float64) float64 {
	return 2 * (-2*xs[0] + xs[1])
}

func secondConstraintForIdealGradSecondMod(xs []float64) float64 {
	if secondConstraintForIdeal(xs) > 0 {
		return secondConstraintForIdealGradSecond(xs)
	}
	return 0
}

func constraintExtForIdealFunc(constraint1 func(xs []float64) float64, constraint2 func(xs []float64) float64) func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return r / float64(2) * (math.Pow(constraint1(xs), 2) + math.Pow(constraint2(xs), 2))
	}
}

func constraintExtFuncForIdealFirstGrad() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return r / float64(2) * (firstConstraintForIdealGradFirstMod(xs) + secondConstraintForIdealGradFirstMod(xs))
	}
}

func constraintExtFuncForIdealSecondGrad() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return r / float64(2) * (firstConstraintForIdealGradSecondMod(xs) + secondConstraintForIdealGradSecondMod(xs))
	}
}

func constraintExtFuncForIdealThirdGrad() func(xs []float64, r float64) float64 {
	return func(xs []float64, r float64) float64 {
		return 0
	}
}
