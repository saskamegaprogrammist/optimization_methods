package main

import "math"

func constraintLagrangeFunc(constraint1 func(xs []float64, r float64, m []float64) float64, constraint2 func(xs []float64, r float64, m []float64) float64, constraint3 func(xs []float64, r float64, m []float64) float64) func(xs []float64, r float64, m []float64) float64 {
	return func(xs []float64, r float64, m []float64) float64 {
		return float64(1) / (r * float64(2)) * (math.Pow(constraint1(xs, r, m), 2-math.Pow(m[0], 2)) + math.Pow(constraint2(xs, r, m), 2-math.Pow(m[1], 2)) + math.Pow(constraint3(xs, r, m), 2) - math.Pow(m[2], 2))
	}
}

func constraintLagrangeFuncFirstGrad() func(xs []float64, r float64, m []float64) float64 {
	return func(xs []float64, r float64, m []float64) float64 {
		return float64(1) / (r * float64(2)) * (firstLagrangeConstraintFirstGradMod(xs, r, m) + secondLagrangeConstraintFirstGradMod(xs, r, m))
	}
}

func constraintLagrangeFuncSecondGrad() func(xs []float64, r float64, m []float64) float64 {
	return func(xs []float64, r float64, m []float64) float64 {
		return float64(1) / (r * float64(2)) * (firstLagrangeConstraintSecondGradMod(xs, r, m) + thirdLagrangeConstraintSecondGradMod(xs, r, m))
	}
}

func firstLagrangeConstraintMod(xs []float64, r float64, m []float64) float64 {
	val := m[0] + r*firstConstraint(xs)
	if val > 0 {
		return val
	}
	return 0
}

func secondLagrangeConstraintMod(xs []float64, r float64, m []float64) float64 {
	val := m[1] + r*secondConstraint(xs)
	if val > 0 {
		return val
	}
	return 0
}

func thirdLagrangeConstraintMod(xs []float64, r float64, m []float64) float64 {
	val := m[2] + r*thirdConstraint(xs)
	if val > 0 {
		return val
	}
	return 0
}

func firstLagrangeConstraintFirstGrad(xs []float64, r float64) float64 {
	return 2 * r * firstConstraint(xs) * 2 * xs[0]
}

func firstLagrangeConstraintFirstGradMod(xs []float64, r float64, m []float64) float64 {
	val := m[0] + r*firstConstraint(xs)
	if val > 0 {
		return firstLagrangeConstraintFirstGrad(xs, r)
	}
	return 0
}

func firstLagrangeConstraintSecondGrad(xs []float64, r float64) float64 {
	return 2 * r * firstConstraint(xs) * 2 * xs[0]
}

func firstLagrangeConstraintSecondGradMod(xs []float64, r float64, m []float64) float64 {
	val := m[0] + r*firstConstraint(xs)
	if val > 0 {
		return firstLagrangeConstraintSecondGrad(xs, r)
	}
	return 0
}

func secondLagrangeConstraintFirstGrad(xs []float64, r float64) float64 {
	return 2 * xs[0] * math.Pow(r, 2)
}

func secondLagrangeConstraintFirstGradMod(xs []float64, r float64, m []float64) float64 {
	val := m[1] + r*secondConstraint(xs)
	if val > 0 {
		return secondLagrangeConstraintFirstGrad(xs, r)
	}
	return 0
}

func thirdLagrangeConstraintSecondGrad(xs []float64, r float64) float64 {
	return 2 * xs[1] * math.Pow(r, 2)
}

func thirdLagrangeConstraintSecondGradMod(xs []float64, r float64, m []float64) float64 {
	val := m[2] + r*thirdConstraint(xs)
	if val > 0 {
		return thirdLagrangeConstraintSecondGrad(xs, r)
	}
	return 0
}
