package constraint_methods

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"github.com/saskamegaprogrammist/optimization_methods/many_dimension_search"
	"math"
)

type PenaltyCombined struct {
	startPoint            []float64
	dimension             int
	targetFunc            func(xs []float64) float64
	targetFuncConstraint  func(xs []float64, r float64) float64
	gradient              []func(xs []float64) float64
	gradientConstraintExt []func(xs []float64, r float64) float64
	gradientConstraintInt []func(xs []float64, r float64) float64
	constraintExt         func(xs []float64, r float64) float64
	constraintInt         func(xs []float64, r float64) float64
	c1                    float64
	c2                    float64
	eps                   float64
	method                string
	methodMap             map[string]func(x []float64, r float64, targetFunc func(xs []float64) float64,
		gradient []func(xs []float64) float64) ([]float64, float64, error)
}

func (pc *PenaltyCombined) Init(startPoint []float64, dimension int,
	targetFunc func(xs []float64) float64,
	gradient []func(xs []float64) float64,
	gradientConstraintExt []func(xs []float64, r float64) float64,
	gradientConstraintInt []func(xs []float64, r float64) float64,
	constraintExt func(xs []float64, r float64) float64,
	constraintInt func(xs []float64, r float64) float64,
	eps float64, c1 float64, c2 float64, method string) {
	pc.startPoint = startPoint
	pc.targetFunc = targetFunc
	pc.dimension = dimension
	pc.gradient = gradient
	pc.gradientConstraintExt = gradientConstraintExt
	pc.gradientConstraintInt = gradientConstraintInt
	pc.constraintExt = constraintExt
	pc.constraintInt = constraintInt
	pc.c1 = c1
	pc.c2 = c2
	pc.eps = eps
	pc.method = method
	pc.methodMap = map[string]func(x []float64, r float64, targetFunc func(xs []float64) float64,
		gradient []func(xs []float64) float64) ([]float64, float64, error){
		"hooke jeeves":            pc.hookeJeevesSearch,
		"fast gradient":           pc.fastGradientDescendSearch,
		"nelder mead":             pc.nelderMeadSearch,
		"fletcher reeves":         pc.fletcherReevesSearch,
		"pollac":                  pc.pollacSearch,
		"davidon fletcher powell": pc.davidonFletcherPowell,
	}
}

func (pc *PenaltyCombined) Solve() ([]float64, float64, error) {
	var err error
	var x la_methods.Vector
	var k int
	var r1, r2 float64
	var xMin []float64
	var yMin, valConstraintMin float64
	var yMinOld, valConstrOld float64
	var constrFunc func(xs []float64) float64
	var targFunc func(xs []float64) float64
	var constrFuncGrad []func(xs []float64) float64
	var targFuncGrad []func(xs []float64) float64
	r1 = 1
	r2 = 1
	err = x.InitWithPoints(pc.dimension, pc.startPoint)
	if err != nil {
		return nil, 0, fmt.Errorf("error during vector initializing: %v", err)
	}
	var constraints []func(xs []float64, r float64) float64
	constraints = append(constraints, pc.constraintExt, pc.constraintInt)
	var constraintsGrads [][]func(xs []float64, r float64) float64
	constraintsGrads = append(constraintsGrads, pc.gradientConstraintExt, pc.gradientConstraintInt)
	constrFunc = pc.addConstraints(constraints, r1)
	valConstrOld = constrFunc(x.Points)
	for {
		constrFunc = pc.addConstraints(constraints, r1)
		constrFuncGrad = pc.addGradientsConstraints(constraintsGrads, r1)

		xMin, valConstraintMin, err = pc.methodMap[pc.method](x.Points, r1, constrFunc, constrFuncGrad)
		//fmt.Println(xMin, valConstraintMin)
		if math.Abs(valConstraintMin-valConstrOld) < pc.eps {
			break
		} else {
			valConstrOld = valConstraintMin
			k++
			r1 *= pc.c1
			err = x.InitWithPoints(pc.dimension, xMin)
			if err != nil {
				return nil, 0, fmt.Errorf("error during vector initializing: %v", err)
			}
		}
	}

	targFunc = pc.addFunctions(pc.targetFunc, constraints, r2)
	yMinOld = targFunc(x.Points)
	for {
		targFunc = pc.addFunctions(pc.targetFunc, constraints, r2)
		targFuncGrad = pc.addGradients(pc.gradient, constraintsGrads, r2)

		xMin, yMin, err = pc.methodMap[pc.method](x.Points, r2, targFunc, targFuncGrad)
		//fmt.Println(xMin, yMin)
		if math.Abs(yMin-yMinOld) < pc.eps {
			fmt.Printf("k value: %d\n", k)
			return xMin, yMin, nil
		} else {
			yMinOld = yMin
			k++
			r2 *= pc.c2
			err = x.InitWithPoints(pc.dimension, xMin)
			if err != nil {
				return nil, 0, fmt.Errorf("error during vector initializing: %v", err)
			}
		}
	}
}

func (pc *PenaltyCombined) hookeJeevesSearch(x []float64, r float64, targetFunc func(xs []float64) float64,
	gradient []func(xs []float64) float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var hjs many_dimension_search.HookeJeevesSearch
	hjs.Init(x, 0.1, pc.dimension, 2, 0.0001, 0.1,
		0.1, targetFunc, "break in two")
	xMin, yMin, err = hjs.Solve()
	if err != nil {
		return nil, 0, fmt.Errorf("error solving hooke jeeves : %v\n", err)
	}
	return xMin, yMin, nil
}

func (pc *PenaltyCombined) nelderMeadSearch(x []float64, r float64, targetFunc func(xs []float64) float64,
	gradient []func(xs []float64) float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var nms many_dimension_search.NelderMeadSearch
	nms.Init(x, 0.1, pc.dimension, pc.eps, targetFunc)
	xMin, yMin, err = nms.Solve()
	if err != nil {
		return nil, 0, fmt.Errorf("error solving nelder mead : %v\n", err)
	}
	return xMin, yMin, nil
}

func (pc *PenaltyCombined) fastGradientDescendSearch(x []float64, r float64, targetFunc func(xs []float64) float64,
	gradient []func(xs []float64) float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var fgd many_dimension_search.FastGradientDescendSearch
	fgd.Init(x, pc.eps, pc.eps, targetFunc, gradient, pc.dimension, pc.eps, pc.eps, "golden ratio")
	xMin, yMin, err = fgd.Solve()
	if err != nil {
		return nil, 0, fmt.Errorf("error solving fast gradient descent method : %v\n", err)
	}
	return xMin, yMin, nil
}

func (pc *PenaltyCombined) pollacSearch(x []float64, r float64, targetFunc func(xs []float64) float64,
	gradient []func(xs []float64) float64) ([]float64, float64, error) {
	return pc.fletcherReevesSearchWithParam(x, r, true, targetFunc, gradient)
}

func (pc *PenaltyCombined) fletcherReevesSearch(x []float64, r float64, targetFunc func(xs []float64) float64,
	gradient []func(xs []float64) float64) ([]float64, float64, error) {
	return pc.fletcherReevesSearchWithParam(x, r, false, targetFunc, gradient)
}

func (pc *PenaltyCombined) fletcherReevesSearchWithParam(x []float64, r float64, pollac bool, targetFunc func(xs []float64) float64,
	gradient []func(xs []float64) float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var frs many_dimension_search.FletcherReevesSearch
	frs.Init(x, 0.001, 2, pc.eps, pc.eps, pc.eps, 0.00011, 100,
		targetFunc, gradient, "break in two", pollac)
	xMin, yMin, err = frs.Solve()
	if err != nil {
		fmt.Printf("error solving pollak : %v\n", err)
		return nil, 0, fmt.Errorf("error solving pollak : %v\n", err)
	}
	return xMin, yMin, nil
}

func (pc *PenaltyCombined) davidonFletcherPowell(x []float64, r float64, targetFunc func(xs []float64) float64,
	gradient []func(xs []float64) float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var dfps many_dimension_search.DavidonFletcherPowellSearch
	dfps.Init(x, pc.eps, pc.dimension, pc.eps, pc.eps, pc.eps, 0.00011, 100, targetFunc,
		gradient, "golden ratio")
	xMin, yMin, err = dfps.Solve()
	if err != nil {
		fmt.Printf("error solving davidon fletcher powell : %v\n", err)
		return nil, 0, fmt.Errorf("error solving davidon fletcher powell : %v\n", err)
	}
	return xMin, yMin, nil
}

func (pc *PenaltyCombined) addFunctions(function func(xs []float64) float64,
	constraintFunctions []func(xs []float64, r float64) float64, r float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		var sum float64
		for _, f := range constraintFunctions {
			sum += f(xs, r)
		}
		return function(xs) + sum
	}
}

func (pc *PenaltyCombined) addGradients(gradient []func(xs []float64) float64,
	gradientConstraints [][]func(xs []float64, r float64) float64, r float64) []func(xs []float64) float64 {
	var newGrad []func(xs []float64) float64
	newGrad = append(newGrad, func(xs []float64) float64 {
		var sum float64
		for _, f := range gradientConstraints {
			sum += f[0](xs, r)
		}
		return gradient[0](xs) + sum
	}, func(xs []float64) float64 {
		var sum float64
		for _, f := range gradientConstraints {
			sum += f[1](xs, r)
		}
		return gradient[1](xs) + sum
	})
	return newGrad
}

func (pc *PenaltyCombined) addConstraints(constraintFunctions []func(xs []float64, r float64) float64, r float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		var sum float64
		for _, f := range constraintFunctions {
			sum += f(xs, r)
		}
		return sum
	}
}

func (pc *PenaltyCombined) addGradientsConstraints(gradientConstraints [][]func(xs []float64, r float64) float64, r float64) []func(xs []float64) float64 {
	var newGrad []func(xs []float64) float64
	newGrad = append(newGrad, func(xs []float64) float64 {
		var sum float64
		for _, f := range gradientConstraints {
			sum += f[0](xs, r)
		}
		return sum
	}, func(xs []float64) float64 {
		var sum float64
		for _, f := range gradientConstraints {
			sum += f[1](xs, r)
		}
		return sum
	})
	return newGrad
}
