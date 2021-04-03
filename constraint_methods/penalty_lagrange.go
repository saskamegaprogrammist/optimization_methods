package constraint_methods

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"github.com/saskamegaprogrammist/optimization_methods/many_dimension_search"
	"math"
)

type PenaltyLagrange struct {
	startPoint           []float64
	dimension            int
	targetFunc           func(xs []float64) float64
	targetFuncConstraint func(xs []float64, r float64) float64
	gradient             []func(xs []float64) float64
	constraintFunctions  []func(xs []float64, r float64, m []float64) float64
	gradientConstraint   []func(xs []float64, r float64, m []float64) float64
	constraint           func(xs []float64, r float64, m []float64) float64
	c                    float64
	m                    []float64
	eps                  float64
	method               string
	methodMap            map[string]func(x []float64, r float64, m []float64) ([]float64, float64, error)
}

func (pl *PenaltyLagrange) Init(startPoint []float64, dimension int,
	targetFunc func(xs []float64) float64,
	gradient []func(xs []float64) float64,
	constraintFunctions []func(xs []float64, r float64, m []float64) float64,
	gradientConstraint []func(xs []float64, r float64, m []float64) float64,
	constraint func(xs []float64, r float64, m []float64) float64,
	m []float64,
	eps float64, c float64, method string) {
	pl.startPoint = startPoint
	pl.targetFunc = targetFunc
	pl.dimension = dimension
	pl.gradient = gradient
	pl.gradientConstraint = gradientConstraint
	pl.constraint = constraint
	pl.constraintFunctions = constraintFunctions
	pl.c = c
	pl.m = m
	pl.eps = eps
	pl.method = method
	pl.methodMap = map[string]func(x []float64, r float64, m []float64) ([]float64, float64, error){
		"hooke jeeves":            pl.hookeJeevesSearch,
		"fast gradient":           pl.fastGradientDescendSearch,
		"nelder mead":             pl.nelderMeadSearch,
		"fletcher reeves":         pl.fletcherReevesSearch,
		"pollac":                  pl.pollacSearch,
		"davidon fletcher powell": pl.davidonFletcherPowell,
	}
}

func (pl *PenaltyLagrange) Solve() ([]float64, float64, error) {
	var err error
	var x la_methods.Vector
	var k int
	var r float64
	var xMin []float64
	var yMin float64
	var m []float64
	r = 4
	err = x.InitWithPoints(pl.dimension, pl.startPoint)
	if err != nil {
		return nil, 0, fmt.Errorf("error during vector initializing: %v", err)
	}
	m = pl.m
	for {
		//fmt.Println(m)
		xMin, yMin, err = pl.methodMap[pl.method](x.Points, r, m)
		//fmt.Println(xMin, yMin, pl.constraint(xMin, r, m))
		if math.Abs(pl.constraint(xMin, r, m)) < pl.eps {
			fmt.Printf("k value: %d\n", k)
			return xMin, yMin, nil
		} else {
			k++
			//r *= pl.c
			err = x.InitWithPoints(pl.dimension, xMin)
			if err != nil {
				return nil, 0, fmt.Errorf("error during vector initializing: %v", err)
			}
			m = pl.calculateM(m, x.Points, r)
		}
	}
}

func (pl *PenaltyLagrange) calculateM(m []float64, x []float64, r float64) []float64 {
	var newM []float64
	for _, f := range pl.constraintFunctions {
		newM = append(newM, f(x, r, m))
	}
	return newM
}

func (pl *PenaltyLagrange) hookeJeevesSearch(x []float64, r float64, m []float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var hjs many_dimension_search.HookeJeevesSearch
	tf := pl.addFunctions(pl.targetFunc, pl.constraint, r, m)
	hjs.Init(x, 0.1, pl.dimension, 2, 0.0001, 0.1,
		0.1, tf, "golden ratio")
	xMin, yMin, err = hjs.Solve()
	if err != nil {
		return nil, 0, fmt.Errorf("error solving hooke jeeves : %v\n", err)
	}
	return xMin, yMin, nil
}

func (pl *PenaltyLagrange) nelderMeadSearch(x []float64, r float64, m []float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var nms many_dimension_search.NelderMeadSearch
	nms.Init(x, 0.1, pl.dimension, pl.eps, pl.addFunctions(pl.targetFunc, pl.constraint, r, m))
	xMin, yMin, err = nms.Solve()
	if err != nil {
		return nil, 0, fmt.Errorf("error solving nelder mead : %v\n", err)
	}
	return xMin, yMin, nil
}

func (pl *PenaltyLagrange) fastGradientDescendSearch(x []float64, r float64, m []float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var fgd many_dimension_search.FastGradientDescendSearch
	fgd.Init(x, pl.eps, pl.eps, pl.addFunctions(pl.targetFunc, pl.constraint, r, m), pl.addGradients(pl.gradient, pl.gradientConstraint, r, m), pl.dimension, pl.eps, pl.eps, "golden ratio")
	xMin, yMin, err = fgd.Solve()
	if err != nil {
		return nil, 0, fmt.Errorf("error solving fast gradient descent method : %v\n", err)
	}
	return xMin, yMin, nil
}

func (pl *PenaltyLagrange) pollacSearch(x []float64, r float64, m []float64) ([]float64, float64, error) {
	return pl.fletcherReevesSearchWithParam(x, r, m, true)
}

func (pl *PenaltyLagrange) fletcherReevesSearch(x []float64, r float64, m []float64) ([]float64, float64, error) {
	return pl.fletcherReevesSearchWithParam(x, r, m, false)
}

func (pl *PenaltyLagrange) fletcherReevesSearchWithParam(x []float64, r float64, m []float64, pollac bool) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var frs many_dimension_search.FletcherReevesSearch
	frs.Init(x, 0.0001, 2, pl.eps, pl.eps, 0.00001, 0.0001, 10,
		pl.addFunctions(pl.targetFunc, pl.constraint, r, m),
		pl.addGradients(pl.gradient, pl.gradientConstraint, r, m), "golden ratio", pollac)
	xMin, yMin, err = frs.Solve()
	if err != nil {
		fmt.Printf("error solving pollak : %v\n", err)
		return nil, 0, fmt.Errorf("error solving pollak : %v\n", err)
	}
	return xMin, yMin, nil
}

func (pl *PenaltyLagrange) davidonFletcherPowell(x []float64, r float64, m []float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var dfps many_dimension_search.DavidonFletcherPowellSearch
	dfps.Init(x, 0.0001, pl.dimension, pl.eps, pl.eps, 0.00001, 0.0001, 10, pl.addFunctions(pl.targetFunc, pl.constraint, r, m),
		pl.addGradients(pl.gradient, pl.gradientConstraint, r, m), "break in two")
	xMin, yMin, err = dfps.Solve()
	if err != nil {
		fmt.Printf("error solving davidon fletcher powell : %v\n", err)
		return nil, 0, fmt.Errorf("error solving davidon fletcher powell : %v\n", err)
	}
	return xMin, yMin, nil
}

func (pl *PenaltyLagrange) addFunctions(function func(xs []float64) float64,
	constraintFunction func(xs []float64, r float64, m []float64) float64, r float64, m []float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		return function(xs) + constraintFunction(xs, r, m)
	}

}

func (pl *PenaltyLagrange) addGradients(gradient []func(xs []float64) float64,
	gradientConstraint []func(xs []float64, r float64, m []float64) float64, r float64, m []float64) []func(xs []float64) float64 {
	var newGrad []func(xs []float64) float64
	newGrad = append(newGrad, func(xs []float64) float64 {
		return gradient[0](xs) + gradientConstraint[0](xs, r, m)
	}, func(xs []float64) float64 {
		return gradient[1](xs) + gradientConstraint[1](xs, r, m)
	})
	return newGrad
}
