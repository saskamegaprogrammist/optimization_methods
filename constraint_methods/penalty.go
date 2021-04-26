package constraint_methods

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/genetic_methods"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"github.com/saskamegaprogrammist/optimization_methods/many_dimension_search"
	"math"
)

type Penalty struct {
	startPoint           []float64
	dimension            int
	targetFunc           func(xs []float64) float64
	targetFuncConstraint func(xs []float64, r float64) float64
	penalties            []func(xs []float64) float64
	gradient             []func(xs []float64) float64
	hessian              func(xs []float64) la_methods.Matrix
	gradientConstraint   []func(xs []float64, r float64) float64
	hessianConstraint    func(xs []float64, r float64) la_methods.Matrix
	constraint           func(xs []float64, r float64) float64
	c                    float64
	eps                  float64
	method               string
	methodMap            map[string]func(x []float64, r float64) ([]float64, float64, error)
}

func (ep *Penalty) Init(startPoint []float64, dimension int,
	targetFunc func(xs []float64) float64, penalties []func(xs []float64) float64,
	gradient []func(xs []float64) float64, hessian func(xs []float64) la_methods.Matrix,
	gradientConstraint []func(xs []float64, r float64) float64,
	hessianConstraint func(xs []float64, r float64) la_methods.Matrix,
	constraint func(xs []float64, r float64) float64,
	eps float64, c float64, method string) {
	ep.startPoint = startPoint
	ep.targetFunc = targetFunc
	ep.dimension = dimension
	ep.penalties = penalties
	ep.gradient = gradient
	ep.hessian = hessian
	ep.gradientConstraint = gradientConstraint
	ep.hessianConstraint = hessianConstraint
	ep.constraint = constraint
	ep.c = c
	ep.eps = eps
	ep.method = method
	ep.methodMap = map[string]func(x []float64, r float64) ([]float64, float64, error){
		"hooke jeeves":            ep.hookeJeevesSearch,
		"fast gradient":           ep.fastGradientDescendSearch,
		"nelder mead":             ep.nelderMeadSearch,
		"fletcher reeves":         ep.fletcherReevesSearch,
		"pollac":                  ep.pollacSearch,
		"davidon fletcher powell": ep.davidonFletcherPowell,
		"levenberg":               ep.levenbergMarkkvadratSearch,
		"genetic":                 ep.geneticAlgorithm,
	}
}

func (ep *Penalty) InitSimple(startPoint []float64, dimension int,
	targetFunc func(xs []float64) float64, penalties []func(xs []float64) float64,
	gradient []func(xs []float64) float64,
	gradientConstraint []func(xs []float64, r float64) float64,
	constraint func(xs []float64, r float64) float64,
	eps float64, c float64, method string) {
	ep.startPoint = startPoint
	ep.targetFunc = targetFunc
	ep.dimension = dimension
	ep.penalties = penalties
	ep.gradient = gradient
	ep.gradientConstraint = gradientConstraint
	ep.constraint = constraint
	ep.c = c
	ep.eps = eps
	ep.method = method
	ep.methodMap = map[string]func(x []float64, r float64) ([]float64, float64, error){
		"hooke jeeves":            ep.hookeJeevesSearch,
		"fast gradient":           ep.fastGradientDescendSearch,
		"nelder mead":             ep.nelderMeadSearch,
		"fletcher reeves":         ep.fletcherReevesSearch,
		"pollac":                  ep.pollacSearch,
		"davidon fletcher powell": ep.davidonFletcherPowell,
		"genetic":                 ep.geneticAlgorithm,
	}
}

func (ep *Penalty) Solve() ([]float64, float64, error) {
	var err error
	var x la_methods.Vector
	var k int
	var r float64
	var xMin []float64
	var yMin float64
	r = 1
	err = x.InitWithPoints(ep.dimension, ep.startPoint)
	if err != nil {
		return nil, 0, fmt.Errorf("error during vector initializing: %v", err)
	}
	for {
		xMin, yMin, err = ep.methodMap[ep.method](x.Points, r)
		//fmt.Println(xMin, yMin, ep.constraint(xMin, r))
		if math.Abs(ep.constraint(xMin, r)) < ep.eps {
			//fmt.Printf("k value: %d\n", k)
			return xMin, yMin, nil
		} else {
			k++
			r *= ep.c
			err = x.InitWithPoints(ep.dimension, xMin)
			if err != nil {
				return nil, 0, fmt.Errorf("error during vector initializing: %v", err)
			}
		}
	}
}

func (ep *Penalty) hookeJeevesSearch(x []float64, r float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var hjs many_dimension_search.HookeJeevesSearch
	hjs.Init(x, 0.1, ep.dimension, 2, 0.0001, 0.1,
		0.1, ep.addFunctions(ep.targetFunc, ep.constraint, r), "fibonacci")
	xMin, yMin, err = hjs.Solve()
	if err != nil {
		return nil, 0, fmt.Errorf("error solving hooke jeeves : %v\n", err)
	}
	return xMin, yMin, nil
}

func (ep *Penalty) nelderMeadSearch(x []float64, r float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var nms many_dimension_search.NelderMeadSearch
	nms.Init(x, 0.1, ep.dimension, ep.eps, ep.addFunctions(ep.targetFunc, ep.constraint, r))
	xMin, yMin, err = nms.Solve()
	if err != nil {
		return nil, 0, fmt.Errorf("error solving nelder mead : %v\n", err)
	}
	return xMin, yMin, nil
}

func (ep *Penalty) geneticAlgorithm(x []float64, r float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var ga genetic_methods.GeneticAlgorithm
	ga.Init(0, 4, 1, 2000, x, ep.dimension, ep.addFunctions(ep.targetFunc, ep.constraint, r), func(xs []float64) float64 {
		return float64(1) / ep.addFunctions(ep.targetFunc, ep.constraint, r)(xs)
	})
	xMin, yMin, err = ga.Solve()
	if err != nil {
		return nil, 0, fmt.Errorf("error solving genetic algorithm : %v\n", err)
	}
	return xMin, yMin, nil
}

func (ep *Penalty) fastGradientDescendSearch(x []float64, r float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var fgd many_dimension_search.FastGradientDescendSearch
	fgd.Init(x, ep.eps, ep.eps, ep.addFunctions(ep.targetFunc, ep.constraint, r), ep.addGradients(ep.gradient, ep.gradientConstraint, r), ep.dimension, ep.eps, ep.eps, "fibonacci")
	xMin, yMin, err = fgd.Solve()
	if err != nil {
		return nil, 0, fmt.Errorf("error solving fast gradient descent method : %v\n", err)
	}
	return xMin, yMin, nil
}

func (ep *Penalty) pollacSearch(x []float64, r float64) ([]float64, float64, error) {
	return ep.fletcherReevesSearchWithParam(x, r, true)
}

func (ep *Penalty) fletcherReevesSearch(x []float64, r float64) ([]float64, float64, error) {
	return ep.fletcherReevesSearchWithParam(x, r, false)
}

func (ep *Penalty) fletcherReevesSearchWithParam(x []float64, r float64, pollac bool) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var frs many_dimension_search.FletcherReevesSearch
	frs.Init(x, 0.001, 3, ep.eps, ep.eps, ep.eps, 0.00011, 100,
		ep.addFunctions(ep.targetFunc, ep.constraint, r),
		ep.addGradients(ep.gradient, ep.gradientConstraint, r), "golden ratio", pollac)
	xMin, yMin, err = frs.Solve()
	if err != nil {
		fmt.Printf("error solving pollak : %v\n", err)
		return nil, 0, fmt.Errorf("error solving pollak : %v\n", err)
	}
	return xMin, yMin, nil
}

func (ep *Penalty) davidonFletcherPowell(x []float64, r float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var dfps many_dimension_search.DavidonFletcherPowellSearch
	dfps.Init(x, ep.eps, ep.dimension, ep.eps, ep.eps, ep.eps, 0.00011, 100, ep.addFunctions(ep.targetFunc, ep.constraint, r),
		ep.addGradients(ep.gradient, ep.gradientConstraint, r), "fibonacci")
	xMin, yMin, err = dfps.Solve()
	if err != nil {
		fmt.Printf("error solving davidon fletcher powell : %v\n", err)
		return nil, 0, fmt.Errorf("error solving davidon fletcher powell : %v\n", err)
	}
	return xMin, yMin, nil
}

func (ep *Penalty) levenbergMarkkvadratSearch(x []float64, r float64) ([]float64, float64, error) {
	var xMin []float64
	var yMin float64
	var err error
	var lms many_dimension_search.LevenbergMarkkvadratSearch
	lms.Init(x, ep.dimension, ep.addFunctions(ep.targetFunc, ep.constraint, r),
		ep.addGradients(ep.gradient, ep.gradientConstraint, r),
		ep.addHessians(ep.hessian, ep.hessianConstraint, r), 1000, 10, 0.00001)
	xMin, yMin, err = lms.Solve()
	if err != nil {
		return nil, 0, fmt.Errorf("error solving levenberg markkvadrat method : %v\n", err)
	}
	return xMin, yMin, nil
}

func (ep *Penalty) addFunctions(function func(xs []float64) float64,
	constraintFunction func(xs []float64, r float64) float64, r float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		return function(xs) + constraintFunction(xs, r)
	}

}

func (ep *Penalty) addHessians(funcHessian func(xs []float64) la_methods.Matrix,
	constraintHessian func(xs []float64, r float64) la_methods.Matrix, r float64) func(xs []float64) la_methods.Matrix {
	return func(xs []float64) la_methods.Matrix {
		matrixFunc := funcHessian(xs)
		newM, _ := matrixFunc.AddM(constraintHessian(xs, r))
		return newM
	}
}

func (ep *Penalty) addGradients(gradient []func(xs []float64) float64,
	gradientConstraint []func(xs []float64, r float64) float64, r float64) []func(xs []float64) float64 {
	var newGrad []func(xs []float64) float64
	newGrad = append(newGrad, func(xs []float64) float64 {
		return gradient[0](xs) + gradientConstraint[0](xs, r)
	}, func(xs []float64) float64 {
		return gradient[1](xs) + gradientConstraint[1](xs, r)
	}, func(xs []float64) float64 {
		return gradient[2](xs) + gradientConstraint[2](xs, r)
	})
	return newGrad
}
