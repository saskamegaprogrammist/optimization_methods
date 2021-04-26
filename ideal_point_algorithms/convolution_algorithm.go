package ideal_point_algorithms

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/constraint_methods"
	"github.com/saskamegaprogrammist/optimization_methods/genetic_methods"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"github.com/saskamegaprogrammist/optimization_methods/many_dimension_search"
	"gonum.org/v1/gonum/mat"
)

type ConvolutionMulticriteria struct {
	startPoint         []float64
	dimension          int
	critPriority       [][]float64
	targetFuncs        []func(xs []float64) float64                                // f1 and f2
	penalties          []func(xs []float64) float64                                // g1 and g2
	helpFunc           func(xs []float64, ideal []float64, ws []float64) float64   // f1 and f2 combined function
	gradient           []func(xs []float64, ideal []float64, ws []float64) float64 // f1 and f2 combined function gradient
	gradientConstraint []func(xs []float64, r float64) float64                     // g1 and g2 combined function gradient
	constraint         func(xs []float64, r float64) float64                       // g1 and g2 combined function

	penaltyAlg constraint_methods.Penalty
	search     many_dimension_search.NelderMeadSearch
	genetic    genetic_methods.GeneticAlgorithm

	useGenetic bool
}

func (cm *ConvolutionMulticriteria) Init(dimension int, startPoint []float64, critPriority [][]float64, targetFuncs []func(xs []float64) float64, penalties []func(xs []float64) float64,
	helpFunc func(xs []float64, ideal []float64, ws []float64) float64, gradient []func(xs []float64, ideal []float64, ws []float64) float64, gradientConstraint []func(xs []float64, r float64) float64, constraint func(xs []float64, r float64) float64, useGenetic bool) {
	cm.startPoint = startPoint
	cm.targetFuncs = targetFuncs
	cm.dimension = dimension
	cm.penalties = penalties
	cm.gradient = gradient
	cm.helpFunc = helpFunc
	cm.critPriority = critPriority
	cm.gradientConstraint = gradientConstraint
	cm.constraint = constraint
	cm.useGenetic = useGenetic
}

func (cm *ConvolutionMulticriteria) Solve() ([][]float64, [][]float64, error) {
	var err error

	var xMin []float64
	var xMins, yMin [][]float64
	var fIdeal = make([]float64, 2)

	var xIdeal1, xIdeal2 []float64
	if cm.useGenetic {
		cm.genetic.Init(0, 1, 1, 2000, cm.startPoint, cm.dimension, cm.targetFuncs[0], func(xs []float64) float64 {
			return float64(1) / cm.targetFuncs[0](xs)
		})
		xIdeal1, fIdeal[0], err = cm.genetic.Solve()
	} else {
		cm.search.Init(cm.startPoint, 0.1, 3, 0.001, cm.targetFuncs[0])
		xIdeal1, fIdeal[0], err = cm.search.Solve()
	}

	if err != nil {
		return nil, nil, fmt.Errorf("error finding first ideal point: %v", err)
	}

	fmt.Println("xIdeal1 ", xIdeal1)

	if cm.useGenetic {
		cm.genetic.Init(0, 4, 1, 3000, cm.startPoint, cm.dimension, cm.targetFuncs[1], func(xs []float64) float64 {
			return float64(1) / cm.targetFuncs[1](xs)
		})
		xIdeal2, fIdeal[1], err = cm.genetic.Solve()
	} else {
		cm.search.Init(cm.startPoint, 0.1, 3, 0.001, cm.targetFuncs[1])
		xIdeal2, fIdeal[1], err = cm.search.Solve()
	}

	if err != nil {
		return nil, nil, fmt.Errorf("error finding second ideal point: %v", err)
	}

	fmt.Println("xIdeal2 ", xIdeal2)

	fmt.Println("ideal values: ", fIdeal)

	var frontLen = len(cm.critPriority)
	for i := 0; i < frontLen; i++ {
		a1 := cm.critPriority[i][0] / cm.critPriority[i][1]
		a2 := float64(1) / a1

		aDense := mat.NewDense(2, 2, []float64{1, a1, a2, 1})
		var eigen mat.Eigen
		eigen.Factorize(aDense, mat.EigenBoth)
		var vectors *mat.CDense = mat.NewCDense(2, 2, make([]complex128, 4))
		var values []complex128 = make([]complex128, 2)
		eigen.VectorsTo(vectors)
		eigen.Values(values)
		var sum float64
		var vec la_methods.Vector
		if real(values[0]) > real(values[1]) {
			err = vec.InitWithPoints(2, []float64{real(vectors.At(0, 0)), real(vectors.At(1, 0))})
			if err != nil {
				return nil, nil, fmt.Errorf("error initing vector: %v", err)
			}
		} else {
			err = vec.InitWithPoints(2, []float64{real(vectors.At(0, 1)), real(vectors.At(1, 1))})
			if err != nil {
				return nil, nil, fmt.Errorf("error initing vector: %v", err)
			}
		}
		for i := 0; i < 2; i++ {
			sum += vec.Points[i]
		}

		var ws = make([]float64, 2)
		ws[0] = vec.Points[0] / sum
		ws[1] = vec.Points[1] / sum

		//fmt.Println(fIdeal)
		fmt.Println(ws)

		var gradient = []func(xs []float64) float64{
			func(xs []float64) float64 { return cm.gradient[0](xs, fIdeal, ws) },
			func(xs []float64) float64 { return cm.gradient[1](xs, fIdeal, ws) },
			func(xs []float64) float64 { return cm.gradient[2](xs, fIdeal, ws) },
		}

		var ep constraint_methods.Penalty
		var mthd string = "nelder mead"
		if cm.useGenetic {
			mthd = "genetic"
		}
		ep.InitSimple(cm.startPoint, cm.dimension, func(xs []float64) float64 { return cm.helpFunc(xs, fIdeal, ws) }, cm.penalties,
			gradient, cm.gradientConstraint, cm.constraint, 0.0001, 1.618, mthd)

		xMin, _, err = ep.Solve()
		if err != nil {
			return nil, nil, fmt.Errorf("error solving external penalty method : %v", err)
		}
		xMins = append(xMins, xMin)
		yMin = append(yMin, []float64{cm.targetFuncs[0](xMin), cm.targetFuncs[1](xMin)})

	}

	return xMins, yMin, nil
}
