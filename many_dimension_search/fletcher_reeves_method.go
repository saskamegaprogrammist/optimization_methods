package many_dimension_search

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/interpolation_search"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"github.com/saskamegaprogrammist/optimization_methods/one_dimension_search"
	"math"
)

type FletcherReevesSearch struct {
	startPoint     []float64
	eps1           float64
	eps2           float64
	delta          float64
	alphaPrecision float64
	dimension      int
	targetFunc     func(xs []float64) float64
	gradient       []func(xs []float64) float64
	method         string
	svennAlgorithm one_dimension_search.Svenn
	oneDStep       float64
	maxIter        int
}

func (frs *FletcherReevesSearch) Init(startPoint []float64, delta float64, dimension int,
	eps1 float64, eps2 float64, alphaPrecision float64,
	oneDStep float64, maxIter int, targetFunc func(xs []float64) float64, gradient []func(xs []float64) float64,
	method string) {
	frs.startPoint = startPoint
	frs.delta = delta
	frs.eps1 = eps1
	frs.eps2 = eps2
	frs.targetFunc = targetFunc
	frs.dimension = dimension
	frs.method = method
	frs.alphaPrecision = alphaPrecision
	frs.oneDStep = oneDStep
	frs.maxIter = maxIter
	frs.gradient = gradient
}

func (frs *FletcherReevesSearch) Solve() ([]float64, float64, error) {
	var err error
	var x, xOld, xSub la_methods.Vector
	var k int
	var grad, gradOld, gradMinus, d, dNew, dInter la_methods.Vector
	var alpha, w float64
	var fib, sqrInterp, lastIter bool
	alpha = frs.alphaPrecision
	var search one_dimension_search.OneDimensionSearchI
	if frs.method == "break in two" {
		search = &one_dimension_search.BreakInTwoSearch{}
	} else if frs.method == "golden ratio" {
		search = &one_dimension_search.GoldenRatioSearch{}
	} else if frs.method == "fibonacci" {
		fib = true
	} else if frs.method == "square interpolation" {
		sqrInterp = true
	} else {
		return []float64{}, 0, fmt.Errorf("wrong one dimensional method: %s", frs.method)
	}
	alpha = frs.alphaPrecision
	err = x.InitWithPoints(frs.dimension, frs.startPoint)
	if err != nil {
		return nil, 0, fmt.Errorf("error during vector initializing: %v", err)
	}
	xOld = x.Copy()
	for {
		grad, err = frs.calculateGradient(x)
		if err != nil {
			return nil, 0, fmt.Errorf("error calculcating gradient: %v", err)
		}
		if grad.Len() < frs.eps1 || k >= frs.maxIter {
			fmt.Printf("k value: %d\n", k)
			return x.Points, frs.targetFunc(x.Points), nil
		}
		if k == 0 {
			d = grad.MulOnValue(-1)
		}
		gradOld, err = frs.calculateGradient(xOld)
		if err != nil {
			return nil, 0, fmt.Errorf("error calculcating gradient: %v", err)
		}
		if k%frs.dimension == 0 {
			w = 0
		} else {
			w = math.Pow(grad.Len(), 2) / math.Pow(gradOld.Len(), 2)
		}
		dInter = d.MulOnValue(w)
		gradMinus = grad.MulOnValue(-1)
		dNew, err = gradMinus.Add(dInter)
		if err != nil {
			return nil, 0, fmt.Errorf("error during vector adding: %v", err)
		}
		if fib {
			alpha, err = frs.oneDimensionFibonacciSearch(x, dNew, alpha)
		} else if sqrInterp {
			alpha, err = frs.oneDimensionSquareInterpolation(x, dNew, alpha)
		} else {
			alpha, err = frs.oneDimensionSearch(x, dNew, alpha, search)
		}
		dInter = dNew.MulOnValue(alpha)
		xOld = x
		x, err = x.Add(dInter)
		if err != nil {
			return nil, 0, fmt.Errorf("error during vector adding: %v", err)
		}
		xSub, err = x.Sub(xOld)
		if err != nil {
			return nil, 0, fmt.Errorf("error during vector substracting: %v", err)
		}
		if xSub.Len() < frs.delta && math.Abs(frs.targetFunc(x.Points)-frs.targetFunc(xOld.Points)) < frs.eps2 {
			if lastIter {
				fmt.Printf("k value: %d\n", k)
				return x.Points, frs.targetFunc(x.Points), nil
			} else {
				lastIter = true
			}
		}
		d = dNew
		//fmt.Println(frs.targetFunc(x.Points))
		k++
	}
}

func (frs *FletcherReevesSearch) calculateGradient(x la_methods.Vector) (la_methods.Vector, error) {
	var gradPoints []float64
	var grad la_methods.Vector
	for i := 0; i < frs.dimension; i++ {
		gradPoints = append(gradPoints, frs.gradient[i](x.Points))
	}
	err := grad.InitWithPoints(frs.dimension, gradPoints)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error during vector initializing: %v", err)
	}
	return grad, nil
}

func (frs *FletcherReevesSearch) findBounds(alpha float64, targetFunc func(x float64) float64) (float64, float64, error) {
	frs.svennAlgorithm.Init(frs.eps1, alpha, targetFunc)
	return frs.svennAlgorithm.Solve()
}

func (frs *FletcherReevesSearch) oneDimensionFibonacciSearch(y la_methods.Vector, d la_methods.Vector, alpha float64) (float64, error) {
	var fib one_dimension_search.FibonacciSearch
	fOneD := frs.getOneDimensionFunc(d, y)
	a, b, err := frs.findBounds(alpha, fOneD)
	if err != nil {
		return 0, fmt.Errorf("error during svenn algorithm: %v", err)
	}
	fib.Init(a, b, frs.oneDStep, frs.oneDStep, fOneD)
	min, _, err := fib.Solve()
	if err != nil {
		return 0, fmt.Errorf("error during fibonacci search: %v", err)
	}
	return min, nil
}

func (frs *FletcherReevesSearch) oneDimensionSquareInterpolation(y la_methods.Vector, d la_methods.Vector, alpha float64) (float64, error) {
	var sqrInter interpolation_search.SquareInterpolation
	fOneD := frs.getOneDimensionFunc(d, y)
	sqrInter.Init(frs.alphaPrecision, frs.oneDStep, frs.eps1, frs.eps1, fOneD)
	min, _ := sqrInter.Solve()
	return min, nil
}

func (frs *FletcherReevesSearch) oneDimensionSearch(y la_methods.Vector, d la_methods.Vector, alpha float64, search one_dimension_search.OneDimensionSearchI) (float64, error) {
	fOneD := frs.getOneDimensionFunc(d, y)
	a, b, err := frs.findBounds(alpha, fOneD)
	if err != nil {
		return 0, fmt.Errorf("error during svenn algorithm: %v", err)
	}
	search.Init(a, b, frs.oneDStep, fOneD)
	min, _ := search.Solve()
	return min, nil
}

func (frs *FletcherReevesSearch) getOneDimensionFunc(d la_methods.Vector, y la_methods.Vector) func(x float64) float64 {
	return func(x float64) float64 {
		aplhaD := d.MulOnValue(x)
		sumYAlphD, _ := y.Add(aplhaD)
		return frs.targetFunc(sumYAlphD.Points)
	}
}
