package many_dimension_search

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/interpolation_search"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"github.com/saskamegaprogrammist/optimization_methods/one_dimension_search"
	"math"
)

type FastGradientDescendSearch struct {
	startPoint     []float64
	eps1           float64
	eps2           float64
	dimension      int
	targetFunc     func(xs []float64) float64
	gradient       []func(xs []float64) float64
	method         string
	precision      float64
	alphaPrecision float64
	svennAlgorithm one_dimension_search.Svenn
}

func (fgd *FastGradientDescendSearch) Init(startPoint []float64, eps1 float64, eps2 float64,
	targetFunc func(xs []float64) float64, gradient []func(xs []float64) float64, dimension int, precision float64, alphaPrecision float64, method string) {
	fgd.startPoint = startPoint
	fgd.eps1 = eps1
	fgd.eps2 = eps2
	fgd.targetFunc = targetFunc
	fgd.gradient = gradient
	fgd.dimension = dimension
	fgd.method = method
	fgd.precision = precision
	fgd.alphaPrecision = alphaPrecision
}

func (fgd *FastGradientDescendSearch) Solve() ([]float64, float64, error) {
	var err error
	var alpha float64
	var grad, d, x, xNew, alphaGrad la_methods.Vector
	var k int
	var fib, sqrInter bool
	var search one_dimension_search.OneDimensionSearchI
	if fgd.method == "break in two" {
		search = &one_dimension_search.BreakInTwoSearch{}
	} else if fgd.method == "golden ratio" {
		search = &one_dimension_search.GoldenRatioSearch{}
	} else if fgd.method == "fibonacci" {
		fib = true
	} else if fgd.method == "square interpolation" {
		sqrInter = true
	} else {
		return []float64{}, 0, fmt.Errorf("wrong one dimensional method: %s", fgd.method)
	}
	alpha = fgd.alphaPrecision
	err = x.InitWithPoints(fgd.dimension, fgd.startPoint)
	if err != nil {
		return nil, 0, fmt.Errorf("error during vector initializing: %v", err)
	}
	for {
		grad, err = fgd.calculateGradient(x)
		if err != nil {
			return nil, 0, fmt.Errorf("error calculcating gradient: %v", err)
		}
		if grad.Len() < fgd.eps1 {
			fmt.Printf("k value: %d\n", k)
			return x.Points, fgd.targetFunc(x.Points), nil
		}
		d = grad.MulOnValue(-1)
		if fib {
			alpha, err = fgd.oneDimensionFibonacciSearch(x, d, alpha)
		} else if sqrInter {
			alpha, err = fgd.oneDimensionSquareInterpolation(x, d, alpha)
		} else {
			alpha, err = fgd.oneDimensionSearch(x, d, alpha, search)
		}
		alphaGrad = grad.MulOnValue(alpha)
		xNew, err = x.Sub(alphaGrad)
		if err != nil {
			return nil, 0, fmt.Errorf("error during vector substracting: %v", err)
		}
		f := fgd.targetFunc(x.Points)
		fNew := fgd.targetFunc(xNew.Points)
		if alphaGrad.Len() < fgd.eps1 && math.Abs(fNew-f) < fgd.eps2 {
			fmt.Printf("k value: %d\n", k)
			return xNew.Points, fNew, nil
		} else {
			x = xNew
			k++
		}
	}
}

func (fgd *FastGradientDescendSearch) calculateGradient(x la_methods.Vector) (la_methods.Vector, error) {
	var gradPoints []float64
	var grad la_methods.Vector
	for i := 0; i < fgd.dimension; i++ {
		gradPoints = append(gradPoints, fgd.gradient[i](x.Points))
	}
	err := grad.InitWithPoints(fgd.dimension, gradPoints)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error during vector initializing: %v", err)
	}
	return grad, nil
}

func (fgd *FastGradientDescendSearch) findBounds(alpha float64, targetFunc func(x float64) float64) (float64, float64, error) {
	fgd.svennAlgorithm.Init(fgd.precision, alpha, targetFunc)
	return fgd.svennAlgorithm.Solve()
}

func (fgd *FastGradientDescendSearch) oneDimensionSquareInterpolation(y la_methods.Vector, d la_methods.Vector, alpha float64) (float64, error) {
	var sqrInter interpolation_search.SquareInterpolation
	fOneD := fgd.getOneDimensionFunc(d, y)
	sqrInter.Init(alpha, 0.001, 0.0001, 0.00001, fOneD)
	min, _ := sqrInter.Solve()
	return min, nil
}

func (fgd *FastGradientDescendSearch) oneDimensionFibonacciSearch(y la_methods.Vector, d la_methods.Vector, alpha float64) (float64, error) {
	var fib one_dimension_search.FibonacciSearch
	fOneD := fgd.getOneDimensionFunc(d, y)
	a, b, err := fgd.findBounds(alpha, fOneD)
	if err != nil {
		return 0, fmt.Errorf("error during svenn algorithm: %v", err)
	}
	fib.Init(a, b, fgd.precision, fgd.precision, fOneD)
	min, _, err := fib.Solve()
	if err != nil {
		return 0, fmt.Errorf("error during fibonacci search: %v", err)
	}
	return min, nil
}

func (fgd *FastGradientDescendSearch) oneDimensionSearch(y la_methods.Vector, d la_methods.Vector, alpha float64, search one_dimension_search.OneDimensionSearchI) (float64, error) {
	fOneD := fgd.getOneDimensionFunc(d, y)
	a, b, err := fgd.findBounds(alpha, fOneD)
	if err != nil {
		return 0, fmt.Errorf("error during svenn algorithm: %v", err)
	}
	search.Init(a, b, fgd.precision, fOneD)
	min, _ := search.Solve()
	return min, nil
}

func (fgd *FastGradientDescendSearch) getOneDimensionFunc(d la_methods.Vector, y la_methods.Vector) func(x float64) float64 {
	return func(x float64) float64 {
		aplhaD := d.MulOnValue(x)
		sumYAlphD, _ := y.Add(aplhaD)
		return fgd.targetFunc(sumYAlphD.Points)
	}
}
