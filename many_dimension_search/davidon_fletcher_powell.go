package many_dimension_search

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/interpolation_search"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"github.com/saskamegaprogrammist/optimization_methods/one_dimension_search"
	"math"
)

type DavidonFletcherPowellSearch struct {
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

func (dfps *DavidonFletcherPowellSearch) Init(startPoint []float64, delta float64, dimension int,
	eps1 float64, eps2 float64, alphaPrecision float64,
	oneDStep float64, maxIter int, targetFunc func(xs []float64) float64, gradient []func(xs []float64) float64,
	method string) {
	dfps.startPoint = startPoint
	dfps.delta = delta
	dfps.eps1 = eps1
	dfps.eps2 = eps2
	dfps.targetFunc = targetFunc
	dfps.dimension = dimension
	dfps.method = method
	dfps.alphaPrecision = alphaPrecision
	dfps.oneDStep = oneDStep
	dfps.maxIter = maxIter
	dfps.gradient = gradient
}

func (dfps *DavidonFletcherPowellSearch) Solve() ([]float64, float64, error) {
	var err error
	var G, GNew la_methods.Matrix
	var x, xOld, xSub la_methods.Vector
	var k int
	var grad, gradOld, gradMinus, d, dInter la_methods.Vector
	var alpha float64
	var fib, sqrInterp, lastIter bool
	alpha = dfps.alphaPrecision
	var search one_dimension_search.OneDimensionSearchI
	if dfps.method == "break in two" {
		search = &one_dimension_search.BreakInTwoSearch{}
	} else if dfps.method == "golden ratio" {
		search = &one_dimension_search.GoldenRatioSearch{}
	} else if dfps.method == "fibonacci" {
		fib = true
	} else if dfps.method == "square interpolation" {
		sqrInterp = true
	} else {
		return []float64{}, 0, fmt.Errorf("wrong one dimensional method: %s", dfps.method)
	}
	alpha = dfps.alphaPrecision
	err = x.InitWithPoints(dfps.dimension, dfps.startPoint)
	if err != nil {
		return nil, 0, fmt.Errorf("error during vector initializing: %v", err)
	}
	xOld = x.Copy()
	G.Init(dfps.dimension, dfps.dimension)
	G.E()
	for {
		grad, err = dfps.calculateGradient(x)
		if err != nil {
			return nil, 0, fmt.Errorf("error calculcating gradient: %v", err)
		}
		if grad.Len() < dfps.eps1 || k >= dfps.maxIter {
			fmt.Printf("k value: %d\n", k)
			return x.Points, dfps.targetFunc(x.Points), nil
		}
		if k > 0 {
			gradOld, err = dfps.calculateGradient(xOld)
			if err != nil {
				return nil, 0, fmt.Errorf("error calculcating gradient: %v", err)
			}
			xSub, err = x.Sub(xOld)
			if err != nil {
				return nil, 0, fmt.Errorf("error during vector substracting: %v", err)
			}
			GNew, err = dfps.calculateG(grad, gradOld, xSub, G)
		} else {
			GNew = G
		}

		gradMinus = grad.MulOnValue(-1)
		d, err = GNew.MulV(gradMinus)
		if err != nil {
			return nil, 0, fmt.Errorf("error during vector and matrix multiplying: %v", err)
		}
		if fib {
			alpha, err = dfps.oneDimensionFibonacciSearch(x, d, alpha)
		} else if sqrInterp {
			alpha, err = dfps.oneDimensionSquareInterpolation(x, d, alpha)
		} else {
			alpha, err = dfps.oneDimensionSearch(x, d, alpha, search)
		}
		dInter = d.MulOnValue(alpha)
		xOld = x
		x, err = x.Add(dInter)
		if err != nil {
			return nil, 0, fmt.Errorf("error during vector adding: %v", err)
		}
		xSub, err = x.Sub(xOld)
		if err != nil {
			return nil, 0, fmt.Errorf("error during vector substracting: %v", err)
		}
		if xSub.Len() < dfps.delta && math.Abs(dfps.targetFunc(x.Points)-dfps.targetFunc(xOld.Points)) < dfps.eps2 {
			if lastIter {
				fmt.Printf("k value: %d\n", k)
				return x.Points, dfps.targetFunc(x.Points), nil
			} else {
				lastIter = true
			}
		}
		//fmt.Println(dfps.targetFunc(x.Points))
		k++
	}
}

func (dfps *DavidonFletcherPowellSearch) calculateG(gradNew la_methods.Vector, gradOld la_methods.Vector,
	deltaX la_methods.Vector, G la_methods.Matrix) (la_methods.Matrix, error) {
	var err error
	var deltaXdeltaGrad, ggF float64
	var deltaXMRow, deltaXMCol, A, B, deltaX2, ggCol, ggRow, ggM, deltaG, GNew la_methods.Matrix
	var gradDelta, gg la_methods.Vector
	gradDelta, err = gradNew.Sub(gradOld)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during vector substracting: %v", err)
	}
	err = deltaXMCol.InitWithVectorColumn(dfps.dimension, dfps.dimension, deltaX)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during vector initalizing: %v", err)
	}
	err = deltaXMRow.InitWithVectorRow(dfps.dimension, dfps.dimension, deltaX)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during vector initalizing: %v", err)
	}
	deltaXdeltaGrad, err = deltaX.Mul(gradDelta)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during vector multiplying: %v", err)
	}
	deltaX2, err = deltaXMCol.MulM(deltaXMRow)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during matrix multiplying: %v", err)
	}
	A = deltaX2.MulVal(float64(1) / deltaXdeltaGrad)

	gg, err = G.MulV(gradDelta)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during matrix and vector multiplying: %v", err)
	}
	err = ggCol.InitWithVectorColumn(dfps.dimension, dfps.dimension, gg)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during vector initalizing: %v", err)
	}
	ggRow, err = ggCol.Transponate()
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during vector transonation: %v", err)
	}
	ggM, err = ggCol.MulM(ggRow)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during matrix multiplying: %v", err)
	}
	ggF, err = gradDelta.Mul(gg)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during vector multiplying: %v", err)
	}
	B = ggM.MulVal(float64(-1) / ggF)
	deltaG, err = A.AddM(B)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during matrix adding: %v", err)
	}
	GNew, err = G.AddM(deltaG)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during matrix adding: %v", err)
	}
	return GNew, nil
}

func (dfps *DavidonFletcherPowellSearch) calculateGradient(x la_methods.Vector) (la_methods.Vector, error) {
	var gradPoints []float64
	var grad la_methods.Vector
	for i := 0; i < dfps.dimension; i++ {
		gradPoints = append(gradPoints, dfps.gradient[i](x.Points))
	}
	err := grad.InitWithPoints(dfps.dimension, gradPoints)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error during vector initializing: %v", err)
	}
	return grad, nil
}

func (dfps *DavidonFletcherPowellSearch) findBounds(alpha float64, targetFunc func(x float64) float64) (float64, float64, error) {
	dfps.svennAlgorithm.Init(dfps.eps1, alpha, targetFunc)
	return dfps.svennAlgorithm.Solve()
}

func (dfps *DavidonFletcherPowellSearch) oneDimensionFibonacciSearch(y la_methods.Vector, d la_methods.Vector, alpha float64) (float64, error) {
	var fib one_dimension_search.FibonacciSearch
	fOneD := dfps.getOneDimensionFunc(d, y)
	a, b, err := dfps.findBounds(alpha, fOneD)
	if err != nil {
		return 0, fmt.Errorf("error during svenn algorithm: %v", err)
	}
	fib.Init(a, b, dfps.oneDStep, dfps.oneDStep, fOneD)
	min, _, err := fib.Solve()
	if err != nil {
		return 0, fmt.Errorf("error during fibonacci search: %v", err)
	}
	return min, nil
}

func (dfps *DavidonFletcherPowellSearch) oneDimensionSquareInterpolation(y la_methods.Vector, d la_methods.Vector, alpha float64) (float64, error) {
	var sqrInter interpolation_search.SquareInterpolation
	fOneD := dfps.getOneDimensionFunc(d, y)
	sqrInter.Init(dfps.alphaPrecision, dfps.oneDStep, dfps.eps1, dfps.eps1, fOneD)
	min, _ := sqrInter.Solve()
	return min, nil
}

func (dfps *DavidonFletcherPowellSearch) oneDimensionSearch(y la_methods.Vector, d la_methods.Vector, alpha float64, search one_dimension_search.OneDimensionSearchI) (float64, error) {
	fOneD := dfps.getOneDimensionFunc(d, y)
	a, b, err := dfps.findBounds(alpha, fOneD)
	if err != nil {
		return 0, fmt.Errorf("error during svenn algorithm: %v", err)
	}
	search.Init(a, b, dfps.oneDStep, fOneD)
	min, _ := search.Solve()
	return min, nil
}

func (dfps *DavidonFletcherPowellSearch) getOneDimensionFunc(d la_methods.Vector, y la_methods.Vector) func(x float64) float64 {
	return func(x float64) float64 {
		aplhaD := d.MulOnValue(x)
		sumYAlphD, _ := y.Add(aplhaD)
		return dfps.targetFunc(sumYAlphD.Points)
	}
}
