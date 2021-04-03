package many_dimension_search

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"github.com/saskamegaprogrammist/optimization_methods/one_dimension_search"
)

type HookeJeevesSearch struct {
	startPoint     []float64
	precision      float64
	alphaPrecision float64
	delta          float64
	lambda         float64
	dimension      int
	targetFunc     func(xs []float64) float64
	method         string
	oneDStep       float64
	svennAlgorithm one_dimension_search.Svenn
}

func (hjs *HookeJeevesSearch) Init(startPoint []float64, delta float64, dimension int,
	lambda float64, precision float64, alphaPrecision float64,
	oneDStep float64, targetFunc func(xs []float64) float64, method string) {
	hjs.startPoint = startPoint
	hjs.delta = delta
	hjs.lambda = lambda
	hjs.precision = precision
	hjs.targetFunc = targetFunc
	hjs.dimension = dimension
	hjs.method = method
	hjs.alphaPrecision = alphaPrecision
	hjs.oneDStep = oneDStep
}

func (hjs *HookeJeevesSearch) Solve() ([]float64, float64, error) {
	var err error
	var y, yPrev la_methods.Vector
	var x la_methods.Vector
	var i, k int
	var delta la_methods.Vector
	var alpha float64
	var stop, fib bool
	alpha = hjs.alphaPrecision
	var search one_dimension_search.OneDimensionSearchI
	if hjs.method == "break in two" {
		search = &one_dimension_search.BreakInTwoSearch{}
	} else if hjs.method == "golden ratio" {
		search = &one_dimension_search.GoldenRatioSearch{}
	} else if hjs.method == "fibonacci" {
		fib = true
	} else {
		return []float64{}, 0, fmt.Errorf("wrong one dimensional method: %s", hjs.method)
	}
	err = x.InitWithPoints(hjs.dimension, hjs.startPoint)
	if err != nil {
		return []float64{}, 0, fmt.Errorf("error initializing vector: %v", err)
	}
	y = x.Copy()
	yPrev = y.Copy()
	k = 1
	delta.InitWithValue(hjs.dimension, hjs.delta)
	for {
		x, err := hjs.research(i, y, delta)
		if err != nil {
			return []float64{}, 0, fmt.Errorf("error researching: %v", err)
		}

		fY := hjs.targetFunc(y.Points)
		fX := hjs.targetFunc(x.Points)
		if fX < fY {
			yPrev = y
			y = x
			yInterm, err := y.Sub(yPrev)
			if err != nil {
				return []float64{}, 0, fmt.Errorf("error substracting: %v", err)
			}
			x, err = y.Add(yInterm.MulOnValue(hjs.lambda))
			if err != nil {
				return []float64{}, 0, fmt.Errorf("error adding: %v", err)
			}
			k++
		}

		d, err := y.Sub(yPrev)
		if err != nil {
			return []float64{}, 0, fmt.Errorf("error substracting: %v", err)
		}
		if fib {
			alpha, err = hjs.oneDimensionFibonacciSearch(y, d, alpha)
		} else {
			alpha, err = hjs.oneDimensionSearch(y, d, alpha, search)
		}
		if err != nil {
			return []float64{}, 0, fmt.Errorf("error during one dimension search: %v", err)
		}
		delta, stop = hjs.checkStop(alpha, delta, hjs.precision)
		//fmt.Println(hjs.targetFunc(y.Points))
		if stop {
			//fmt.Printf("k value: %d\n", k)
			return y.Points, hjs.targetFunc(y.Points), nil
		}
		aplhaD := d.MulOnValue(alpha)
		y, err = y.Add(aplhaD)
		yPrev = y
		if err != nil {
			return []float64{}, 0, fmt.Errorf("error adding: %v", err)
		}
		k++
	}
}

func (hjs *HookeJeevesSearch) findBounds(alpha float64, targetFunc func(x float64) float64) (float64, float64, error) {
	hjs.svennAlgorithm.Init(hjs.precision, alpha, targetFunc)
	return hjs.svennAlgorithm.Solve()
}

func (hjs *HookeJeevesSearch) oneDimensionFibonacciSearch(y la_methods.Vector, d la_methods.Vector, alpha float64) (float64, error) {
	var fib one_dimension_search.FibonacciSearch
	fOneD := hjs.getOneDimensionFunc(d, y)
	a, b, err := hjs.findBounds(alpha, fOneD)
	if err != nil {
		return 0, fmt.Errorf("error during svenn algorithm: %v", err)
	}
	fib.Init(a, b, hjs.oneDStep, hjs.oneDStep, fOneD)
	min, _, err := fib.Solve()
	if err != nil {
		return 0, fmt.Errorf("error during fibonacci search: %v", err)
	}
	return min, nil
}

func (hjs *HookeJeevesSearch) oneDimensionSearch(y la_methods.Vector, d la_methods.Vector, alpha float64, search one_dimension_search.OneDimensionSearchI) (float64, error) {
	fOneD := hjs.getOneDimensionFunc(d, y)
	a, b, err := hjs.findBounds(alpha, fOneD)
	if err != nil {
		return 0, fmt.Errorf("error during svenn algorithm: %v", err)
	}
	search.Init(a, b, hjs.oneDStep, fOneD)
	min, _ := search.Solve()
	return min, nil
}

func (hjs *HookeJeevesSearch) getOneDimensionFunc(d la_methods.Vector, y la_methods.Vector) func(x float64) float64 {
	return func(x float64) float64 {
		aplhaD := d.MulOnValue(x)
		sumYAlphD, _ := y.Add(aplhaD)
		return hjs.targetFunc(sumYAlphD.Points)
	}
}

func (hjs *HookeJeevesSearch) checkStop(alpha float64, delta la_methods.Vector, eps float64) (la_methods.Vector, bool) {
	//fmt.Println(alpha)
	stop := alpha < eps
	if !stop {
		for i, d := range delta.Points {
			if d > eps {
				delta.Points[i] = d / 2
			}
		}
	}
	return delta, stop
}

func (hjs *HookeJeevesSearch) research(iInit int, x la_methods.Vector, delta la_methods.Vector) (la_methods.Vector, error) {
	for i := iInit; i < hjs.dimension; i++ {
		f := hjs.targetFunc(x.Points)
		xAdd, err := x.AddKOnIndex(delta.Points[i], i)
		if err != nil {
			return x, fmt.Errorf("error adding delta value: %v", err)
		}
		f1 := hjs.targetFunc(xAdd.Points)

		if f1 < f {
			x = xAdd
		} else {
			xSub, err := x.SubKOnIndex(delta.Points[i], i)
			if err != nil {
				return x, fmt.Errorf("error substracting delta value: %v", err)
			}
			f2 := hjs.targetFunc(xSub.Points)
			if f2 < f {
				x = xSub
			}
		}
	}
	return x, nil
}
