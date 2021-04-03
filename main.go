package main

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/constraint_methods"
	"github.com/saskamegaprogrammist/optimization_methods/interpolation_search"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"github.com/saskamegaprogrammist/optimization_methods/many_dimension_search"
	"github.com/saskamegaprogrammist/optimization_methods/one_dimension_search"
	"math"
	"time"
)

func rozenbrokeOneDimensionFunction(a float64, b float64, f float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		var result float64
		n := len(xs)
		for i, x := range xs {
			if i != n-1 {
				result += a*math.Pow(math.Pow(x, 2)-xs[i+1], 2) + b*math.Pow(x-1, 2)
			}
		}
		result += f
		return result
	}
}

func rozenbrokeOneDimensionFunction2() func(xs []float64) float64 {
	var a, b, f float64
	a = 158
	b = 2
	f = 40
	return func(xs []float64) float64 {
		var result float64
		n := len(xs)
		for i, x := range xs {
			if i != n-1 {
				result += a*math.Pow(math.Pow(x, 2)-xs[i+1], 2) + b*math.Pow(x-1, 2)
			}
		}
		result += f
		return result
	}
}

func rozenbrokeFirstGrad(a, b, f float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		return a*2*2*xs[0]*(math.Pow(xs[0], 2)-xs[1]) + 2*b*(xs[0]-1)
	}
}

func rozenbrokeFirstGrad2() func(xs []float64) float64 {
	var a, b float64
	a = 158
	b = 2
	return func(xs []float64) float64 {
		return a*2*2*xs[0]*(math.Pow(xs[0], 2)-xs[1]) + 2*b*(xs[0]-1)
	}
}

func rozenbrokeFirstGradFirst(a float64, b float64, f float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		return a*2*2*(3*math.Pow(xs[0], 2)-xs[1]) + 2*b
	}
}

func rozenbrokeFirstGradFirst2() func(xs []float64) float64 {
	var a, b float64
	a = 158
	b = 2
	return func(xs []float64) float64 {
		return a*2*2*(3*math.Pow(xs[0], 2)-xs[1]) + 2*b
	}
}

func rozenbrokeFirstGradSecond(a float64, b float64, f float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		return (-1) * a * 2 * 2 * xs[0]
	}
}

func rozenbrokeFirstGradSecond2() func(xs []float64) float64 {
	var a float64 = 158
	return func(xs []float64) float64 {
		return (-1) * a * 2 * 2 * xs[0]
	}
}

func rozenbrokeSecondGrad(a float64, b float64, f float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		return -a*2*(math.Pow(xs[0], 2)-xs[1]) + 2*b*(xs[1]-1) + 2*2*a*xs[1]*(math.Pow(xs[1], 2)-xs[2])
	}
}

func rozenbrokeSecondGradFirst(a float64, b float64, f float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		return -2 * 2 * a * xs[0]
	}
}

func rozenbrokeSecondGradSecond(a float64, b float64, f float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		return 12*a*math.Pow(xs[1], 2) - 4*a*xs[2] + 2*a + 2*b
	}
}

func rozenbrokeSecondGradThird(a float64, b float64, f float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		return -2 * 2 * a * xs[1]
	}
}

func rozenbrokeThirdGrad(a float64, b float64, f float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		return -a * 2 * (math.Pow(xs[1], 2) - xs[2])
	}
}

func rozenbrokeThirdGradSecond(a float64, b float64, f float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		return -2 * 2 * a * xs[1]
	}
}

func rozenbrokeThirdGradThird(a float64, b float64, f float64) func(xs []float64) float64 {
	return func(xs []float64) float64 {
		return 2 * a
	}
}

func rozenbroke2SecondGrad() func(xs []float64) float64 {
	var a float64 = 158
	return func(xs []float64) float64 {
		return -a * 2 * (math.Pow(xs[0], 2) - xs[1])
	}
}

func rozenbroke2SecondGradFirst() func(xs []float64) float64 {
	var a float64 = 158
	return func(xs []float64) float64 {
		return -2 * 2 * a * xs[0]
	}
}

func rozenbroke2SecondGradSecond() func(xs []float64) float64 {
	var a float64 = 158
	return func(xs []float64) float64 {
		return 2 * a
	}
}

func hess(a float64, b float64, f float64) func(xs []float64) la_methods.Matrix {
	return func(xs []float64) la_methods.Matrix {
		return hessian(xs, a, b, f)
	}
}

func hess2() func(xs []float64) la_methods.Matrix {
	return func(xs []float64) la_methods.Matrix {
		return hessian2(xs)
	}
}

func hessian(xs []float64, a float64, b float64, f float64) la_methods.Matrix {
	var hess la_methods.Matrix
	dim := len(xs)
	hess.Init(dim, dim)
	hess.Points[0][0] = rozenbrokeFirstGradFirst(a, b, f)(xs)
	hess.Points[0][1] = rozenbrokeFirstGradSecond(a, b, f)(xs)
	hess.Points[1][0] = rozenbrokeSecondGradFirst(a, b, f)(xs)
	hess.Points[1][1] = rozenbrokeSecondGradSecond(a, b, f)(xs)
	hess.Points[1][2] = rozenbrokeSecondGradThird(a, b, f)(xs)
	hess.Points[2][1] = rozenbrokeThirdGradSecond(a, b, f)(xs)
	hess.Points[2][1] = rozenbrokeThirdGradThird(a, b, f)(xs)
	return hess
}

func hessian2(xs []float64) la_methods.Matrix {
	var hess la_methods.Matrix
	dim := len(xs)
	hess.Init(dim, dim)
	hess.Points[0][0] = rozenbrokeFirstGradFirst2()(xs)
	hess.Points[0][1] = rozenbrokeFirstGradSecond2()(xs)
	hess.Points[1][0] = rozenbroke2SecondGradFirst()(xs)
	hess.Points[1][1] = rozenbroke2SecondGradSecond()(xs)
	return hess
}

func targetFunction(x float64) float64 {
	return 100*math.Pow(math.Pow(x, 2)-2, 3) + math.Pow(x-1, 2) - math.Abs(10+x)
}

func targetFunctionDerivative(x float64) float64 {
	return 600*math.Pow(math.Pow(x, 2)-2, 2)*x + 2*(x-1) - (10+x)/math.Abs(10+x)
}

func first() {
	var timeStart, timeEnd time.Time
	var svenn one_dimension_search.Svenn
	var a, b float64
	var err error
	precision := 0.000001
	timeStart = time.Now()
	svenn.Init(0.005, -2, targetFunction)
	a, b, err = svenn.Solve()
	timeEnd = time.Now()
	if err != nil {
		fmt.Printf("error finding uncertainty interval: %v\n", err)
		return
	}
	fmt.Printf("interval: %f, %f\n", a, b)
	fmt.Printf("svenn algorithm took : %v\n", timeEnd.Sub(timeStart))

	var bit one_dimension_search.BreakInTwoSearch
	var xMin, fMin float64
	timeStart = time.Now()
	bit.Init(a, b, precision, targetFunction)
	xMin, fMin = bit.Solve()
	timeEnd = time.Now()
	fmt.Printf("minimum: %f, %f\n", xMin, fMin)
	bitConvergence, err := bit.CountConvergence()
	if err != nil {
		fmt.Printf("error counting convergence: %v\n", err)
		return
	}
	fmt.Printf("break in two convergence: %f\n", bitConvergence)
	fmt.Printf("break in two search algorithm took : %v\n", timeEnd.Sub(timeStart))

	var gr one_dimension_search.GoldenRatioSearch
	timeStart = time.Now()
	gr.Init(a, b, precision, targetFunction)
	xMin, fMin = gr.Solve()
	timeEnd = time.Now()
	fmt.Printf("minimum: %f, %f\n", xMin, fMin)
	grConvergence, err := gr.CountConvergence()
	if err != nil {
		fmt.Printf("error counting convergence: %v\n", err)
		return
	}
	fmt.Printf("golden ratio convergence: %f\n", grConvergence)
	fmt.Printf("golden ratio search algorithm took : %v\n", timeEnd.Sub(timeStart))

	var fs one_dimension_search.FibonacciSearch
	timeStart = time.Now()
	fs.Init(a, b, precision, precision, targetFunction)
	xMin, fMin, err = fs.Solve()
	timeEnd = time.Now()
	if err != nil {
		fmt.Printf("error solving fibonacci: %v\n", err)
		return
	}
	fmt.Printf("minimum: %f, %f\n", xMin, fMin)
	fsConvergence, err := fs.CountConvergence()
	if err != nil {
		fmt.Printf("error counting convergence: %v\n", err)
		return
	}
	fmt.Printf("fibonacci convergence: %f\n", fsConvergence)
	fmt.Printf("fibonacci search algorithm took : %v\n", timeEnd.Sub(timeStart))
}

func second() {
	var timeStart, timeEnd time.Time
	var xMin, fMin float64
	precision := 0.000001

	var sqrInt interpolation_search.SquareInterpolation
	timeStart = time.Now()
	sqrInt.Init(-2, 0.005, precision, precision, targetFunction)
	xMin, fMin = sqrInt.Solve()
	timeEnd = time.Now()

	fmt.Printf("minimum: %f, %f\n", xMin, fMin)
	fmt.Printf("square interpolation search algorithm took : %v\n", timeEnd.Sub(timeStart))

	var cInt interpolation_search.CubicInterpolation
	timeStart = time.Now()
	cInt.Init(-2, 0.005, precision, precision, targetFunction, targetFunctionDerivative)
	xMin, fMin = cInt.Solve()
	timeEnd = time.Now()

	fmt.Printf("minimum: %f, %f\n", xMin, fMin)
	fmt.Printf("cubic interpolation search algorithm took : %v\n", timeEnd.Sub(timeStart))
}

func third() {
	var err error
	var timeStart, timeEnd time.Time
	var xMin []float64
	var yMin float64
	precision := 0.0001
	alphaPrecision := 0.0000001
	oneDStep := 0.1
	rFunc := rozenbrokeOneDimensionFunction(100, 2, 45)

	var hjs many_dimension_search.HookeJeevesSearch
	timeStart = time.Now()
	hjs.Init([]float64{0, 0, 0}, precision*10, 3, 2, precision, alphaPrecision,
		oneDStep, rFunc, "break in two")
	xMin, yMin, err = hjs.Solve()
	if err != nil {
		fmt.Printf("error solving hooke jeeves: %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("hooke jeeves search break in two algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	hjs.Init([]float64{0, 0, 0}, precision*10, 3, 2, precision, alphaPrecision,
		oneDStep, rFunc, "golden ratio")
	xMin, yMin, err = hjs.Solve()
	if err != nil {
		fmt.Printf("error solving hooke jeeves: %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("hooke jeeves search golden ratio algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	hjs.Init([]float64{0, 0, 0}, precision*10, 3, 2, precision, alphaPrecision,
		oneDStep, rFunc, "fibonacci")
	xMin, yMin, err = hjs.Solve()
	if err != nil {
		fmt.Printf("error solving hooke jeeves: %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("hooke jeeves search fibonacci algorithm took : %v\n", timeEnd.Sub(timeStart))

	var nms many_dimension_search.NelderMeadSearch
	timeStart = time.Now()
	nms.Init([]float64{0, 0, 0}, 0.1, 3, precision, rFunc)
	xMin, yMin, err = nms.Solve()
	if err != nil {
		fmt.Printf("error solving nelder mead: %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("nelder mead search algorithm took : %v\n", timeEnd.Sub(timeStart))
}

func fourth() {
	var err error
	var timeStart, timeEnd time.Time
	var xMin []float64
	var yMin float64
	precision := 0.000001
	alphaPrecision := 0.000001
	oneDStepFib := 0.001
	precisionInterp := 0.001
	oneDStepInterp := 0.01
	maxIter := 5000

	rFunc := rozenbrokeOneDimensionFunction(100, 2, 45)
	gradFunctions := []func(xs []float64) float64{rozenbrokeFirstGrad(100, 2, 45),
		rozenbrokeSecondGrad(100, 2, 45), rozenbrokeThirdGrad(100, 2, 45)}
	hessian := hess(100, 2, 45)
	var fgd many_dimension_search.FastGradientDescendSearch
	var frs many_dimension_search.FletcherReevesSearch
	var dfps many_dimension_search.DavidonFletcherPowellSearch
	var lms many_dimension_search.LevenbergMarkkvadratSearch

	timeStart = time.Now()
	fgd.Init([]float64{0, 0, 0}, precision, precision, rFunc, gradFunctions, 3, precision, alphaPrecision, "break in two")
	xMin, yMin, err = fgd.Solve()
	if err != nil {
		fmt.Printf("error solving fast gradient descent: %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("fast gradient descent search break in two algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	fgd.Init([]float64{0, 0, 0}, precision, precision, rFunc, gradFunctions, 3, precision, alphaPrecision, "golden ratio")
	xMin, yMin, err = fgd.Solve()
	if err != nil {
		fmt.Printf("error solving fast gradient descent: %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("fast gradient descent golden ratio algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	fgd.Init([]float64{0, 0, 0}, precision, precision, rFunc, gradFunctions, 3, precision, alphaPrecision, "square interpolation")
	xMin, yMin, err = fgd.Solve()
	if err != nil {
		fmt.Printf("error solving fast gradient descent: %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("fast gradient descent square interpolation algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	fgd.Init([]float64{0, 0, 0}, precision, precision, rFunc, gradFunctions, 3, precision, alphaPrecision, "fibonacci")
	xMin, yMin, err = fgd.Solve()
	if err != nil {
		fmt.Printf("error solving fast gradient descent: %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("fast gradient descent fibonacci algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	frs.Init([]float64{0, 0, 0}, alphaPrecision, 3, precision, precision, alphaPrecision, oneDStepFib, maxIter, rFunc, gradFunctions, "break in two", true)
	xMin, yMin, err = frs.Solve()
	if err != nil {
		fmt.Printf("error solving pollak : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("pollak break in two algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	frs.Init([]float64{0, 0, 0}, 0.0001, 3, precision, precision, alphaPrecision, 0.01, maxIter, rFunc, gradFunctions, "golden ratio", true)
	xMin, yMin, err = frs.Solve()
	if err != nil {
		fmt.Printf("error solving pollak : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("pollak golden ratio algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	frs.Init([]float64{0, 0, 0}, alphaPrecision, 3, precision, precision, alphaPrecision, oneDStepFib, maxIter, rFunc, gradFunctions, "fibonacci", true)
	xMin, yMin, err = frs.Solve()
	if err != nil {
		fmt.Printf("error solving pollak : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("pollak fibonacci algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	frs.Init([]float64{0, 0, 0}, alphaPrecision, 3, precision, precision, alphaPrecision, oneDStepFib, maxIter, rFunc, gradFunctions, "break in two", false)
	xMin, yMin, err = frs.Solve()
	if err != nil {
		fmt.Printf("error solving fletcher reeves : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("fletcher reeves break in two algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	frs.Init([]float64{0, 0, 0}, 0.0001, 3, precision, precision, alphaPrecision, 0.01, maxIter, rFunc, gradFunctions, "golden ratio", false)
	xMin, yMin, err = frs.Solve()
	if err != nil {
		fmt.Printf("error solving fletcher reeves : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("fletcher reeves golden ratio algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	frs.Init([]float64{0, 0, 0}, alphaPrecision, 3, precision, precision, alphaPrecision, oneDStepFib, maxIter, rFunc, gradFunctions, "fibonacci", false)
	xMin, yMin, err = frs.Solve()
	if err != nil {
		fmt.Printf("error solving fletcher reeves : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("fletcher reeves fibonacci algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	dfps.Init([]float64{0, 0, 0}, alphaPrecision, 3, precisionInterp, precisionInterp, precisionInterp, oneDStepInterp, maxIter, rFunc, gradFunctions, "break in two")
	xMin, yMin, err = dfps.Solve()
	if err != nil {
		fmt.Printf("error solving davidon fletcher powell : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("davidon fletcher powell break in two algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	dfps.Init([]float64{0, 0, 0}, alphaPrecision, 3, precisionInterp, precisionInterp, precisionInterp, oneDStepInterp, maxIter, rFunc, gradFunctions, "golden ratio")
	xMin, yMin, err = dfps.Solve()
	if err != nil {
		fmt.Printf("error solving davidon fletcher powell : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("davidon fletcher powell golden ratio algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	dfps.Init([]float64{0, 0, 0}, alphaPrecision, 3, precisionInterp, precisionInterp, precisionInterp, oneDStepInterp, maxIter, rFunc, gradFunctions, "fibonacci")
	xMin, yMin, err = dfps.Solve()
	if err != nil {
		fmt.Printf("error solving davidon fletcher powell : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("davidon fletcher powell fibonacci algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	dfps.Init([]float64{0, 0, 0}, alphaPrecision, 3, precisionInterp, precisionInterp, precisionInterp, oneDStepInterp, maxIter, rFunc, gradFunctions, "square interpolation")
	xMin, yMin, err = dfps.Solve()
	if err != nil {
		fmt.Printf("error solving davidon fletcher powell : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("davidon fletcher powell square interpolation algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	lms.Init([]float64{0, 0, 0}, 3, rFunc, gradFunctions, hessian, 10000, 100000, 0.001)
	xMin, yMin, err = lms.Solve()
	if err != nil {
		fmt.Printf("error solving levenberg markkvadrat method : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("levenberg markkvadrat algorithm took : %v\n", timeEnd.Sub(timeStart))
}

func testLMS() {
	var err error
	var timeStart, timeEnd time.Time
	var xMin []float64
	var yMin float64

	rFunc := rozenbrokeOneDimensionFunction(158, 2, 40)
	gradFunctions := []func(xs []float64) float64{rozenbrokeFirstGrad(158, 2, 40),
		rozenbrokeSecondGrad(158, 2, 40), rozenbrokeThirdGrad(158, 2, 40)}
	hessian := hess(158, 2, 40)
	var lms many_dimension_search.LevenbergMarkkvadratSearch
	timeStart = time.Now()
	lms.Init([]float64{0, 0, 0}, 3, rFunc, gradFunctions, hessian, 10000, 100000, 0.001)
	xMin, yMin, err = lms.Solve()
	if err != nil {
		fmt.Printf("error solving levenberg markkvadrat method : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("levenberg markkvadrat algorithm took : %v\n", timeEnd.Sub(timeStart))
}

func testInverted() {
	var m la_methods.Matrix
	_ = m.InitWithPoints(2, 2, [][]float64{{1, 1}, {5, 2}})
	inv, err := m.Inverted()
	if err != nil {
		fmt.Printf("%v", err)
	}
	inv.Print()
}

func fifth() {
	var err error
	var timeStart, timeEnd time.Time
	var xMin []float64
	var yMin float64
	precision := 0.001

	rFunc := rozenbrokeOneDimensionFunction2()
	gradFunctions := []func(xs []float64) float64{rozenbrokeFirstGrad2(),
		rozenbroke2SecondGrad()}
	hessian := hess2()

	constraintExtFunc := constraintExtFunc(firstExtConstraintMod, secondExtConstraintMod, thirdExtConstraintMod)
	constraintExtGradFunctions := []func(xs []float64, r float64) float64{constraintExtFuncFirstGrad(),
		constraintExtFuncSecondGrad()}
	hessianExtConstraint := hessConstraintExt()

	constraintInt1Func := constraintInt1Func(firstConstraint, secondConstraint, thirdConstraint)
	constraintInt1GradFunctions := []func(xs []float64, r float64) float64{constraintInt1FuncFirstGrad(),
		constraintInt1FuncSecondGrad()}
	hessianInt1Constraint := hessConstraintInt1()

	constraintInt2Func := constraintInt2Func(firstConstraint, secondConstraint, thirdConstraint)
	constraintInt2GradFunctions := []func(xs []float64, r float64) float64{constraintInt2FuncFirstGrad(),
		constraintInt2FuncSecondGrad()}
	hessianInt2Constraint := hessConstraintInt2()

	constraintLagrangeFunc := constraintLagrangeFunc(firstLagrangeConstraintMod, secondLagrangeConstraintMod, thirdLagrangeConstraintMod)
	constraintLagrangeGradFunctions := []func(xs []float64, r float64, m []float64) float64{constraintLagrangeFuncFirstGrad(),
		constraintLagrangeFuncSecondGrad()}

	var ep constraint_methods.Penalty
	var pc constraint_methods.PenaltyCombined
	var pl constraint_methods.PenaltyLagrange
	//var gm constraint_methods.GradientMethod

	timeStart = time.Now()
	ep.Init([]float64{-4, -4}, 2, rFunc, []func(xs []float64) float64{firstExtConstraintMod, secondExtConstraintMod, thirdExtConstraintMod},
		gradFunctions, hessian, constraintExtGradFunctions, hessianExtConstraint, constraintExtFunc, precision, 1.618, "hooke jeeves")

	xMin, yMin, err = ep.Solve()
	if err != nil {
		fmt.Printf("error solving external penalty method : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("external penalty algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	ep.Init([]float64{-4, -4}, 2, rFunc, []func(xs []float64) float64{firstConstraint, secondConstraint, thirdConstraint},
		gradFunctions, hessian, constraintInt1GradFunctions, hessianInt1Constraint, constraintInt1Func, precision, float64(1)/float64(10), "hooke jeeves")

	xMin, yMin, err = ep.Solve()
	if err != nil {
		fmt.Printf("error solving internal penalty method : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("internal penalty algorithm took : %v\n", timeEnd.Sub(timeStart))

	//fmt.Println(constraintInt2Func([]float64{2, 2}, 1))

	timeStart = time.Now()
	ep.Init([]float64{0.5, 0.5}, 2, rFunc, []func(xs []float64) float64{firstConstraint, secondConstraint, thirdConstraint},
		gradFunctions, hessian, constraintInt2GradFunctions, hessianInt2Constraint, constraintInt2Func, precision, float64(1)/float64(10), "nelder mead")

	xMin, yMin, err = ep.Solve()
	if err != nil {
		fmt.Printf("error solving internal penalty method : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("internal penalty algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	pc.Init([]float64{-4, -4}, 2, rFunc,
		gradFunctions, constraintExtGradFunctions, constraintInt1GradFunctions, constraintExtFunc, constraintInt1Func, precision, float64(1)/float64(10), float64(1)/float64(10), "hooke jeeves")

	xMin, yMin, err = pc.Solve()
	if err != nil {
		fmt.Printf("error solving combined penalty method : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("combined penalty algorithm took : %v\n", timeEnd.Sub(timeStart))

	timeStart = time.Now()
	pl.Init([]float64{-5, -5}, 2, rFunc,
		gradFunctions, []func(xs []float64, r float64, m []float64) float64{firstLagrangeConstraintMod, secondLagrangeConstraintMod, thirdLagrangeConstraintMod}, constraintLagrangeGradFunctions, constraintLagrangeFunc, []float64{2, 2, 2}, precision, 1.6, "hooke jeeves")

	xMin, yMin, err = pl.Solve()
	if err != nil {
		fmt.Printf("error solving lagrange penalty method : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("lagrange penalty algorithm took : %v\n", timeEnd.Sub(timeStart))
}

func sixth() {
	var err error
	var timeStart, timeEnd time.Time
	var xMin []float64
	var yMin float64
	precision := 0.001

	rFunc := rozenbrokeOneDimensionFunction2()
	gradFunctions := []func(xs []float64) float64{rozenbrokeFirstGrad2(),
		rozenbroke2SecondGrad()}

	var gm constraint_methods.GradientMethod

	timeStart = time.Now()
	gm.Init([]float64{2, 1}, 2, rFunc,
		[]func(xs []float64) float64{firstConstraint, secondConstraint, thirdConstraint}, gradFunctions, A(3, 2), -10, precision, 30, "break in two")

	xMin, yMin, err = gm.Solve()
	if err != nil {
		fmt.Printf("error solving gradient penalty method : %v\n", err)
		return
	}
	timeEnd = time.Now()

	fmt.Printf("minimum: %f\n", yMin)
	for _, p := range xMin {
		fmt.Printf("minimum point: %f ", p)

	}
	fmt.Println()
	fmt.Printf("gradient penalty algorithm took : %v\n", timeEnd.Sub(timeStart))
}

func main() {
	sixth()

}
