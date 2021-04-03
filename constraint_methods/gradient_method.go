package constraint_methods

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"github.com/saskamegaprogrammist/optimization_methods/one_dimension_search"
	"gonum.org/v1/gonum/mat"
	"math"
)

type GradientMethod struct {
	startPoint     []float64
	dimension      int
	targetFunc     func(xs []float64) float64
	penalties      []func(xs []float64) float64
	gradient       []func(xs []float64) float64
	A              func(xs []float64) la_methods.Matrix
	eps2           float64
	eps1           float64
	maxIter        int
	method         string
	svennAlgorithm one_dimension_search.Svenn
}

func (gm *GradientMethod) Init(startPoint []float64, dimension int,
	targetFunc func(xs []float64) float64, penalties []func(xs []float64) float64,
	gradient []func(xs []float64) float64, A func(xs []float64) la_methods.Matrix,
	eps1 float64, eps2 float64, M int, method string) {
	gm.startPoint = startPoint
	gm.targetFunc = targetFunc
	gm.dimension = dimension
	gm.penalties = penalties
	gm.gradient = gradient
	gm.A = A
	gm.maxIter = M
	gm.eps2 = eps2
	gm.eps1 = eps1
	gm.method = method
}

func (gm *GradientMethod) Solve() ([]float64, float64, error) {
	var err error
	var x, deltaX, gradV, lambda la_methods.Vector
	var A, AClean la_methods.Matrix
	var k, minIndex int
	var stop, has, zero bool
	var alphMin float64
	var val, min float64
	var grad []float64
	var excluded []int
	var grVal float64
	var lastAlpha float64
	var hasExcl, fib bool
	var search one_dimension_search.OneDimensionSearchI
	if gm.method == "break in two" {
		search = &one_dimension_search.BreakInTwoSearch{}
	} else if gm.method == "golden ratio" {
		search = &one_dimension_search.GoldenRatioSearch{}
	} else if gm.method == "fibonacci" {
		fib = true
	}
	err = x.InitWithPoints(gm.dimension, gm.startPoint)
	if err != nil {
		return nil, 0, fmt.Errorf("error during vector initializing: %v", err)
	}
	lastAlpha = 0.1
	for {
		AClean, err = gm.getA(x, []int{})
		if err != nil {
			return nil, 0, fmt.Errorf("error getting A clean: %v", err)
		}
		A, err = gm.getA(x, excluded)
		if err != nil {
			return nil, 0, fmt.Errorf("error getting A: %v", err)
		}
		if k >= gm.maxIter {
			fmt.Printf("k value: %d\n", k)
			return x.Points, gm.targetFunc(x.Points), nil
			//goto NINE
		}
		has = false
		for _, g := range gm.penalties {
			val = g(x.Points)
			//fmt.Println(val)
			if val <= 0 && val >= gm.eps1 {
				has = true
				break
			}
		}
		if !has {
			x, err = gm.getNewX(x, AClean)
			if err != nil {
				return nil, 0, fmt.Errorf("error getting new x: %v", err)
			}
		}
		zero = true
		grad = []float64{}
		for _, grF := range gm.gradient {
			grVal = grF(x.Points)
			grad = append(grad, grVal)
			zero = zero && (grVal == 0)
		}
		err = gradV.InitWithPoints(gm.dimension, grad)
		if err != nil {
			return nil, 0, fmt.Errorf("error initializing vector: %v", err)
		}
		if zero {
			if k == 0 && !(gm.penalties[0](x.Points) <= 0 &&
				gm.penalties[1](x.Points) <= 0 && gm.penalties[2](x.Points) <= 0) {
				return nil, 0, fmt.Errorf("select another starting point")
			}
			goto NINE
		}
		deltaX, err = gm.getDeltaX(A, gradV)
		if err != nil {
			return nil, 0, fmt.Errorf("error getting delta x: %v", err)
		}
		if !(deltaX.Len() <= gm.eps2) {
			goto TEN
		}
	NINE:
		lambda, err = gm.getLambda(A, gradV)
		if err != nil {
			return nil, 0, fmt.Errorf("error getting lambda: %v", err)
		}
		stop = true
		minIndex = 0
		min = lambda.Points[0]
		for i, l := range lambda.Points {
			stop = stop && (l >= 0)
			if l < min {
				minIndex = i
				min = l
			}
		}
		if stop {
			fmt.Printf("k value: %d\n", k)
			return x.Points, gm.targetFunc(x.Points), nil
		} else {
			hasExcl = gm.excludeConstraints(minIndex, &excluded)
		}
	TEN:
		tF := gm.getOneDimensionFunc(x, deltaX, gm.targetFunc)
		if fib {
			alphMin, err = gm.oneDimensionFibonacciSearch(lastAlpha, tF)
		} else if search != nil {
			alphMin, err = gm.oneDimensionSearch(lastAlpha, search, tF)
		} else {
			return nil, 0, fmt.Errorf("wrong search method")
		}
		if err != nil {
			return nil, 0, fmt.Errorf("error finding minimum: %v", err)
		}
		//fmt.Println(alphMin)
		if !hasExcl {
			minIndex = len(lambda.Points)
		}
		alphasPass := gm.findAlphaPass(excluded, minIndex, x, deltaX)
		//fmt.Println(excluded)
		//fmt.Println(alphasPass)

		alphasPass = append(alphasPass, alphMin)
		minAlpha := gm.findMinAlpha(lastAlpha, alphasPass)
		//fmt.Println(minAlpha)
		lastAlpha = minAlpha
		k++
		aplhaM := deltaX.MulOnValue(minAlpha)
		sumXAlphM, err := x.Add(aplhaM)
		if err != nil {
			return nil, 0, fmt.Errorf("error adding vectors: %v", err)
		}
		err = x.InitWithPoints(gm.dimension, sumXAlphM.Points)
		//fmt.Println(gm.targetFunc(sumXAlphM.Points))
		if err != nil {
			return nil, 0, fmt.Errorf("error during vector initializing: %v", err)
		}
	}
}
func (gm *GradientMethod) findMinAlpha(currMin float64, alphas []float64) float64 {
	var minimum = currMin
	for _, alpha := range alphas {
		if alpha <= -100 {
			continue
		}
		if alpha < minimum {
			minimum = alpha
		}
	}
	return minimum
}

func (gm *GradientMethod) findAlphaPass(excluded []int, lastExcluded int, x la_methods.Vector, deltaX la_methods.Vector) []float64 {
	var alphas []float64
	var alph float64
	for _, e := range excluded {
		if e == lastExcluded {
			continue
		}
		switch e {
		case 0:
			alph = firstConstraintLambda(gm.getOneDimensionFunc(x, deltaX, gm.penalties[0]))
			if !(math.Abs(gm.getOneDimensionFunc(x, deltaX, gm.penalties[0])(alph)) < 0.000001) {
				continue
			}
		case 1:
			alph = secondConstraintLambda(x, deltaX)
		case 2:
			alph = thirdConstraintLambda(x, deltaX)
		}
		alphas = append(alphas, alph)
	}
	return alphas
}

func firstConstraintLambda(
	g func(alpha float64) float64) float64 {
	return findBinaryAlpha(g, -math.Pow(10, 20), math.Pow(10, 20))
}

func findBinaryAlpha(g func(alpha float64) float64, start float64, end float64) float64 {
	var val, point float64
	var k int
	max := 100
	for {
		point = (math.Abs(end) - math.Abs(start)) / 2
		val = g(point)
		if k >= max {
			return point
		}
		if math.Abs(val) < 0.000001 {
		} else if val >= 0.000001 {
			end = point
		} else {
			start = point
		}
		k++
	}
}

func secondConstraintLambda(x la_methods.Vector, deltaX la_methods.Vector) float64 {
	return -x.Points[0] / deltaX.Points[0]
}

func thirdConstraintLambda(x la_methods.Vector, deltaX la_methods.Vector) float64 {
	return -x.Points[1] / deltaX.Points[1]
}

func (gm *GradientMethod) findBounds(alpha float64, targetFunc func(x float64) float64) (float64, float64, error) {
	gm.svennAlgorithm.Init(0.0001, alpha, targetFunc)
	return gm.svennAlgorithm.Solve()
}

func (gm *GradientMethod) oneDimensionFibonacciSearch(alpha float64, tF func(alpha float64) float64) (float64, error) {
	var fib one_dimension_search.FibonacciSearch
	a, b, err := gm.findBounds(alpha, tF)
	if err != nil {
		return 0, fmt.Errorf("error during svenn algorithm: %v", err)
	}
	fib.Init(a, b, 0.00001, 0.00001, tF)
	min, _, err := fib.Solve()
	if err != nil {
		return 0, fmt.Errorf("error during fibonacci search: %v", err)
	}
	return min, nil
}

func (gm *GradientMethod) oneDimensionSearch(
	alpha float64, search one_dimension_search.OneDimensionSearchI, tF func(alpha float64) float64) (float64, error) {
	a, b, err := gm.findBounds(alpha, tF)
	if err != nil {
		return 0, fmt.Errorf("error during svenn algorithm: %v", err)
	}
	search.Init(a, b, 0.00001, tF)
	min, _ := search.Solve()
	return min, nil
}

func (gm *GradientMethod) excludeConstraints(minIndex int, excluded *[]int) bool {
	for _, v := range *excluded {
		if minIndex == v {
			return true
		}
	}
	*excluded = append(*excluded, minIndex)
	return false
}

func (gm *GradientMethod) getOneDimensionFunc(x la_methods.Vector, deltaX la_methods.Vector,
	tF func(xs []float64) float64) func(alpha float64) float64 {
	return func(alpha float64) float64 {
		aplhaM := deltaX.MulOnValue(alpha)
		sumXAlphM, _ := x.Add(aplhaM)
		return tF(sumXAlphM.Points)
	}
}

func (gm *GradientMethod) getA(x la_methods.Vector, excluded []int) (la_methods.Matrix, error) {
	var A, AExcl la_methods.Matrix
	var err error
	A = gm.A(x.Points)
	AExcl, err = A.RemoveRows(excluded)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during rows removing: %v", err)
	}
	return AExcl, nil
}

func (gm *GradientMethod) getNewX(x la_methods.Vector, A la_methods.Matrix) (la_methods.Vector, error) {
	var Atransp, Amul, AInv, AMul2 la_methods.Matrix
	var tVec, tMulVec, xNewVec la_methods.Vector
	var t []float64
	var err error
	Atransp, err = A.Transponate()
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error transponating matrix: %v", err)
	}
	Amul, err = A.MulM(Atransp)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error multiplying matrix: %v", err)
	}

	var a []float64
	for _, p := range Amul.Points {
		a = append(a, p...)
	}

	aMulN := mat.NewDense(Amul.DimensionRows, Amul.DimensionColumns, a)
	aMulNInverted := mat.NewDense(Amul.DimensionRows, Amul.DimensionColumns, nil)
	_ = aMulNInverted.Inverse(aMulN)
	//if err != nil {
	//	return la_methods.Vector{}, fmt.Errorf("error inverting matrix: %v", err)
	//}
	var aI [][]float64
	rows, cols := aMulNInverted.Dims()
	for i := 0; i < rows; i++ {
		aI = append(aI, aMulNInverted.RawRowView(i))
	}
	err = AInv.InitWithPoints(rows, cols, aI)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error initing matrix: %v", err)
	}

	AMul2, err = Atransp.MulM(AInv)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error multiplying matrix: %v", err)
	}
	for _, g := range gm.penalties {
		t = append(t, g(x.Points))
	}
	err = tVec.InitWithPoints(len(gm.penalties), t)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error initializing vector: %v", err)
	}
	tMulVec, err = AMul2.MulV(tVec)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error multiplying matrix and vector: %v", err)
	}
	xNewVec, err = x.Add(tMulVec)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error adding vectors: %v", err)
	}
	return xNewVec, nil
}

func (gm *GradientMethod) getDeltaX(A la_methods.Matrix, gradient la_methods.Vector) (la_methods.Vector, error) {
	var Atransp, Amul, AInv, AMul2, AMul3, AE, AESub la_methods.Matrix
	var xDelta la_methods.Vector
	var err error
	Atransp, err = A.Transponate()
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error transponating matrix: %v", err)
	}
	Amul, err = A.MulM(Atransp)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error multiplying matrix: %v", err)
	}
	var a []float64
	for _, p := range Amul.Points {
		a = append(a, p...)
	}
	aMulN := mat.NewDense(Amul.DimensionRows, Amul.DimensionColumns, a)
	aMulNInverted := mat.NewDense(Amul.DimensionRows, Amul.DimensionColumns, nil)
	_ = aMulNInverted.Inverse(aMulN)
	//if err != nil {
	//	return la_methods.Vector{}, fmt.Errorf("error inverting matrix: %v", err)
	//}
	var aI [][]float64
	rows, cols := aMulNInverted.Dims()
	for i := 0; i < rows; i++ {
		aI = append(aI, aMulNInverted.RawRowView(i))
	}
	err = AInv.InitWithPoints(rows, cols, aI)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error initing matrix: %v", err)
	}
	AMul2, err = Atransp.MulM(AInv)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error multiplying matrix: %v", err)
	}
	AMul3, err = AMul2.MulM(A)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error multiplying matrix: %v", err)
	}
	AE.Init(AMul3.DimensionRows, AMul3.DimensionColumns)
	AE.E()
	AMul3 = AMul3.MulVal(-1)
	AESub, err = AE.AddM(AMul3)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error substracting matrices: %v", err)
	}
	AESub = AESub.MulVal(-1)
	xDelta, err = AESub.MulV(gradient)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error multiplying matrix and vector: %v", err)
	}
	return xDelta, nil
}

func (gm *GradientMethod) getLambda(A la_methods.Matrix, gradient la_methods.Vector) (la_methods.Vector, error) {
	var Atransp, Amul, AInv, AInvMinus, AMul1 la_methods.Matrix
	var lambda la_methods.Vector
	var err error
	Atransp, err = A.Transponate()
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error transponating matrix: %v", err)
	}
	Amul, err = A.MulM(Atransp)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error multiplying matrix: %v", err)
	}
	var a []float64
	for _, p := range Amul.Points {
		a = append(a, p...)
	}
	aMulN := mat.NewDense(Amul.DimensionRows, Amul.DimensionColumns, a)
	aMulNInverted := mat.NewDense(Amul.DimensionRows, Amul.DimensionColumns, nil)
	_ = aMulNInverted.Inverse(aMulN)
	//if err != nil {
	//	return la_methods.Vector{}, fmt.Errorf("error inverting matrix: %v", err)
	//}
	var aI [][]float64
	rows, cols := aMulNInverted.Dims()
	for i := 0; i < rows; i++ {
		aI = append(aI, aMulNInverted.RawRowView(i))
	}
	err = AInv.InitWithPoints(rows, cols, aI)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error initing matrix: %v", err)
	}
	AInvMinus = AInv.MulVal(-1)
	AMul1, err = AInvMinus.MulM(A)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error multiplying matrix: %v", err)
	}
	lambda, err = AMul1.MulV(gradient)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error multiplying matrix and vector: %v", err)
	}
	return lambda, nil
}
