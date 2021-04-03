package many_dimension_search

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
)

type LevenbergMarkkvadratSearch struct {
	startPoint    []float64
	dimension     int
	targetFunc    func(xs []float64) float64
	gradient      []func(xs []float64) float64
	hessian       func(xs []float64) la_methods.Matrix
	m             float64
	eps           float64
	maxIterations int
}

func (lms *LevenbergMarkkvadratSearch) Init(startPoint []float64, dimension int,
	targetFunc func(xs []float64) float64, gradient []func(xs []float64) float64,
	hessian func(xs []float64) la_methods.Matrix, m float64, maxIterations int, eps float64) {
	lms.startPoint = startPoint
	lms.targetFunc = targetFunc
	lms.dimension = dimension
	lms.gradient = gradient
	lms.m = m
	lms.hessian = hessian
	lms.maxIterations = maxIterations
	lms.eps = eps
}

func (lms *LevenbergMarkkvadratSearch) Solve() ([]float64, float64, error) {
	var err error
	var Hess, mM, HessInter, HessInterInv la_methods.Matrix
	var x, xOld la_methods.Vector
	var k int
	var m float64
	var grad, d la_methods.Vector
	m = lms.m
	err = x.InitWithPoints(lms.dimension, lms.startPoint)
	if err != nil {
		return nil, 0, fmt.Errorf("error during vector initializing: %v", err)
	}
	xOld = x.Copy()
OUTER:
	for {
		grad, err = lms.calculateGradient(x)
		if err != nil {
			return nil, 0, fmt.Errorf("error calculcating gradient: %v", err)
		}
		if grad.Len() < lms.eps {
			//fmt.Printf("k value: %d\n", k)
			return x.Points, lms.targetFunc(x.Points), nil
		}
		Hess = lms.hessian(x.Points)
		for {
			if k >= lms.maxIterations {
				//fmt.Printf("k value: %d\n", k)
				return x.Points, lms.targetFunc(x.Points), nil
			}
			mM.Init(lms.dimension, lms.dimension)
			mM.E()
			mM = mM.MulVal(m)
			HessInter, err = Hess.AddM(mM)
			if err != nil {
				return nil, 0, fmt.Errorf("error adding matrices: %v", err)
			}
			HessInterInv, err = HessInter.Inverted()
			if err != nil {
				return nil, 0, fmt.Errorf("error inverting matrix: %v", err)
			}

			HessInterInv.Print()

			//gradMinus = grad.MulOnValue(-1)
			d, err = HessInterInv.MulV(grad)
			if err != nil {
				return nil, 0, fmt.Errorf("error during vector and matrix multiplying: %v", err)
			}
			xOld = x
			x, err = x.Sub(d)
			if err != nil {
				return nil, 0, fmt.Errorf("error during vector adding: %v", err)
			}
			//fmt.Println(lms.targetFunc(x.Points))
			if lms.targetFunc(x.Points) < lms.targetFunc(xOld.Points) {
				m /= 2
				k++
				continue OUTER
			} else {
				m *= 2
				k++
			}
		}
	}
}

func (lms *LevenbergMarkkvadratSearch) calculateG(gradNew la_methods.Vector, gradOld la_methods.Vector,
	deltaX la_methods.Vector, G la_methods.Matrix) (la_methods.Matrix, error) {
	var err error
	var deltaXdeltaGrad, ggF float64
	var deltaXMRow, deltaXMCol, A, B, deltaX2, ggCol, ggRow, ggM, deltaG, GNew la_methods.Matrix
	var gradDelta, gg la_methods.Vector
	gradDelta, err = gradNew.Sub(gradOld)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during vector substracting: %v", err)
	}
	err = deltaXMCol.InitWithVectorColumn(lms.dimension, lms.dimension, deltaX)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error during vector initalizing: %v", err)
	}
	err = deltaXMRow.InitWithVectorRow(lms.dimension, lms.dimension, deltaX)
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
	err = ggCol.InitWithVectorColumn(lms.dimension, lms.dimension, gg)
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

func (lms *LevenbergMarkkvadratSearch) calculateGradient(x la_methods.Vector) (la_methods.Vector, error) {
	var gradPoints []float64
	var grad la_methods.Vector
	for i := 0; i < lms.dimension; i++ {
		gradPoints = append(gradPoints, lms.gradient[i](x.Points))
	}
	err := grad.InitWithPoints(lms.dimension, gradPoints)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error during vector initializing: %v", err)
	}
	return grad, nil
}
