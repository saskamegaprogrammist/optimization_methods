package simplex_methods

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"gonum.org/v1/gonum/mat"
	"sort"
)

const FIRST = 1
const SECOND = 2

type SimplexMethod struct {
	n           int
	m           int
	fr          int
	constraints [][]float64
	f           []float64
	firstPhase  int
	eps         float64
}

func (sm *SimplexMethod) Init(n int, m int, constraints [][]float64, f []float64, firstPhase int) error {
	if len(constraints) != m {
		return fmt.Errorf("wrong constraints dimension:%d != %d", len(constraints), m)
	}
	sm.m = m
	sm.constraints = constraints
	if len(f) != n {
		return fmt.Errorf("wrong f dimension:%d != %d", len(f), n)
	}
	sm.n = n
	sm.f = f
	sm.fr = n - m
	if firstPhase != FIRST && firstPhase != SECOND {
		return fmt.Errorf("wrong phase const")
	}
	sm.firstPhase = firstPhase

	sm.eps = -0.01
	return nil
}

func (sm *SimplexMethod) Solve() ([]float64, float64, error) {
	var baseIndexA, freeIndexA []int
	var bMatrix la_methods.Matrix
	var err error
	var system la_methods.Matrix
	if sm.firstPhase == FIRST {
		baseIndexA, freeIndexA, bMatrix, system, err = sm.firstPhaseGauss()
		if err != nil {
			return nil, 0, fmt.Errorf("error during first phase: %v", err)
		}
	} else if sm.firstPhase == SECOND {
		baseIndexA, freeIndexA, bMatrix, system, err = sm.firstPhaseSimplex()
		if err != nil {
			return nil, 0, fmt.Errorf("error during first phase: %v", err)
		}
	}
	_, _, b, min, err := simplex(sm.m, sm.fr, sm.n, bMatrix, freeIndexA, baseIndexA, sm.f, system, sm.eps)
	return b, min, err
}

func simplex(m int, fr int, n int, bMatrix la_methods.Matrix, freeIndexA []int, baseIndexA []int, f []float64, system la_methods.Matrix, eps float64) (la_methods.Matrix, []int, []float64, float64, error) {
	var err error
	var bPoints = make([][]float64, m)
	for i := 0; i < m; i++ {
		bPoints[i] = make([]float64, m)
	}
	var bInverted, cT, piMatrix la_methods.Matrix
	var piVector la_methods.Vector
	var bPointsFlat = make([]float64, m*m)
	var cTPoints = make([]float64, m)
	var ds = make([]float64, fr)
	var has, has2 bool
	var dMinI, baMinI int
	var aVecNew la_methods.Vector
	var bOld = make([]float64, m)
	var bOldVec, bNewVec la_methods.Vector

	if err != nil {
		return la_methods.Matrix{}, []int{}, nil, 0, fmt.Errorf("error initing matrix: %v", err)
	}

	for {
		//bMatrix.Print()
		//system.Print()
		for i := 0; i < m; i++ {
			for j := 0; j < m; j++ {
				bPointsFlat[i*m+j] = bMatrix.Points[i][j]
			}
		}

		bDense := mat.NewDense(m, m, bPointsFlat)
		bDenseInv := mat.NewDense(m, m, nil)
		_ = bDenseInv.Inverse(bDense)

		var bI [][]float64
		rows, cols := bDenseInv.Dims()
		for i := 0; i < rows; i++ {
			bI = append(bI, bDenseInv.RawRowView(i))
		}
		err = bInverted.InitWithPoints(rows, cols, bI)
		if err != nil {
			return la_methods.Matrix{}, []int{}, nil, 0, fmt.Errorf("error initing matrix: %v", err)
		}

		//fmt.Println(baseIndexA)
		//fmt.Println(freeIndexA)

		for i := 0; i < m; i++ {
			cTPoints[i] = f[baseIndexA[i]]
		}

		err = cT.InitWithPoints(1, m, [][]float64{cTPoints})
		if err != nil {
			return la_methods.Matrix{}, []int{}, nil, 0, fmt.Errorf("error initing matrix: %v", err)
		}

		for i := 0; i < m; i++ {
			bOld[i] = system.Points[i][n]
		}
		err = bOldVec.InitWithPoints(m, bOld)
		if err != nil {
			return la_methods.Matrix{}, []int{}, nil, 0, fmt.Errorf("error initing vector: %v", err)
		}
		bNewVec, err = bInverted.MulV(bOldVec)
		if err != nil {
			return la_methods.Matrix{}, []int{}, nil, 0, fmt.Errorf("error multiplying: %v", err)
		}

		piMatrix, err = cT.MulM(bInverted)
		if err != nil {
			return la_methods.Matrix{}, []int{}, nil, 0, fmt.Errorf("error multiplying matrix: %v", err)
		}
		err = piVector.InitWithPoints(m, piMatrix.Points[0])
		if err != nil {
			return la_methods.Matrix{}, []int{}, nil, 0, fmt.Errorf("error initing vector: %v", err)
		}
		var a = make([]float64, m)
		var aVec la_methods.Vector
		var d float64
		for i := 0; i < fr; i++ {
			for j := 0; j < m; j++ {
				a[j] = system.Points[j][freeIndexA[i]]
			}
			err = aVec.InitWithPoints(m, a)
			if err != nil {
				return la_methods.Matrix{}, []int{}, nil, 0, fmt.Errorf("error initing vector: %v", err)
			}
			d, err = piVector.Mul(aVec)
			if err != nil {
				return la_methods.Matrix{}, []int{}, nil, 0, fmt.Errorf("error multiplying vectors: %v", err)
			}
			d -= f[freeIndexA[i]]
			ds[i] = d
		}
		has, _, dMinI = findMinNegative(ds, eps)
		if !has {
			var val float64
			for i := 0; i < m; i++ {
				val += f[baseIndexA[i]] * bNewVec.Points[i]
			}

			var xVec la_methods.Vector
			xVec.Init(n)
			for i := 0; i < m; i++ {
				xVec.Points[baseIndexA[i]] = bNewVec.Points[i]
			}

			var newMatrix la_methods.Matrix
			var newMatrixPoints = make([][]float64, m)
			for i := 0; i < m; i++ {
				newMatrixPoints[i] = make([]float64, m)
				for j := 0; j < m; j++ {
					newMatrixPoints[i][j] = system.Points[i][baseIndexA[j]]
				}
			}
			err = newMatrix.InitWithPoints(m, m, newMatrixPoints)
			if err != nil {
				return la_methods.Matrix{}, []int{}, nil, 0, fmt.Errorf("error initing matrix:%v", err)
			}

			return newMatrix, baseIndexA, xVec.Points, val, nil
		}
		for j := 0; j < m; j++ {
			a[j] = system.Points[j][freeIndexA[dMinI]]
		}
		err = aVec.InitWithPoints(m, a)
		if err != nil {
			return la_methods.Matrix{}, []int{}, nil, 0, fmt.Errorf("error initing vector: %v", err)
		}
		aVecNew, err = bInverted.MulV(aVec)
		if err != nil {
			return la_methods.Matrix{}, []int{}, nil, 0, fmt.Errorf("error multiplying: %v", err)
		}

		has2, _, baMinI = findMin(bNewVec, aVecNew, m)
		if !has2 {
			return la_methods.Matrix{}, []int{}, nil, 0, fmt.Errorf("function is limitless")
		}

		for i := 0; i < m; i++ {
			bMatrix.Points[i][baMinI] = aVec.Points[i]
		}
		sw := freeIndexA[dMinI]
		freeIndexA[dMinI] = baseIndexA[baMinI]
		baseIndexA[baMinI] = sw
	}
}

func findMin(bVecNew la_methods.Vector, aVecNew la_methods.Vector, m int) (bool, float64, int) {
	var min = float64(100000000000)
	var has bool
	var minI int
	var val float64
	for i := 0; i < m; i++ {
		if aVecNew.Points[i] < 0 {
			continue
		}
		val = bVecNew.Points[i] / aVecNew.Points[i]
		if val < min {
			has = true
			min = val
			minI = i
		}
	}
	return has, min, minI
}

func findMinNegative(ds []float64, eps float64) (bool, float64, int) {
	var has bool
	var min float64 = eps
	var minI int
	for i, d := range ds {
		if d < eps {
			has = true
			if d < min {
				min = d
				minI = i
			}
		}
	}
	return has, min, minI
}

func (sm *SimplexMethod) firstPhaseGauss() ([]int, []int, la_methods.Matrix, la_methods.Matrix, error) {
	var matrix la_methods.Matrix
	var err error
	err = matrix.InitWithPoints(sm.m, sm.n+1, sm.constraints)
	if err != nil {
		return nil, nil, la_methods.Matrix{}, la_methods.Matrix{}, fmt.Errorf("error initializing matrix: %v", err)
	}
	matrix = matrix.MakeE()
	var freeIndexA = make([]int, sm.fr)
	var baseIndexA = make([]int, sm.m)
	for i := 0; i < sm.fr; i++ {
		freeIndexA[i] = i + sm.m
	}
	for i := 0; i < sm.m; i++ {
		baseIndexA[i] = i
	}
	var bMatrix la_methods.Matrix
	bMatrix.Init(sm.m, sm.m)
	bMatrix.E()
	return baseIndexA, freeIndexA, bMatrix, matrix, nil
}

func (sm *SimplexMethod) firstPhaseSimplex() ([]int, []int, la_methods.Matrix, la_methods.Matrix, error) {
	var matrix, bFirst la_methods.Matrix
	var baseAFirst, freeAFirst []int
	var mn = sm.m + sm.n
	var sMatrixPoints = make([][]float64, sm.m)
	for i := 0; i < sm.m; i++ {
		sMatrixPoints[i] = make([]float64, mn+1)
	}
	var err error
	matrix.Init(sm.m, sm.n+1)
	bFirst.Init(sm.m, sm.m)
	bFirst.E()
	for i := 0; i < sm.m; i++ {
		for j := 0; j < sm.n; j++ {
			sMatrixPoints[i][j] = sm.constraints[i][j]
			matrix.Points[i][j] = sm.constraints[i][j]
		}
		sMatrixPoints[i][mn] = sm.constraints[i][sm.n]
		matrix.Points[i][sm.n] = sm.constraints[i][sm.n]
	}
	for i := 0; i < sm.m; i++ {
		sMatrixPoints[i][i+sm.n] = 1
	}
	for i := 0; i < sm.n; i++ {
		freeAFirst = append(freeAFirst, i)
	}
	for i := sm.n; i < mn; i++ {
		baseAFirst = append(baseAFirst, i)
	}
	var system la_methods.Matrix
	err = system.InitWithPoints(sm.m, mn+1, sMatrixPoints)
	if err != nil {
		return nil, nil, la_methods.Matrix{}, la_methods.Matrix{}, fmt.Errorf("error initing matrix: %v", err)
	}
	var f = make([]float64, mn)
	for i := 0; i < sm.m; i++ {
		f[mn-1-i] = -1
	}
	var baseIndex, freeIndex []int
	var newMatrix la_methods.Matrix
	newMatrix, baseIndex, _, _, err = simplex(sm.m, sm.n, mn, bFirst, freeAFirst, baseAFirst, f, system, sm.eps)
	if err != nil {
		return nil, nil, la_methods.Matrix{}, la_methods.Matrix{}, fmt.Errorf("error initing matrix: %v", err)
	}
	var baseIndexHelp []int = make([]int, len(baseIndex))
	copy(baseIndexHelp, baseIndex)
	sort.Ints(baseIndexHelp)
	j := 0
	for i := 0; i < sm.n; i++ {
		if j < len(baseIndexHelp) && baseIndexHelp[j] == i {
			j++
			continue
		}
		freeIndex = append(freeIndex, i)
	}

	return baseIndex, freeIndex, newMatrix, matrix, nil
}
