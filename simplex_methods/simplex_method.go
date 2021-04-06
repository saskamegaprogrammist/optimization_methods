package simplex_methods

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"gonum.org/v1/gonum/mat"
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
	var baseIndex, freeIndex int
	var err error
	var system la_methods.Matrix
	if sm.firstPhase == FIRST {
		baseIndex = 0
		freeIndex = sm.m
		system, err = sm.firstPhaseGauss()
		if err != nil {
			return nil, 0, fmt.Errorf("error during first phase: %v", err)
		}
	}
	var bPoints = make([][]float64, sm.m)
	for i := 0; i < sm.m; i++ {
		bPoints[i] = make([]float64, sm.m)
	}
	var bMatrix, bInverted, cT, piMatrix la_methods.Matrix
	var piVector la_methods.Vector
	var bPointsFlat = make([]float64, sm.m*sm.m)
	var cTPoints = make([]float64, sm.m)
	var ds = make([]float64, sm.fr)
	var has, has2 bool
	var dMinI, baMinI int
	//var dMin, baMin float64
	var aVecNew la_methods.Vector
	//var ePVec, beP la_methods.Vector
	//var aSub la_methods.Vector
	//var aSubM, bKM, epM la_methods.Matrix
	var bOld = make([]float64, sm.m)
	var bOldVec, bNewVec la_methods.Vector

	var f = sm.f

	bMatrix.Init(sm.m, sm.m)
	bMatrix.E()

	var freeIndexA = make([]int, sm.n-sm.m)
	var baseIndexA = make([]int, sm.m)
	for i := 0; i < sm.n-sm.m; i++ {
		freeIndexA[i] = i + freeIndex
	}
	for i := 0; i < sm.n-sm.m; i++ {
		baseIndexA[i] = i + baseIndex
	}

	for {
		for i := 0; i < sm.m; i++ {
			for j := 0; j < sm.m; j++ {
				bPointsFlat[i*sm.m+j] = bMatrix.Points[i][j]
			}
		}

		bDense := mat.NewDense(sm.m, sm.m, bPointsFlat)
		bDenseInv := mat.NewDense(sm.m, sm.m, nil)
		_ = bDenseInv.Inverse(bDense)

		var bI [][]float64
		rows, cols := bDenseInv.Dims()
		for i := 0; i < rows; i++ {
			bI = append(bI, bDenseInv.RawRowView(i))
		}
		err = bInverted.InitWithPoints(rows, cols, bI)
		if err != nil {
			return nil, 0, fmt.Errorf("error initing matrix: %v", err)
		}

		for i := 0; i < sm.m; i++ {
			cTPoints[i] = f[baseIndexA[i]]
		}

		err = cT.InitWithPoints(1, sm.m, [][]float64{cTPoints})
		if err != nil {
			return nil, 0, fmt.Errorf("error initing matrix: %v", err)
		}

		for i := 0; i < sm.m; i++ {
			bOld[i] = system.Points[i][sm.n]
		}
		err = bOldVec.InitWithPoints(sm.m, bOld)
		if err != nil {
			return nil, 0, fmt.Errorf("error initing vector: %v", err)
		}
		bNewVec, err = bInverted.MulV(bOldVec)
		if err != nil {
			return nil, 0, fmt.Errorf("error multiplying: %v", err)
		}
		bNewVec.Print()

		//for j := 0; j < sm.m; j++ {
		//	system.Points[j][sm.n] = bNewVec.Points[j]
		//}

		system.Print()

		piMatrix, err = cT.MulM(bInverted)
		if err != nil {
			return nil, 0, fmt.Errorf("error multiplying matrix: %v", err)
		}
		err = piVector.InitWithPoints(sm.m, piMatrix.Points[0])
		if err != nil {
			return nil, 0, fmt.Errorf("error initing vector: %v", err)
		}
		var a = make([]float64, sm.m)
		var aVec la_methods.Vector
		var d float64
		for i := 0; i < sm.fr; i++ {
			for j := 0; j < sm.m; j++ {
				a[j] = system.Points[j][freeIndexA[i]]
			}
			err = aVec.InitWithPoints(sm.m, a)
			if err != nil {
				return nil, 0, fmt.Errorf("error initing vector: %v", err)
			}
			d, err = piVector.Mul(aVec)
			if err != nil {
				return nil, 0, fmt.Errorf("error multiplying vectors: %v", err)
			}
			d -= sm.f[freeIndexA[i]]
			ds[i] = d
		}

		has, _, dMinI = sm.findMinNegative(ds)
		var val float64
		for i := 0; i < sm.m; i++ {
			val += f[baseIndexA[i]] * bNewVec.Points[i]
		}
		fmt.Println(val)
		if !has {
			var val float64
			for i := 0; i < sm.m; i++ {
				val += f[baseIndexA[i]] * bNewVec.Points[i]
			}
			fmt.Println(val)
			return nil, val, nil
		}
		for j := 0; j < sm.m; j++ {
			a[j] = system.Points[j][freeIndexA[dMinI]]
		}
		err = aVec.InitWithPoints(sm.m, a)
		if err != nil {
			return nil, 0, fmt.Errorf("error initing vector: %v", err)
		}
		aVecNew, err = bInverted.MulV(aVec)
		if err != nil {
			return nil, 0, fmt.Errorf("error multiplying: %v", err)
		}

		has2, _, baMinI = sm.findMin(bNewVec, aVecNew)
		if !has2 {
			return nil, 0, fmt.Errorf("function is limitless")
		}
		//var eP = make([]float64, sm.m)
		//eP[baMinI] = 1
		//err = ePVec.InitWithPoints(sm.m, eP)
		//if err != nil {
		//	return nil, 0, fmt.Errorf("error initing vector: %v", err)
		//}
		//beP, err = bMatrix.MulV(ePVec)
		//if err != nil {
		//	return nil, 0, fmt.Errorf("error initing vector: %v", err)
		//}
		//aSub, err = aVec.Sub(beP)
		//if err != nil {
		//	return nil, 0, fmt.Errorf("error substracting vectors: %v", err)
		//}
		//
		//err = aSubM.InitWithVectorRow(sm.m, sm.m, aSub)
		//if err != nil {
		//	return nil, 0, fmt.Errorf("error initing matrix: %v", err)
		//}
		//err = epM.InitWithVectorColumn(sm.m, sm.m, ePVec)
		//if err != nil {
		//	return nil, 0, fmt.Errorf("error initing matrix: %v", err)
		//}
		//bKM, err = aSubM.MulM(epM)
		//if err != nil {
		//	return nil, 0, fmt.Errorf("error multiplying: %v", err)
		//}
		//bMatrix, err = bMatrix.AddM(bKM)
		for i := 0; i < sm.m; i++ {
			bMatrix.Points[i][baMinI] = aVec.Points[i]
		}
		bMatrix.Print()
		fmt.Println(dMinI, baMinI)
		fmt.Println(freeIndexA[dMinI], baseIndexA[baMinI])
		sw := freeIndexA[dMinI]
		freeIndexA[dMinI] = baseIndexA[baMinI]
		baseIndexA[baMinI] = sw
		//
		//sw := f[baMinI]
		//f[baMinI] = f[dMinI+freeIndex]
		//f[dMinI+freeIndex] = sw

		fmt.Println(baseIndexA, freeIndexA)

	}
}

func (sm *SimplexMethod) findMin(bVecNew la_methods.Vector, aVecNew la_methods.Vector) (bool, float64, int) {
	var min = float64(1000000000001234)
	fmt.Print("test2")
	var has bool
	var minI int
	var val float64
	for i := 0; i < sm.m; i++ {
		if aVecNew.Points[i] < 0 {
			continue
		}
		val = bVecNew.Points[i] / aVecNew.Points[i]
		fmt.Println(val)
		if val < min {
			has = true
			min = val
			minI = i
		}
	}
	return has, min, minI
}

func (sm *SimplexMethod) findMinNegative(ds []float64) (bool, float64, int) {
	var has bool
	var min float64
	var minI int
	for i, d := range ds {
		if d < sm.eps {
			has = true
			if d < min {
				min = d
				minI = i
			}
		}
	}
	return has, min, minI
}

func (sm *SimplexMethod) firstPhaseGauss() (la_methods.Matrix, error) {
	var matrix la_methods.Matrix
	var err error
	err = matrix.InitWithPoints(sm.m, sm.n+1, sm.constraints)
	if err != nil {
		return la_methods.Matrix{}, fmt.Errorf("error initializing matrix: %v", err)
	}
	return matrix.MakeE(), nil
}
