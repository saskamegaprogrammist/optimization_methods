package simplex_methods

import (
	"fmt"
	"math"
)

type SimplexMethodReal struct {
	n           int
	m           int
	fr          int
	constraints [][]float64
	f           []float64
	firstPhase  int
	eps         float64
	bestF       float64
	bestPoints  []Point
}

func (smr *SimplexMethodReal) Init(n int, m int, constraints [][]float64, f []float64, firstPhase int) error {
	if len(constraints) != m {
		return fmt.Errorf("wrong constraints dimension:%d != %d", len(constraints), m)
	}
	smr.m = m
	smr.constraints = constraints
	if len(f) != n {
		return fmt.Errorf("wrong f dimension:%d != %d", len(f), n)
	}
	smr.n = n
	smr.f = f
	smr.fr = n - m
	if firstPhase != FIRST && firstPhase != SECOND {
		return fmt.Errorf("wrong phase const")
	}
	smr.firstPhase = firstPhase

	smr.eps = -0.01
	return nil
}

func (smr *SimplexMethodReal) SolveReal() ([]float64, float64, error) {
	var sm SimplexMethod
	var x []float64
	var fVal float64
	var err error
	err = sm.Init(smr.n, smr.m, smr.constraints, smr.f, smr.firstPhase)
	if err != nil {
		return nil, 0, fmt.Errorf("error initing simplex method: %v", err)
	}
	x, fVal, err = sm.Solve()
	var hasFr bool
	var optimizeIndex int
	var realPart float64
	for i := 0; i < sm.n; i++ {
		realPart = math.Floor(x[i])
		if x[i]-realPart != 0 {
			hasFr = true
			optimizeIndex = i
			//break
		}
	}
	if !hasFr {
		return x, fVal, nil
	}
	var n, m int
	m = sm.m + 1
	n = sm.n + 1
	var constraints = make([][]float64, m)
	var f []float64
	f = sm.f
	f = append(f, 0)
	for i := 0; i < m-1; i++ {
		constraints[i] = sm.constraints[i]
		constraints[i] = append(constraints[i], constraints[i][sm.n])
		constraints[i][sm.n] = 0
	}
	var J1Constraints [][]float64
	var J []PointHelp
	var firstNewConstraint = make([]float64, n+1)
	var secondNewConstraint = make([]float64, n+1)
	firstNewConstraint[optimizeIndex] = 1
	secondNewConstraint[optimizeIndex] = -1
	firstNewConstraint[n-1] = 1
	secondNewConstraint[n-1] = 1
	firstNewConstraint[n] = realPart
	secondNewConstraint[n] = -(realPart + 1)
	J1Constraints = append(J1Constraints, firstNewConstraint, secondNewConstraint)
	//fmt.Println("CONSTRAINT")
	//fmt.Println(firstNewConstraint)
	//fmt.Println(secondNewConstraint)

	point := Point{x: x, f: fVal}
	J = append(J, PointHelp{point: point, ind: 1}, PointHelp{point: point, ind: 2})
	var wasSetBestF bool
	for {
		if len(J) == 0 {
			break
		}
		var x1 []float64
		var f1 float64
		//var currPoint PointHelp = J[0]
		J = J[1:]
		constraints[m-1] = J1Constraints[0]
		J1Constraints = J1Constraints[1:]
		err = sm.Init(n, m, constraints, f, smr.firstPhase)
		if err != nil {
			return nil, 0, fmt.Errorf("error initing simplex method: %v", err)
		}
		x1, f1, err = sm.Solve()
		//fmt.Println(f1)
		//fmt.Println(x1)
		realPart = math.Floor(x1[optimizeIndex])
		if x1[optimizeIndex]-realPart != 0 {
			if wasSetBestF {
				if smr.bestF >= f1 {
					continue
				} else {
					//??
				}
			} else {
				var firstNewConstraint1 = make([]float64, n+1)
				var secondNewConstraint1 = make([]float64, n+1)
				copy(firstNewConstraint1, firstNewConstraint)
				copy(secondNewConstraint1, secondNewConstraint)
				firstNewConstraint1[n] = realPart
				secondNewConstraint1[n] = -(realPart + 1)
				J1Constraints = append(J1Constraints, firstNewConstraint1, secondNewConstraint1)
				//fmt.Println("CONSTRAINT")
				//fmt.Println(firstNewConstraint1)
				//fmt.Println(secondNewConstraint1)
				point := Point{x: x1, f: f1}
				J = append(J, PointHelp{point: point, ind: 1}, PointHelp{point: point, ind: 2})
			}
		} else {
			if wasSetBestF {
				if f1 > smr.bestF {
					smr.bestPoints = append(smr.bestPoints, Point{x: x1, f: f1})
				}
			} else {
				smr.bestF = f1
				smr.bestPoints = append(smr.bestPoints, Point{x: x1, f: f1})
				var newJ []PointHelp
				var newJConstraints [][]float64
				for i, j := range J {
					if j.point.f > f1 {
						newJ = append(newJ, j)
						newJConstraints = append(newJConstraints, J1Constraints[i])
					}
				}
				J1Constraints = newJConstraints
				J = newJ
			}
		}
	}
	if len(smr.bestPoints) == 0 {
		return x, fVal, nil
	}
	var max = smr.bestPoints[0].f
	var bP = smr.bestPoints[0]
	for i := 0; i < len(smr.bestPoints); i++ {
		if smr.bestPoints[i].f > max {
			bP = smr.bestPoints[i]
		}
	}
	return bP.x, bP.f, nil
}
