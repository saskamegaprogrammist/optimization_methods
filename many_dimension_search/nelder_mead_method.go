package many_dimension_search

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"math"
)

type NelderMeadSearch struct {
	startPoint []float64
	s          float64
	precision  float64
	alpha      float64
	beta       float64
	gamma      float64
	m          float64
	teta       float64
	dimension  int
	targetFunc func(xs []float64) float64
}

func (nms *NelderMeadSearch) Init(startPoint []float64, s float64, dimension int,
	precision float64, targetFunc func(xs []float64) float64) {
	nms.startPoint = startPoint
	nms.precision = precision
	nms.targetFunc = targetFunc
	nms.dimension = dimension
	nms.s = s

	nms.alpha = 1
	nms.beta = 0.5
	nms.gamma = 2
	nms.m = 0.5
	nms.teta = 0.01
}

func getL1(s float64, n float64) float64 {
	return s * (math.Sqrt(n+1) + n - 1) / (n * math.Sqrt(2))
}

func getL2(s float64, n float64) float64 {
	return s * (math.Sqrt(n+1) - 1) / (n * math.Sqrt(2))
}

func (nms *NelderMeadSearch) getLVector(i int, s float64) la_methods.Vector {
	var lVector la_methods.Vector
	lPoints := make([]float64, nms.dimension)
	l1 := getL1(s, float64(nms.dimension))
	l2 := getL2(s, float64(nms.dimension))
	for j := 0; j < nms.dimension; j++ {
		if j != i {
			lPoints[j] = l2
		} else {
			lPoints[j] = l1
			lPoints[j] += s
		}
	}
	_ = lVector.InitWithPoints(nms.dimension, lPoints)

	return lVector
}

func (nms *NelderMeadSearch) Solve() ([]float64, float64, error) {
	var err error
	var k int
	var s float64
	var xStart la_methods.Vector
	var minVOld la_methods.Vector
	err = xStart.InitWithPoints(nms.dimension, nms.startPoint)
	if err != nil {
		return []float64{}, 0, fmt.Errorf("error initializing vector: %v", err)
	}
	minVOld = xStart.Copy()
	s = nms.s
	vectors := make([]la_methods.Vector, nms.dimension+1)
	for i := 0; i < nms.dimension; i++ {
		lVec := nms.getLVector(i, s)
		vectors[i], err = xStart.Add(lVec)
		if err != nil {
			return []float64{}, 0, fmt.Errorf("error adding vector: %v", err)
		}
	}
	vectors[nms.dimension] = xStart
	for {

		maxV, max, maxI, minV, min, minI := nms.findMaxAndMin(vectors)
		xC, _, xCMax, _, err := nms.findMassCenter(vectors, maxI)

		if err != nil {
			return []float64{}, 0, fmt.Errorf("error finding mass center: %v", err)
		}
		xF, err := nms.getTestPoint(xC, maxV)
		if err != nil {
			return []float64{}, 0, fmt.Errorf("error getting test point: %v", err)
		}
		fF := nms.targetFunc(xF.Points)

		var xN la_methods.Vector
		var fN float64
		if fF < min {
			xInter, err := xC.Sub(maxV)
			if err != nil {
				return []float64{}, 0, fmt.Errorf("error substracting vectors: %v", err)
			}
			xInter = xInter.MulOnValue(nms.gamma)
			xN, err = xC.Add(xInter)
			if err != nil {
				return []float64{}, 0, fmt.Errorf("error adding vectors: %v", err)
			}
			fN = nms.targetFunc(xN.Points)

			if fN < fF {
				maxV = xN
				vectors[maxI] = xN
				max = nms.targetFunc(xN.Points)
			} else {
				maxV = xF
				vectors[maxI] = xF
				max = nms.targetFunc(xF.Points)
			}
		} else if min <= fF && fF <= xCMax {
			maxV = xF
			vectors[maxI] = xF
			max = nms.targetFunc(xF.Points)
		} else if xCMax < fF && fF <= max {
			xInter, err := xF.Sub(xC)
			if err != nil {
				return []float64{}, 0, fmt.Errorf("error substracting vectors: %v", err)
			}
			xInter = xInter.MulOnValue(nms.beta)
			xN, err = xC.Add(xInter)
			if err != nil {
				return []float64{}, 0, fmt.Errorf("error adding vectors: %v", err)
			}
			fN = nms.targetFunc(xN.Points)

			if fN < fF {
				maxV = xN
				vectors[maxI] = xN
				max = nms.targetFunc(xN.Points)
			} else {
				vectors, err = nms.reduction(vectors, minV, min, minI)
				if err != nil {
					return []float64{}, 0, fmt.Errorf("error during reduction: %v", err)
				}
			}
		} else if fF > max {
			xInter, err := maxV.Sub(xC)
			if err != nil {
				return []float64{}, 0, fmt.Errorf("error substracting vectors: %v", err)
			}
			xInter = xInter.MulOnValue(nms.beta)
			xN, err = xC.Add(xInter)
			if err != nil {
				return []float64{}, 0, fmt.Errorf("error adding vectors: %v", err)
			}
			fN = nms.targetFunc(xN.Points)

			if fN < max {
				maxV = xN
				vectors[maxI] = xN
				max = nms.targetFunc(xN.Points)
			} else {
				vectors, err = nms.reduction(vectors, minV, min, minI)
				if err != nil {
					return []float64{}, 0, fmt.Errorf("error during reduction: %v", err)
				}
			}
		}
		stopFirst, err := nms.checkStopFirst(minV, minVOld)
		if err != nil {
			return []float64{}, 0, fmt.Errorf("error during checking first condition: %v", err)
		}
		if stopFirst && nms.checkStopSecond(vectors, min) {
			fmt.Printf("k value: %d\n", k)
			return minV.Points, nms.targetFunc(minV.Points), nil
		}
		if k%10 == 0 {
			tetaCur, err := nms.findTetaCurrent(vectors)
			if err != nil {
				return []float64{}, 0, fmt.Errorf("error during finding teta current: %v", err)
			}
			if tetaCur <= nms.teta {
				minIR := nms.findMin(vectors, minI)
				xInterI, err := minV.Sub(minIR)
				if err != nil {
					return []float64{}, 0, fmt.Errorf("error substracting vectors: %v", err)
				}
				xStart = minV
				s = xInterI.Len()

				for i := 0; i < nms.dimension; i++ {
					lVec := nms.getLVector(i, s)
					vectors[i], err = xStart.Add(lVec)
					if err != nil {
						return []float64{}, 0, fmt.Errorf("error adding vector: %v", err)
					}
				}
			}
		}
		k++
		minVOld = minV
	}
}

func (nms *NelderMeadSearch) findMin(vectors []la_methods.Vector, minI int) la_methods.Vector {
	var minINew int
	var minNew float64
	var testI int
	if minI == 0 {
		testI = 1
	}
	minNew = nms.targetFunc(vectors[testI].Points)
	minINew = testI
	for i, vec := range vectors[1:] {
		if (i + 1) != minI {
			v := nms.targetFunc(vec.Points)
			if v < minNew {
				minNew = v
			}
			minINew = i + 1
		}
	}
	return vectors[minINew]
}

func (nms *NelderMeadSearch) findTetaCurrent(vectors []la_methods.Vector) (float64, error) {
	min := math.Pi * float64(2)
	for i := 0; i < nms.dimension; i++ {
		for j := 0; j < nms.dimension; j++ {
			for k := 0; k < nms.dimension; k++ {
				var v float64
				if i != k && j != k && j != i {
					xInterI, err := vectors[i].Sub(vectors[k])
					if err != nil {
						return 0, fmt.Errorf("error substracting vectors: %v", err)
					}
					xInterJ, err := vectors[j].Sub(vectors[k])
					if err != nil {
						return 0, fmt.Errorf("error substracting vectors: %v", err)
					}
					mulF, err := xInterI.Mul(xInterJ)
					if err != nil {
						return 0, fmt.Errorf("error multiplicating vectors: %v", err)
					}
					mulS := xInterJ.Len() * xInterI.Len()
					v = math.Acos(math.Abs(mulF / mulS))
					if v < min {
						min = v
					}
				}
			}
		}
	}
	return min, nil
}

func (nms *NelderMeadSearch) checkStopFirst(minV la_methods.Vector, minVOld la_methods.Vector) (bool, error) {
	xInter, err := minV.Sub(minVOld)
	if err != nil {
		return false, fmt.Errorf("error substracting vectors: %v", err)
	}
	return xInter.Len() <= nms.precision, nil
}

func (nms *NelderMeadSearch) checkStopSecond(vectors []la_methods.Vector, min float64) bool {
	var sum float64
	for _, vec := range vectors {
		sum += math.Pow(nms.targetFunc(vec.Points)-min, 2)
	}
	sum = math.Sqrt(sum)
	sum /= float64(nms.dimension + 1)
	return sum <= nms.precision
}

func (nms *NelderMeadSearch) reduction(vectors []la_methods.Vector, minV la_methods.Vector, min float64, minI int) ([]la_methods.Vector, error) {
	for i, vec := range vectors {
		if i != minI {
			xInter, err := vec.Sub(minV)
			if err != nil {
				return vectors, fmt.Errorf("error substracting vectors: %v", err)
			}
			xInter = xInter.MulOnValue(nms.m)
			xInter, err = minV.Add(xInter)
			if err != nil {
				return vectors, fmt.Errorf("error adding vectors: %v", err)
			}
			vectors[i] = xInter
		}
	}
	return vectors, nil
}

func (nms *NelderMeadSearch) getTestPoint(xC la_methods.Vector, xH la_methods.Vector) (la_methods.Vector, error) {
	var err error
	xInter, err := xC.Sub(xH)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error substracting vectors: %v", err)
	}
	xInter = xInter.MulOnValue(nms.alpha)
	xInter, err = xInter.Add(xC)
	if err != nil {
		return la_methods.Vector{}, fmt.Errorf("error adding vectors: %v", err)
	}
	return xInter, nil
}

func (nms *NelderMeadSearch) findMassCenter(vectors []la_methods.Vector, maxI int) (la_methods.Vector, la_methods.Vector, float64, int, error) {
	var sum la_methods.Vector
	var err error
	var maxCI int
	var max float64
	if maxI == 0 && nms.dimension > 1 {
		maxCI = 1
		max = nms.targetFunc(vectors[1].Points)
	} else {
		max = nms.targetFunc(vectors[0].Points)
	}
	sum.Init(nms.dimension)

	for i, vec := range vectors {
		if i != maxI {
			sum, err = sum.Add(vec)
			if err != nil {
				return la_methods.Vector{}, vectors[maxCI], nms.targetFunc(vectors[maxCI].Points), maxCI, fmt.Errorf("error adding vectors: %v", err)
			}
			fInter := nms.targetFunc(vec.Points)
			if fInter > max {
				max = fInter
				maxCI = i
			}
		}
	}
	sum = sum.MulOnValue(float64(1) / float64(nms.dimension))
	return sum, vectors[maxCI], nms.targetFunc(vectors[maxCI].Points), maxCI, err
}

func (nms *NelderMeadSearch) findMaxAndMin(vectors []la_methods.Vector) (la_methods.Vector, float64, int, la_methods.Vector, float64, int) {
	var maxI, minI int
	var max, min float64
	max = nms.targetFunc(vectors[0].Points)
	min = max
	for i, vec := range vectors[1:] {
		f := nms.targetFunc(vec.Points)
		if f > max {
			max = f
			maxI = i + 1
		} else if f < min {
			min = f
			minI = i + 1
		}
	}
	return vectors[maxI], max, maxI, vectors[minI], min, minI
}
