package one_dimension_search

import (
	"fmt"
)

type FibonacciSearch struct {
	aStart     float64
	bStart     float64
	abLen      float64
	precision  float64
	diffConst  float64
	targetFunc func(x float64) float64
	k          int
	r          float64
	fibNumbers []float64
	fibN       float64
}

func (fs *FibonacciSearch) Init(aStart float64, bStart float64, precision float64, diffConst float64, targetFunc func(x float64) float64) {
	fs.aStart = aStart
	fs.bStart = bStart
	fs.abLen = bStart - aStart
	fs.precision = precision
	fs.diffConst = diffConst
	fs.targetFunc = targetFunc
	fs.fibNumbers = make([]float64, 0)
}

func (fs *FibonacciSearch) Solve() (float64, float64, error) {
	var err error
	var k int
	var y, z float64
	var fY, fZ float64
	var a, b float64
	n, fibN := fs.GetFibonacciNumber()
	a = fs.aStart
	b = fs.bStart
	fibN2, err := fs.GetFibonacci(n - 2)
	if err != nil {
		return 0, 0, fmt.Errorf("error counting fibonacci number: %v", err)
	}
	fibN1, err := fs.GetFibonacci(n - 1)
	if err != nil {
		return 0, 0, fmt.Errorf("error counting fibonacci number: %v", err)
	}
	y = a + (b-a)*fibN2/fibN
	z = a + (b-a)*fibN1/fibN
	for {
		fY = fs.targetFunc(y)
		fZ = fs.targetFunc(z)
		if checkDecreasing(fZ, fY) {
			a = y
			y = z
			fibNk1, err := fs.GetFibonacci(n - k - 1)
			if err != nil {
				return 0, 0, fmt.Errorf("error counting fibonacci number: %v", err)
			}
			fibNk2, err := fs.GetFibonacci(n - k - 2)
			if err != nil {
				return 0, 0, fmt.Errorf("error counting fibonacci number: %v", err)
			}
			z = a + (b-a)*fibNk2/fibNk1
		} else {
			b = z
			z = y
			fibNk1, err := fs.GetFibonacci(n - k - 1)
			if err != nil {
				return 0, 0, fmt.Errorf("error counting fibonacci number: %v", err)
			}
			fibNk3, err := fs.GetFibonacci(n - k - 3)
			if err != nil {
				return 0, 0, fmt.Errorf("error counting fibonacci number: %v", err)
			}
			y = a + (b-a)*fibNk3/fibNk1
		}
		if k >= n-3 {
			y = z
			z = y + fs.diffConst
			fY = fs.targetFunc(y)
			fZ = fs.targetFunc(z)
			if checkDecreasing(fZ, fY) {
				a = y
			} else {
				b = z
			}
			//fmt.Printf("k value: %v\n", k)
			fs.k = k
			xMid := (a + b) / 2
			return xMid, fs.targetFunc(xMid), nil
		}
		k++
	}
}

func (fs *FibonacciSearch) GetFibonacci(i int) (float64, error) {
	if len(fs.fibNumbers) == 0 {
		return 0, fmt.Errorf("fibonacci numbers aren't yet computed")
	}
	if i < 0 {
		return 1, nil
	}
	return fs.fibNumbers[i], nil
}

func (fs *FibonacciSearch) GetFibonacciNumber() (int, float64) {
	var i int
	var fibN, fibN1, fibN2 float64
	stopCond := fs.abLen / fs.precision
	for {
		if i == 0 {
			fibN1 = 1
			fibN = 1
		} else if i == 1 {
			fibN2 = 1
			fibN1 = 1
			fibN = 1
		} else {
			fibN = fibN2 + fibN1
			fibN2 = fibN1
			fibN1 = fibN
		}
		fs.fibNumbers = append(fs.fibNumbers, fibN)
		if fibN >= stopCond {
			fs.fibN = fibN
			return i, fibN
		}
		i++
	}
}

func (fs *FibonacciSearch) CountConvergence() (float64, error) {
	if fs.k != 0 {
		return float64(1) / fs.fibN, nil
	}
	return 0, fmt.Errorf("the algorithm hasn't been run")
}
