package la_methods

import (
	"fmt"
)

type ZeidelMethod struct {
	precision float64
}

func (zm *ZeidelMethod) Init(precision float64) {
	zm.precision = precision
}

func (zm *ZeidelMethod) checkConvergence(matrix Matrix) bool {
	return matrix.CheckDiagonallyDominant()
}

func (zm *ZeidelMethod) Solve(matrix Matrix, vec Vector) (Vector, error) {
	//if !zm.checkConvergence(matrix) {
	//	return Vector{}, fmt.Errorf("method doens't converge")
	//}
	return zm.solve(matrix, vec)
}

func (zm *ZeidelMethod) solve(matrix Matrix, vec Vector) (Vector, error) {
	var dimension, k int
	var err error
	var stop bool
	var initial, nextIter Vector
	if matrix.DimensionRows != matrix.DimensionColumns {
		return Vector{}, fmt.Errorf("matrix dimensions don't match")
	}
	dimension = matrix.DimensionRows
	if dimension != vec.Dimension {
		return Vector{}, fmt.Errorf("matrix and vector dimensions don't match")
	}
	initial, err = zm.setInitial(matrix, vec)
	if err != nil {
		return Vector{}, fmt.Errorf("error during vector initializing: %v", err)
	}
	nextIter.Init(dimension)
	for {
		k++
		//nextIter.Print()
		for i := 0; i < dimension; i++ {
			numerator := vec.Points[i]
			for j := 0; j < i; j++ {
				numerator -= matrix.Points[i][j] * nextIter.Points[j]
			}
			for j := i + 1; j < dimension; j++ {
				numerator -= matrix.Points[i][j] * initial.Points[j]
			}
			nextIter.Points[i] = numerator / matrix.Points[i][i]
		}
		stop, err = zm.checkStop(initial, nextIter)
		if err != nil {
			return Vector{}, fmt.Errorf("error during checking condition: %v", err)
		}
		if stop {
			break
		}
		initial = nextIter.Copy()
	}
	return nextIter, nil
}

func (zm *ZeidelMethod) checkStop(initial Vector, nextIter Vector) (bool, error) {
	diff, err := nextIter.Sub(initial)
	if err != nil {
		return false, fmt.Errorf("error substractiong vectors: %v", err)
	}
	return diff.Len() < zm.precision, nil
}

func (zm *ZeidelMethod) setInitial(matrix Matrix, vec Vector) (Vector, error) {
	var initial Vector
	initial.Init(vec.Dimension)
	for i := 0; i < matrix.DimensionColumns; i++ {
		if matrix.Points[i][i] == 0 {
			return Vector{}, fmt.Errorf("main matrix has zero on diagonal")
		}
		initial.Points[i] = vec.Points[i] / matrix.Points[i][i]
	}
	return initial, nil
}
