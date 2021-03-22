package la_methods

import (
	"fmt"
)

type GaussMethod struct {
}

func (gm *GaussMethod) Solve(matrix Matrix, vec Vector) (Vector, error) {
	var system, systemTriang Matrix
	system = gm.makeSystem(matrix, vec)
	systemTriang = system.MakeTriangular()
	return gm.backwards(systemTriang)
}

func (gm *GaussMethod) makeSystem(matrix Matrix, vec Vector) Matrix {
	var newM Matrix
	newM.Init(matrix.DimensionRows, matrix.DimensionColumns+1)
	for i := 0; i < matrix.DimensionRows; i++ {
		for j := 0; j < matrix.DimensionColumns; j++ {
			newM.Points[i][j] = matrix.Points[i][j]
		}
		newM.Points[i][matrix.DimensionColumns] = vec.Points[i]
	}
	return newM
}

func (gm *GaussMethod) backwards(system Matrix) (Vector, error) {
	var answer, answerCopy, row Vector
	var err error
	answer.Init(system.DimensionRows)
	answerCopy.Init(system.DimensionRows)
	for i := system.DimensionRows - 1; i >= 0; i-- {
		answerCopy = answer.Copy()
		err = row.InitWithPoints(system.DimensionRows, system.Points[i][:system.DimensionRows])
		if err != nil {
			return Vector{}, fmt.Errorf("error during vector initializing: %v", err)
		}
		answerCopy, err = answerCopy.MulVFrom(row, i)
		if err != nil {
			return Vector{}, fmt.Errorf("error during vector multiplying: %v", err)
		}
		elementSum := answerCopy.ElemsSum(i)
		answer.Points[i] = (system.Points[i][system.DimensionRows] - elementSum) / system.Points[i][i]
	}
	return answer, nil
}
