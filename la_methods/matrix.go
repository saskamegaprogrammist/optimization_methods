package la_methods

import (
	"fmt"
	"math"
)

type Matrix struct {
	DimensionRows    int
	DimensionColumns int
	Points           [][]float64
}

func (m *Matrix) Init(dimensionRows int, dimensionColumns int) {
	m.DimensionRows = dimensionRows
	m.DimensionColumns = dimensionColumns
	m.Points = make([][]float64, dimensionRows)
	for i := range m.Points {
		m.Points[i] = make([]float64, dimensionColumns)
	}
}

func (m *Matrix) InitWithPoints(dimensionRows int, dimensionColumns int, points [][]float64) error {
	if len(points) != dimensionRows {
		return fmt.Errorf("dimension rows doesn't match points length: %d", len(points))
	}
	if len(points[0]) != dimensionColumns {
		return fmt.Errorf("dimension columns doesn't match points[0] length: %d", len(points))
	}
	m.DimensionRows = dimensionRows
	m.DimensionColumns = dimensionColumns
	m.Points = points
	return nil
}

func (m *Matrix) InitWithVectorRow(dimensionRows int, dimensionColumns int, vector Vector) error {
	if vector.Dimension != dimensionColumns {
		return fmt.Errorf("dimension vector doesn't match dimension columns: %d", vector.Dimension)
	}
	m.DimensionRows = dimensionRows
	m.DimensionColumns = dimensionColumns
	m.Points = make([][]float64, dimensionRows)
	for i := range m.Points {
		m.Points[i] = make([]float64, dimensionColumns)
	}
	for i, point := range vector.Points {
		m.Points[0][i] = point
	}
	return nil
}

func (m *Matrix) InitWithVectorColumn(dimensionRows int, dimensionColumns int, vector Vector) error {
	if vector.Dimension != dimensionRows {
		return fmt.Errorf("dimension vector doesn't match dimension rows: %d", vector.Dimension)
	}
	m.DimensionRows = dimensionRows
	m.DimensionColumns = dimensionColumns
	m.Points = make([][]float64, dimensionRows)
	for i := range m.Points {
		m.Points[i] = make([]float64, dimensionColumns)
	}
	for i, point := range vector.Points {
		m.Points[i][0] = point
	}
	return nil
}

func (m *Matrix) E() {
	for i := 0; i < m.DimensionRows; i++ {
		for j := 0; j < m.DimensionColumns; j++ {
			if i == j {
				m.Points[i][j] = 1
			} else {
				m.Points[i][j] = 0
			}
		}
	}
}

func (m *Matrix) MulV(vector Vector) (Vector, error) {
	var mulVector Vector
	var sum float64
	if m.DimensionColumns != vector.Dimension {
		return mulVector, fmt.Errorf("vector dimension doesn't match matrix columns dimension: %d", vector.Dimension)
	}
	mulVector.Init(m.DimensionRows)
	for i := 0; i < m.DimensionRows; i++ {
		sum = 0
		for j := 0; j < m.DimensionColumns; j++ {
			sum += m.Points[i][j] * vector.Points[j]
		}
		mulVector.Points[i] = sum
	}
	return mulVector, nil
}

func (m *Matrix) MulM(matrix Matrix) (Matrix, error) {
	var mulM Matrix
	var sum float64
	if m.DimensionColumns != matrix.DimensionRows {
		return mulM, fmt.Errorf("matrix dimensions doesn't match")
	}
	mulM.Init(m.DimensionRows, matrix.DimensionColumns)
	for i := 0; i < m.DimensionRows; i++ {
		for k := 0; k < matrix.DimensionColumns; k++ {
			sum = 0
			for j := 0; j < m.DimensionColumns; j++ {
				sum += m.Points[i][j] * matrix.Points[j][k]
			}
			mulM.Points[i][k] = sum
		}

	}
	return mulM, nil
}

func (m *Matrix) MulVal(val float64) Matrix {
	var mulM Matrix
	_ = mulM.InitWithPoints(m.DimensionRows, m.DimensionColumns, m.Points)

	for i := 0; i < m.DimensionRows; i++ {
		for j := 0; j < m.DimensionColumns; j++ {
			mulM.Points[i][j] *= val
		}
	}
	return mulM
}

func (m *Matrix) AddM(matrix Matrix) (Matrix, error) {
	var newM Matrix
	if m.DimensionColumns != matrix.DimensionColumns || m.DimensionRows != matrix.DimensionRows {
		return Matrix{}, fmt.Errorf("matrix dimensions doesn't match: %d, %d", matrix.DimensionRows, matrix.DimensionColumns)
	}
	_ = newM.InitWithPoints(m.DimensionRows, m.DimensionColumns, m.Points)
	for i := 0; i < m.DimensionRows; i++ {
		for j := 0; j < m.DimensionColumns; j++ {
			newM.Points[i][j] += matrix.Points[i][j]
		}
	}
	return newM, nil
}

func (m *Matrix) Copy() Matrix {
	var copyM Matrix
	copyM.Init(m.DimensionRows, m.DimensionColumns)
	for i := 0; i < m.DimensionRows; i++ {
		for j := 0; j < m.DimensionColumns; j++ {
			copyM.Points[i][j] = m.Points[i][j]
		}
	}
	return copyM
}

func (m *Matrix) Transponate() (Matrix, error) {
	var trM Matrix
	//if m.DimensionColumns != m.DimensionRows {
	//	return Matrix{}, fmt.Errorf("matrix is not square")
	//}
	trM.Init(m.DimensionColumns, m.DimensionRows)
	for i := 0; i < m.DimensionColumns; i++ {
		for j := 0; j < m.DimensionRows; j++ {
			trM.Points[i][j] = m.Points[j][i]
		}
	}
	return trM, nil
}

func (m *Matrix) CheckDiagonallyDominant() bool {
	var sum float64
	for i := 0; i < m.DimensionRows; i++ {
		sum = 0
		for j := 0; j < m.DimensionColumns; j++ {
			if i != j {
				sum += math.Abs(m.Points[i][j])
			}
		}
		if math.Abs(m.Points[i][i]) <= sum {
			return false
		}
	}
	return true
}

func (m *Matrix) MakeTriangular() Matrix {
	var newM Matrix
	newM.InitWithPoints(m.DimensionRows, m.DimensionColumns, m.Points)
	for i := 0; i < m.DimensionRows; i++ {
		m.MakeNullColumn(i, &newM)
	}
	return newM
}

func (m *Matrix) MakeE() Matrix {
	var newM Matrix
	newM.InitWithPoints(m.DimensionRows, m.DimensionColumns, m.Points)
	for i := 0; i < m.DimensionRows; i++ {
		m.MakeNullColumnE(i, &newM)
	}
	for i := m.DimensionRows - 1; i >= 0; i-- {
		m.MakeNullColumnEBackwards(i, &newM)
	}
	return newM
}

func (m *Matrix) MakeNullColumn(diagonalIndex int, newMatrix *Matrix) {
	var newValue, koeff, newKoeff float64
	dimension := m.DimensionColumns
	koeff = m.Points[diagonalIndex][diagonalIndex]
	normLine := make([]float64, dimension)
	// добавить проверку нулевой строки
	for j := 0; j < dimension; j++ {
		normLine[j] = m.Points[diagonalIndex][j] / koeff
	}
	for k := diagonalIndex + 1; k < m.DimensionRows; k++ {
		newKoeff = m.Points[k][diagonalIndex]
		if newKoeff != 0 {
			for t := 0; t < dimension; t++ {
				newValue = m.Points[k][t] - newKoeff*normLine[t]
				newMatrix.Points[k][t] = newValue
			}
		}
	}
}

func (m *Matrix) MakeNullColumnE(diagonalIndex int, newMatrix *Matrix) {
	var newValue, koeff, newKoeff float64
	dimension := m.DimensionColumns
	koeff = m.Points[diagonalIndex][diagonalIndex]
	if koeff != 0 {
		for j := 0; j < dimension; j++ {
			m.Points[diagonalIndex][j] /= koeff
		}
	}
	for k := diagonalIndex + 1; k < m.DimensionRows; k++ {
		newKoeff = m.Points[k][diagonalIndex]
		if newKoeff != 0 {
			for t := 0; t < dimension; t++ {
				newValue = m.Points[k][t] - newKoeff*m.Points[diagonalIndex][t]
				newMatrix.Points[k][t] = newValue
			}
		}
	}
}

func (m *Matrix) MakeNullColumnEBackwards(diagonalIndex int, newMatrix *Matrix) {
	var newValue, koeff, newKoeff float64
	dimension := m.DimensionColumns
	koeff = m.Points[diagonalIndex][diagonalIndex]
	if koeff != 0 {
		for j := 0; j < dimension; j++ {
			m.Points[diagonalIndex][j] /= koeff
		}
	}
	for k := diagonalIndex - 1; k >= 0; k-- {
		newKoeff = m.Points[k][diagonalIndex]
		if newKoeff != 0 {
			for t := 0; t < dimension; t++ {
				newValue = m.Points[k][t] - newKoeff*m.Points[diagonalIndex][t]
				newMatrix.Points[k][t] = newValue
			}
		}
	}
}

func (m *Matrix) LU() (Matrix, Matrix, error) {
	if m.DimensionRows != m.DimensionColumns {
		return Matrix{}, Matrix{}, fmt.Errorf("matrix dimensions don't match")
	}
	var size int
	size = m.DimensionColumns
	var uMatrix, lMatrix Matrix
	uMatrix.Init(m.DimensionRows, m.DimensionColumns)
	lMatrix.Init(m.DimensionRows, m.DimensionColumns)
	var sum float64

	for i := 0; i < size; i++ {
		for j := 0; j < size; j++ {
			uMatrix.Points[i][j] = 0
			lMatrix.Points[i][j] = 0
			lMatrix.Points[i][i] = 1
		}
	}

	for i := 0; i < size; i++ {
		for j := 0; j < size; j++ {
			if i <= j {
				sum = 0
				for k := 0; k <= i-1; k++ {
					sum += lMatrix.Points[i][k] * uMatrix.Points[k][j]
				}
				uMatrix.Points[i][j] = m.Points[i][j] - sum

			}
			if i > j {
				sum = 0
				for k := 0; k <= j-1; k++ {
					sum += lMatrix.Points[i][k] * uMatrix.Points[k][j]
				}
				lMatrix.Points[i][j] = (m.Points[i][j] - sum) / uMatrix.Points[j][j]
			}
		}
	}
	return lMatrix, uMatrix, nil
}

func (m *Matrix) Inverted() (Matrix, error) {
	var inverted, lM, uM Matrix
	var vecOne, gV, zV Vector
	var err error
	var gaussMethod GaussMethod
	var zeidelMethod ZeidelMethod
	zeidelMethod.Init(0.0001)
	lM, uM, err = m.LU()
	if err != nil {
		return Matrix{}, fmt.Errorf("error during LU: %v", err)
	}
	//checkM, _ := lM.MulM(uM)
	//checkM.Print()
	inverted.Init(m.DimensionRows, m.DimensionColumns)
	for i := 0; i < m.DimensionRows; i++ {
		vecOne.Init(m.DimensionRows)
		vecOne.Points[i] = 1
		gV, err = gaussMethod.Solve(lM, vecOne)
		//gV.Print()
		if err != nil {
			return Matrix{}, fmt.Errorf("error during gauss solving: %v", err)
		}

		zV, err = zeidelMethod.Solve(uM, gV)
		//zV.Print()
		if err != nil {
			return Matrix{}, fmt.Errorf("error during zeidel solving: %v", err)
		}
		for j := 0; j < m.DimensionRows; j++ {
			inverted.Points[j][i] = zV.Points[j]
		}
	}
	return inverted, nil
}

func (m *Matrix) Print() {
	for i := 0; i < m.DimensionRows; i++ {
		for j := 0; j < m.DimensionColumns; j++ {
			fmt.Printf("%f ", m.Points[i][j])
		}
		fmt.Println()
	}
	fmt.Println()
}

func (m *Matrix) RemoveRows(rows []int) (Matrix, error) {
	var err error
	var newPoints [][]float64
	for i := 0; i < m.DimensionRows; i++ {
		found := false
		for _, r := range rows {
			if i == r {
				found = true
				break
			}
		}
		if !found {
			newPoints = append(newPoints, m.Points[i])
		}
	}
	var newMatrix Matrix
	err = newMatrix.InitWithPoints(len(newPoints), m.DimensionColumns, newPoints)
	if err != nil {
		return Matrix{}, fmt.Errorf("error initializing vector: %v", err)
	}
	return newMatrix, nil
}
