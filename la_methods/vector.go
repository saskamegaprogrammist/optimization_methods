package la_methods

import (
	"fmt"
	"math"
)

type Vector struct {
	Points    []float64
	Dimension int
}

func (v *Vector) Init(dimension int) {
	v.Dimension = dimension
	v.Points = make([]float64, dimension)
}

func (v *Vector) InitWithPoints(dimension int, points []float64) error {
	if len(points) != dimension {
		return fmt.Errorf("dimension doesn't match points length: %d", len(points))
	}
	v.Dimension = dimension
	v.Points = points
	return nil
}

func (v *Vector) InitWithValue(dimension int, value float64) {
	v.Dimension = dimension
	v.Points = make([]float64, dimension)
	for i, _ := range v.Points {
		v.Points[i] = value
	}
}

func (v *Vector) MulOnValue(val float64) Vector {
	var mulVector Vector
	mulVector.Init(v.Dimension)
	for i, point := range v.Points {
		mulVector.Points[i] = point * val
	}
	return mulVector
}

func (v *Vector) Add(vector Vector) (Vector, error) {
	var sumVector Vector
	if vector.Dimension != v.Dimension {
		return sumVector, fmt.Errorf("vector dimensions doesn't match: %d", vector.Dimension)
	}
	sumVector.Init(v.Dimension)
	for i, point := range v.Points {
		sumVector.Points[i] = point + vector.Points[i]
	}
	return sumVector, nil
}

func (v *Vector) AddKOnIndex(k float64, i int) (Vector, error) {
	var sumVector Vector
	if i >= v.Dimension {
		return sumVector, fmt.Errorf("vector dimensions doesn't match: %d", i)
	}
	sumVector = v.Copy()
	sumVector.Points[i] += k
	return sumVector, nil
}

func (v *Vector) AddK(k float64) Vector {
	var sumVector Vector
	sumVector.Init(v.Dimension)
	for i, point := range v.Points {
		sumVector.Points[i] = point + k
	}
	return sumVector
}

func (v *Vector) Sub(vector Vector) (Vector, error) {
	var subVector Vector
	if vector.Dimension != v.Dimension {
		return subVector, fmt.Errorf("vector dimensions doesn't match: %d", vector.Dimension)
	}
	subVector.Init(v.Dimension)
	for i, point := range v.Points {
		subVector.Points[i] = point - vector.Points[i]
	}
	return subVector, nil
}

func (v *Vector) SubKOnIndex(k float64, i int) (Vector, error) {
	var subVector Vector
	if i >= v.Dimension {
		return subVector, fmt.Errorf("vector dimensions doesn't match: %d", i)
	}
	subVector = v.Copy()
	subVector.Points[i] -= k
	return subVector, nil
}

func (v *Vector) SubK(k float64) Vector {
	var subVector Vector
	subVector.Init(v.Dimension)
	for i, point := range v.Points {
		subVector.Points[i] = point - k
	}
	return subVector
}

func (v *Vector) Mul(vec Vector) (float64, error) {
	if vec.Dimension != v.Dimension {
		return 0, fmt.Errorf("vectors dimensions are different: %d", vec.Dimension)
	}
	var sum float64
	for i, point := range v.Points {
		sum += point * vec.Points[i]
	}
	return sum, nil
}

func (v *Vector) MulV(vec Vector) (Vector, error) {
	return v.MulVFrom(vec, 0)
}

func (v *Vector) MulVFrom(vec Vector, ind int) (Vector, error) {
	var mulVector Vector
	mulVector.Init(vec.Dimension)
	if vec.Dimension != v.Dimension {
		return mulVector, fmt.Errorf("vectors dimensions are different: %d", vec.Dimension)
	}
	for i := ind; i < v.Dimension; i++ {
		mulVector.Points[i] = v.Points[i] * vec.Points[i]
	}
	return mulVector, nil
}

func (v *Vector) Copy() Vector {
	var copyVector Vector
	copyVector.Init(v.Dimension)
	for i, point := range v.Points {
		copyVector.Points[i] = point
	}
	return copyVector
}

func (v *Vector) Len() float64 {
	var sum float64
	for _, point := range v.Points {
		sum += math.Pow(point, 2)
	}
	return math.Sqrt(sum)
}

func (v *Vector) ElemsSum(ind int) float64 {
	var sum float64
	for i := ind; i < v.Dimension; i++ {
		sum += v.Points[i]
	}
	return sum
}

func (v *Vector) Print() {
	for i, p := range v.Points {
		fmt.Printf("point %d: %f\n", i, p)
	}
}
