package random_points_gen

import "math"

type LinGen struct {
	startPoint []float64
	a          []float64
	c          []float64
	end        float64
	n          int
}

func (lg *LinGen) Init(startPoint []float64, a []float64, c []float64,
	end float64, n int) {
	lg.startPoint = startPoint
	lg.a = a
	lg.c = c
	lg.end = end
	lg.n = n
}

func (lg *LinGen) Generate() [][]float64 {
	dimension := len(lg.startPoint)
	points := make([][]float64, lg.n)
	points[0] = lg.startPoint
	for j := 1; j < lg.n; j++ {
		points[j] = make([]float64, dimension)
	}
	for i := 0; i < dimension; i++ {
		for j := 1; j < lg.n; j++ {
			points[j][i] = gen(lg.a[i], lg.c[i], points[j-1][i], lg.end)
		}
	}
	return points
}

func gen(a float64, c float64, last float64, m float64) float64 {
	return math.Mod(a*last+c, m)
}
