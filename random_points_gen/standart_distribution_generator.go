package random_points_gen

import (
	"math/rand"
	"time"
)

type StDistributionGen struct {
	start     float64
	end       float64
	n         int
	dimension int
	shift     float64
}

func (sdg *StDistributionGen) Init(start float64,
	end float64, n int, dimension int) {
	sdg.start = start
	sdg.end = end
	sdg.n = n
	sdg.dimension = dimension
	sdg.shift = end - start
}

func (sdg *StDistributionGen) Generate() [][]float64 {

	randomSource := rand.NewSource(time.Now().UnixNano())
	random := rand.New(randomSource)

	points := make([][]float64, sdg.n)
	for i := 0; i < sdg.n; i++ {
		points[i] = make([]float64, sdg.dimension)
	}
	//var val float64
	for i := 0; i < sdg.n; i++ {
		for j := 0; j < sdg.dimension; j++ {
			points[i][j] = random.Float64()*float64(sdg.n)*sdg.shift/(float64(sdg.n)-1) + sdg.start
		}
	}
	return points
}
