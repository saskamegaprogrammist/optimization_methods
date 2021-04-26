package genetic_methods

import (
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"github.com/saskamegaprogrammist/optimization_methods/random_points_gen"
	"math/rand"
	"time"
)

var (
	Pm = 0.5
	Pc = 0.5
)

type HelpPoint struct {
	index int
	point []float64
}

type GeneticAlgorithm struct {
	Np          int // max populations size
	Mp          int // population size
	dimension   int
	targetFunc  func(xs []float64) float64
	fitnessFunc func(xs []float64) float64
	generator   random_points_gen.StDistributionGen
	//search      many_dimension_search.NelderMeadSearch
	startPoint []float64
	alpha      float64
	beta       float64
}

func (ga *GeneticAlgorithm) Init(alpha float64, beta float64, Np int, Mp int, startPoint []float64, dimension int,
	targetFunc func(xs []float64) float64, fitnessFunc func(xs []float64) float64) {
	ga.targetFunc = targetFunc
	ga.fitnessFunc = fitnessFunc
	ga.startPoint = startPoint
	ga.Np = Np
	ga.Mp = Mp
	ga.beta = beta
	ga.alpha = alpha
	ga.dimension = dimension
	ga.generator.Init(alpha, beta, Mp, dimension)
}

func (ga *GeneticAlgorithm) Solve() ([]float64, float64, error) {
	//var err error
	var t = 1
	var k = 1

	// generate initial population

	pointsInit := ga.generator.Generate()
	var vectorsInit = make([]la_methods.Vector, ga.Mp)
	var fitnessVals = make([]float64, ga.Mp)

	for i, p := range pointsInit {
		vectorsInit[i] = la_methods.Vector{
			Points:    p,
			Dimension: ga.dimension,
		}
		fitnessVals[i] = ga.fitnessFunc(p)
		//fmt.Println(p, ga.fitnessFunc(p))
	}

	var max float64
	var maxPoint []float64
	for ; t <= ga.Np; t++ {
		points := pointsInit

		var fitnessVals = make([]float64, ga.Mp)
		var fitnessValsSum float64
		var minFitness = float64(100000000)
		var minFitnessIndex int
		for ; k <= ga.Mp; k++ {
			//fmt.Println(k)
			// calculate cumulative probability

			var qs = make([]float64, ga.Mp)
			var sum float64
			var fV float64
			for i, p := range points {
				fV = ga.fitnessFunc(p)
				if fV < minFitness {
					minFitness = ga.fitnessFunc(p)
					minFitnessIndex = i
				}
				fitnessValsSum += fV
				sum += fV
				qs[i] = sum
			}

			//fmt.Println(qs, fitnessValsSum)

			// selection

			randomSource := rand.NewSource(time.Now().UnixNano())
			random := rand.New(randomSource)

			var indexMap = make(map[int]bool, ga.Mp)

			var newPoints [][]float64
			var r float64
			for i := 0; i < ga.Mp; i++ {
				r = random.Float64() * fitnessValsSum
				var ind = 1
				for ; ind < ga.Mp; ind++ {
					//fmt.Println(r, ind, qs[ind])
					if r > qs[ind-1] && r <= qs[ind] {
						break
					}
				}
				if ind == ga.Mp {
					continue
				}
				if !indexMap[ind-1] {
					newPoints = append(newPoints, points[ind-1])
				}
				indexMap[ind-1] = true
			}

			//fmt.Println(newPoints)

			var selectedLen = len(newPoints)

			// crossingover

			var newIndexMap = make(map[int]bool, selectedLen)

			var crossNumber = float64(selectedLen) * Pc
			var parentPoints []HelpPoint
			for i := 0; i < int(crossNumber); i++ {
				r = random.Float64()
				if r < Pc {
					if !newIndexMap[i] {
						parentPoints = append(parentPoints, HelpPoint{point: newPoints[i], index: i})
					}
					newIndexMap[i] = true
				}
			}
			//fmt.Println(parentPoints)

			var parentPairsLen = len(parentPoints) / 2
			var parentsFirst = make([]HelpPoint, parentPairsLen)
			var parentsSecond = make([]HelpPoint, parentPairsLen)
			for i := 0; i < parentPairsLen; i++ {
				parentsFirst[i] = parentPoints[i]
				parentsSecond[i] = parentPoints[parentPairsLen+i]
			}

			var x = make([]float64, ga.dimension)
			var y = make([]float64, ga.dimension)

			c := random.Float64()
			for i := 0; i < parentPairsLen; i++ {
				var acceptableX = true
				var acceptableY = true
				for j := 0; j < ga.dimension; j++ {
					x[j] = c*parentsFirst[i].point[j] + (1-c)*parentsSecond[i].point[j]
					y[j] = (1-c)*parentsFirst[i].point[j] + c*parentsSecond[i].point[j]

					acceptableX = acceptableX && x[j] <= ga.beta && x[j] >= ga.alpha
					acceptableY = acceptableY && y[j] <= ga.beta && y[j] >= ga.alpha
				}
				if acceptableX && acceptableY {
					newPoints[parentsFirst[i].index] = x
					newPoints[parentsSecond[i].index] = y
				}
			}

			// mutation

			var mutationNumber = float64(selectedLen) * Pm

			var mutationParentPoints []HelpPoint

			newIndexMap = make(map[int]bool, selectedLen)

			for i := 0; i < int(mutationNumber); i++ {
				r = random.Float64()
				if r < Pm {
					if !newIndexMap[i] {
						mutationParentPoints = append(mutationParentPoints, HelpPoint{point: newPoints[i], index: i})
					}
					newIndexMap[i] = true
				}
			}
			//fmt.Println(mutationParentPoints)

			var mutationParentPairsLen = len(mutationParentPoints)

			var rInt int
			for i := 0; i < mutationParentPairsLen; i++ {
				rInt = random.Intn(ga.dimension)
				r = rand.Float64()*(ga.beta-ga.alpha) + ga.alpha
				mutationParentPoints[i].point[rInt] = r
			}

			// choosing mutation point

			if mutationParentPairsLen != 0 {
				rInt = random.Intn(mutationParentPairsLen)
				points[minFitnessIndex] = mutationParentPoints[rInt].point
				fitnessVals[minFitnessIndex] = ga.fitnessFunc(mutationParentPoints[rInt].point)
			}

			//for i := 0; i < len(points); i++ {
			//	fmt.Println(points[i], ga.targetFunc(points[i]))
			//}
		}
		var maxI int
		var maxFitness = float64(0)
		for i := 0; i < ga.Mp; i++ {
			if fitnessVals[i] > maxFitness {
				maxFitness = fitnessVals[i]
				maxI = i
			}
		}
		max = ga.targetFunc(points[maxI])
		maxPoint = points[maxI]
	}

	return maxPoint, max, nil
}
