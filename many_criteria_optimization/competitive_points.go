package many_criteria_optimization

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"github.com/saskamegaprogrammist/optimization_methods/many_dimension_search"
	"github.com/saskamegaprogrammist/optimization_methods/random_points_gen"
)

//var nPow = 5

type CompetitivePointsMultistart struct {
	alpha      float64
	beta       float64
	dimension  int
	targetFunc func(xs []float64) float64
	generator  random_points_gen.StDistributionGen
	search     many_dimension_search.NelderMeadSearch
}

func (cpm *CompetitivePointsMultistart) Init(dimension int,
	alpha float64, beta float64, targetFunc func(xs []float64) float64) {
	cpm.targetFunc = targetFunc
	cpm.dimension = dimension
	cpm.alpha = alpha
	cpm.beta = beta
	cpm.generator.Init(alpha, beta, n, dimension)
}

func (cpm *CompetitivePointsMultistart) Solve() ([]float64, float64, error) {
	var err error
	var clustersLen = n
	points := cpm.generator.Generate()

	var vectors = make([]la_methods.Vector, clustersLen)

	for i, p := range points {
		vectors[i] = la_methods.Vector{
			Points:    p,
			Dimension: cpm.dimension,
		}
	}

	var clusterMin = make([]float64, clustersLen)
	var clusterMinV = make([]la_methods.Vector, clustersLen)
	var clusters = make([][]la_methods.Vector, clustersLen)
	for i := 0; i < clustersLen; i++ {
		clusters[i] = append(clusters[i], vectors[i])
		clusterMinV[i] = vectors[i]
		clusterMin[i] = cpm.targetFunc(vectors[i].Points)
	}

	for clustersLen != 1 {
		var xMin []float64
		var yMin float64

		for i := 0; i < clustersLen; i++ {
			cpm.search.Init(clusters[i][0].Points, 0.1, cpm.dimension, 0.001, cpm.targetFunc)
			xMin, yMin, err = cpm.search.Solve()
			//fmt.Println(yMin)
			if err != nil {
				return nil, 0, fmt.Errorf("error solver nelder mead: %v", err)
			}
			newV := la_methods.Vector{
				Points:    xMin,
				Dimension: cpm.dimension,
			}
			clusters[i][0] = newV
			clusterMin[i] = yMin
			clusterMinV[i] = newV
		}

		minDist := float64(1000000000)
		var dist float64
		var minI1, minI2 int
		for i := 0; i < clustersLen; i++ {
			for j := i + 1; j < clustersLen; j++ {

				len1 := len(clusters[i])
				len2 := len(clusters[j])
				minDistInner := float64(1000000000)
				for k := 0; k < len1; k++ {
					for r := 0; r < len2; r++ {
						dist, err = clusters[i][k].EqDist(clusters[j][r])
						if err != nil {
							return nil, 0, fmt.Errorf("error caluclating dist: %v", err)
						}
						if dist < minDistInner {
							minDistInner = dist
						}
					}
				}

				if dist < minDist {
					minDist = dist
					minI1 = i
					minI2 = j
				}
			}
		}
		var newClusters [][]la_methods.Vector
		var newClustersMin []float64
		var newClustersMinV []la_methods.Vector
		clustersLen--
		for i, c := range clusters {
			if i != minI2 && i != minI1 {
				newClusters = append(newClusters, c)
				newClustersMin = append(newClustersMin, clusterMin[i])
				newClustersMinV = append(newClustersMinV, clusterMinV[i])
			}
		}
		if clusterMin[minI1] < clusterMin[minI2] {
			newClustersMin = append(newClustersMin, clusterMin[minI1])
			newClustersMinV = append(newClustersMinV, clusterMinV[minI1])
			newClusters = append(newClusters, clusters[minI1])

		} else {
			newClustersMin = append(newClustersMin, clusterMin[minI2])
			newClustersMinV = append(newClustersMinV, clusterMinV[minI2])
			newClusters = append(newClusters, clusters[minI2])
		}
		clusters = newClusters
		clusterMinV = newClustersMinV
		clusterMin = newClustersMin
	}
	return clusterMinV[0].Points, clusterMin[0], nil
}
