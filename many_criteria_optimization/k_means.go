package many_criteria_optimization

import (
	"fmt"
	"github.com/saskamegaprogrammist/optimization_methods/la_methods"
	"github.com/saskamegaprogrammist/optimization_methods/random_points_gen"
)

var n = 20

type KMeansMultistart struct {
	startPoint []float64
	alpha      float64
	beta       float64
	delta      float64
	dimension  int
	targetFunc func(xs []float64) float64
	generator  random_points_gen.LinGen
}

func (kmm *KMeansMultistart) Init(startPoint []float64, delta float64, dimension int,
	alpha float64, beta float64, targetFunc func(xs []float64) float64) {
	kmm.startPoint = startPoint
	kmm.delta = delta
	kmm.targetFunc = targetFunc
	kmm.dimension = dimension
	kmm.alpha = alpha
	kmm.beta = beta
	kmm.generator.Init(startPoint, []float64{1, 5, 3}, []float64{9, 2, 1}, beta, n)
}

func (kmm *KMeansMultistart) Solve() ([]float64, float64, error) {
	var err error
	var clustersLen = n
	var vectors = make([]la_methods.Vector, clustersLen)
	points := kmm.generator.Generate()
	for i, p := range points {
		vectors[i] = la_methods.Vector{
			Points:    p,
			Dimension: kmm.dimension,
		}
	}

	var clusterMin = make([]float64, clustersLen)
	var clusterMinV = make([]la_methods.Vector, clustersLen)
	var clusters = make([][]la_methods.Vector, clustersLen)
	for i := 0; i < clustersLen; i++ {
		clusters[i] = append(clusters[i], vectors[i])
		clusterMinV[i] = vectors[i]
		clusterMin[i] = kmm.targetFunc(vectors[i].Points)
	}

	//fmt.Println(clusterMin)

	for clustersLen != 1 {
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
		if minDist < kmm.delta {
			clustersLen--
			for i, c := range clusters {
				if i != minI2 && i != minI1 {
					newClusters = append(newClusters, c)
					newClustersMin = append(newClustersMin, clusterMin[i])
					newClustersMinV = append(newClustersMinV, clusterMinV[i])
				}
			}
			var newCluster []la_methods.Vector
			newCluster = clusters[minI1]
			newCluster = append(newCluster, clusters[minI2]...)
			if clusterMin[minI1] < clusterMin[minI2] {
				newClustersMin = append(newClustersMin, clusterMin[minI1])
				newClustersMinV = append(newClustersMinV, clusterMinV[minI1])
			} else {
				newClustersMin = append(newClustersMin, clusterMin[minI2])
				newClustersMinV = append(newClustersMinV, clusterMinV[minI2])
			}
			newClusters = append(newClusters, newCluster)
		}
		clusters = newClusters
		clusterMinV = newClustersMinV
		clusterMin = newClustersMin
	}

	return clusterMinV[0].Points, clusterMin[0], nil
}
