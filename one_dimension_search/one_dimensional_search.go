package one_dimension_search

type OneDimensionSearchI interface {
	Init(aStart float64, bStart float64, precision float64, targetFunc func(x float64) float64)
	Solve() (float64, float64)
	CountConvergence() (float64, error)
}
