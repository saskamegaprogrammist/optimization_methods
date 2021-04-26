package simplex_methods

type Point struct {
	x         []float64
	f         float64
	dimension int
}

type PointHelp struct {
	point Point
	ind   int
}
