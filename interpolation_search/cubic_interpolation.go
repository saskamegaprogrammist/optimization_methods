package interpolation_search

import (
	"math"
)

type CubicInterpolation struct {
	startPoint           float64
	step                 float64
	precisionDelta       float64
	precisionEps         float64
	targetFunc           func(x float64) float64
	targetFuncDerivative func(x float64) float64
}

func (cInt *CubicInterpolation) Init(startPoint float64, step float64, precisionDelta float64, precisionEps float64,
	targetFunc func(x float64) float64, targetFuncDerivative func(x float64) float64) {
	cInt.startPoint = startPoint
	cInt.step = step
	cInt.precisionDelta = precisionDelta
	cInt.precisionEps = precisionEps
	cInt.targetFunc = targetFunc
	cInt.targetFuncDerivative = targetFuncDerivative
}

func (cInt *CubicInterpolation) Solve() (float64, float64) {
	var alph0, alphM, alphM_1 float64
	var der0 float64
	var sign, i float64
	var alph1, alph2 float64
	var f1, f2, der1, der2 float64
	var alph, f float64
	alph0 = cInt.startPoint
	der0 = cInt.targetFuncDerivative(alph0)
	if der0 < 0 {
		sign = 1
	} else {
		sign = -1
	}
	alphM_1 = alph0
	for {
		alphM = alphM_1 + sign*math.Pow(2, i)*cInt.step
		if cInt.targetFuncDerivative(alphM_1)*cInt.targetFuncDerivative(alphM) <= 0 {
			break
		}
		i++
	}
	alph1 = alphM_1
	alph2 = alphM
	for {
		f1 = cInt.targetFunc(alph1)
		f2 = cInt.targetFunc(alph2)
		der1 = cInt.targetFuncDerivative(alph1)
		der2 = cInt.targetFuncDerivative(alph2)
		alph, f = cInt.interpolation(alph1, alph2, f1, f2, der1, der2)
		for {
			if f <= f1 {
				break
			}
			alph -= (alph - alph1) / 2
			f = cInt.targetFunc(alph)
		}
		der := cInt.targetFuncDerivative(alph)
		if cInt.checkFirstCondition(der) && cInt.checkSecondCondition(alph, alph1) {
			return alph, cInt.targetFunc(alph)
		}
		if der*der1 < 0 {
			alph2 = alph1
			alph1 = alph
		}
		if der*der2 < 0 {
			alph1 = alph
		}
	}
}

func (cInt *CubicInterpolation) interpolation(alph1 float64, alph2 float64,
	f1 float64, f2 float64, der1 float64, der2 float64) (float64, float64) {
	var w, z, m, alph float64
	z = 3*(f1-f2)/(alph2-alph1) + der1 + der2
	w = math.Sqrt(math.Pow(z, 2) - der1*der2)
	if alph1 > alph2 {
		w = -w
	}
	m = (der2 + w - z) / (der2 - der1 + 2*w)
	if m < 0 {
		alph = alph2
	} else if m >= 0 && m <= 1 {
		alph = alph2 - m*(alph2-alph1)
	} else if m > 1 {
		alph = alph1
	}
	return alph, cInt.targetFunc(alph)
}

func (cInt *CubicInterpolation) checkFirstCondition(der float64) bool {
	return math.Abs(der) <= cInt.precisionEps
}

func (cInt *CubicInterpolation) checkSecondCondition(alph float64, alph1 float64) bool {
	return math.Abs((alph-alph1)/alph) <= cInt.precisionDelta
}
