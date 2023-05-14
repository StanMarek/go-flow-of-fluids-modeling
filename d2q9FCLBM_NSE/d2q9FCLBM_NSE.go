package d2q9FCLBM_NSE

import (
	"math"
	"math/rand"
)

const (
	LX = 32 * 2
	LY = 32

	// F d2q9
	F    = 9
	csSq = 1. / 3.
	wC   = 4. / 9.
	wS   = 1. / 9.
	wL   = 1. / 36.

	// FT d2q5
	FT    = 5
	csSqT = 1. / 3.
	wCT   = 1. / 3.
	wNCT  = 1. / 6.

	// FLUID boundary conditions
	FLUID       = 0
	SOLID       = 1
	BcIn        = 0
	BcOut       = 2
	BcHot       = 0
	BcCold      = 4
	BcAdiabatic = 8
	BcTCond     = 16

	// LES turbulence model
	CLes = 0.19

	RandMax = 2147483647
)

// d2q9 enum
const (
	_C  int = 0
	_N  int = 1
	_S  int = 2
	_W  int = 3
	_E  int = 4
	_NE int = 5
	_NW int = 6
	_SE int = 7
	_SW int = 8
)

// d2q5 enum
const (
	_CT int = 0
	_NT int = 1
	_ST int = 2
	_WT int = 3
	_ET int = 4
)

var (
	cx = [F]float64{0, 0, 0, -1, 1, 1, -1, 1, -1}
	cy = [F]float64{0, 1, -1, 0, 0, 1, 1, -1, -1}
	w  = [F]float64{wC, wS, wS, wS, wS, wL, wL, wL, wL}

	cxt = [FT]float64{0, 0, 0, -1, 1}
	cyt = [FT]float64{0, 1, -1, 0, 0}
	wt  = [FT]float64{wCT, wNCT, wNCT, wNCT, wNCT}

	cells    = [LX * LY * F]float64{}
	tmpCells = [LX * LY * F]float64{}

	mapArray = [LX * LY]int{}
	Rho      = [LX * LY]float64{}
	Ux       = [LX * LY]float64{}
	Uy       = [LX * LY]float64{}
	Temp     = [LX * LY]float64{}

	tau        float64
	uinLB      float64
	vinLB      float64
	betagxLB   float64
	betagyLB   float64
	taut       float64
	tauts      float64
	tempHotLB  = 1.0
	tempColdLB = 0.0
	t0LB       = 0.0
	VelMax     = 0.0
	TempMax    = 0.0
)

func AR(x, y, f int) int {
	return f + F*(x) + F*LX*(y)
}

func ART(x, y, f int) int {
	return f + FT*(x) + FT*LX*(y)
}

func ARM(x, y int) int {
	return x + LX*y
}

func feq(i int, rho, u, v float64) float64 {
	cu := (cx[i]*u + cy[i]*v) / csSq
	return w[i] * rho * (1 + cu + 0.5*cu*cu - 0.5*(u*u+v*v)/csSq)
}

func feqT(i int, temp, u, v float64) float64 {
	cu := (cxt[i]*u + cyt[i]*v) / csSqT
	return wt[i] * temp * (1 + cu + 0.5*cu*cu - 0.5*(u*u+v*v)/csSqT)
}

func fi(i int, tau, u, v, Fx, Fy float64) float64 {
	return w[i] * (cx[i]*Fx + cy[i]*Fy) / csSqT
}

func collideStream(c, tc, cT, tcT []float64) {
	//var rho, ux, uy, temp float64
	//var fstar, fstart float64[FT]{}

	velMax := 0.0
	tempMax := 0.0

	for y := 0; y < LY; y++ {
		for x := 0; x < LX; x++ {
			ctype := mapArray[ARM(x, y)]
			if ctype == FLUID {
				rho := 0.0
				ux := 0.0
				uy := 0.0
				temp := 0.0

				for i := 0; i < F; i++ {
					rho += c[AR(x, y, i)]
					ux += cx[i] * c[AR(x, y, i)]
					uy += cy[i] * c[AR(x, y, i)]
				}
				for i := 0; i < FT; i++ {
					temp += cT[ART(x, y, i)]
				}

				temp += 0.000001 * (.5 - float64(rand.Int()/RandMax))
				Temp[ARM(x, y)] = temp

				ux /= rho
				uy /= rho

				// Force
				Fx := betagxLB * (temp - t0LB)
				Fy := betagyLB * (temp - t0LB)
				ux += 0.5 * Fx / rho
				uy += 0.5 * Fy / rho
				Rho[ARM(x, y)] = rho
				Ux[ARM(x, y)] = ux
				Uy[ARM(x, y)] = uy
				if math.Sqrt(ux*ux+uy*uy) > VelMax {
					velMax = math.Sqrt(ux*ux + uy*uy)
				} else {
					velMax = VelMax
				}
				if temp > TempMax {
					tempMax = temp
				} else {
					tempMax = TempMax
				}

				// LES turbulence model
				// compute non-equilibrium part of stress tensor PI
			}
		}
	}
}
