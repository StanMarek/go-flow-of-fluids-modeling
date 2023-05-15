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
	var fstar, fstart [F]float64

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
				if math.Sqrt(ux*ux+uy*uy) > velMax {
					VelMax = math.Sqrt(ux*ux + uy*uy)
				} else {
					VelMax = velMax
				}
				if temp > tempMax {
					TempMax = temp
				} else {
					TempMax = tempMax
				}

				// LES turbulence model
				// compute non-equilibrium part of stress tensor PI
				mxx := c[AR(x, y, _E)] - feq(_E, rho, ux, uy) + c[AR(x, y, _W)] - feq(_W, rho, ux, uy) + c[AR(x, y, _NE)] - feq(_NE, rho, ux, uy) + c[AR(x, y, _NW)] - feq(_NW, rho, ux, uy) + c[AR(x, y, _SE)] - feq(_SE, rho, ux, uy) + c[AR(x, y, _SW)] - feq(_E, rho, ux, uy)
				myy := c[AR(x, y, _S)] - feq(_S, rho, ux, uy) + c[AR(x, y, _N)] - feq(_N, rho, ux, uy) + c[AR(x, y, _NE)] - feq(_NE, rho, ux, uy) + c[AR(x, y, _NW)] - feq(_NW, rho, ux, uy) + c[AR(x, y, _SE)] - feq(_SE, rho, ux, uy) + c[AR(x, y, _SW)] - feq(_E, rho, ux, uy)
				mxy := c[AR(x, y, _NE)] - feq(_NE, rho, ux, uy) - (c[AR(x, y, _NW)] - feq(_NW, rho, ux, uy)) + c[AR(x, y, _SE)] - feq(_SE, rho, ux, uy) - (c[AR(x, y, _SW)] - feq(_SW, rho, ux, uy))
				PI := math.Sqrt(2 * (mxx*mxx + 2*mxy*mxy + myy*myy))

				// Smagorinsky LES
				tauLes := 0.5 * (tau + math.Sqrt(tau*tau+2*CLes*CLes*PI/(csSq*csSq*rho)))
				tautLes := 0.5 * (taut + math.Sqrt(taut*taut+2*CLes*CLes*PI/(csSqT*csSqT*rho)))

				fc := c[AR(x, y, _C)]
				fe := c[AR(x, y, _E)]
				fn := c[AR(x, y, _N)]
				fw := c[AR(x, y, _W)]
				fs := c[AR(x, y, _S)]
				fne := c[AR(x, y, _NE)]
				fnw := c[AR(x, y, _NW)]
				fsw := c[AR(x, y, _SW)]
				fse := c[AR(x, y, _SE)]

				k3 := 1. / 12. * (rho*(ux*ux+uy*uy) - fe - fn - fs - fw - 2.*(fse+fsw+fne+fnw-rho/3.))
				k4 := .25 / tauLes * (fn + fs - fe - fw + rho*(ux*ux-uy*uy))
				k5 := .25 / tauLes * (fne + fsw - fnw - fse - ux*uy*rho)

				kxxyyAt := (6*k3 + 2*k4 - rho*ux*ux + fe + fw + fne + fnw + fse + fsw) *
					(6*k3 - 2*k4 - rho*uy*uy + fn + fs + fne + fnw + fse + fsw)

				k6 := -((fse+fsw-fne-fnw-2*ux*ux*uy*rho+uy*(rho-fn-fs-fc))*.25 + .5*ux*(fne-fnw-fse+fsw) - .5*uy*(-3*k3-k4) - 2*ux*k5)
				k7 := -((fsw+fnw-fse-fne-2*uy*uy*ux*rho+ux*(rho-fw-fe-fc))*.25 + .5*uy*(fne+fsw-fse-fnw) - .5*ux*(-3*k3+k4) - 2*uy*k5)
				k8 := .25*(kxxyyAt-fne-fnw-fse-fsw+
					2*(ux*(fne-fnw+fse-fsw)+uy*(fne+fnw-fse-fsw))+
					4*ux*uy*(fnw-fne+fse-fsw)-
					ux*ux*(fn+fne+fnw+fs+fse+fsw)+
					uy*uy*(3*ux*ux*rho-fe-fne-fnw-fse-fsw-fw)) -
					2*k3 - 2*ux*k7 - 2*uy*k6 + 4*ux*uy*k5 +
					(-1.5*k3+.5*k4)*ux*ux +
					(-1.5*k3-.5*k4)*uy*uy

				fstar[_C] = fc + 4*(-k3+k8) + fi(_C, tauLes, ux, uy, Fx, Fy)
				fstar[_W] = fw - k3 - 2*k8 + k4 - 2*k7 + fi(_W, tauLes, ux, uy, Fx, Fy)
				fstar[_E] = fe - k3 - 2*k8 + k4 + 2*k7 + fi(_E, tauLes, ux, uy, Fx, Fy)
				fstar[_N] = fn - k3 - 2*k8 - k4 + 2*k6 + fi(_N, tauLes, ux, uy, Fx, Fy)
				fstar[_S] = fs - k3 - 2*k8 - k4 - 2*k6 + fi(_S, tauLes, ux, uy, Fx, Fy)
				fstar[_NW] = fnw + 2*k3 + k8 + k5 - k6 + k7 + fi(_NW, tauLes, ux, uy, Fx, Fy)
				fstar[_NE] = fne + 2*k3 + k8 - k5 - k6 - k7 + fi(_NE, tauLes, ux, uy, Fx, Fy)
				fstar[_SW] = fsw + 2*k3 + k8 - k5 + k6 + k7 + fi(_SW, tauLes, ux, uy, Fx, Fy)
				fstar[_SE] = fse + 2*k3 + k8 + k5 + k6 - k7 + fi(_SE, tauLes, ux, uy, Fx, Fy)

				for i := 0; i < FT; i++ {
					fstart[i] = (1.-1./tautLes)*cT[ART(x, y, i)] + feqT(i, temp, ux, uy)/tautLes
				}
			} else {
				// SOLID CELLS
				// Bounce-back for NSE
				fstar[_C] = c[AR(x, y, _C)]
				fstar[_W] = c[AR(x, y, _E)]
				fstar[_E] = c[AR(x, y, _W)]
				fstar[_N] = c[AR(x, y, _S)]
				fstar[_S] = c[AR(x, y, _N)]
				fstar[_NW] = c[AR(x, y, _SE)]
				fstar[_NE] = c[AR(x, y, _SW)]
				fstar[_SW] = c[AR(x, y, _NE)]
				fstar[_SE] = c[AR(x, y, _NW)]

				temp := 0.0

				for i := 0; i < FT; i++ {
					temp += cT[ART(x, y, i)]
				}

				Temp[ARM(x, y)] = temp
				switch ctype & 0x1c {
				case 0:
					for i := 0; i < FT; i++ {
						fstart[i] = feqT(i, tempHotLB, 0, 0)
					}
					Temp[ARM(x, y)] = tempHotLB
					break
				case 4:
					for i := 0; i < FT; i++ {
						fstart[i] = feqT(i, tempColdLB, 0, 0)
					}
					Temp[ARM(x, y)] = tempColdLB
					break
				case 8:
					fstart[_WT] = cT[ART(x, y, _ET)]
					fstart[_ET] = cT[ART(x, y, _WT)]
					fstart[_NT] = cT[ART(x, y, _ST)]
					fstart[_ST] = cT[ART(x, y, _NT)]
					break
				case 16:
					for i := 0; i < FT; i++ {
						// fstart[i] = (1-1/tauts)*cT[ART(x,y,i)] + feqT(i, temp, 0, 0)/tauts;
					}
					break
				}
			}

			// PERIODIC BC
			var (
				//xp := (x == LX - 1) ? 0 : x + 1;
				//yp := (y == LY - 1) ? 0 : y + 1;
				//xm := (x == 0) ? LX - 1 : x - 1;
				//ym := (y == 0) ? LY - 1 : y - 1;
				xp, yp, xm, ym int
			)

			if x == LX-1 {
				xp = 0
			} else {
				xp = x + 1
			}

			if y == LY-1 {
				yp = 0
			} else {
				yp = y + 1
			}

			if x == 0 {
				xm = LX - 1
			} else {
				xm = x - 1
			}

			if y == 0 {
				ym = LY - 1
			} else {
				ym = y - 1
			}

			// STREAMING
			tc[AR(x, y, _C)] = fstar[_C]
			tc[AR(x, yp, _N)] = fstar[_N]
			tc[AR(x, ym, _S)] = fstar[_S]
			tc[AR(xm, y, _W)] = fstar[_W]
			tc[AR(xp, y, _E)] = fstar[_E]
			tc[AR(xp, yp, _NE)] = fstar[_NE]
			tc[AR(xp, ym, _SE)] = fstar[_SE]
			tc[AR(xm, yp, _NW)] = fstar[_NW]
			tc[AR(xm, ym, _SW)] = fstar[_SW]

			tcT[ART(x, y, _CT)] = fstart[_CT]
			tcT[ART(x, yp, _NT)] = fstart[_NT]
			tcT[ART(x, ym, _ST)] = fstart[_ST]
			tcT[ART(xm, y, _WT)] = fstart[_WT]
			tcT[ART(xp, y, _ET)] = fstart[_ET]
		}
	}
}
