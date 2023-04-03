package d2q9NSE

import (
	"fmt"
	"math"
	"mppwc/util"
	"os"
)

const (
	LX = 128 * 2
	LY = 128

	F    = 9
	csSq = 1. / 3.
	wC   = 4. / 9.
	wS   = 1. / 9.
	wL   = 1. / 36.
)

var (
	cx = [F]float64{0, 0, 0, -1, 1, 1, -1, 1, -1}
	cy = [F]float64{0, 1, -1, 0, 0, 1, 1, -1, -1}
	w  = [F]float64{wC, wS, wS, wS, wS, wL, wL, wL, wL}

	cells    = [LX * LY * F]float64{}
	tmpCells = [LX * LY * F]float64{}

	mapArray = [LX * LY]int{}
	Rho      = [LX * LY]float64{}
	Ux       = [LX * LY]float64{}
	Uy       = [LX * LY]float64{}

	tau float64
)

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

var AID = util.GetFromIndex

func AR(x, y, f int, getFromIndex func(x, y, f, F, LX int) int) int {
	return getFromIndex(x, y, f, F, LX)
}

func ARM(x, y int) int {
	return x + LX*y
}

func feq(i int, rho, u, v float64) float64 {
	cu := (cx[i]*u + cy[i]*v) / csSq
	return w[i] * rho * (1 + cu + 0.5*cu*cu - 0.5*(u*u+v*v)/csSq)
}

func collideStream(c, tc []float64) {
	var rho, ux, uy float64
	var fStar = [F]float64{}
	var xp, xm, ym, yp int

	for y := 0; y < LY; y++ {
		for x := 0; x < LX; x++ {
			rho = 0
			ux = 0
			uy = 0
			for i := 0; i < F; i++ {
				rho += c[AR(x, y, i, AID)]
				ux += cx[i] * c[AR(x, y, i, AID)]
				uy += cy[i] * c[AR(x, y, i, AID)]
			}
			ux /= rho
			uy /= rho
			Rho[ARM(x, y)] = rho
			Ux[ARM(x, y)] = ux
			Uy[ARM(x, y)] = uy

			if mapArray[ARM(x, y)] != 0 {
				fStar[_C] = c[AR(x, y, _C, AID)]
				fStar[_W] = c[AR(x, y, _E, AID)]
				fStar[_E] = c[AR(x, y, _W, AID)]
				fStar[_N] = c[AR(x, y, _S, AID)]
				fStar[_S] = c[AR(x, y, _N, AID)]
				fStar[_NW] = c[AR(x, y, _SE, AID)]
				fStar[_NE] = c[AR(x, y, _SW, AID)]
				fStar[_SW] = c[AR(x, y, _NE, AID)]
				fStar[_SE] = c[AR(x, y, _NW, AID)]
			} else {
				for i := 0; i < F; i++ {
					fStar[i] = (1-1/tau)*c[AR(x, y, i, AID)] + feq(i, rho, ux, uy)/tau
				}
			}
			// BC
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
			tc[AR(x, y, _C, AID)] = fStar[_C]
			tc[AR(x, yp, _N, AID)] = fStar[_N]
			tc[AR(x, ym, _S, AID)] = fStar[_S]
			tc[AR(xm, y, _W, AID)] = fStar[_W]
			tc[AR(xp, y, _E, AID)] = fStar[_E]
			tc[AR(xp, yp, _NE, AID)] = fStar[_NE]
			tc[AR(xp, ym, _SE, AID)] = fStar[_SE]
			tc[AR(xm, yp, _NW, AID)] = fStar[_NW]
			tc[AR(xm, ym, _SW, AID)] = fStar[_SW]
		}
	}
}

func setIC(rho0, u0, v0 float64) {
	for y := 0; y < LY; y++ {
		for x := 0; x < LX; x++ {
			for i := 0; i < F; i++ {
				cells[AR(x, y, i, AID)] = feq(i, rho0, u0, v0)
				Rho[ARM(x, y)] = rho0
				Ux[ARM(x, y)] = u0
				Uy[ARM(x, y)] = v0
				mapArray[ARM(x, y)] = 0
			}
		}
	}

	for x := LX/2 - 8 - 32; x < LX/2+8-32; x++ {
		for y := LY/2 - 8 - 32; y < LY/2+8-32; y++ {
			mapArray[ARM(x, y)] = 1.
		}
	}

	for x := LX/2 - 8; x < LX/2+8; x++ {
		for y := LY/2 - 8; y < LY/2+8; y++ {
			mapArray[ARM(x, y)] = 1.
		}
	}
}

func dumpStateVTK(filename string) {
	file, err := os.Create(filename)
	util.Check(err)

	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			panic(err)
		}
	}(file)

	appendFile := util.AppendFile
	intToString := util.IntToString
	floatToStringENotation := util.FloatToStringENotation

	appendFile(filename, "# vtk DataFile Version 2.0\n")
	appendFile(filename, "2D-ADE data file \n")
	appendFile(filename, "ASCII\n")
	appendFile(filename, "DATASET RECTILINEAR_GRID\n")
	appendFile(filename, "DIMENSIONS "+intToString(LX)+" "+intToString(LY)+" 1\n")
	appendFile(filename, "X_COORDINATES "+intToString(LX)+"\n")
	for x := 0; x < LX; x++ {
		appendFile(filename, intToString(x)+" ")
	}
	appendFile(filename, "\n")
	appendFile(filename, "Y_COORDINATES "+intToString(LY)+"\n")
	for x := 0; x < LY; x++ {
		appendFile(filename, intToString(x)+" ")
	}
	appendFile(filename, "\n")
	appendFile(filename, "Z_COORDINATES 1 int\n")
	appendFile(filename, "0\n")
	appendFile(filename, "POINT_DATA "+intToString(LX*LY)+"\n")
	appendFile(filename, "SCALARS density double 1\n")
	appendFile(filename, "LOOKUP_TABLE default\n")
	for y := 0; y < LY; y++ {
		for x := 0; x < LX; x++ {
			appendFile(filename, floatToStringENotation(Rho[ARM(x, y)])+"\n")
		}
	}
	appendFile(filename, "VECTORS velocity double\n")
	for y := 0; y < LY; y++ {
		for x := 0; x < LX; x++ {
			appendFile(filename, floatToStringENotation(Ux[ARM(x, y)])+" "+floatToStringENotation(Uy[ARM(x, y)])+" 0.0\n")
		}
	}
	appendFile(filename, "\n")
}

func Handle() {
	iter := 0
	var MaxIter int = 1e4
	var rho0LB, u0LB, v0LB, nuLB float64
	var Re, Ma, N float64

	Ma = 0.1  // Ma = U/c_s
	Re = 1000 // Re = UL/nu
	N = LX

	// IC setup
	rho0LB = 1.
	u0LB = Ma * math.Sqrt(csSq)
	v0LB = u0LB

	setIC(rho0LB, u0LB, v0LB)
	dumpStateVTK("state0.vtk")

	// PHYS setup
	nuLB = u0LB * N / Re

	// LBM setup
	tau = nuLB/csSq + 0.5

	fmt.Printf(" nu_LB: %f\n tau: %f\n u0LB: %f\n", nuLB, tau, u0LB)

	for {
		if iter%2 != 0 {
			collideStream(tmpCells[:], cells[:])
		} else {
			collideStream(cells[:], tmpCells[:])
		}
		if iter%1000 == 0 {
			fmt.Printf("#\n")
		}

		iter++
		if iter >= MaxIter {
			break
		}
	}

	dumpStateVTK("state.vtk")

}
