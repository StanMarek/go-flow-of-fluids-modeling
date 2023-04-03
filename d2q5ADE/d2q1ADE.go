package d2q5ADE

import (
	"mppwc/util"
	"os"
)

const (
	LX   = 128
	LY   = 128
	F    = 5
	wC   = 1. / 3.
	wNC  = 1. / 6.
	csSq = 1. / 3.
)

var (
	tau float64
	uLB float64
	vLB float64

	cx = [F]float64{0., 0., 0., -1., 1.}
	cy = [F]float64{0., 1., -1., 0., 0.}
	w  = [F]float64{wC, wNC, wNC, wNC, wNC}

	cells    = [LX * LY * F]float64{}
	tmpCells = [LX * LY * F]float64{}
)

const (
	_C int = 0
	_N int = 1
	_S int = 2
	_W int = 3
	_E int = 4
)

func AR(x, y, f int, getFromIndex func(x, y, f, F, LX int) int) int {
	return getFromIndex(x, y, f, F, LX)
}

var AID = util.GetFromIndex

func feq(i int, phi, u, v float64) float64 {
	return w[i] * phi * (1 + (cx[i]*u+cy[i]*v)/csSq)
}

func setIC(phi0, u0, v0 float64) {
	for y := 0; y < LY; y++ {
		for x := 0; x < LX; x++ {
			for f := 0; f < F; f++ {
				cells[AR(x, y, f, AID)] = feq(f, phi0, u0, v0)
			}
		}
	}

	for i := 0; i < F; i++ {
		cells[AR(LX/2, LY/2, i, AID)] = feq(i, 1., u0, v0)
	}
}

func collideStream(c, tc []float64) {
	var phi float64
	var fStar = [F]float64{}
	var i, x, y, xp, xm, yp, ym int

	for y = 0; y < LY; y++ {
		for x = 0; x < LX; x++ {
			phi = 0
			for i = 0; i < F; i++ {
				phi += c[AR(x, y, i, AID)]
			}

			for i = 0; i < F; i++ {
				fStar[i] = (1-1/tau)*c[AR(x, y, i, AID)] + feq(i, phi, uLB, vLB)
			}

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

			tc[AR(x, y, _C, AID)] = fStar[_C]
			tc[AR(x, yp, _N, AID)] = fStar[_N]
			tc[AR(x, ym, _S, AID)] = fStar[_S]
			tc[AR(xm, y, _W, AID)] = fStar[_W]
			tc[AR(xp, y, _E, AID)] = fStar[_E]
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

	appendFile(filename, "# vtk DataFile Version 3.0\n")
	appendFile(filename, "2D-ADE data file \n")
	appendFile(filename, "ASCII\n")
	appendFile(filename, "DATASET RECTILINEAR_GRID\n")
	appendFile(filename, "DIMENSIONS "+intToString(LX)+" "+intToString(LY)+" 1\n")
	appendFile(filename, "X_COORDINATES "+intToString(LX)+" int\n")
	for i := 0; i < LX; i++ {
		appendFile(filename, intToString(i)+" ")
	}
	appendFile(filename, "\n")
	appendFile(filename, "Y_COORDINATES "+intToString(LY)+" int\n")
	for i := 0; i < LY; i++ {
		appendFile(filename, intToString(i)+" ")
	}
	appendFile(filename, "\n")
	appendFile(filename, "Z_COORDINATES 1 int\n")
	appendFile(filename, "0\n")
	appendFile(filename, "POINT_DATA "+intToString(LX*LY)+"\n")
	appendFile(filename, "SCALARS temperature double 1\n")
	appendFile(filename, "LOOKUP_TABLE default\n")
	for y := 0; y < LY; y++ {
		for x := 0; x < LX; x++ {
			rho := 0.
			var f int
			for f = 0; f < F; f++ {
				rho += cells[AR(x, y, f, AID)]
			}
			appendFile(filename, floatToStringENotation(rho)+" \n")
			util.Check(err)
		}
		appendFile(filename, "\n")
	}
}

func Handle() {
	iter := 0
	var MaxIter int = 1e3
	phi0LB := 0.
	u0LB := 0.
	v0LB := 0.
	alphaLB := 0.01

	setIC(phi0LB, u0LB, v0LB)
	dumpStateVTK("state0.vtk")

	tau = alphaLB/csSq + 0.5

	for {
		if iter%2 == 0 {
			collideStream(cells[:], tmpCells[:])
		} else {
			collideStream(tmpCells[:], cells[:])
		}
		iter++
		if iter >= MaxIter {
			break
		}
	}

	dumpStateVTK("state.vtk")
}
