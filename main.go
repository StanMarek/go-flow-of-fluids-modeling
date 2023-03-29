package main

import (
	"fmt"
	"os"
	"strconv"
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

func getFromIndex(x, y, f int) int {
	return f + F*x + F*LX*y
}

func feq(i int, phi, u, v float64) float64 {
	return w[i] * phi * (1 + (cx[i]*u+cy[i]*v)/csSq)
}

func setIC(phi0, u0, v0 float64) {
	for y := 0; y < LY; y++ {
		for x := 0; x < LX; x++ {
			for f := 0; f < F; f++ {
				cells[getFromIndex(x, y, f)] = feq(f, phi0, u0, v0)
			}
		}
	}

	for i := 0; i < F; i++ {
		cells[getFromIndex(LX/2, LY/2, i)] = feq(i, 1., u0, v0)
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
				phi += c[getFromIndex(x, y, i)]
			}

			for i = 0; i < F; i++ {
				fStar[i] = (1-1/tau)*c[getFromIndex(x, y, i)] + feq(i, phi, uLB, vLB)
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

			tc[getFromIndex(x, y, _C)] = fStar[_C]
			tc[getFromIndex(x, yp, _N)] = fStar[_N]
			tc[getFromIndex(x, ym, _S)] = fStar[_S]
			tc[getFromIndex(xm, y, _W)] = fStar[_W]
			tc[getFromIndex(xp, y, _E)] = fStar[_E]
		}
	}
}

func check(err error) {
	if err != nil {
		panic(err)
	}
}

func intToString(value interface{}) string {
	return strconv.FormatInt(int64(value.(int)), 10)
}

func floatToStringENotation(value interface{}) string {
	return fmt.Sprintf("%e", value.(float64))
}

func appendFile(filename string, data string) {
	file, err := os.OpenFile(filename, os.O_APPEND|os.O_WRONLY, 0600)
	check(err)

	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			panic(err)
		}
	}(file)

	_, err = file.WriteString(data)
	check(err)
}

func dumpStateVTK(filename string) {
	file, err := os.Create(filename)
	check(err)

	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			panic(err)
		}
	}(file)

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
				rho += cells[getFromIndex(x, y, f)]
			}
			appendFile(filename, floatToStringENotation(rho)+" \n")
			check(err)
		}
		appendFile(filename, "\n")
	}
}

func main() {
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
