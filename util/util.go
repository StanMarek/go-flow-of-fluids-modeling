package util

import (
	"fmt"
	"os"
	"strconv"
)

func GetFromIndex(x, y, f, F, LX int) int {
	return f + F*x + F*LX*y
}

func Check(err error) {
	if err != nil {
		panic(err)
	}
}

func IntToString(value interface{}) string {
	return strconv.FormatInt(int64(value.(int)), 10)
}

func FloatToStringENotation(value interface{}) string {
	return fmt.Sprintf("%e", value.(float64))
}

func AppendFile(filename string, data string) {
	file, err := os.OpenFile(filename, os.O_APPEND|os.O_WRONLY, 0600)
	Check(err)

	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			panic(err)
		}
	}(file)

	_, err = file.WriteString(data)
	Check(err)
}
