package covs

import (
	"github.com/mingzhi/gomath/stat/desc"
)

type CMatrix struct {
	Size   int
	Length int
	Matrix [][]int
}

func NewCMatrix(size, length int, matrix [][]int) (cmatrix *CMatrix) {
	return &CMatrix{Size: size, Length: length, Matrix: matrix}
}

// calculate the KS and VarD.
func (cm *CMatrix) D() (m float64, v float64) {
	mean := desc.NewMean()
	variance := desc.NewVariance()
	for i := 0; i < len(cm.Matrix); i++ {
		d := float64(len(cm.Matrix[i])) / float64(cm.Length)
		mean.Increment(d)
		variance.Increment(d)
	}

	m = mean.GetResult()
	v = variance.GetResult()
	return
}
