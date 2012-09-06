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

func (cm *CMatrix) KS() (ks float64) {
	mean := desc.NewMean()
	for i := 0; i < len(cm.Matrix); i++ {
		m := float64(len(cm.Matrix[i])) / float64(cm.Length)
		mean.Increment(m)
	}

	ks = mean.GetResult()
	return
}

func (cm *CMatrix) D() (m float64, v float64) {
	mean := desc.NewMean()
	variance := desc.NewVariance()
	for i := 0; i < len(cm.Matrix); i++ {
		m := float64(len(cm.Matrix[i])) / float64(cm.Length)
		mean.Increment(m)
		variance.Increment(m)
	}

	m = mean.GetResult()
	v = variance.GetResult()
	return
}
