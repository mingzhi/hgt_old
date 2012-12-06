package covs

import (
	"math/rand"
	"testing"
)

func TestCircleCov(t *testing.T) {
	// create testing distance sequences
	size := 10
	leng := 1000
	seqs := make([][]int, size)
	dmatrix := make([][]int, size)
	for i := 0; i < len(seqs); i++ {
		seqs[i] = make([]int, leng)
		dmatrix[i] = []int{}
		for j := 0; j < len(seqs[i]); j++ {
			d := rand.Intn(2)
			seqs[i][j] = d
			if d == 1 {
				dmatrix[i] = append(dmatrix[i], j)
			}
		}
	}

	maxl := 100
	n_xy := make([]int, maxl)
	for l := 0; l < maxl; l++ {
		for i := 0; i < len(seqs); i++ {
			for j := 0; j < len(seqs[i]); j++ {
				if j+l < len(seqs[i]) {
					if seqs[i][j] == 1 && seqs[i][j+l] == 1 {
						n_xy[l]++
					}
				} else {
					if seqs[i][j] == 1 && seqs[i][j+l-len(seqs[i])] == 1 {
						n_xy[l]++
					}
				}

			}
		}
	}

	p_xy := make([]float64, maxl)
	for i := 0; i < maxl; i++ {
		p_xy[i] = float64(n_xy[i]) / float64(size*leng)
	}
	cmatrix := NewCMatrix(size, leng, dmatrix)
	_, _, xy, _, _ := cmatrix.CovCircle(maxl)
	for i := 0; i < maxl; i++ {
		if p_xy[i] != xy[i] {
			t.Errorf("%d: p_xy = %g, xy = %g", i, p_xy[i], xy[i])
		}
	}

}
