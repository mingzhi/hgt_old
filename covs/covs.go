package covs

import (
	"github.com/mingzhi/gomath/stat/desc"
	"math"
	"runtime"
	"sort"
)

type CMatrix struct {
	Size   int
	Length int
	Matrix [][]int
}

func NewCMatrix(size, length int, matrix [][]int) (cmatrix *CMatrix) {
	return &CMatrix{Size: size, Length: length, Matrix: matrix}
}

// calculate KS (math formula)
func CalculateKS(n int, u float64, r float64, l int, a int) float64 {
	ut := 2.0 * (u + r*float64(l)) // probability of total events
	ratio := float64(l) * r / u    // ratio between transfer and mutation rate
	ustar := 1.0 - math.Exp(-ut*(1.0+1.0/(float64(a-1)*(1.0+ratio))))
	return float64(n) * ustar / (ratio + float64(a)/float64(a-1)*(float64(n)*ustar+1.0-ustar))
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

// Calculate the covs.
func (cm *CMatrix) Cov(maxL int) (scovs, rcovs, xyPL, xsysPL, smXYPL []float64) {
	// calculate xs and xy
	xs := make([]int, cm.Length) // total counts of each column.
	xy := make([]int, maxL)      // total counts of cocurrence.
	// number of cpu
	ncpu := runtime.GOMAXPROCS(0)
	ch := make(chan result)
	for i := 0; i < ncpu; i++ {
		begin := i * cm.Size / ncpu
		end := (i + 1) * cm.Size / ncpu
		go cm.covPart(begin, end, maxL, ch)
	}
	// collect results
	for i := 0; i < ncpu; i++ {
		r := <-ch
		for j := 0; j < cm.Length; j++ {
			xs[j] += r.xs[j]
		}
		for j := 0; j < maxL; j++ {
			xy[j] += r.xy[j]
		}
	}

	// calculate xsP (frequency of xs)
	xsP := make([]float64, len(xs))
	for i := 0; i < len(xs); i++ {
		xsP[i] = float64(xs[i]) / float64(cm.Size)
	}

	// calculate xyP (frequence of xy)
	xyP := make([]float64, len(xy))
	for i := 0; i < len(xy); i++ {
		xyP[i] = float64(xy[i]) / float64(cm.Size*cm.Length)
	}

	xsysP := make([]float64, maxL) // <X.Y>
	smXsP := make([]float64, maxL) // sum(X)
	smYsP := make([]float64, maxL) // sum(Y)
	for l := 0; l < maxL; l++ {
		for x := 0; x < cm.Length; x++ {
			if x+l >= cm.Length {
				xsysP[l] += xsP[x] * xsP[l-(cm.Length-x)]
				smXsP[l] += xsP[x]
				smYsP[l] += xsP[l-(cm.Length-x)]
			} else {
				xsysP[l] += xsP[x] * xsP[x+l]
				smXsP[l] += xsP[x]
				smYsP[l] += xsP[x+l]
			}
		}
		xsysP[l] /= float64(cm.Length)
		smXsP[l] /= float64(cm.Length)
		smYsP[l] /= float64(cm.Length)
	}

	scovs = make([]float64, maxL)
	rcovs = make([]float64, maxL)
	xyPL = make([]float64, maxL)
	xsysPL = make([]float64, maxL)
	smXYPL = make([]float64, maxL)
	for l := 0; l < maxL; l++ {
		scovs[l] = xyP[l] - xsysP[l]
		rcovs[l] = xsysP[l] - smXsP[l]*smYsP[l]
		xyPL[l] = xyP[l]
		xsysPL[l] = xsysP[l]
		smXYPL[l] = smXsP[l] * smYsP[l]
	}

	return
}

type result struct {
	xs []int
	xy []int
}

func (cm *CMatrix) covPart(begin, end, maxL int, ch chan result) {
	// calculate xs and xy
	xs := make([]int, cm.Length) // total counts of each column.
	xy := make([]int, maxL)      // total counts of cocurrence.
	// sort each row
	for i := begin; i < end; i++ {
		sort.Ints(cm.Matrix[i])
	}
	// loop
	for i := begin; i < end; i++ {
		row := cm.Matrix[i]
		for xi, xv := range row {
			// add xs
			xs[xv]++
			// add xy
			for yi := 0; yi < len(row); {
				yv := row[yi]
				if yv < xv { // deal with circle genome
					if yv+cm.Length-xv < maxL {
						xy[yv+cm.Length-xv]++
						yi++
					} else {
						yi = xi
					}
				} else {
					if yv-xv < maxL {
						xy[yv-xv]++
						yi++
					} else {
						break
					}
				}
			}
		}
	}

	ch <- result{xs: xs, xy: xy}
}
