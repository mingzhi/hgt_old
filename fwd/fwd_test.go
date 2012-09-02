package fwd

import (
	"github.com/mingzhi/gomath/random"
	"github.com/mingzhi/hgt/covs"
	"math/rand"
	"runtime"
	"testing"
)

func TestEvolve(t *testing.T) {
	runtime.GOMAXPROCS(runtime.NumCPU())
	size := 10000
	length := 10000
	mutation := 0.00001
	transfer := 0.0
	fragment := 0
	time := 1000
	src := random.NewLockedSource(rand.NewSource(1))
	sampleSize := 100

	pop := NewSeqPop(size, length, mutation, transfer, fragment, src)

	pop.Evolve(time)

	samples := pop.Sample(sampleSize)
	dmatrix := GenerateDistanceMatrix(samples)
	cmatrix := covs.NewCMatrix(sampleSize, length, dmatrix)
	ks := cmatrix.KS()
	t.Error(ks)
}
