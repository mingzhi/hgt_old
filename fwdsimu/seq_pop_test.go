package fwdsimu

import (
	"github.com/mingzhi/gomath/stat/desc"
	"github.com/mingzhi/hgt/covs"
	"math"
	"testing"
)

func TestEvolve(t *testing.T) {
	// population parameter
	numofgens := 10000
	samplesize := 100
	size := 1000
	length := 1000
	mutation := 1e-4
	fragment := 100
	// no transfer
	transfer := 0.0

	// construct population
	pop := NewSeqPop(size, length, mutation, transfer, fragment)
	ksmean := desc.NewMean()
	ksvar := desc.NewVariance()
	for i := 0; i < numofgens; i++ {
		pop.Evolve()
		if pop.NumOfGens > 1000 {
			sample := Sample(pop.Genomes, samplesize)
			dmatrix := GenerateDistanceMatrix(sample)
			cmatrix := covs.NewCMatrix(samplesize, length, dmatrix)
			ks, _ := cmatrix.D()
			ksmean.Increment(ks)
			ksvar.Increment(ks)
		}
	}
	expectKS := covs.CalculatKS(size, 4, mutation, transfer, fragment)
	if math.Abs(ksmean.GetResult()-expectKS) > ksvar.GetResult() {
		t.Errorf("ks = %g, but expect %g", ksmean.GetResult(), expectKS)
	}
}

func BenchmarkEvolve(b *testing.B) {
	b.StopTimer()
	size := 1000
	length := 1000
	mutation := 1e-4
	transfer := 1e-4
	fragment := 100
	pop := NewSeqPop(size, length, mutation, transfer, fragment)
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		pop.Evolve()
	}
}

func BenchmarkGeneration(b *testing.B) {
	b.StopTimer()
	size := 1000
	length := 1000
	mutation := 1e-4
	transfer := 1e-4
	fragment := 100
	pop := NewSeqPop(size, length, mutation, transfer, fragment)
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		pop.reproduce()
	}
}

func BenchmarkManipulation(b *testing.B) {
	b.StopTimer()
	size := 1000
	length := 1000
	mutation := 1e-3
	transfer := 1e-3
	fragment := 100
	pop := NewSeqPop(size, length, mutation, transfer, fragment)
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		pop.manipulate()
	}
}
