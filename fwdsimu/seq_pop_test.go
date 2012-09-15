package fwdsimu

import (
	"fmt"
	"github.com/mingzhi/gomath/stat/desc"
	"github.com/mingzhi/hgt/covs"
	"math"
	"testing"
)

func TestEvolve(t *testing.T) {
	// population parameter
	var (
		numofgens, samplesize  int
		size, length, fragment int
		mutation, transfer     float64
		pop                    *SeqPop
		ksmean                 *desc.Mean
		ksvar                  *desc.Variance
		expectedKS, tolerance  float64
	)
	samplesize = 100
	size = 1000
	length = 1000
	mutation = 1e-4
	fragment = 100
	fmt.Println("--Population parameters:")
	fmt.Printf("---Size: %d\n", size)
	fmt.Printf("---Genome Length: %d\n", length)
	fmt.Printf("---Mutation rate: %g\n", mutation)
	fmt.Printf("---Transferred fragment: %d\n", fragment)

	fmt.Println("--Test no horizontal transfer")
	transfer = 0
	pop = NewSeqPop(size, length, mutation, transfer, fragment)
	ksmean = desc.NewMean()
	ksvar = desc.NewVariance()
	numofgens = 100
	for {
		for i := 0; i < numofgens; i++ {
			pop.Evolve()
			sample := Sample(pop.Genomes, samplesize)
			dmatrix := GenerateDistanceMatrix(sample)
			cmatrix := covs.NewCMatrix(samplesize, length, dmatrix)
			ks, _ := cmatrix.D()
			ksmean.Increment(ks)
			ksvar.Increment(ks)
		}
		expectedKS = covs.CalculatKS(size, 4, mutation, transfer, fragment)
		tolerance = math.Sqrt(ksvar.GetResult())
		if math.Abs(ksmean.GetResult()-expectedKS) > tolerance {
			if numofgens < 10000 {
				numofgens *= 10
				fmt.Printf("try %d generation again ...\n", numofgens)
			} else {
				t.Errorf("ks = %g, with tolerance = %g, expect %g\n", ksmean.GetResult(), tolerance, expectedKS)
			}
		} else {
			break
		}
	}

	fmt.Println("--Test horizontal transfer with same rate as mutation")
	transfer = mutation
	pop = NewSeqPop(size, length, mutation, transfer, fragment)
	ksmean = desc.NewMean()
	ksvar = desc.NewVariance()
	numofgens = 100
	for {
		for i := 0; i < numofgens; i++ {
			pop.Evolve()
			sample := Sample(pop.Genomes, samplesize)
			dmatrix := GenerateDistanceMatrix(sample)
			cmatrix := covs.NewCMatrix(samplesize, length, dmatrix)
			ks, _ := cmatrix.D()
			ksmean.Increment(ks)
			ksvar.Increment(ks)
		}
		expectedKS = covs.CalculatKS(size, 4, mutation, transfer, fragment)
		tolerance = math.Sqrt(ksvar.GetResult())
		if math.Abs(ksmean.GetResult()-expectedKS) > tolerance {
			if numofgens < 10000 {
				numofgens *= 10
				fmt.Printf("try %d generation again ...\n", numofgens)
			} else {
				t.Errorf("ks = %g, with tolerance = %g, expect %g\n", ksmean.GetResult(), tolerance, expectedKS)
			}
		} else {
			break
		}
	}

	fmt.Println("--Test horizontal transfer with 100 times of mutation rate")
	transfer = mutation * 100.0
	pop = NewSeqPop(size, length, mutation, transfer, fragment)
	ksmean = desc.NewMean()
	ksvar = desc.NewVariance()
	numofgens = 100
	for {
		for i := 0; i < numofgens; i++ {
			pop.Evolve()
			sample := Sample(pop.Genomes, samplesize)
			dmatrix := GenerateDistanceMatrix(sample)
			cmatrix := covs.NewCMatrix(samplesize, length, dmatrix)
			ks, _ := cmatrix.D()
			ksmean.Increment(ks)
			ksvar.Increment(ks)
		}
		expectedKS = covs.CalculatKS(size, 4, mutation, transfer, fragment)
		tolerance = math.Sqrt(ksvar.GetResult())
		if math.Abs(ksmean.GetResult()-expectedKS) > tolerance {
			if numofgens < 10000 {
				numofgens *= 10
				fmt.Printf("try %d generation again ...\n", numofgens)
			} else {
				t.Errorf("ks = %g, with tolerance = %g, expect %g\n", ksmean.GetResult(), tolerance, expectedKS)
			}
		} else {
			break
		}
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
