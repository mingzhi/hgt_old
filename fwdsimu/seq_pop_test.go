package fwdsimu

import (
	"testing"
)

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
