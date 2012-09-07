package fwd2

import (
	"testing"
)

func BenchmarkGeneration(b *testing.B) {
	b.StopTimer()
	size := 1000
	length := 1000
	mutation := 0.00001
	transfer := 0.0
	fragment := 0
	pop := NewSeqPop(size, length, mutation, transfer, fragment)
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		pop.generation()
	}
}
