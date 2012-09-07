package fwd1

import (
	"github.com/mingzhi/gomath/random"
	"math/rand"
	"testing"
)

func BenchmarkGeneration(b *testing.B) {
	b.StopTimer()
	size := 1000
	length := 1000
	mutation := 0.00001
	transfer := 0.0
	fragment := 0
	src := random.NewLockedSource(rand.NewSource(1))
	pop := NewSeqPop(size, length, mutation, transfer, fragment, src)
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		pop.generation()
	}
}
