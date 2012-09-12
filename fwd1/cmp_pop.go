// This file implements CmpPop, which is a population of compacted sequence genomes.
// To use it:
// 1. create a CmpPop object using NewSeqPop by providing population Size, genome Length, Mutation rate, Transfer rate, Fragment length, and a LockedSource.
// 2. do the evoluation using Evolve
package fwd1

import (
	"github.com/mingzhi/gomath/random"
	"math/rand"
	"time"
)

// CmpPop: represents a population with compacted geneomes
type CmpPop struct {
	Size     int           // Population size
	Length   int           // Genome length
	Fragment int           // Transferred fragment length
	Mutation float64       // Mutation rate
	Transfer float64       // Transfer rate
	Genomes  []CmpSequence // genomes containing compact sequences.
	NumOfGen int           // Number of generation

	src   rand.Source     // random source
	rng   *rand.Rand      // random number generator
	mPois *random.Poisson // poisson random generator (mutator)
	tPois *random.Poisson // poisson random generator (transfer)
}

// MutatedSite: represent a mutated site in a sequence
type MutatedSite struct {
	Position int // Position of this mutation in the sequence
	Time     int // Time when this mutation happens
}

// CmpSequence: represent a compacted sequence
// which contains mutations
type CmpSequence struct {
	MutatedSites []MutatedSite
}

// NewCmpPop: return a CmpPop
func NewCmpPop(size, length int, mutation, transfer float64, fragment int) *CmpPop {
	// verify size and length
	if size <= 0 || length <= 0 {
		panic("Population size or length <= 0!")
	}

	// create a CmpPop
	pop := &CmpPop{
		Size:     size,
		Length:   length,
		Mutation: mutation,
		Transfer: transfer,
		Fragment: fragment,
	}

	// create random sources
	seed := time.Now().UnixNano()                          // use time now as a seed.
	pop.src = random.NewLockedSource(rand.NewSource(seed)) // create a locked source, because we want concurrency.
	pop.rng = rand.New(pop.src)                            // create a random number generator with locked source.
	mMean := pop.Mutation * float64(pop.Length)            // mutation rate of a genome
	pop.mPois = random.NewPoisson(mMean, pop.src)          // poisson random generator for mutation process
	tMean := pop.Transfer * float64(pop.Length)            // transfer rate of a genome
	pop.tPois = random.NewPoisson(tMean, pop.src)          // poisson random generator for transfer process

	// create genomes
	pop.Genomes = make([]CmpSequence, pop.Size)

	return pop
}

// generation: wright-fisher generate a fixed number of population.
func (pop *CmpPop) generation() {
	newGenomes := make([]CmpSequence, pop.Size) // new genomes
	for i := 0; i < pop.Size; i++ {
		o := pop.rng.Intn(pop.Size)    // randomly determine a parent.
		newGenomes[i] = pop.Genomes[o] // assign this parent to the child, do not need hard coping.
	}
}

// mutation: K2P model.
func (pop *CmpPop) mutation() {
	ch := make(chan bool) // create a commutation channel
	for i := 0; i < pop.Size; i++ {
		go pop.mutateGenome(i, ch)
	}

	// waiting for workers finishing their job
	for i := 0; i < pop.Size; i++ {
		<-ch
	}
}

// mutateGenome: mutate one genome
func (pop *CmpPop) mutateGenome(idx int, ch chan bool) {
	count := pop.mPois.Int() // determine the number of mutations happening in the genome
	for i := 0; i < count; i++ {
		p := pop.rng.Intn(pop.Length)                                            // randomly select a position
		t := pop.NumOfGen                                                        // time when this mutation happens
		m := MutatedSite{Position: p, Time: t}                                   // create a mutated site
		pop.Genomes[idx].MutatedSites = append(pop.Genomes[idx].MutatedSites, m) // add the new mutated site to the genome
	}

	ch <- true // tell others that I am done.
}
