// This file implements SeqPop, which is a population of full sequence genomes.
// To use it:
// 1. create a SeqPop object using NewSeqPop by providing population Size, genome Length, Mutation rate, Transfer rate, Fragment length, and a LockedSource.
// 2. do the evoluation using Evolve
package fwd1

import (
	"github.com/mingzhi/gomath/random"
	"math/rand"
	"time"
)

// SeqPop: represents a population with full sequence genomes
type SeqPop struct {
	Size     int        // Population size
	Length   int        // Genome length
	Mutation float64    // Mutation rate
	Transfer float64    // Transfer rate
	Fragment int        // Transferred fragment
	NumOfGen int        // Number of generation
	Genomes  []Sequence // Store genomes

	// private variables
	src    rand.Source     // random sources (locked source)
	rng    *rand.Rand      // random number generator with locked source
	mPois  *random.Poisson // poisson random generator (mutator)
	tPois  *random.Poisson // poisson random generator (transfer)
	states string          // nucleotides
}

// Genome sequence
type Sequence []byte

// NewSeqPop: returns a SeqPop.
func NewSeqPop(size, length int, mutation, transfer float64, fragment int) (pop *SeqPop) {
	// verify size and length
	if size <= 0 || length <= 0 {
		panic("Population size or length should be positive values!")
	}

	// construct a SeqPop instance
	pop = &SeqPop{
		Size:     size,
		Length:   length,
		Mutation: mutation,
		Transfer: transfer,
		Fragment: fragment,
	}

	// create random sources
	seed := time.Now().UnixNano()
	pop.src = random.NewLockedSource(rand.NewSource(seed))
	pop.rng = rand.New(pop.src)
	// mutation poisson distribution: mean = mutation rate * genome length
	pop.mPois = random.NewPoisson(pop.Mutation*float64(pop.Length), pop.src)
	// transfer poisson distribution: mean = transfer rate * genome length
	pop.tPois = random.NewPoisson(pop.Transfer*float64(pop.Length), pop.src)

	// nucleotides
	pop.states = "ATGC"
	// create a random reference sequence
	refseq := make(Sequence, pop.Length)
	for i := 0; i < pop.Length; i++ {
		s := pop.states[pop.rng.Intn(len(pop.states))]
		refseq[i] = s
	}
	// create a population size of genomes by coping the reference sequence
	pop.Genomes = make([]Sequence, pop.Size)
	for i := 0; i < pop.Size; i++ {
		genome := make(Sequence, pop.Length)
		copy(genome, refseq)
		pop.Genomes[i] = genome
	}

	return
}

// Evolve: perform a step of evolution, including:
// 1. generation (wright-fisher model)
// 2. mutation (K2P model, nucleotide frequencies equal)
// 3. homologous recombination with provided rate and fragment length.
func (pop *SeqPop) Evolve() {
	pop.generation()
	pop.mutation()
	pop.transfer()
	pop.NumOfGen++
}

// generation: wright-fisher generate a fixed number of population.
func (pop *SeqPop) generation() {
	// new generation
	newGenomes := make([]Sequence, pop.Size)

	// wright-fisher sampling
	used := make([]bool, pop.Size) // indicates whether the pointer of a old genome has been used
	for i := 0; i < pop.Size; i++ {
		o := pop.rng.Intn(pop.Size) // sample with replacement
		if used[o] {
			newGenomes[i] = make(Sequence, pop.Length)
			copy(newGenomes[i], pop.Genomes[o])
		} else {
			newGenomes[i] = pop.Genomes[o]
			used[o] = true
		}
	}

	// update genomes
	pop.Genomes = newGenomes
}

// mutation: K2P model.
// Each genome is independent with each other so that concurrency is implemented.
func (pop *SeqPop) mutation() {
	// in case that the mutation rate is changed
	mean := pop.Mutation * float64(pop.Length)
	if pop.mPois.Mean != mean {
		pop.mPois.Mean = mean
	}

	// for each genome, do mutation
	for i := 0; i < pop.Size; i++ {
		pop.mutateGenome(i)
	}
}

// mutate one genome (atomic job)
func (pop *SeqPop) mutateGenome(idx int) {
	genome := pop.Genomes[idx]
	// how many mutations in this genome (following a poisson distribution)
	count := pop.mPois.Int()
	// for each mutation, randomly choose a position and replace it with another nucleotide.
	for i := 0; i < count; i++ {
		p := pop.rng.Intn(pop.Length)
		s := pop.states[pop.rng.Intn(len(pop.states))]
		for s == genome[p] { // loop until a random nucleotide that is different to the original one
			s = pop.states[pop.rng.Intn(len(pop.states))]
		}
		genome[p] = s
	}
}

// transfer
// Warning: transaction (read and write) is not atomic and not safe.
func (pop *SeqPop) transfer() {
	// in case that the transfer rate is changed
	mean := pop.Transfer * float64(pop.Length)
	if pop.tPois.Mean != mean {
		pop.tPois.Mean = mean
	}

	// for each genome, do transfer
	for i := 0; i < pop.Size; i++ {
		pop.transferGenome(i)
	}
}

// transfer one genome (receiver)
func (pop *SeqPop) transferGenome(idx int) {
	receiver := pop.Genomes[idx]

	// how many times of homologous recombination happen.
	count := pop.tPois.Int()
	// for each transfer, a donor genome is randomly chose, then a fragment is randomly determined.
	for c := 0; c < count; c++ {
		// a different donor genome is randomly selected.
		donor_id := pop.rng.Intn(pop.Size)
		for donor_id == idx {
			donor_id = pop.rng.Intn(pop.Size)
		}
		donor := pop.Genomes[donor_id]

		// a fragment is randomly determined
		fragL := pop.rng.Intn(pop.Length)
		fragR := fragL + pop.Fragment
		if fragR < pop.Length {
			for i := fragL; i < fragR; i++ {
				receiver[i] = donor[i]
			}
		} else {
			for i := fragL; i < pop.Length; i++ {
				receiver[i] = donor[i]
			}
			for i := 0; i < fragR-pop.Length; i++ {
				receiver[i] = donor[i]
			}
		}
	}
}
