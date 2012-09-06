// This file implements SeqPop, which is a population of full sequence genomes.
// To use it:
// 1. create a SeqPop object using NewSeqPop by providing population Size, genome Length, Mutation rate, Transfer rate, Fragment length, and a LockedSource.
// 2. do the evoluation using Evolve
package fwd1

import (
	"github.com/mingzhi/gomath/random"
	"math/rand"
)

// SeqPop: represents a population with full sequence
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

type Sequence []byte

// NewSeqPop: returns a SeqPop.
func NewSeqPop(size, length int, mutation, transfer float64, fragment int, src rand.Source) (pop *SeqPop) {
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
	pop.src = src
	pop.rng = rand.New(src)
	// mutation poisson distribution: mean = mutation rate * genome length
	pop.mPois = random.NewPoisson(pop.Mutation*float64(pop.Length), src)
	// transfer poisson distribution: mean = transfer rate * genome length
	pop.tPois = random.NewPoisson(pop.Transfer*float64(pop.Length), src)

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
	ch := make(chan bool)          // channel to communicate with workers
	n := 0                         // number of genomes for hard copying
	used := make([]bool, pop.Size) // indicates whether the pointer of a old genome has been used
	for i := 0; i < pop.Size; i++ {
		o := pop.rng.Intn(pop.Size) // sample with replacement
		if used[o] {
			go pop.copyGenome(i, o, newGenomes, ch)
			n++
		} else {
			newGenomes[i] = pop.Genomes[o]
			used[i] = true
		}
	}

	// waiting for workers finishing hard genome copy
	for i := 0; i < n; i++ {
		<-ch
	}

	// update genomes
	pop.Genomes = newGenomes
}

// hard copy genome (worker)
func (pop *SeqPop) copyGenome(n, o int, newGenomes []Sequence, ch chan bool) {
	newGenomes[n] = make(Sequence, pop.Length)
	copy(newGenomes[n], pop.Genomes[o])
	// tell others that I am done
	ch <- true
}

// mutation: K2P model.
// Each genome is independent with each other so that concurrency is implemented.
func (pop *SeqPop) mutation() {
	// in case that the mutation rate is changed
	mean := pop.Mutation * float64(pop.Length)
	if pop.mPois.Mean != mean {
		pop.mPois.Mean = mean
	}
	ch := make(chan bool) // channel for communicating with workers
	// for each genome, do mutation
	for i := 0; i < pop.Size; i++ {
		go pop.mutateGenome(i, ch)
	}

	// waiting for worker finishing doing mutation
	for i := 0; i < pop.Size; i++ {
		<-ch
	}
}

// mutate one genome (atomic job)
func (pop *SeqPop) mutateGenome(idx int, ch chan bool) {
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
	// tell others that I am done
	ch <- true
}

// transfer
// Warning: transaction (read and write) is not atomic and not safe.
func (pop *SeqPop) transfer() {
	// in case that the transfer rate is changed
	mean := pop.Transfer * float64(pop.Length)
	if pop.tPois.Mean != mean {
		pop.tPois.Mean = mean
	}

	ch := make(chan bool) // channel for communicating with other workers
	// for each genome, do transfer
	for i := 0; i < pop.Size; i++ {
		go pop.transferGenome(i, ch)
	}

	// waiting for worker finishing their jobs
	for i := 0; i < pop.Size; i++ {
		<-ch
	}
}

// transfer one genome (receiver)
func (pop *SeqPop) transferGenome(idx int, ch chan bool) {
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

	ch <- true
}
