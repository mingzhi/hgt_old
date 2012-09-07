package fwd2

import (
	"github.com/mingzhi/gsl-cgo/randist"
)

type SeqPop struct {
	Size     int
	Length   int
	Fragment int
	Mutation float64
	Transfer float64
	NumOfGen int

	Genomes []Sequence

	// private variables
	states string
	rng    *randist.RNG
}

type Sequence []byte

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

	// construct random source
	pop.rng = randist.NewRNG(randist.MT19937)

	// construct genomes
	pop.states = "ATGC"
	refseq := make(Sequence, pop.Length)
	for i := 0; i < pop.Length; i++ {
		s := pop.states[randist.UniformRandomInt(pop.rng, len(pop.states))]
		refseq[i] = s
	}
	pop.Genomes = make([]Sequence, pop.Size)
	for i := 0; i < pop.Size; i++ {
		pop.Genomes[i] = make(Sequence, pop.Length)
		copy(pop.Genomes[i], refseq)
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

func (pop *SeqPop) generation() {
	// new generation
	newGenomes := make([]Sequence, pop.Size)

	// wright-fisher sampling
	ch := make(chan bool)          // channel to communicate with workers
	n := 0                         // number of genomes for hard copying
	used := make([]bool, pop.Size) // indicates whether the pointer of a old genome has been used
	for i := 0; i < pop.Size; i++ {
		o := randist.UniformRandomInt(pop.rng, pop.Size) // sample with replacement
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
	count := randist.PoissonRandomInt(pop.rng, pop.Mutation*float64(pop.Length))
	// for each mutation, randomly choose a position and replace it with another nucleotide.
	for i := 0; i < count; i++ {
		p := randist.UniformRandomInt(pop.rng, pop.Length)
		s := pop.states[randist.UniformRandomInt(pop.rng, len(pop.states))]
		for s == genome[p] { // loop until a random nucleotide that is different to the original one
			s = pop.states[randist.UniformRandomInt(pop.rng, len(pop.states))]
		}
		genome[p] = s
	}
	// tell others that I am done
	ch <- true
}

// transfer
// Warning: transaction (read and write) is not atomic and not safe.
func (pop *SeqPop) transfer() {
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
	count := randist.PoissonRandomInt(pop.rng, pop.Transfer*float64(pop.Length))
	// for each transfer, a donor genome is randomly chose, then a fragment is randomly determined.
	for c := 0; c < count; c++ {
		// a different donor genome is randomly selected.
		donor_id := randist.UniformRandomInt(pop.rng, pop.Size)
		for donor_id == idx {
			donor_id = randist.UniformRandomInt(pop.rng, pop.Size)
		}
		donor := pop.Genomes[donor_id]

		// a fragment is randomly determined
		fragL := randist.UniformRandomInt(pop.rng, pop.Length)
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
