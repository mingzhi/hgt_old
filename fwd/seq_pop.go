package fwd

import (
	"bitbucket.org/mingzhi/gsl/randist"
	"encoding/json"
	"log"
)

// A population contains full sequence genomes
type SeqPop struct {
	Size     int     // population size
	Length   int     // genome length
	Mutation float64 // mutation rate
	Transfer float64 // transfer rate
	Fragment int     // transferred fragment length

	Genomes []Sequence // genomes

	NumOfGens int // number of generations

	states int // nucleotides

	// random source
	rng *randist.RNG
}

// full sequence genome
type Sequence []byte

// NewSeqPop returns a population with full sequence genomes
func NewSeqPop(size, length int, mutation, transfer float64, fragment int) *SeqPop {
	// verify population size and genome length
	if size <= 0 || length <= 0 {
		log.Panic("Population size or genome length should be positive!")
	}
	// verify transfer rate and transferred fragment
	if transfer > 0 && fragment <= 0 {
		log.Panic("Transferred fragment should be positive when transfer rate > 0!")
	}

	// construct a population with parameters
	pop := SeqPop{
		Size:     size,
		Length:   length,
		Mutation: mutation,
		Transfer: transfer,
		Fragment: fragment,
	}

	// create random sources
	pop.rng = randist.NewRNG(randist.MT19937)

	// create genomes
	// randomly create a parent sequence
	pop.states = 4 // typical nucleotides, "ATGC"
	refseq := make(Sequence, pop.Length)
	for i := 0; i < len(refseq); i++ {
		refseq[i] = byte(randist.UniformRandomInt(pop.rng, pop.states)) // randomly assign a character
	}
	// copy the parent sequence to every genomes
	pop.Genomes = make([]Sequence, pop.Size)
	for i := 0; i < len(pop.Genomes); i++ {
		pop.Genomes[i] = make(Sequence, pop.Length)
		copy(pop.Genomes[i], refseq)
	}

	return &pop
}

// set seed
func (s *SeqPop) Seed(seed int) {
	s.rng.SetSeed(seed)
}

// Evolve: do one generation of population evolution, which includes:
// 1. Wright-Fisher reproduction
// 2. Mutation
// 3. Horizontal gene transfer
func (pop *SeqPop) Evolve() {
	pop.reproduce()
	pop.manipulate()
	pop.NumOfGens++
}

// GetGenomes: return genome sequences
func (pop *SeqPop) GetGenomes() []Sequence {
	return pop.Genomes
}

// GetLength: return genome length
func (pop *SeqPop) GetLength() int {
	return pop.Length
}

// GetTime: return evolved time
func (pop *SeqPop) GetTime() int {
	return pop.NumOfGens
}

// Json: return the entire population in JSON format
func (pop *SeqPop) Json() []byte {
	b, err := json.Marshal(pop)
	if err != nil {
		log.Panic(err)
	}

	return b
}

// Wright-Fisher generation
func (pop *SeqPop) reproduce() {
	newGenomes := make([]Sequence, pop.Size)
	used := make([]bool, pop.Size) // records used genomes
	for n := 0; n < pop.Size; n++ {
		o := randist.UniformRandomInt(pop.rng, pop.Size) // randomly select a parent
		if used[o] {
			// do hard copy
			newGenomes[n] = make(Sequence, pop.Length)
			copy(newGenomes[n], pop.Genomes[o])
		} else {
			// just pass the pointer of the genome
			newGenomes[n] = pop.Genomes[o]
			used[o] = true
		}
	}

	// update genomes
	pop.Genomes = newGenomes
}

// evolutionary operators: mutation and transfer
func (pop *SeqPop) manipulate() {
	// calculate the number of events, which is poisson distribution
	lambda := float64(pop.Size) * float64(pop.Length) * (pop.Mutation + pop.Transfer)
	k := randist.PoissonRandomInt(pop.rng, lambda)
	// given k, the number of mutations or transfers are following Binormial distribution.
	// therefore, we can do k times of Bernoulli test
	for i := 0; i < k; i++ {
		// determine whether this event is a mutation or transfer
		ratio := pop.Mutation / (pop.Mutation + pop.Transfer) // the ratio of mutation
		r := randist.UniformRandomFloat64(pop.rng)            // randomly produce a probability
		if r <= ratio {                                       // determin whether is a mutation or transfer
			pop.mutate()
		} else {
			pop.transfer()
		}
	}
}

// mutation operator
// parameters: g is genome idx.
func (pop *SeqPop) mutate() {
	g := randist.UniformRandomInt(pop.rng, pop.Size)   // randomly choose a genome
	l := randist.UniformRandomInt(pop.rng, pop.Length) // randomly determine a position
	s := randist.UniformRandomInt(pop.rng, pop.states) // randomly find a idx of the mutated character in pop.states
	for byte(s) == pop.Genomes[g][l] {
		s = randist.UniformRandomInt(pop.rng, pop.states)
	}

	// update the character.
	pop.Genomes[g][l] = byte(s)
}

// transfer operator
// parameter: g is genome idx.
func (pop *SeqPop) transfer() {
	g := randist.UniformRandomInt(pop.rng, pop.Size) // randomly choose a receiver
	d := randist.UniformRandomInt(pop.rng, pop.Size) // randomly choose a donor
	for g == d {
		d = randist.UniformRandomInt(pop.rng, pop.Size)
	}

	l := randist.UniformRandomInt(pop.rng, pop.Length) // randomly choose a left index
	r := l + pop.Fragment                              // fragment right index

	if r < pop.Length {
		copy(pop.Genomes[g][l:r], pop.Genomes[d][l:r])
	} else {
		copy(pop.Genomes[g][l:pop.Length], pop.Genomes[d][l:pop.Length])
		copy(pop.Genomes[g][0:r-pop.Length], pop.Genomes[d][0:r-pop.Length])
	}
}
