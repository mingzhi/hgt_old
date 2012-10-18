package fwd

import (
	"encoding/json"
	"github.com/mingzhi/gomath/random"
	"log"
	"math/rand"
	"time"
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

	states string // nucleotides

	// random variables
	src  rand.Source     // random source
	rng  *rand.Rand      // random number generator
	pois *random.Poisson // poisson sampling
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
	// use math/rand source
	seed := time.Now().UnixNano() // use time now as a seed
	pop.src = rand.NewSource(seed)
	pop.rng = rand.New(pop.src)
	// determin poisson lambda (expected mean)
	// mutation events and transfer events are independent poisson distribution.
	// therefore, the total number of events are also following poission distrubtion
	// with mean = (u + r) * L * N
	mean := float64(pop.Size*pop.Length) * (pop.Mutation + pop.Transfer)
	pop.pois = random.NewPoisson(mean, pop.src)

	// create genomes
	// randomly create a parent sequence
	pop.states = "atgc" // typical nucleotides, "ATGC"
	refseq := make(Sequence, pop.Length)
	for i := 0; i < len(refseq); i++ {
		refseq[i] = pop.states[pop.rng.Intn(len(pop.states))] // randomly assign a character
	}
	// copy the parent sequence to every genomes
	pop.Genomes = make([]Sequence, pop.Size)
	for i := 0; i < len(pop.Genomes); i++ {
		pop.Genomes[i] = make(Sequence, pop.Length)
		copy(pop.Genomes[i], refseq)
	}

	return &pop
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
		o := pop.rng.Intn(pop.Size) // randomly select a parent
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
	k := pop.pois.Int() // total number of events
	// given k, the number of mutations or transfers are following Binormial distribution.
	// therefore, we can do k times of Bernoulli test
	for i := 0; i < k; i++ {
		// determine a genome
		g := pop.rng.Intn(pop.Size)
		// determine whether this event is a mutation or transfer by flipping a coin
		ratio := pop.Mutation / (pop.Mutation + pop.Transfer) // the ratio of mutation
		r := pop.rng.Float64()                                // randomly produce a probability
		if r <= ratio {                                       // determin whether is a mutation or transfer
			pop.mutate(g)
		} else {
			pop.transfer(g)
		}
	}
}

// mutation operator
// parameters: g is genome idx.
func (pop *SeqPop) mutate(g int) {
	l := pop.rng.Intn(pop.Length)            // randomly determine a position
	ri := pop.rng.Intn(len(pop.states))      // randomly find a idx of the mutated character in pop.states
	if pop.states[ri] == pop.Genomes[g][l] { // if the mutated character is same as the old one
		ri = ri + pop.rng.Intn(len(pop.states)-1) + 1 // then shuffle 3 positions, so that we don't need for loop.
		if ri >= len(pop.states) {
			ri = ri - len(pop.states)
		}
	}

	// update the character.
	pop.Genomes[g][l] = pop.states[ri]
}

// transfer operator
// parameter: g is genome idx.
func (pop *SeqPop) transfer(g int) {
	// randomly choose a donor
	d := pop.rng.Intn(pop.Size)
	for d == g {
		d = pop.rng.Intn(pop.Size)
	}

	// randomly decide the boundaris of the transferred fragment.
	left := pop.rng.Intn(pop.Length)
	right := left + pop.Fragment
	if right < pop.Length {
		for i := left; i < right; i++ {
			pop.Genomes[g][i] = pop.Genomes[d][i]
		}
	} else { // microorganizm genomes are circle.
		for i := left; i < pop.Length; i++ {
			pop.Genomes[g][i] = pop.Genomes[d][i]
		}
		for i := 0; i < right-pop.Length; i++ {
			pop.Genomes[g][i] = pop.Genomes[d][i]
		}
	}
}
