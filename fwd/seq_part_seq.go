package fwd

import (
	"encoding/json"
	"github.com/mingzhi/gomath/random"
	"log"
	"math/rand"
	"time"
)

// a SeqPartPop: a population with partial sequence genomes
type SeqPartPop struct {
	Size      int        // population size
	Length    int        // partial genome length
	Fragment  int        // length of transferred fragment
	NumOfGens int        // number of generations
	Mutation  float64    // mutation rate
	Transfer  float64    // transfer rate
	Genomes   []Sequence // genomes

	states string // nucleotide characters

	// random variables
	src  rand.Source     // random source
	rng  *rand.Rand      // random number generator
	pois *random.Poisson // poisson sampling
}

// NewSeqPartPop: returns a new created partial genome population.
func NewSeqPartPop(size, length int, mutation, transfer float64, fragment int) *SeqPartPop {
	// verify population size and genome length.
	if size <= 0 || length <= 0 {
		log.Panic("Population size or genome length or partial length should be positive!")
	}

	// construct a population with parameters
	pop := SeqPartPop{
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
	// with mean = (u + r) * (P+F) * N, where P is the partial length of a genome.
	// mutation events: u * l * N, transfer events: r * (l + f) * N
	mean := float64(pop.Size) * (pop.Mutation*float64(pop.Length) + pop.Transfer*float64(pop.Length+pop.Fragment))
	pop.pois = random.NewPoisson(mean, pop.src)

	// create genomes (partial length)
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
func (pop *SeqPartPop) Evolve() {
	pop.reproduce()
	pop.manipulate()
	pop.NumOfGens++
}

// Json: return the entire population in JSON format
func (pop *SeqPartPop) Json() []byte {
	b, err := json.Marshal(pop)
	if err != nil {
		log.Panic(err)
	}

	return b
}

// Wright-Fisher generation
func (pop *SeqPartPop) reproduce() {
	newGenomes := make([]Sequence, pop.Size)
	used := make([]bool, pop.Size)  // records used genomes
	for n := 0; n < pop.Size; n++ { // n is new genome index
		o := pop.rng.Intn(pop.Size) // randomly generate a old genome index
		if used[o] {                // this old genome has been reassigned
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

// mutation and transfer
// apply evolutionary operators on the population
func (pop *SeqPartPop) manipulate() {
	k := pop.pois.Int() // total number of events
	// given k, the number of mutations or transfers are following Binormial distribution.
	// therefore, we can do k times of Bernoulli test
	for i := 0; i < k; i++ {
		// determine whether this event is a mutation or transfer by flipping a coin
		// lambda = (mu * length) / ((mu * length) * (tr * (length + fragment)))
		pmu := pop.Mutation * float64(pop.Length)              // proportion of mutation events
		ptr := pop.Transfer * float64(pop.Length+pop.Fragment) // proportion of transfer events
		lambda := pmu / (pmu + ptr)                            // ratio of mutation and transfer events
		r := pop.rng.Float64()
		if r < lambda {
			pop.mutate()
		} else {
			pop.transfer()
		}
	}
}

// mutation operator.
func (pop *SeqPartPop) mutate() {
	g := pop.rng.Intn(pop.Size)               // choose a genome
	l := pop.rng.Intn(pop.Length)             // choose a position
	ri := pop.rng.Intn(len(pop.states))       // randomly find a idx of the mutated character in pop.states
	for pop.states[ri] == pop.Genomes[g][l] { // make sure to change to a different state
		ri = pop.rng.Intn(len(pop.states))
	}

	// update the character.
	pop.Genomes[g][l] = pop.states[ri]
}

// transfer operator
func (pop *SeqPartPop) transfer() {
	g := pop.rng.Intn(pop.Size) // choose a receiver
	d := pop.rng.Intn(pop.Size) // choose a donor
	for d == g {
		d = pop.rng.Intn(pop.Size)
	}

	// choose a start position
	l := pop.rng.Intn(pop.Length+pop.Fragment) - pop.Fragment

	// determin the boundaries
	// begin of the transferred fragment
	begin := 0
	if l < 0 {
		begin = 0
	} else {
		begin = l
	}

	// end of the transferred fragment
	end := pop.Length
	if l+pop.Fragment < pop.Length {
		end = l + pop.Fragment
	} else {
		end = pop.Length
	}

	// do the transfer
	for i := begin; i < end; i++ {
		pop.Genomes[g][i] = pop.Genomes[d][i]
	}
}
