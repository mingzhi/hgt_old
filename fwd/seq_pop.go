package fwd

import (
	"github.com/mingzhi/gomath/random"
	"math/rand"
)

type SeqPop struct {
	Size   int // population size
	Length int // genome length

	Mutation float64 // mutation rate

	Transfer float64 // transfer rate
	Fragment int     // transferred fragment length

	Genomes []Sequence // store genomes

	NumOfGeneration int // number of generations

	states string

	src      rand.Source     // random source
	rng      *rand.Rand      // random number generator
	mpoisson *random.Poisson // mutation poisson distribution generator
	tpoisson *random.Poisson // transfer poisson distribution generator
}

func NewSeqPop(size, length int, mutation, transfer float64, fragment int, src rand.Source) *SeqPop {
	// verify
	if size <= 0 || length <= 0 {
		panic("Population size or length is <= 0!")
	}

	// construct population
	pop := SeqPop{
		Size:     size,
		Length:   length,
		Mutation: mutation,
		Transfer: transfer,
		Fragment: fragment,
	}

	// random source
	pop.src = src
	pop.rng = rand.New(src)
	pop.mpoisson = random.NewPoisson(pop.Mutation*float64(pop.Length), src)
	pop.tpoisson = random.NewPoisson(pop.Transfer*float64(pop.Length), src)

	// initialize genomes
	pop.states = "ATGC"
	refseq := make(Sequence, pop.Length)
	for i := 0; i < pop.Length; i++ {
		s := pop.states[pop.rng.Intn(len(pop.states))]
		refseq[i] = s
	}
	pop.Genomes = make([]Sequence, pop.Size)
	for i := 0; i < pop.Size; i++ {
		pop.Genomes[i] = make(Sequence, pop.Length)
		copy(pop.Genomes[i], refseq)
	}

	return &pop
}

func (pop *SeqPop) Evolve(time int) {
	for i := 0; i < time; i++ {
		pop.generate()
		pop.mutate()
		pop.transfer()
		pop.NumOfGeneration++
	}
}

func (pop *SeqPop) Sample(size int) (genomes []Sequence) {
	for i := 0; i < pop.Size; i++ {
		j := pop.rng.Intn(pop.Size)
		pop.Genomes[i], pop.Genomes[j] = pop.Genomes[j], pop.Genomes[i]
	}

	genomes = make([]Sequence, size)
	for i := 0; i < size; i++ {
		genomes[i] = pop.Genomes[i]
	}

	return
}

func (pop *SeqPop) generate() {
	newGenomes := make([]Sequence, pop.Size)
	oldGenomes := pop.Genomes

	used := make([]bool, pop.Size)
	ch := make(chan bool)
	n := 0

	for i := 0; i < pop.Size; i++ {
		o := pop.rng.Intn(pop.Size)
		if used[o] {
			go copyGenome(i, o, newGenomes, oldGenomes, ch)
			n++
		} else {
			newGenomes[i] = oldGenomes[o]
			used[o] = true
		}
	}

	for i := 0; i < n; i++ {
		<-ch
	}

	pop.Genomes = newGenomes
}

func copyGenome(n, o int, newGenomes, oldGenomes []Sequence, ch chan bool) {
	newGenomes[n] = make(Sequence, len(oldGenomes[o]))
	copy(newGenomes[n], oldGenomes[o])
	ch <- true
}

func (pop *SeqPop) mutate() {
	mean := pop.Mutation * float64(pop.Length)
	if pop.mpoisson.Mean != mean {
		pop.mpoisson.Mean = mean
	}

	ch := make(chan bool)
	for i := 0; i < pop.Size; i++ {
		go pop.mutateGenome(i, ch)
	}

	for i := 0; i < pop.Size; i++ {
		<-ch
	}
}

func (pop *SeqPop) mutateGenome(idx int, ch chan bool) {
	genome := pop.Genomes[idx]

	count := pop.mpoisson.Int()

	for i := 0; i < count; i++ {
		pos := pop.rng.Intn(pop.Length)

		s := pop.states[pop.rng.Intn(len(pop.states))]
		for genome[pos] == s {
			s = pop.states[pop.rng.Intn(len(pop.states))]
		}

		genome[pos] = s
	}

	ch <- true
}

func (pop *SeqPop) transfer() {
	mean := pop.Transfer * float64(pop.Length)
	if pop.tpoisson.Mean != mean {
		pop.tpoisson.Mean = mean
	}

	ch := make(chan bool)
	for i := 0; i < pop.Size; i++ {
		go pop.transferGenome(i, ch)
	}

	for i := 0; i < pop.Size; i++ {
		<-ch
	}
}

func (pop *SeqPop) transferGenome(idx int, ch chan bool) {
	receiver := pop.Genomes[idx]

	count := pop.tpoisson.Int()
	for c := 0; c < count; c++ {
		donor_id := pop.rng.Intn(pop.Size)
		for donor_id == idx {
			donor_id = pop.rng.Intn(pop.Size)
		}
		donor := pop.Genomes[donor_id]

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
