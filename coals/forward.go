package coals

import (
	"bitbucket.org/mingzhi/gsl/randist"
)

const (
	NucleicAcids = "ATGC"
)

func (w *WFPopulation) Fortrace() (seqMap map[int][]byte) {
	seqMap = make(map[int][]byte)
	for ei := len(w.history.Events) - 1; ei >= 0; ei-- {
		event := w.history.Events[ei]
		if event.Type == CoalescenceEvent {
			a := event.Participants[0] // just one participant, the ancestor
			seq, yes := seqMap[a]
			if !yes {
				seq = randomGenerateSequence(w.GenomeLength, w.rng)
			}
			b := w.history.Tree[a].Children[0]
			c := w.history.Tree[a].Children[1]
			seqB := make([]byte, w.GenomeLength)
			seqC := make([]byte, w.GenomeLength)
			copy(seqB, seq)
			copy(seqC, seq)
			delete(seqMap, a)
			seqMap[b] = seqB
			seqMap[c] = seqC
		} else {
			seqC := make([]byte, w.GenomeLength)
			for _, a := range event.Participants {
				seqA := seqMap[a]
				for _, frag := range w.history.Tree[a].Genome {
					for i := frag.Begin; i <= frag.End; i++ {
						seqC[i] = seqA[i]
					}
				}
				delete(seqMap, a)
			}
			c := w.history.Tree[event.Participants[0]].Children[0]
			seqMap[c] = seqC
		}
		lambda := event.Time * w.MutationRate * float64(w.Size) * float64(w.GenomeLength)
		mutateAll(seqMap, lambda, w.GenomeLength, w.rng)
	}

	return
}

func mutateAll(seqMap map[int][]byte, lambda float64, l int, rng *randist.RNG) {
	for _, seq := range seqMap {
		count := randist.PoissonRandomInt(rng, lambda)
		for i := 0; i < count; i++ {
			idx := randist.UniformRandomInt(rng, l)
			a := NucleicAcids[randist.UniformRandomInt(rng, len(NucleicAcids))]
			for a == seq[idx] {
				a = NucleicAcids[randist.UniformRandomInt(rng, len(NucleicAcids))]
			}
			seq[idx] = a
		}
	}
}

func randomGenerateSequence(length int, rng *randist.RNG) (seq []byte) {
	seq = make([]byte, length)
	for i, _ := range seq {
		seq[i] = NucleicAcids[randist.UniformRandomInt(rng, len(NucleicAcids))]
	}
	return
}
