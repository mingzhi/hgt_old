package coals

import (
	"github.com/mingzhi/gomath/random"
	"math/rand"
)

const (
	NucleicAcids = "ATGC"
)

func (w *WFPopulation) Fortrace() (seqMap map[int][]byte) {
	seqMap = make(map[int][]byte)
	for _, event := range w.history.Events {
		if event.Type == CoalescenceEvent {
			a := event.Participants[0] // just one participant, the ancestor
			seq, yes := seqMap[a]
			if !yes {
				seq = randomGenerateSequence(w.GenomeLength)
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
		mutateAll(seqMap, event.Time*float64(w.Size), w.MutationRate, w.GenomeLength, w.seed)
	}

	return
}

func mutateAll(seqMap map[int][]byte, t, u float64, l int, seed int64) {
	lambda := t * u * float64(l)
	poisson := random.NewPoisson(lambda, rand.NewSource(seed))
	for _, seq := range seqMap {
		count := poisson.Int()
		for i := 0; i < count; i++ {
			idx := rand.Intn(l)
			a := NucleicAcids[rand.Intn(len(NucleicAcids))]
			for a == seq[idx] {
				a = NucleicAcids[rand.Intn(len(NucleicAcids))]
			}
			seq[idx] = a
		}
	}
}

func randomGenerateSequence(length int) (seq []byte) {
	seq = make([]byte, length)
	for i, _ := range seq {
		seq[i] = NucleicAcids[rand.Intn(len(NucleicAcids))]
	}
	return
}
