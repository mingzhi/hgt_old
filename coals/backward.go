package coals

import (
	"bitbucket.org/mingzhi/gsl/randist"
)

func (w *WFPopulation) Backtrace() {
	for {
		if len(w.history.CurrentPool) <= 1 {
			break
		}
		etype, etime := w.nextEvent()
		// based on the type of event, backward
		if etype == CoalescenceEvent {
			// do coalescence
			ancestor := w.coalescent()
			// create event node
			event := EventNode{Type: etype, Time: etime, Participants: []int{ancestor}}
			w.history.Events = append(w.history.Events, event)
		} else {
			ancestors := w.transfer()
			// create event node
			event := EventNode{Type: etype, Time: etime, Participants: ancestors}
			w.history.Events = append(w.history.Events, event)
		}
	}
}

// determine the next event
func (w *WFPopulation) nextEvent() (eventType int, eventTime float64) {
	k := float64(len(w.history.CurrentPool))
	p := 2.0 * w.TransferRate * float64(w.Size)
	l := (k*p + k*(k-1.0)) / 2.0
	v := randist.ExponentialRandomFloat64(w.rng, 1.0/l)
	eventTime = v
	r := randist.UniformRandomFloat64(w.rng)
	if r < p/(k-1.0+p) {
		eventType = TransferEvent
	} else {
		eventType = CoalescenceEvent
	}
	return
}

// do coalescent
func (w *WFPopulation) coalescent() int {
	// randomly choose two nodes
	a := w.history.CurrentPool[randist.UniformRandomInt(w.rng, len(w.history.CurrentPool))]
	b := w.history.CurrentPool[randist.UniformRandomInt(w.rng, len(w.history.CurrentPool))]
	for a == b {
		b = w.history.CurrentPool[randist.UniformRandomInt(w.rng, len(w.history.CurrentPool))]
	}

	aTreeNode := w.history.Tree[a]
	bTreeNode := w.history.Tree[b]

	// create their ancestor tree node and add it into the tree
	genome := Merge(aTreeNode.Genome, bTreeNode.Genome)
	children := []int{a, b}
	w.history.Tree = append(w.history.Tree, TreeNode{Genome: genome, Children: children})

	// update the current pool
	pool := []int{}
	for _, pos := range w.history.CurrentPool {
		if pos != a && pos != b {
			pool = append(pool, pos)
		}
	}
	ancestorID := len(w.history.Tree) - 1
	pool = append(pool, ancestorID)
	w.history.CurrentPool = pool

	return ancestorID
}

// do transfer
func (w *WFPopulation) transfer() (parents []int) {
	// randomly choose a node
	c := w.history.CurrentPool[randist.UniformRandomInt(w.rng, len(w.history.CurrentPool))]
	// tranferring fragment
	begin := randist.UniformRandomInt(w.rng, w.GenomeLength)
	end := begin + w.TransferLength
	if end < w.GenomeLength {
		amA, amB := Split(w.history.Tree[c].Genome, begin, end)
		if len(amA) != 0 {
			genome := amA
			children := []int{c}
			w.history.Tree = append(w.history.Tree, TreeNode{Genome: genome, Children: children})
			parents = append(parents, len(w.history.Tree)-1)
		}

		if len(amB) != 0 {
			genome := amB
			children := []int{c}
			w.history.Tree = append(w.history.Tree, TreeNode{Genome: genome, Children: children})
			parents = append(parents, len(w.history.Tree)-1)
		}
	} else {
		amB, amA := Split(w.history.Tree[c].Genome, end-w.GenomeLength+1, begin-1)
		if len(amA) != 0 {
			genome := amA
			children := []int{c}
			w.history.Tree = append(w.history.Tree, TreeNode{Genome: genome, Children: children})
			parents = append(parents, len(w.history.Tree)-1)
		}

		if len(amB) != 0 {
			genome := amB
			children := []int{c}
			w.history.Tree = append(w.history.Tree, TreeNode{Genome: genome, Children: children})
			parents = append(parents, len(w.history.Tree)-1)
		}
	}

	// update the current pool
	pool := []int{}
	for _, pos := range w.history.CurrentPool {
		if pos != c {
			pool = append(pool, pos)
		}
	}
	for _, pos := range parents {
		pool = append(pool, pos)
	}
	w.history.CurrentPool = pool

	return
}
