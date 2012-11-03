package coals

import (
	"math/rand"
	"sort"
)

const (
	CoalescenceEvent = 0
	TransferEvent    = 1
)

type WFPopulation struct {
	Size           int     // population size
	SampleSize     int     // sample size
	GenomeLength   int     // genome length
	MutationRate   float64 // mutation rate
	TransferRate   float64 // transfer rate
	TransferLength int     // transfer fragment length
	history        *EvolutionHistory
}

func NewWFPopulation(size, sample, length int, mutation, transfer float64, tract int) *WFPopulation {
	wf := &WFPopulation{
		Size:           size,
		SampleSize:     sample,
		GenomeLength:   length,
		MutationRate:   mutation,
		TransferRate:   transfer,
		TransferLength: tract,
	}
	wf.history = NewEvolutionHistory(size, length)
	return wf
}

func (w *WFPopulation) Backtrace() {
	w.backward()
}

func (w *WFPopulation) backward() {
	// determine event type and time
	etype, etime := w.nextEvent()
	// based on the type of event, backward
	if etype == CoalescenceEvent {
		// do coalescence
		ancestor := w.coalescent()
		// iterate until only one node left in the pool
		if len(w.history.CurrentPool) > 1 {
			w.backward()
		}
		// create event node
		event := EventNode{Type: etype, Time: etime, Participants: []int{ancestor}}
		w.history.Events = append(w.history.Events, event)
	} else {
		ancestors := w.transfer()
		// iterate until only one node left in the pool
		if len(w.history.CurrentPool) > 1 {
			w.backward()
		}
		// create event node
		event := EventNode{Type: etype, Time: etime, Participants: ancestors}
		w.history.Events = append(w.history.Events, event)
	}
}

// determine the next event
func (w *WFPopulation) nextEvent() (eventType int, eventTime float64) {
	k := float64(len(w.history.CurrentPool))
	p := w.TransferRate * float64(w.GenomeLength*w.Size)
	l := (k*p + k*(k-1.0)) / 2.0
	v := rand.ExpFloat64() * l
	eventTime = v * float64(w.Size)
	r := rand.Float64()
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
	a := w.history.CurrentPool[rand.Intn(len(w.history.CurrentPool))]
	b := w.history.CurrentPool[rand.Intn(len(w.history.CurrentPool))]
	for a == b {
		b = w.history.CurrentPool[rand.Intn(len(w.history.CurrentPool))]
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
	c := w.history.CurrentPool[rand.Intn(len(w.history.CurrentPool))]
	// tranferring fragment
	begin := rand.Intn(w.GenomeLength)
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

type Fragment struct {
	Begin    int    // begin position of this fragment
	End      int    // end position of this fragment
	Sequence []byte // fragment sequence
}

type Assembly []Fragment

// interface functions for sort package
func (a Assembly) Len() int           { return len(a) }
func (a Assembly) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a Assembly) Less(i, j int) bool { return a[i].Begin < a[j].Begin }

func Merge(a, b Assembly) (c Assembly) {
	temp := Assembly{}
	for _, frag := range a {
		temp = append(temp, frag)
	}
	for _, frag := range b {
		temp = append(temp, frag)
	}
	sort.Sort(temp)
	begin := temp[0].Begin
	end := temp[0].End
	for _, frag := range temp {
		if frag.Begin <= end {
			if frag.End > end {
				end = frag.End
			}
		} else {
			c = append(c, Fragment{Begin: begin, End: end})
			begin = frag.Begin
			end = frag.End
		}
	}
	c = append(c, Fragment{Begin: begin, End: end})
	return
}

func Split(a Assembly, begin, end int) (b, c Assembly) {
	sort.Sort(a)
	// find the first fragment, which begin index is in
	idxL := sort.Search(len(a), func(i int) bool { return a[i].End >= begin })
	idxR := sort.Search(len(a), func(i int) bool { return a[i].Begin > end })
	for i, _ := range a {
		if i >= idxL && i < idxR {
			if begin <= a[i].Begin {
				begin = a[i].Begin
			} else {
				b = append(b, Fragment{Begin: a[i].Begin, End: begin - 1})
			}
			if end >= a[i].End {
				c = append(c, Fragment{Begin: begin, End: a[i].End})
			} else {
				c = append(c, Fragment{Begin: begin, End: end})
				b = append(b, Fragment{Begin: end + 1, End: a[i].End})
			}
		} else {
			b = append(b, a[i])
		}
	}

	return
}

type EvolutionHistory struct {
	CurrentPool []int
	Tree        []TreeNode
	Events      []EventNode
}

func NewEvolutionHistory(size, length int) *EvolutionHistory {
	history := EvolutionHistory{}
	// initilize the current pool and tree leave
	for i := 0; i < size; i++ {
		history.CurrentPool = append(history.CurrentPool, i)
		history.Tree = append(history.Tree, TreeNode{Genome: Assembly{Fragment{Begin: 0, End: length - 1}}})
	}

	return &history
}

type TreeNode struct {
	Genome   Assembly
	Children []int
}

type EventNode struct {
	Time         float64
	Type         int
	Participants []int
}
