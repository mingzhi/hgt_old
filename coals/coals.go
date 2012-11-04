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
	TransferRate   float64 // transfer rate per genome
	TransferLength int     // transfer fragment length
	history        *EvolutionHistory

	seed int64 // random source seed
}

func NewWFPopulation(size, sample, length int, mutation, transfer float64, tract int) *WFPopulation {
	wf := &WFPopulation{
		Size:           size,
		SampleSize:     sample,
		GenomeLength:   length,
		MutationRate:   mutation,
		TransferRate:   transfer,
		TransferLength: tract,
		seed:           1,
	}
	wf.history = NewEvolutionHistory(sample, length)
	return wf
}

func (w *WFPopulation) GetHistory() *EvolutionHistory {
	return w.history
}

func (w *WFPopulation) Seed(seed int64) {
	w.seed = seed
	rand.Seed(seed)
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

	if len(temp) > 0 {
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
	}

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
