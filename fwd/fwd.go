package fwd

type Sequence []byte

func GenerateDistanceMatrix(genomes []Sequence) (dmatrix [][]int) {
	ch := make(chan []int)
	n := 0
	for i := 0; i < len(genomes); i++ {
		a := genomes[i]
		for j := i + 1; j < len(genomes); j++ {
			b := genomes[j]
			go calculateDistance(a, b, ch)
			n++
		}
	}

	dmatrix = make([][]int, n)
	for i := 0; i < n; i++ {
		dmatrix[i] = <-ch
	}

	return
}

func calculateDistance(a, b Sequence, ch chan []int) {
	ds := []int{}
	for i := 0; i < len(a); i++ {
		if a[i] != b[i] {
			ds = append(ds, i)
		}
	}

	ch <- ds
}
