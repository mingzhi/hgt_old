package fwdsimu

import (
	"math/rand"
)

// Sample: returns a sample of genomes.
func Sample(genomes []Sequence, size int) []Sequence {
	// shuffle the idx
	idxes := rand.Perm(len(genomes))

	// get the sample according the shuffled idexes.
	sample := make([]Sequence, size)
	for i := 0; i < size; i++ {
		sample[i] = genomes[idxes[i]]
	}

	return sample
}

// GenerateDistanceMatrix: generate distance matrix given a set of genomes.
func GenerateDistanceMatrix(genomes []Sequence) [][]int {
	dmatirx := [][]int{}
	for i := 0; i < len(genomes); i++ {
		for j := i + 1; j < len(genomes); j++ {
			a := genomes[i] // sequence a
			b := genomes[j] // sequence b
			ds := []int{}   // array for mismatch positions
			for k := 0; k < len(a); k++ {
				if a[k] != b[k] {
					ds = append(ds, k)
				}
			}
			dmatirx = append(dmatirx, ds)
		}
	}

	return dmatirx
}
