// This scripts finds all orthologs using Reciptcal Best Hits method.
// It requires blast+ installed and added to the PATH.
// findortholog -evalue=1e-5 -out=ortholog genome_paths.csv

package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"runtime"
	"strings"
)

type species struct {
	name    string   // species name
	genomes []string // genome locations
}

var evalue string
var out string

func init() {
	flag.StringVar(&evalue, "evalue", "1e-5", "blast evalue cutoff")
	flag.StringVar(&out, out, ".", "output dir")

	flag.Parse()
}

func main() {
	// read species map
	speciesmap := readSpeciesMap(flag.Arg(0))
	specieslist := []species{}
	for name, genomes := range speciesmap {
		s := species{name, genomes}
		specieslist = append(specieslist, s)
	}

	// parallel find orthologs
	ncpu := runtime.NumCPU()
	ch := make(chan bool, ncpu)
	for i := 0; i < ncpu; i++ {
		b := i * len(specieslist) / ncpu
		e := (i + 1) * len(specieslist) / ncpu
		go func(begin, end int, list []species, queen chan bool) {
			for j := begin; j < end; j++ {
				s := list[j]
				findOrthologs(s)
			}
			queen <- true
		}(b, e, specieslist, ch)
	}
	for i := 0; i < ncpu; i++ {
		<-ch
	}
}

func findOrthologs(s species) {
	outfile, err := os.Create(out + "/" + s.name + "_orthologs.txt")
	if err != nil {
		log.Fatal(err)
	}
	defer outfile.Close()

	for i := 0; i < len(s.genomes); i++ {
		for j := i + 1; j < len(s.genomes); j++ {
			refgenome := s.genomes[i]
			taggenome := s.genomes[j]
			orthologs := rbhOrthologs(refgenome, taggenome)
			for _, pair := range orthologs {
				outfile.WriteString(pair + "\n")
			}
		}
	}

	log.Printf("Done finding orthologs: %s\n", s.name)
}

func rbhOrthologs(genomeA, genomeB string) []string {
	// recipctal blast
	pair1 := blastp(genomeA, genomeB)
	log.Printf("%s vs %s: %d\n", genomeA, genomeB, len(pair1))
	pair2 := blastp(genomeB, genomeA)
	log.Printf("%s vs %s: %d\n", genomeB, genomeA, len(pair2))

	pairs := []string{}
	for gene1, gene2 := range pair1 {
		if pair2[gene2] == gene1 {
			pairs = append(pairs, fmt.Sprintf("%s %s", gene1, gene2))
		}
	}

	return pairs
}

// blastp to find the best hit
func blastp(db, query string) map[string]string {
	cmd := exec.Command("blastp", "-db", db, "-query", query, "-outfmt", "6 qseqid sseqid evalue", "-evalue", evalue)
	stdout, err := cmd.StdoutPipe()
	if err != nil {
		log.Fatal(err)
	}
	defer stdout.Close()

	err = cmd.Start()
	if err != nil {
		log.Fatal(err)
	}

	pair := make(map[string]string)
	reader := bufio.NewReader(stdout)
	line, err := reader.ReadString('\n')
	for err == nil {
		terms := strings.Split(line[:len(line)-1], "\t")
		_, found := pair[terms[0]]
		if !found {
			pair[terms[0]] = terms[1]
		}
		line, err = reader.ReadString('\n')
	}

	if err != io.EOF {
		log.Fatal(err)
	}

	return pair

}

func readSpeciesMap(filename string) map[string][]string {
	file, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	reader := bufio.NewReader(file)

	speciesmap := make(map[string][]string)
	line, err := reader.ReadString('\n')
	for err == nil {
		parts := strings.Split(line[:len(line)-1], ",")
		name := parts[0]
		genome := parts[1]
		genomes, found := speciesmap[name]
		if !found {
			genomes = []string{}
		}
		genomes = append(genomes, genome)
		speciesmap[name] = genomes

		line, err = reader.ReadString('\n')
	}

	if err != io.EOF {
		log.Fatal(err)
	}

	return speciesmap
}
