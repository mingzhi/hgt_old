// This script first finds species with more than 4 strains.
// Then, it makes blast dbs of the chose genomes.
// It requires blast+ has been installed and added to the PATH.
// The input is a director of all genomes.

package main

import (
	"flag"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"runtime"
	"strings"
)

var cutoff int

func init() {
	flag.IntVar(&cutoff, "cutoff", 4, "number of strains of a chosen species")

	flag.Parse()
}

func main() {
	genomedir := flag.Arg(0)

	sfiles, err := ioutil.ReadDir(genomedir)
	if err != nil {
		log.Fatal(err)
	}

	genomemap := make(map[string][]string)
	for _, file := range sfiles {
		if file.IsDir() {
			parts := strings.Split(file.Name(), "_")
			if len(parts) > 2 {
				name := strings.Join(parts[:2], "_") // species name
				// find the chromosome
				// assumed that only one chromosome
				// others are plasims with smaller sizes
				cfiles, err := ioutil.ReadDir(genomedir + "/" + file.Name())
				if err != nil {
					log.Fatal(err)
				}
				var maxsize int64
				var maxname string
				for _, cfile := range cfiles {
					if cfile.Size() > maxsize {
						maxsize = cfile.Size()
						maxname = cfile.Name()
					}
				}
				chrpath := genomedir + "/" + file.Name() + "/" + maxname

				// add the chrpath into genomemap
				genomes, found := genomemap[name]
				if !found {
					genomes = []string{}
				}
				genomes = append(genomes, chrpath)
				genomemap[name] = genomes
			}
		}
	}

	// filter, only obtain species with more than $cutoff$ strains.
	// output to "genome_paths.csv"
	outfile, err := os.Create("genome_paths.csv")
	if err != nil {
		log.Fatal(err)
	}
	defer outfile.Close()

	totalgenomes := []string{}
	num := 0
	for name, genomes := range genomemap {
		if len(genomes) >= cutoff {
			for _, genome := range genomes {
				outfile.WriteString(name + "," + genome + "\n")
				totalgenomes = append(totalgenomes, genome)
			}
			num++
		}
	}
	log.Printf("Total %d species were chose, writed into genome_paths.csv\n", num)

	// parallel make blast db
	ncpu := runtime.NumCPU()
	ch := make(chan bool, ncpu)
	for i := 0; i < ncpu; i++ {
		b := i * len(totalgenomes) / ncpu
		e := (i + 1) * len(totalgenomes) / ncpu
		go func(begin, end int, paths []string, queen chan bool) {
			for j := begin; j < end; j++ {
				genomepath := paths[j]
				makeblastdb(genomepath)
			}
			queen <- true
		}(b, e, totalgenomes, ch)
	}

	for i := 0; i < ncpu; i++ {
		<-ch
	}
}

func makeblastdb(genome string) {
	// get mask information
	asn := genome + ".asbn"
	maskcmd := exec.Command("segmasker", "-in", genome, "-out", asn, "-outfmt", "maskinfo_asn1_bin", "-parse_seqids")
	maskerr := maskcmd.Run()
	if maskerr != nil {
		log.Fatal(maskerr)
	}

	// make blast db
	out := genome
	mkdbcmd := exec.Command("makeblastdb", "-in", genome, "-out", out, "-dbtype", "prot", "-mask_data", asn, "-parse_seqids")
	mkdberr := mkdbcmd.Run()
	if mkdberr != nil {
		log.Fatal(mkdberr)
	}

	log.Printf("Finished make blast db: %s\n", genome)
}
