// gomakeblastdb: make blast dbs in batch
// argv --
// in: genome folder
// out: blast db folder
// procedure -- 
// 1. looking for .ffn (nucl) or .ffa (prot) in the genome folder.
// 2. blast+ dustmasker (nucl) or segmasker (prot) to get mask information, and store it to .asbn file.
// 3. blast+ makeblastdb to create blast db according the mask data.
package main

import (
	"flag"
	"io/ioutil"
	"log"
	"os/exec"
	"strings"
)

var (
	indir  string // genome folder
	outdir string // blast db folder
)

func init() {
	flag.StringVar(&indir, "in", ".", "genomes folder")
	flag.StringVar(&outdir, "out", ".", "blast db folder")

	flag.Parse()
}

func main() {
	files, err := ioutil.ReadDir(indir)
	if err != nil {
		log.Fatal(err)
	}

	for _, fi := range files {
		terms := strings.Split(fi.Name(), ".")
		if terms[len(terms)-1] == "ffn" {
			log.Printf("Process %s ...\n", fi.Name())
			in := indir + "/" + fi.Name()
			out := outdir + "/" + fi.Name()
			makedb(in, out, "nucl")
			log.Println("Finish!")
		} else if terms[len(terms)-1] == "ffa" {
			log.Printf("Process %s ...\n", fi.Name())
			in := indir + "/" + fi.Name()
			out := outdir + "/" + fi.Name()
			makedb(in, out, "prot")
			log.Println("Finish!")
		}
	}
}

func makedb(in, out, seqtype string) {
	if seqtype == "nucl" {
		dustmasker(in, out+".asbn", "fasta", "maskinfo_asn1_bin")
		makeblastdb(in, out, "fasta", "nucl", out+".asbn")
	} else {
		segmasker(in, out+".asbn", "fasta", "maskinfo_asn1_bin")
		makeblastdb(in, out, "fasta", "prot", out+".asbn")
	}

}

func segmasker(in, out, infmt, outfmt string) {
	cmd := exec.Command("segmasker",
		"-in", in,
		"-out", out,
		"-infmt", infmt,
		"-outfmt", outfmt)
	err := cmd.Run()
	if err != nil {
		log.Fatal(err)
	}
	log.Println("Done segmasker")
}

func dustmasker(in, out, infmt, outfmt string) {
	cmd := exec.Command("dustmasker",
		"-in", in,
		"-out", out,
		"-infmt", infmt,
		"-outfmt", outfmt)
	err := cmd.Run()
	if err != nil {
		log.Fatal(err)
	}
	log.Println("Done dustmasker")
}

func makeblastdb(in, out, input_type, dbtype, mask_data string) {
	cmd := exec.Command("makeblastdb",
		"-in", in,
		"-out", out,
		"-input_type", input_type,
		"-dbtype", dbtype,
		"-mask_data", mask_data)
	err := cmd.Run()
	if err != nil {
		log.Fatal(err)
	}
	log.Println("Done makeblastdb")
}
