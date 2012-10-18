package fwd

type Population interface {
	Evolve()
	GetGenomes() []Sequence
	GetLength() int
	GetTime() int
	Json() []byte
}
