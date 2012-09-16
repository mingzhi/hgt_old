package covs

import (
	"fmt"
	"github.com/mingzhi/gomath/stat/desc"
	"math"
	"math/rand"
	"runtime"
	"testing"
)

func TestCov(t *testing.T) {
	a := 0.9
	b := 0.3
	length := 1000
	size := 10000
	c := simuMarkovMatrix(a, b, size, length)
	maxL := 100
	scovs, rcovs, _, _, _ := c.Cov(maxL)
	sterr := 0.0
	for i := 0; i < maxL; i++ {
		sc := scovs[i]
		f := (1.0 - a) * (1.0 - b) / ((2.0 - a - b) * (2.0 - a - b)) * math.Pow((a+b-1.0), float64(i))
		sterr += (sc - f) * (sc - f)
	}
	if sterr > 1e-6 {
		t.Errorf("The standard error is too big: %g\n", sterr)
	}

	rmean := desc.NewMean()
	for i := 0; i < maxL; i++ {
		rmean.Increment(rcovs[i])
	}
	prec := 1e-6
	if math.Abs(rmean.GetResult()) > prec {
		t.Errorf("Expect average of rcovs as 0.0, but got %g", rmean.GetResult())
	}
}

func TestCalculateKS(t *testing.T) {
	fmt.Println("--Population parameters:")
	var (
		n         int     // population size
		u         float64 // mutation rate
		r         float64 // transfer rate
		a         int     // number of states
		l         int     // fragment length
		ks        float64
		expected  float64
		tolerance float64
	)
	n = 1000
	u = 1e-4
	r = 1e-4
	a = 4
	l = 100
	tolerance = 1e-5
	fmt.Printf("---Population size: %d\n", n)
	fmt.Printf("---Mutation rate: %f\n", u)
	fmt.Printf("---Transfer rate: %f\n", r)
	fmt.Printf("---Character number: %d\n", a)
	fmt.Printf("---Transferred fragment: %d\n", l)
	ks = CalculateKS(n, u, r, l, a)
	expected = 0.15667
	if math.Abs(ks-expected) > tolerance {
		t.Errorf("ks = %g, but expected %g", ks, expected)
	}

	// no transfer
	r = 0.0
	ks = CalculateKS(n, u, r, l, a)
	expected = 0.157911
	if math.Abs(ks-expected) > tolerance {
		t.Errorf("ks = %g, but expected %g", ks, expected)
	}
}

func BenchmarkCov(bench *testing.B) {
	bench.StopTimer()
	runtime.GOMAXPROCS(runtime.NumCPU())
	a := 0.9
	b := 0.3
	length := 1000
	size := 10000
	c := simuMarkovMatrix(a, b, size, length)
	maxL := 100
	bench.StartTimer()
	for i := 0; i < bench.N; i++ {
		c.Cov(maxL)
	}
}

func simuMarkovMatrix(a, b float64, size, length int) (c *CMatrix) {
	distmatrix := make([][]int, size)
	for i := 0; i < size; i++ {
		distmatrix[i] = simuMarkovProcess(a, b, length)
	}

	c = NewCMatrix(size, length, distmatrix)
	return
}

func simuMarkovProcess(a, b float64, length int) (v []int) {
	prev := 0
	v = []int{}
	for i := 0; i < length; i++ {
		if prev == 0 {
			if rand.Float64() < a {
				prev = 0
			} else {
				v = append(v, i)
				prev = 1
			}
		} else {
			if rand.Float64() < b {
				v = append(v, i)
				prev = 1
			} else {
				prev = 0
			}
		}
	}
	return
}
