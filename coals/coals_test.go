package coals

import "testing"

func TestMerge(t *testing.T) {
	frag1 := Fragment{Begin: 0, End: 100}
	frag2 := Fragment{Begin: 44, End: 90}
	a1 := Assembly{frag1}
	a2 := Assembly{frag2}
	a3 := Merge(a1, a2)
	if len(a3) != 1 {
		t.Errorf("Expected to get 1 fragment, but got %d", len(a3))
	}
	if a3[0].Begin != 0 {
		t.Errorf("Begin should be 0, but got %d", a3[0].Begin)
	}
	if a3[0].End != 100 {
		t.Errorf("End should be 100, but got %d", a3[0].End)
	}

	frag3 := Fragment{Begin: 120, End: 130}
	frag4 := Fragment{Begin: 90, End: 110}
	a4 := Assembly{frag3, frag4}
	a5 := Merge(a3, a4)
	if len(a5) != 2 {
		t.Errorf("Expected to get 2 fragments, but got %d", len(a5))
	}
	if a5[0].Begin != 0 {
		t.Errorf("Begin should be 0, but got %d", a5[0].Begin)
	}
	if a5[0].End != 110 {
		t.Errorf("End should be 110, but got %d", a5[0].End)
	}
	if a5[1].Begin != 120 {
		t.Errorf("Begin should be 120, but got %d", a5[1].Begin)
	}
	if a5[1].End != 130 {
		t.Errorf("End should be 130, but got %d", a5[1].End)
	}
}

func TestSplit(t *testing.T) {
	frag1 := Fragment{Begin: 0, End: 30}
	frag2 := Fragment{Begin: 44, End: 90}
	frag3 := Fragment{Begin: 120, End: 150}
	a := Assembly{frag1, frag2, frag3}
	b, c := Split(a, 50, 130)
	if len(b) != 3 {
		t.Errorf("len(b) = 3, but got %d", len(b))
	}
	if c[0].Begin != 50 {
		t.Errorf("Expected to get 50, but got %d", c[0].Begin)
	}
	if c[0].End != 90 {
		t.Errorf("Expected to get 90, but got %d", c[0].End)
	}
	if c[1].Begin != 120 {
		t.Errorf("Expected to get 120, but got %d", c[1].Begin)
	}
	if c[1].End != 130 {
		t.Errorf("Expected to get 130, but got %d", c[1].End)
	}
	if b[1].Begin != 44 {
		t.Errorf("Expected to get 44, but got %d", b[1].Begin)
	}
	if b[1].End != 49 {
		t.Errorf("Expected to get 49, but got %d", b[1].End)
	}
}
