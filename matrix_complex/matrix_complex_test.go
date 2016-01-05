package matrix_complex

import (
	"testing"
	"fmt"
)


func TestNewComplexMatrix(t *testing.T) {
	_, err := NewComplexMatrix(3,3, nil)
	if err!=nil {
		t.Error(err)
	}
	
}

func TestGetDimensions(t *testing.T) {
	numrows, numcolumns := 3,1
	m, _ := NewComplexMatrix(numrows, numcolumns, nil)
	r, c := m.GetDimensions()
	if r!=numrows || c!=numcolumns {
		t.Error("wrong dimensions")
	}
}

func TestGetAt(t *testing.T) {
	identity := NewIdentity(3)
	c := identity.GetAt(2,2)
	if c!=1+0i {
		t.Error("expected complexOne, but not found")
	}
	c = identity.GetAt(2,1)
	if c!=0i {
		t.Error("expected complexZero, but not found")
	}
}

func TestSetAt(t *testing.T) {
	identity := NewIdentity(2)
	identity.SetAt(1,1, 2+0i)
	if identity.GetAt(1,1)!=2+0i {
		t.Error("expeceted 2+0i but found %v", identity.GetAt(1,1))
	}
}

func TestSetColumn(t *testing.T) {
	m, _ := NewComplexMatrix(3,3, nil)
	m2, _ := NewComplexMatrix(3,1, []complex128{1+0i, 0i, 0i})
	err := m.SetColumn(1, m2)
	if err != nil {
		t.Error("%v", err)
	}
}

func TestEquals(t *testing.T) {
	identity := NewIdentity(3)
	identity2 := NewIdentity(3)
	eq := identity.Equals(identity2)
	if !eq {
		t.Error("matrixes should be equal")
	}
	identity3 := NewIdentity(4)
	eq = identity.Equals(identity3)
	if eq {
		t.Error("matrixes shouldn't be equal")
	}
	m, err := NewComplexMatrix(3, 3, []complex128{1+0i, 0i, 0i, 0i, 1+0i, 0i, 0i, 0i, 2+0i})
	if err!=nil {
		t.Error(err)
	}
	eq = identity.Equals(m)
	if eq {
		t.Error("matrixes shouldn't be equal")
	}
}

func TestDeterminant(t *testing.T) {
	d, err := NewIdentity(3).GetDeterminant()
	if err!=nil {
		t.Error(err)
	}
	if d!=1+0i {
		t.Error("determinant wrong")
	}
	e, err := NewComplexMatrix(2,2, []complex128{2+1i, 1+1i, 1+1i, 0+1i})

	// det = 2+i * 0+i - 1+i * 1+i
	// = -1 + 2i   -  2i
	// = -1
	if err !=nil {
		t.Error(err)
	}
	det, err := e.GetDeterminant()
	if err != nil { t.Error(err) }
	if det !=  -1+0i {
		t.Error("determinant wrong: %v", det)
	}
	m, _ := NewComplexMatrix(2,2, []complex128{5+0i, 2+0i, 7+0i, 3+0i})
	det, _ = m.GetDeterminant()
	if det!=1+0i {
		t.Error("determinant wrong: %v", det)
	}
}

func TestGetRow(t *testing.T) {
	r, err := NewIdentity(3).GetRow(4)
	if err==nil {
		t.Error("error should not be nil")
	}
	r, err = NewIdentity(3).GetRow(2)
	if r.GetAt(0,2)!=1+0i {
		t.Error("getrow failed")
	}

}

func TestGetColumn(t *testing.T) {
	r, err := NewIdentity(3).GetColumn(4)
	if err==nil {
		t.Error("error should not be nil")
	}
	r, err = NewIdentity(3).GetColumn(2)
	if r.GetAt(2,0)!=1+0i {
		t.Error("getrow failed")
	}
}

func TestGetSubMatrix(t *testing.T) {
	s, err := NewIdentity(3).GetSubMatrix(1,1)
	if err!=nil {
		t.Error(err)
	}
	if !s.Equals(NewIdentity(2)) {
		t.Error("submatrix should be identity2")
	}
}

func TestTranspose(t *testing.T) {
	m, _ := NewComplexMatrix(2,2, []complex128{5+0i, 2+0i, 7+0i, 3+0i})
	mt, _ := NewComplexMatrix(2,2, []complex128{5+0i, 7+0i, 2+0i, 3+0i})
	if !m.Transpose().Equals(mt) {
		t.Error(fmt.Sprintf("T(%v) should be %v not %v", m, mt, m.Transpose()))
	}
	if !NewIdentity(3).Transpose().Equals(NewIdentity(3)) {
		t.Error("identity should be its own transpose")
	}
}

func TestGetCofactor(t *testing.T) {
	m, _ := NewComplexMatrix(2,2, []complex128{5+0i, 2+0i, 7+0i, 3+0i})
	mc, _ := NewComplexMatrix(2,2, []complex128{3+0i, -7-0i, -2-0i, 5+0i})
	mc2, _ := m.GetCofactor()
	if !mc2.Equals(mc) {
		t.Error(fmt.Sprintf("cofactor of %v should be %v not %v", m, mc, mc2))
	}
}

func TestInverse(t *testing.T) {
	inv, err := NewIdentity(3).Inverse()
	if err!=nil {
		t.Error(err)
	}
	if !inv.Equals(NewIdentity(3)) {
		t.Error("identity should be its own inverse")
	}
	m, _ := NewComplexMatrix(2,2, []complex128{5+0i, 2+0i, 7+0i, 3+0i})
	minv, _ := m.Inverse()
	minv2, _ := NewComplexMatrix(2,2, []complex128{3+0i, -2-0i, -7-0i, 5+0i})
	if !minv2.Equals(minv) {
		t.Error(fmt.Sprintf("%v should be %v", minv, minv2))
	}
}

func TestTimes(t *testing.T) {
	res, err := NewIdentity(3).Times(NewIdentity(3))
	if err!= nil {
		t.Error(err)
	}
	if !res.Equals(NewIdentity(3)) {
		t.Error("I * I should equal I")
	}
	cm1, _ := NewComplexMatrix(1,3, []complex128{1+0i, 2+0i, 3+0i})
	cm2, _ := NewComplexMatrix(3,1, []complex128{4+0i, 5+0i, 6+0i})
	p, err := cm1.Times(cm2)
	if err != nil { t.Error(err) }
	if nr, nc := p.GetDimensions(); nr!=1 || nc!=1 { t.Error("wrong dimensions") }
	if p.GetAt(0,0) != 32+0i {
		t.Error("wrong value: %v", p.GetAt(0,0))
	}
}
