package matrix_complex

import (
	"testing"
	"log"
)

func TestGetAt(t *testing.T) {
	identity := NewIdentity(3)
	c, err := identity.GetAt(3,2)
	if err==nil {
		t.Error("expected error but none was thrown")
	}
	c, err = identity.GetAt(2,2)
	if c!=1+0i {
		t.Error("expected complexOne, but not found")
	}
	c = identity.getAt(2,1)
	if c!=0i {
		t.Error("expected complexZero, but not found")
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
}

func TestGetRow(t *testing.T) {
	r, err := NewIdentity(3).GetRow(4)
	if err==nil {
		t.Error("error should not be nil")
	}
	r, err = NewIdentity(3).GetRow(2)
	if r.getAt(0,2)!=1+0i {
		t.Error("getrow failed")
	}

}

func TestGetColumn(t *testing.T) {
	r, err := NewIdentity(3).GetColumn(4)
	if err==nil {
		t.Error("error should not be nil")
	}
	r, err = NewIdentity(3).GetColumn(2)
	if r.getAt(2,0)!=1+0i {
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
	log.Print(NewIdentity(3).Transpose())
	if !NewIdentity(3).Transpose().Equals(NewIdentity(3)) {
		t.Error("identity should be its own transpose")
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
}

func TestTimes(t *testing.T) {
	res, err := NewIdentity(3).Times(NewIdentity(3))
	if err!= nil {
		t.Error(err)
	}
	if !res.Equals(NewIdentity(3)) {
		t.Error("I * I should equal I")
	}
}
