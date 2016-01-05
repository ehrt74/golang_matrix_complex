package matrix_complex

import (
	"errors"
	"fmt"
	"strings"
)

var WrongValueListLengthError error = errors.New("wrong value list length")
var IndexOutOfBoundsError error = errors.New("index out of bounds")
var NotSquareError error = errors.New("matrix not square")
var DeterminantZeroError error = errors.New("determinant zero")
var IncompatibleMatrixSizesError error = errors.New("incompatible matrix sizes")

type ComplexMatrix struct {
	numrows, numcolumns int
	vals []complex128
}

func (cm *ComplexMatrix)GetDimensions() (int, int) {
	return cm.numrows, cm.numcolumns
}

func (cm *ComplexMatrix)String() string {
	var strs []string
	for _, val := range cm.vals {
		strs=append(strs, fmt.Sprintf("%.3f", val))
	}
	return strings.Join(strs, ",")
}

// zero based. i is row, j is column
func (cm *ComplexMatrix)GetAt(i, j int) complex128 {
	return cm.vals[i*cm.numcolumns + j]
}

func (cm *ComplexMatrix)SetAt(i,j int, c complex128) {
	cm.vals[i*cm.numcolumns +j] = c
}


//NewComplexMatrix returns a complex matrix and an error.
//If vals is nil, the matrix will be filled with 0i
//If vals is not nil and has an incorrect length an error will be returned
func NewComplexMatrix(numrows, numcolumns int, vals []complex128) (*ComplexMatrix, error) {
	if vals==nil {
		vals = make([]complex128, numrows*numcolumns)
	}
	if numrows*numcolumns != len(vals) { return nil, WrongValueListLengthError }
	return &ComplexMatrix{numrows:numrows, numcolumns:numcolumns, vals:vals}, nil
}

func NewIdentity(numrows int) *ComplexMatrix {
	vals := make([]complex128, numrows*numrows)
	for i:=0; i<numrows; i++ {
		for j:=0; j<numrows; j++ {
			if i==j {
				vals[i*numrows+j]=1+0i
			} else {
				vals[i*numrows+j]=0i
			}
		}
	}
	n, _ :=NewComplexMatrix(numrows, numrows, vals)
	return n
}

func (cm *ComplexMatrix)checkBounds(i, j int) error {
	if i>=cm.numrows || j>=cm.numcolumns || i<0 || j<0 { return IndexOutOfBoundsError }
	return nil
}

func (cm *ComplexMatrix)SetRow(i int, cm2 *ComplexMatrix) error {
	if cm2.numrows != 1 || cm2.numcolumns!=cm.numcolumns { return IncompatibleMatrixSizesError }
	for j:=0; j<cm.numcolumns; j++ {
		cm.SetAt(i,j, cm2.GetAt(0,j))
	}
	return nil
}

func (cm *ComplexMatrix)GetRow(i int) (*ComplexMatrix, error) {
	if err:=cm.checkBounds(i, 0); err!=nil { return nil, err }
	return NewComplexMatrix(1, cm.numcolumns, cm.vals[i*cm.numcolumns:(i+1)*cm.numcolumns])
}

func (cm *ComplexMatrix)SetColumn(j int, cm2 *ComplexMatrix) error {
	if cm2.numcolumns!=1 || cm2.numrows!=cm.numrows { return IncompatibleMatrixSizesError }
	for i:=0; i<cm.numrows; i++ {
		cm.SetAt(i,j, cm2.GetAt(i,0))
	}
	return nil
}

func (cm *ComplexMatrix)GetColumn(j int) (*ComplexMatrix, error) {
	if err:=cm.checkBounds(0, j); err!=nil { return nil, err }
	vals := make([]complex128, cm.numrows)
	for i:=0; i<cm.numrows; i++ {
		vals[i] = cm.GetAt(i, j)
	}
	return NewComplexMatrix(cm.numrows, 1, vals)
}

func (cm *ComplexMatrix)GetSubMatrix(not_i, not_j int) (*ComplexMatrix, error) {
	if err:=cm.checkBounds(not_i,not_j); err!= nil { return nil, err }
	if cm.numrows<2 || cm.numcolumns<2 { return nil, errors.New("matrix too small") }
	vals := make([]complex128, 0, (cm.numrows-1)*(cm.numcolumns-1))
	for i:=0; i<cm.numrows; i++ {
		for j:=0; j<cm.numcolumns; j++ {
			if i==not_i || j==not_j { continue }
			vals = append(vals, cm.GetAt(i, j))
		}
	}
	return NewComplexMatrix(cm.numrows-1, cm.numcolumns-1, vals)
}

func (cm *ComplexMatrix)isSquare() bool {
	return cm.numrows == cm.numcolumns
}

func (cm *ComplexMatrix)Transpose() *ComplexMatrix {
	vals := make([]complex128, 0, len(cm.vals))
	for j:=0; j<cm.numcolumns; j++ {
		for i:=0; i<cm.numrows; i++ {
			vals=append(vals, cm.GetAt(i, j))
		}
	}
	n, _ := NewComplexMatrix(cm.numcolumns, cm.numrows, vals)
	return n
}

func (cm *ComplexMatrix)GetDeterminant() (complex128, error) {
	if !cm.isSquare() { return 0i, NotSquareError }
	if cm.numrows == 1 { return cm.GetAt(0,0), nil }
	d := 0i
	f := 1+0i
	for i:=0; i<cm.numrows; i++ {
		sm, err := cm.GetSubMatrix(i,0)
		if err != nil { return 0i, err }
		d2, err := sm.GetDeterminant()
		if err != nil { return 0i, err }
		d += f * cm.GetAt(i, 0) * d2
		f = complex(-real(f), 0)
	}
	return d, nil
}

func (cm *ComplexMatrix)GetCofactor() (*ComplexMatrix, error) {
	if !cm.isSquare() { return nil, NotSquareError }
	vals := make([]complex128, cm.numrows * cm.numcolumns)
	f := 1+0i
	for i := 0; i<cm.numrows; i++ {
		if i%2==1 { f = complex(-1, 0) }
		for j:=0; j<cm.numcolumns; j++ {
			sm, err := cm.GetSubMatrix(i,j)
			if err!=nil { return nil, err }
			d, err := sm.GetDeterminant()
			if err!=nil { return nil, err }
			vals[i*cm.numcolumns + j] = f * d //= f * cm.GetAt(i,j) * d
			f = complex(-real(f), 0)
		}
	}
	return NewComplexMatrix(cm.numrows, cm.numcolumns, vals)
}

func (cm *ComplexMatrix)GetAdjugate() (*ComplexMatrix, error) {
	c, err := cm.GetCofactor()
	if err!=nil { return nil, err }
	return c.Transpose(), nil
}

func (cm *ComplexMatrix)Scale(s complex128) *ComplexMatrix {
	vals := make([]complex128, len(cm.vals))
	for i, v := range cm.vals {
		vals[i] = v*s
	}
	n, _ := NewComplexMatrix(cm.numrows, cm.numcolumns, vals)
	return n
}

func (cm *ComplexMatrix)Inverse() (*ComplexMatrix, error) {
	if !cm.isSquare() { return nil, NotSquareError }
	if cm.numrows==1 { return NewComplexMatrix(1,1, []complex128{(1+0i)/cm.GetAt(0,0)}) }
	adj, err := cm.GetAdjugate()
	if err != nil { return nil, err }
	d, err := cm.GetDeterminant()
	if err != nil { return nil, err }
	if d==0i { return nil, DeterminantZeroError }
	return adj.Scale(complex(1,0)/d), nil
}

func (cm *ComplexMatrix)Equals(o *ComplexMatrix) bool {
	if cm.numrows != o.numrows { return false }
	if cm.numcolumns != o.numcolumns { return false }
	for i:=0; i<len(cm.vals); i++ {
		if cm.vals[i]!=o.vals[i] { return false }
	}
	return true
}

func (cm *ComplexMatrix)Times(o *ComplexMatrix) (*ComplexMatrix, error) {
	if cm.numcolumns != o.numrows { return nil, IncompatibleMatrixSizesError }
	retRows := cm.numrows
	retColumns := o.numcolumns
	vals := make([]complex128, retRows*retColumns)
	for cm_i:=0; cm_i<cm.numrows; cm_i++ {
		cm_row, err := cm.GetRow(cm_i)
		if err != nil { return nil, err }
		for o_j:=0; o_j<o.numcolumns; o_j++ {
			o_column, err := o.GetColumn(o_j)
			if err!=nil { return nil, err }
			res := 0i
			for n:=0; n<cm.numcolumns; n++ {
				res += cm_row.GetAt(0, n) * o_column.GetAt(n, 0)
			}
			vals[cm_i*retRows + o_j]=res
		}
	}
	return NewComplexMatrix(retRows, retColumns, vals)
}
