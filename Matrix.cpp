// avoid double includes
#ifndef _MATRIX_CPP_
#define _MATRIX_CPP_

#include <iostream>
#include <cassert>

#include "Matrix.h"

template <class T>
Matrix<T>::Matrix(): countRows(0), countCols(0)
{
	
}

template <class T>
Matrix<T>::Matrix(int n, int m): countRows(n), countCols(m), elements(new T[n*m])
{
	
}

template <class T>
Matrix<T>::Matrix(int n, int m, T defaultValue): countRows(n), countCols(m), elements(new T[n*m])
{
	for (int i = 0; i < countRows; i++){
		for (int j = 0; j < countCols; j++){
			elements[i*countCols+j] = defaultValue;
		}
	}
}

// Copy constructor
template <class T>
Matrix<T>::Matrix(const Matrix<T>& M):
	countRows(M.getDimension()[0]),
	countCols(M.getDimension()[1]),
	elements(new T[M.getDimension()[0]*M.getDimension()[1]])
{
	std::copy(M.elements, M.elements+(countRows*countCols), elements);
}

// Assignment operator
template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& M)
{
	if (M.elements == elements) return *this;
	if (elements != nullptr) delete[] elements;
	countRows = M.countRows;
	countCols = M.countCols;
	elements = new T[countRows*countCols];
	std::copy(M.elements, M.elements+(countRows*countCols), elements);
	return *this;
}

//destructor
template <class T>
Matrix<T>::~Matrix()
{
	delete[] elements;
}

// simple output in the console
template <class T>
void Matrix<T>::display()
{
	for (int i = 0; i < countRows; i++){
		for (int j = 0; j < countCols; j++){
			//std::cout << elements[i*countCols+j] << ' ';
			std::cout << getElement(i,j) << "  ";
			//printf("%6.2f  ",round(getElement(i,j)*100.0)/100.0);
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

// set one single element of the matrix
template <class T>
void Matrix<T>::setElement(int i, int j, T value)
{
	elements[i*countCols+j] = value;
}


/**
	read a single element of the matrix
	@param i row
	@param j column
	@return value of the element at i,j
 */
template <class T>
T Matrix<T>::getElement(int i, int j) const
{
	return elements[i*countCols+j];
}

// gets a block as matrix from the matrix
template <class T>
Matrix<T> Matrix<T>::getBlock(int k, int l, int rows, int cols) const
{
	assert(k+rows <= countRows && l+cols <= countCols);

	Matrix *B = new Matrix<T>(rows,cols);
	for (int i = k; i < k+rows; i++){
		for (int j = l; j < l+cols; j++){
			B->setElement(i-k,j-l,(*this).getElement(i,j));
		}
	}
	return *B;
}

// get amount of rows and cols
// returns a array with the two values
template <class T>
int* Matrix<T>::getDimension() const
{
	int* dimensions = new int[2];
	dimensions[0] = countRows;
	dimensions[1] = countCols;
	return dimensions;
}


// transposes the matrix
// the amount of rows becomes the amount of cols and cols to rows
template <class T>
void Matrix<T>::transpose()
{
	// create a copy of the matrix elements
	T* elementsCopy = new T[countCols*countRows];
	std::copy(elements, elements+(countRows*countCols), elementsCopy);

	// flip row- to column-elements and cols to rows
	for (int i = 0; i < countRows; i++){
		for (int j = 0; j < countCols; j++){
			elements[j*countRows+i] = elementsCopy[i*countCols+j];
		}
	}

	int temp = countRows;
	countRows = countCols;
	countCols = temp;
}

// removes one column at position k
template <class T>
void Matrix<T>::removeCol(const int k)
{
	// create a copy of the matrix elements
	T* elementsCopy = new T[countRows*countCols];
	std::copy(elements, elements+(countRows*countCols), elementsCopy);

	elements = new T[countRows*(countCols-1)]; // less space needed, because one column will be deleted

	// copy the old elements into the new elements space (without col k)
	for (int i = 0; i < countRows; i++){
		for (int j = 0; j < countCols; j++){
			if (j < k){
				elements[i*(countCols-1)+j] = elementsCopy[i*countCols+j];
			}
			if (j > k){
				elements[i*(countCols-1)+(j-1)] = elementsCopy[i*countCols+j];
			}
		}
	}

	countCols = countCols-1;
}


// calculate the determinant of the matrix
template <class T>
T Matrix<T>::determinant() const
{
	assert(countRows == countCols); // determinant is only defined for n x n

	if (countRows == 1){
		return getElement(0,0);
	}

	if (countRows == 2){
		return getElement(0,0)*getElement(1,1) - getElement(0,1)*getElement(1,0);
	}

	// always develop the first row
	Matrix<T> B = getBlock(1,0,countRows-1,countCols);
	T det = 0;
	for (int i = 0; i < countCols; i++){
		Matrix A = B;
		A.removeCol(i);
		int sign = 1-(i%2)*2;
		det = det + sign*getElement(0,i)*A.determinant();
	}

	return det;
}


// sum of all elements
template <class T>
T Matrix<T>::sum() const
{
	T sum = elements[0];
	for (int i = 1; i < countRows*countCols; i++){
		sum += elements[i];
	}
	return sum;
}

// calculate simmple average: sum of all elements divided with size
template <class T>
T Matrix<T>::average() const
{
	T avg = sum();
	avg *= 1.0/(double)(countRows*countCols);
	return avg;
}

// addition
template <class T>
Matrix<T>& operator+= (Matrix<T>& A, Matrix<T> B)
{
	int* dimA = A.getDimension();
	int* dimB = B.getDimension();
	assert(dimA[0] == dimB[0] && dimA[1] == dimB[1]); // matrix dimensions must be the same!

	for (int i = 0; i < dimA[0]; i++){
		for (int j = 0; j < dimA[1]; j++){
			A.setElement(i,j,A.getElement(i,j)+B.getElement(i,j));
		}
	}

	return A;
}

template <class T>
Matrix<T> operator+ (Matrix<T> A, Matrix<T> B){
	Matrix<T> C = A;
	C += B;
	return C;
}

// difference
template <class T>
Matrix<T>& operator-= (Matrix<T>& A, Matrix<T> B)
{
	int* dimA = A.getDimension();
	int* dimB = B.getDimension();
	assert(dimA[0] == dimB[0] && dimA[1] == dimB[1]); // matrix dimensions must be the same!

	for (int i = 0; i < dimA[0]; i++){
		for (int j = 0; j < dimA[1]; j++){
			A.setElement(i,j,A.getElement(i,j) - B.getElement(i,j));
		}
	}

	return A;
}

template <class T>
Matrix<T> operator- (Matrix<T> A, Matrix<T> B){
	Matrix<T> C = A;
	C -= B;
	return C;
}

// scalar multipilication
template <class T>
Matrix<T>& operator*= (Matrix<T>& A, T scalar){
	int* dimA = A.getDimension();
	for (int i = 0; i < dimA[0]; i++){
		for (int j = 0; j < dimA[1]; j++){
			A.setElement(i,j,A.getElement(i,j)*scalar);
		}
	}
	return A;
}
template <class T>
Matrix<T> operator* (Matrix<T> A, T scalar){
	Matrix<T> B = A;
	B *= scalar;
	return B;
}

// matrix multiplication
template <class T>
Matrix<T> operator* (Matrix<T> A, Matrix<T> B){
	int* dimA = A.getDimension();
	int* dimB = B.getDimension();
	assert(dimA[1] == dimB[0]); // amount of rows of A must be the same as amount of cols of B

	Matrix<T> C(dimA[0],dimB[1],0);
	for (int i = 0; i < dimA[0]; i++){
		for (int j = 0; j < dimB[1]; j++){
			T value = 0;
			for (int k = 0; k < dimA[1]; k++){
				value += A.getElement(i,k) * B.getElement(k,j);
			}
			C.setElement(i,j,value);
		}
	}

	return C;
}


#endif