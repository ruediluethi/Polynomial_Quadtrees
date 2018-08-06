// avoid double includes
#ifndef _MATRIX_H_
#define _MATRIX_H_

template <class T>
class Matrix{
public:
	Matrix();
	Matrix(int, int);
	Matrix(int, int, T);

	// rule of three
	Matrix(const Matrix<T>&); // copy
	Matrix<T>& operator=(const Matrix<T>&); // assignment
	~Matrix(); // destructor
	// TODO: Move-Konstruktor und Move-Zuweisung

	void display();
	void setElement(int, int, T);
	T getElement(int, int) const;
	Matrix<T> getBlock(int, int, int, int) const;
	int* getDimension() const;
	void transpose();
	void removeCol(const int);
	T determinant() const;
	T sum() const;
	T average() const;

protected:
	int countRows;
	int countCols;
	T* elements;
};


// addition
template <class T>
Matrix<T>& operator+= (Matrix<T>&, Matrix<T>);
template <class T>
Matrix<T> operator+ (Matrix<T>, Matrix<T>);

// difference
template <class T>
Matrix<T>& operator+= (Matrix<T>&, Matrix<T>);
template <class T>
Matrix<T> operator+ (Matrix<T>, Matrix<T>);

// multiplication with scalar
template <class T>
Matrix<T>& operator*= (Matrix<T>&, T);
template <class T>
Matrix<T> operator* (Matrix<T>, T);

// matrix multiplication
template <class T>
Matrix<T> operator* (Matrix<T>, Matrix<T>);

#endif