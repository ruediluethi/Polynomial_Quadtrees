// avoid double includes
#ifndef _MATRIXGENERATOR_H_
#define _MATRIXGENERATOR_H_

#include "Matrix.h"

class MatrixGenerator{
public:
	static Matrix<double> eye(int);
	static Matrix<double>* eyePoint(int);
	static Matrix<double> random(int, int, int);
	static Matrix<double> swapTwoRows(int, int, int);
	static Matrix<double> addOneToAnother(int, int, int, double);
	static Matrix<double> multiplyRowWithScalar(int, int, double);
	static Matrix<double> givensRotation(int, int, Matrix<double>);
};

#endif