#include "Matrix.h"
#include "MatrixGenerator.h"

class Solver{
public:
	static Matrix<double> QRsolve(Matrix<double>, Matrix<double>);
	static Matrix<double> polynomialApprox(Matrix<double>, int);
	static Matrix<double> interpolate(Matrix<double>, Matrix<double>, double);
};