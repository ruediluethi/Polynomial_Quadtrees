// avoid double includes
#ifndef _COLORMATRIX_H_
#define _COLORMATRIX_H_

#include "Matrix.h"
#include "Matrix.cpp" // to avoid linking errors
#include "Color.h"

class ColorMatrix : public Matrix<Color>
{
public:
	ColorMatrix();
	ColorMatrix(Matrix<Color>);

	Matrix<double> distance(Color);
	Matrix<double> distance(ColorMatrix);

	ColorMatrix generateQuadtree(double,int,int,int);
	Matrix<double> split();
	Matrix<double> split(bool);
	ColorMatrix polynomialApprox(Matrix<double>,int,bool);
	ColorMatrix diffuse(double);

private:
	
};

#endif