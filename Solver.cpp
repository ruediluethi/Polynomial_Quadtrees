#include <iostream>
#include <cassert>
#include <cmath>

#include "Matrix.cpp" // to avoid linking errors
#include "Solver.h"

// Ax = (QR)x = b 
// Q'*Q = 1  =>  Rx = Q'b
Matrix<double> Solver::QRsolve(Matrix<double> A, Matrix<double> b)
{
	int n = A.getDimension()[0];

	Matrix<double> R = A;
	Matrix<double> Q = MatrixGenerator::eye(n);

	// generate QR decomposition
	for (int i = 1; i < n; i++){
		for (int j = 0; j < i; j++){

			Matrix<double> G = MatrixGenerator::givensRotation(i,j,R);
			R = G*R;
			G.transpose();
			Q = Q*G;

		}
	}

	Matrix<double> Qt = Q;
	Qt.transpose();

	b = Qt*b; // => Rx = Q'b

	// reverse substitution
	Matrix<double> x(n,1,0);
	for (int i = n-1; i >= 0; i--){
		double xi = b.getElement(i,0);
		for (int j = n-1; j >= i; j--){
			xi -= R.getElement(i,j)*x.getElement(j,0);
		}
		xi = xi/R.getElement(i,i);
		x.setElement(i,0,xi);
	}

	/*
	// DEBUG
	std::cout << "Q = " << std::endl;
	Q.display();

	std::cout << "R = " << std::endl;
	R.display();

	A = Q*R;
	std::cout << "A = QR = " << std::endl;
	A.display();

	std::cout << "Q'Q = " << std::endl;
	(Qt*Q).display();

	std::cout << "Q'b = " << std::endl;
	(Qt*b).display();
	*/

	return x;
}


/**
 approximates the vector b with a polynomial function of grade g
 @param b a nx1 matrix (vector)
 @param g the order of the polynomial function
 @return the approximated vector
 */
Matrix<double> Solver::polynomialApprox(Matrix<double> b, int g)
{
	int n = b.getDimension()[0];
	Matrix<double> polynom(n,1,0);

	// generate matrix A (needed for the system of linear equations)
	Matrix<double> A(n,g+1,1);
	for (int i = 0; i < n; i++){
		for (int j = 1; j < g+1; j++){
			double v = i+1;
			for (int p = 0; p < j-1; p++){
				v *= i+1;
			}
			A.setElement(i,g-j,v);
		}
	}

	// transpose A
	Matrix<double> At = A;
	At.transpose();

	// least squares approximation
	A = At*A; 
	b = At*b;
	Matrix<double> coeff = Solver::QRsolve(A,b); // coefficients of the polynom

	// calc the polynom
	for (int x = 0; x < n; x++){
		
		// y = f(x) = c_1*x^(n-1) + c_2*x^(n-2) + ... + c_(n-1)*x + c_n
		double y = 0;
		for (int i = 0; i < g+1; i++){
			double c = coeff.getElement(i,0);
			double powX = x;
			for (int p = 0; p < g-i-1; p++){
				powX *= x;
			}
			if (i == g){ // exception for x^0 = 1
				y += c;
			}else{
				y += c*powX;
			}
		}

		polynom.setElement(x,0,y);
	}

	return polynom;
}

/**
 interpolates difference between vector a and vector b
 @param a nx1 matrix (vector) of length n
 @param b nx1 matrix (vector) of the same length as a
 @param t value between 0 and 1 (0 => output = a, 1 => output = b)
 */
Matrix<double> Solver::interpolate(Matrix<double> a, Matrix<double> b, double t)
{
	int n = a.getDimension()[0];
	Matrix<double> output(n,1,0);
	for (int i = 0; i < n; i++){
		double d = b.getElement(i,0) - a.getElement(i,0);
		output.setElement(i,0,a.getElement(i,0) + d*t);
	}
	return output;
}