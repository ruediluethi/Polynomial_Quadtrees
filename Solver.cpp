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