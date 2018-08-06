#include <iostream>
#include <cassert>
#include <cmath>

#include "MatrixGenerator.h"

/**
	generates the identity n x n matrix
	@param n size of the matrix
	@return the identity matrix as Matrix object
*/
Matrix<double> MatrixGenerator::eye(int n)
{
	Matrix<double> M(n,n,0);
	for (int i = 0; i < n; i++){
		M.setElement(i,i,1);
	}
	return M;
}
Matrix<double>* MatrixGenerator::eyePoint(int n) // TODO
{
	Matrix<double>* M = new Matrix<double>(n,n,0);
	for (int i = 0; i < n; i++){
		M->setElement(i,i,1);
	}
	return M;
}

/**
	generates a n x m matrix with random double entries
	@param n amount of rows
	@param m amount of cols
	@return random matrix
*/
Matrix<double> MatrixGenerator::random(int n, int m, int range)
{
	Matrix<double> M(n,m,0);
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
			M.setElement(i,j,rand()%range);
		}
	}
	return M;
}

/**
	generates a matrix that swaps two rows
	@param n size of the matrix
	@param rowA position of the first row
	@param rowB position of the second row
	@return the matrix that swaps the two rows
*/
Matrix<double> MatrixGenerator::swapTwoRows(int n, int rowA, int rowB)
{
	Matrix<double> M(n,n,0);
	for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == rowA) {
                if (j == rowB) {
                    M.setElement(i,j,1);
                }
            }else if (i == rowB){
                if (j == rowA) {
                    M.setElement(i,j,1);
                }
            }else if (i == j){
                M.setElement(i,j,1);
            }
        }
    }
    return M;
}

/**
	rowA' = rowA + scalar*rowB
	@param n size of the matrix
	@param rowA position of rowA
	@param rowB position of rowB
	@param scalar multiplies with that value
	@return the matrix that does the trick
*/
Matrix<double> MatrixGenerator::addOneToAnother(int n, int rowA, int rowB, double scalar)
{
	Matrix<double> M = MatrixGenerator::eye(n);
	M.setElement(rowA,rowB,scalar);
	return M;
}

/**
	rowA' = scalar*rowA
	@param n size of the matrix
	@param row position of the row
	@param scalar multiplies with that value
	@return the matrix that does the trick
*/
Matrix<double> MatrixGenerator::multiplyRowWithScalar(int n, int row, double scalar)
{
	Matrix<double> M = MatrixGenerator::eye(n);
	M.setElement(row,row,scalar);
	return M;
}


/**
	givens rotation, sets the element at i,j to zero
	@param p row
	@param q col
	@param A matrix
	@return the matrix that does the trick
*/
Matrix<double> MatrixGenerator::givensRotation(int p, int q, Matrix<double> A)
{
	int n = A.getDimension()[0];
	Matrix<double> M = MatrixGenerator::eye(n);

	double a_pp = A.getElement(q,q);
	double a_pq = A.getElement(q,p);
	double a_qp = A.getElement(p,q);
	double a_qq = A.getElement(p,p);

	double omega = sqrt(a_pp*a_pp + a_qp*a_qp);
	double sin_phi = - a_qp / omega;
	double cos_phi =   a_pp / omega;

	M.setElement(p,p, cos_phi);
	M.setElement(p,q, sin_phi);
	M.setElement(q,p,-sin_phi);
	M.setElement(q,q, cos_phi);

	return M;
}
