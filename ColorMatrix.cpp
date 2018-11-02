#include <iostream>
#include <cassert>
#include <cmath>
#include <float.h>

#include "ColorMatrix.h"
#include "Color.h"
#include "Solver.h"

ColorMatrix::ColorMatrix()
{

}

ColorMatrix::ColorMatrix(Matrix<Color> M)
{
	countRows = M.getDimension()[0];
	countCols = M.getDimension()[1];
	elements = new Color[countRows*countCols];

	for (int i = 0; i < countRows; i++){
		for (int j = 0; j < countCols; j++){
			setElement(i,j,M.getElement(i,j));
		}
	}
}

/**
 calculates the distance for every pixel to a specific color
 @param color the color for comparison
 @return a matrix with the distance for every pixel as double value
 */
Matrix<double> ColorMatrix::distance(Color color)
{
	Matrix<double> M(countRows,countCols,0);

	for (int i = 0; i < countRows; i++){
		for (int j = 0; j < countCols; j++){
			M.setElement(i,j,getElement(i,j).distance(color));
		}
	}

	return M;
}

/**
 calculates the distance for every pixel to another matrix
 @param A the matrix for comparison
 @return a matrix with the distance for every pixel as double value
 */
Matrix<double> ColorMatrix::distance(ColorMatrix A)
{
	Matrix<double> M(countRows,countCols,0);

	for (int i = 0; i < countRows; i++){
		for (int j = 0; j < countCols; j++){
			M.setElement(i,j,getElement(i,j).distance(A.getElement(i,j)));
		}
	}

	return M;
}

ColorMatrix ColorMatrix::generateQuadtree(double levelOfDetails, int maxSize, int minSize, int zeros)
{
	ColorMatrix M(Matrix<Color>(countRows,countCols,Color(0,0,0)));

	std::cout << countRows << "x" << countCols << " : " << std::endl;

	bool doSplit = false;
	if (countRows > maxSize || countCols > maxSize){
		doSplit = true;
	}

	if (!doSplit){
		Color avgColor = average();

		ColorMatrix Pcol = ColorMatrix(Matrix<Color>(countRows,countCols,avgColor));
		ColorMatrix Prow = ColorMatrix(Matrix<Color>(countRows,countCols,avgColor));
		if (zeros > -1){ // do polynomial split only if zeros > -1
			Pcol = polynomialApprox(split(true), zeros, true);
			Prow = polynomialApprox(split(false), zeros, false);
		}
		// calc distance
		Matrix<double> distToAvgCol = Pcol.distance(*this);
		Matrix<double> distToAvgRow = Prow.distance(*this);
		double avgCol = distToAvgCol.average();
		double avgRow = distToAvgRow.average();
		double avg = avgCol;
		if (avgRow < avgCol){
			avg = avgRow;
			std::cout << "row-"; // DEBUG
		}else{
			std::cout << "col-"; // DEBUG
		}
		std::cout << "avg = " << avg; // DEBUG
		if (avg < levelOfDetails || countRows <= minSize || countCols <= minSize){
			std::cout << " => paint" << std::endl;
			if (avgRow < avgCol){
				M = Prow;
			}else{
				M = Pcol;
			}
		}else{
			doSplit = true;
		}
	}

	
	if (doSplit){
		std::cout << " => split" << std::endl;
		ColorMatrix UL = getBlock(0,0,countRows/2,countCols/2); // upper left
		ColorMatrix UR = getBlock(0,countCols/2,countRows/2,countCols/2); // upper right
		ColorMatrix LL = getBlock(countRows/2,0,countRows/2,countCols/2); // lower left
		ColorMatrix LR = getBlock(countRows/2,countCols/2,countRows/2,countCols/2); // lower right
		// do the same with the leafs
		UL = UL.generateQuadtree(levelOfDetails, maxSize, minSize, zeros);
		UR = UR.generateQuadtree(levelOfDetails, maxSize, minSize, zeros);
		LL = LL.generateQuadtree(levelOfDetails, maxSize, minSize, zeros);
		LR = LR.generateQuadtree(levelOfDetails, maxSize, minSize, zeros);
		// copy all together
		// upper left
		for (int i = 0; i < countRows/2; i++){
			for (int j = 0; j < countCols/2; j++){
				M.setElement(i,j,UL.getElement(i,j));
			}
		}
		// upper right
		for (int i = 0; i < countRows/2; i++){
			for (int j = 0; j < countCols/2; j++){
				M.setElement(i,countCols/2+j,UR.getElement(i,j));
			}
		}
		// lower left
		for (int i = 0; i < countRows/2; i++){
			for (int j = 0; j < countCols/2; j++){
				M.setElement(countRows/2+i,j,LL.getElement(i,j));
			}
		}
		// lower right
		for (int i = 0; i < countRows/2; i++){
			for (int j = 0; j < countCols/2; j++){
				M.setElement(countRows/2+i,countCols/2+j,LR.getElement(i,j));
			}
		}
	}

	return M;
}

/**
 splits the matrix into two areas
 @param columnWise true: split column wise, false: split row wise
 @return a vector with the split position
 */
Matrix<double> ColorMatrix::split(bool columnWise)
{
	int xAxis = countRows;
	int yAxis = countCols;
	if (columnWise){
		std::cout << "column wise: ";
		xAxis = countCols;
		yAxis = countRows;
	}else{
		std::cout << "row wise: ";
	}

	int progress = 0;
	Matrix<double> b(xAxis,1,0); // this vector will be used later for the solver
	for (int k = 0; k < xAxis; k++){
		Matrix<Color> oneLine(yAxis,1,Color(0,0,0));
		if (columnWise){
			oneLine = getBlock(0,k,countRows,1);
		}else{
			oneLine = getBlock(k,0,1,countCols);
			oneLine.transpose();
		}
		double minSum = DBL_MAX;
		int minCutPos = -1;

		// brute force search for the best split
		for (int l = 0; l < yAxis; l++){
			// split
			ColorMatrix partLeft(oneLine.getBlock(0,0,l,1));
            ColorMatrix partRight(oneLine.getBlock(l,0,yAxis-l,1));

            // calc left difference to average
            Matrix<double> diffLeft = partLeft.distance(partLeft.average());
            // then calc the right side
            Matrix<double> diffRight = partRight.distance(partRight.average());

            // sum up
            double sumLeft = diffLeft.sum();
            double sumRight = diffRight.sum();
            if (sumLeft+sumRight < minSum){ // compare with other results
            	minSum = sumLeft + sumRight;
            	minCutPos = l;
            }
		}

		// save the best split in the vector b
		b.setElement(k,0,(double)minCutPos);

		// show progress in console (in steps of 10%)
		int crntProgress = round((double)k/(double)countCols*10);
		if (crntProgress > progress){
			progress = crntProgress;
			std::cout << (progress*10) << "% " << std::flush;
		}
	}
	std::cout << std::endl;

	return b;
}


/* DEPRECATED */
Matrix<double> ColorMatrix::split()
{
	int progress = 0;
	Matrix<double> b(countCols,1,0); // this vector will be used later for the solver
	for (int k = 0; k < countCols; k++){ // loop column wise through the image
		Matrix<Color> oneLine = getBlock(0,k,countRows,1);
		double minSum = DBL_MAX;
		int minCutPos = -1;

		/*
		TODO: this part does not work ???
		Color avgLeft = oneLine.average();
		Color avgRight(0.0,0.0,0.0);
		for (int l = 1; l < countRows; l++){
			
			// avg can be calced directly (no loop needed)
			double p = (double)(countRows-l);
			double q = (double)(l);
			Color ap = oneLine.getElement(p,k);
			avgLeft = (avgLeft*(p+1) - ap)*(1/p);
			avgRight = (avgRight*(q-1) + ap)*(1/q);

			// split			
			ColorMatrix partLeft(oneLine.getBlock(0,0,(int)p,1));
            ColorMatrix partRight(oneLine.getBlock((int)p,0,(int)q,1));

            // DEBUG
            // std::cout << q << ": [left] " << avgLeft << std::endl;
            // std::cout << "avg left: " << partLeft.average() << std::endl << std::endl;
            // std::cout << q << ": [right] " << avgRight << std::endl;
            // std::cout << "avg right: " << partRight.average() << std::endl << std::endl;
			
            // TODO: calc this part directly (not sure if its possible)
            // calc left difference to average
            Matrix<double> diffLeft = partLeft.distance(avgLeft);
            // then calc the right side
            Matrix<double> diffRight = partRight.distance(avgRight);

            // sum up
            double sumLeft = diffLeft.sum();
            double sumRight = diffRight.sum();
            if (sumLeft+sumRight < minSum){ // compare with other results
            	minSum = sumLeft + sumRight;
            	minCutPos = l;
            }
		}
		*/

		// brute force search for the best split
		for (int l = 0; l < countRows+1; l++){
			// split
			ColorMatrix partLeft(oneLine.getBlock(0,0,l,1));
            ColorMatrix partRight(oneLine.getBlock(l,0,countRows-l,1));

            // calc left difference to average
            Matrix<double> diffLeft = partLeft.distance(partLeft.average());
            // then calc the right side
            Matrix<double> diffRight = partRight.distance(partRight.average());

            // sum up
            double sumLeft = diffLeft.sum();
            double sumRight = diffRight.sum();
            if (sumLeft+sumRight < minSum){ // compare with other results
            	minSum = sumLeft + sumRight;
            	minCutPos = l;
            }
		}

		// save the best split in the vector b
		b.setElement(k,0,(double)minCutPos);

		// show progress in console (in steps of 10%)
		int crntProgress = round((double)k/(double)countCols*10);
		if (crntProgress > progress){
			progress = crntProgress;
			std::cout << (progress*10) << "% " << std::flush;
		}
	}
	std::cout << std::endl;

	return b;
}

/**
 @param b vector with the positions of the best split
 @param zeros order of the polynomial function 
              -1 => only average color is used, 0 => constant, 1 => line, 2 => squared function, ...
 @return image splitted by a polynomial function into two areas with average color
 */
ColorMatrix ColorMatrix::polynomialApprox(Matrix<double> b, int zeros, bool columnWise)
{
	int xAxis = countRows;
	int yAxis = countCols;
	if (columnWise){
		xAxis = countCols;
		yAxis = countRows;
	}

	Matrix<int> output(yAxis,xAxis,0);

	// generate matrix A (needed for the system of linear equations)
	Matrix<double> A(xAxis,zeros+1,1);
	for (int i = 0; i < xAxis; i++){
		for (int j = 1; j < zeros+1; j++){
			double v = i+1;
			for (int p = 0; p < j-1; p++){
				v *= i+1;
			}
			A.setElement(i,zeros-j,v);
		}
	}

	// transpose A
	Matrix<double> At = A;
	At.transpose();

	// least squares approximation
	A = At*A; 
	b = At*b;
	Matrix<double> coeff = Solver::QRsolve(A,b); // coefficients of the polynom

	// paint the polynom
	int countPixelsInA = 0;
	for (int x = 0; x < xAxis; x++){
		
		// f(x) = c_1*x^(n-1) + c_2*x^(n-2) + ... + c_(n-1)*x + c_n
		double fx = 0;
		for (int i = 0; i < zeros+1; i++){
			double c = coeff.getElement(i,0);
			double powX = x;
			for (int p = 0; p < zeros-i-1; p++){
				powX *= x;
			}
			if (i == zeros){ // exception for x^0 = 1
				fx += c;
			}else{
				fx += c*powX;
			}
		}

		for (int y = 0; y < fx && y < yAxis; y++){
			if (columnWise){
				output.setElement(y,x,1);
			}else{
				output.setElement(x,y,1);
			}
			countPixelsInA++;
		}
	}

	// calc average colors for the two areas A/B
	Matrix<Color> areaA(countPixelsInA,1,Color(0,0,0));
	Matrix<Color> areaB(yAxis*xAxis-countPixelsInA,1,Color(0,0,0));
	int aCount = 0;
	int bCount = 0;
	for (int i = 0; i < yAxis; i++){
		for (int j = 0; j < xAxis; j++){
			if (output.getElement(i,j) == 1){
				areaA.setElement(aCount,0,getElement(i,j)); // save all pixels in area A
				aCount++;
			}else{
				areaB.setElement(bCount,0,getElement(i,j)); // save all pixels in area B
				bCount++;
			}
		}
	}
	Color avgColorA = areaA.average();
	Color avgColorB = areaB.average();
	// paint color for output
	Matrix<Color> coloredOutput(countRows,countCols,Color(1,1,1));
	for (int i = 0; i < yAxis; i++){
		for (int j = 0; j < xAxis; j++){
			if (output.getElement(i,j) == 1){
				coloredOutput.setElement(i,j,avgColorA);
			}else{
				coloredOutput.setElement(i,j,avgColorB);
			}
		}
	}

	return ColorMatrix(coloredOutput);
}

ColorMatrix ColorMatrix::diffuse(double factor)
{
	ColorMatrix M(*this);

	for (int i = 1; i < countRows-1; i++){
		for (int j = 1; j < countCols-1; j++){

			int l = j-1;
			int r = j+1;
			int u = i-1;
			int d = i+1;

			Color pixelValue = getElement(i,j);

			pixelValue += (getElement(i,l) - getElement(i,j))*factor;
			pixelValue += (getElement(i,r) - getElement(i,j))*factor;
			pixelValue += (getElement(u,j) - getElement(i,j))*factor;
			pixelValue += (getElement(d,j) - getElement(i,j))*factor;

			double* rgb = pixelValue.getRGB();
			if (rgb[0] > 1.0){
				rgb[0] = 1.0;
			}
			if (rgb[1] > 1.0){
				rgb[1] = 1.0;
			}
			if (rgb[2] > 1.0){
				rgb[2] = 1.0;
			}
			pixelValue.setRGB(rgb[0],rgb[1],rgb[2]);

			M.setElement(i,j,pixelValue);
		}
	}
	return M;
}