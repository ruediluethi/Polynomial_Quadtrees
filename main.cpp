#include <iostream>
#include <time.h>
#include <cmath>
#include <float.h>
#include <sstream>

#include "Color.h"
#include "Matrix.h"
#include "Matrix.cpp" // to avoid linking errors
#include "ColorMatrix.h"
#include "MatrixGenerator.h"
#include "IOcolorMatrix.h"
#include "Solver.h"


Matrix<double> calcColorDistance(Matrix<Color> A, Matrix<Color> B)
{
	int* dimA = A.getDimension();
	int* dimB = B.getDimension();
	assert(dimA[0] == dimB[0] && dimA[1] == dimB[1]); // matrix dimensions must be the same!

	Matrix<double> D(dimA[0],dimA[1],-1);

	for (int i = 0; i < dimA[0]; i++){
		for (int j = 0; j < dimA[0]; j++){
			D.setElement(i,j,A.getElement(i,j).distance(B.getElement(i,j)));
		}
	}

	return D;
}


int main(int argc, char* argv[])
{
	if (argc < 2){
		std::cout << "ERROR: the first two parameters are needed! first one is the input path, second is output path" << std::endl;
		return 0;
	}

	char* inputPath = argv[1];
	char* outputPath = argv[2];

	// DEBUG
	// for (int i = 0; i < argc; i++){
	// 	std::cout << i << ": " << argv[i] << std::endl;
	// }

	double levelOfDetails = 0.3;
	if (argc > 3){
		levelOfDetails = atof(argv[3]);
	}

	int zeros = 2;
	if (argc > 4){
		zeros = atoi(argv[4]);
	}

	// load image from file
	ColorMatrix sourceImg(IOcolorMatrix::readFromFile(inputPath));

	int maxSize = sourceImg.getDimension()[0]/4;
	if (argc > 5){
		maxSize = atoi(argv[5]);
	}

	int minSize = sourceImg.getDimension()[0]/16;;
	if (argc > 6){
		minSize = atoi(argv[6]);
	}

	std::cout << "input path:       " << inputPath << std::endl;
	std::cout << "output path:      " << outputPath << std::endl;
	std::cout << "level of details: " << levelOfDetails << std::endl;
	std::cout << "polynomial order: " << zeros << std::endl;
	std::cout << "max size:         " << maxSize << std::endl;
	std::cout << "min size:         " << minSize << std::endl;

	// start the process
	ColorMatrix quadTreeImg = sourceImg.generateQuadtree(levelOfDetails, maxSize, minSize, zeros);

	// write to disk
	IOcolorMatrix::writeToFile(quadTreeImg, outputPath);


	// for polynomial interpolate animation
	/*
	Matrix<double> split = sourceImg.split(false);
	std::cout << std::endl;

	ColorMatrix debugImg = sourceImg.polynomialApprox(split,2,false);

	ColorMatrix debugImg(Matrix<Color>(sourceImg.getDimension()[0],sourceImg.getDimension()[1],Color(0,0,0)));
	for (int j = 0; j < sourceImg.getDimension()[1]; j++){
		for (int i = 0; i < (int)round(split.getElement(j,0)); i++){
			debugImg.setElement(i,j,Color(1,1,1));
		}
	}

	IOcolorMatrix::writeToFile(debugImg, "output/test_row.bmp");

	// for (int i = 0; i < 3; i++){
	// 	debug = debug.diffuse(0.1);
	// }

	// write to file
	std::ostringstream debugFile;
	debugFile << "split.bmp";
	std::string debugString = outputPath + debugFile.str();
	char* debugChars = new char[debugString.length()];
	strcpy(debugChars, debugString.c_str());
	std::cout << "write: [" << debugChars << "]" << std::endl;
	IOcolorMatrix::writeToFile(debugImg, debugChars);

	int fileCount = 0;
	for (int p = 0; p < 20; p++){
		
		std::cout << std::endl << "p = " << p;
		Matrix<double> polyA = Solver::polynomialApprox(split, p);
		Matrix<double> polyB = Solver::polynomialApprox(split, p+1);
		std::cout << std::endl;
		
		int interpolateLength = 12;
		for (int t = 0; t < interpolateLength; t++){
			Matrix<double> inbetween = Solver::interpolate(polyA, polyB, (double)t/(double)interpolateLength);

			// ColorMatrix splittedImg(Matrix<Color>(sourceImg.getDimension()[0],sourceImg.getDimension()[1],Color(1,1,1)));
			ColorMatrix splittedImg(Matrix<Color>(sourceImg.getDimension()[0],sourceImg.getDimension()[1],Color(87.0/256.0,96.0/256.,58.0/256.0)));
			for (int i = 0; i < sourceImg.getDimension()[0]; i++){
				for (int j = 0; j < sourceImg.getDimension()[1]; j++){
					if (i < (int)round(inbetween.getElement(j,0))){
						splittedImg.setElement(i,j,Color(1,1,1));
					}
				}
			}

			// for (int i = 0; i < 3; i++){
			// 	splittedImg = splittedImg.diffuse(0.1);
			// }

			// save frame
			fileCount++;
			std::ostringstream fileString;
			fileString << fileCount;
			std::string outputString = outputPath + fileString.str() + ".bmp";
			char* outputChars = new char[outputString.length()];
			strcpy(outputChars, outputString.c_str());
			std::cout << "write: [" << outputChars << "]" << std::endl;
			IOcolorMatrix::writeToFile(splittedImg, outputChars);
		}
	}
	*/


	return 0;


	/*
	// std::string fileName = "berg_pano";
	// std::string fileName = "eigermoenchjungfrau128";
	// std::string fileName = "eigermoenchjungfrau_pano";
	std::string fileName = "abendrot";
	// std::string fileName = "mindelheim";
	// std::string fileName = "eigermoenchjungfrau512";


	// input
	std::string inputPath = "input/" + fileName + ".bmp";
	char* inputChars = new char[inputPath.length()];
	strcpy(inputChars, inputPath.c_str());
	std::cout << "input: [" << inputChars << "]" << std::endl;
	// load image from file
	ColorMatrix sourceImg(IOcolorMatrix::readFromFile(inputChars));

	// calc
	// ColorMatrix polynomImg = sourceImg.polynomialSplit(10);
	ColorMatrix quadTreeImg = sourceImg.generateQuadtree(0.1, sourceImg.getDimension()[0]/2, 8, 2);

	// output
	std::string outputPath = "output/" + fileName + ".bmp";
	char* outputChars = new char[outputPath.length()];
	strcpy(outputChars, outputPath.c_str());
	std::cout << "output: [" << outputChars << "]" << std::endl;
	// write to file
	IOcolorMatrix::writeToFile(quadTreeImg, outputChars);
	*/

	/*
	// TEST: for I/O

	// read color matrix from file
	Matrix<Color> sourceImg = IOcolorMatrix::readFromFile("input/eigermoenchjungfrau.bmp");

	// calc average of image
	Matrix<Color> avgImg(sourceImg.getDimension()[0],sourceImg.getDimension()[1],sourceImg.average());

	// calc difference
	Matrix<Color> diffImg = sourceImg - avgImg;

	// save color matrix to file
	IOcolorMatrix::writeToFile(avgImg, "output/diffemj.bmp");
	*/
	

	/*
	// TEST: for row manipulation
	srand(42); // set the random seed
	int n = 5;

	Matrix<double> A = MatrixGenerator::random(n, 90);
	A += Matrix<double>(n,n,10);
	A.display();
	std::cout << "------------------" << std::endl;

	Matrix<double> S = MatrixGenerator::swapTwoRows(n,1,3);
	A = S*A;
	A.display();
	std::cout << "------------------" << std::endl;

	Matrix<double> T = MatrixGenerator::addOneToAnother(n,2,1,-(R.getElement(2,0)/R.getElement(1,0)));
	A = T*A;
	A.display();
	std::cout << "------------------" << std::endl;

	Matrix<double> M = MatrixGenerator::multiplyRowWithScalar(n,2,(1/R.getElement(2,2)));
	A = M*A;
	A.display();
	*/


	/*
	// TEST: for solver with givens rotation
	// Ax = (QR)x = b 
	// Q'*Q = 1  =>  Rx = Q'b

	srand(42); // set the random seed
	int n = 10;

	Matrix<double> A = MatrixGenerator::random(n, n, 90);
	A += Matrix<double>(n,n,10);
	std::cout << "A = " << std::endl;
	A.display();

	Matrix<double> b = MatrixGenerator::random(n, 1, 90);
	std::cout << "b = " << std::endl;
	b.display();

	Matrix<double> x = Solver::QRsolve(A,b);
	std::cout << "x = " << std::endl;
	x.display();
	*/

	
	/*
	// more TESTs
	int n = 5;
	int r = 10;

	Matrix<double> *A = new Matrix<double>(n,n,-1);

	int* dimA = A->getDimension();
	for (int i = 0; i < dimA[0]; i++){
		for (int j = 0; j < dimA[1]; j++){
			//A->setElement(i,j,i*dimA[1]+j);
			// A->setElement(i,j,rand()%(r*2)-r);
			double x = ((double)(i*dimA[1]+j))/((double)(n*n));
			A->setElement(i,j,x);
		}
	}

	A->display();
	
	std::cout << "sum(A) = " << A->sum() << std::endl;
	std::cout << "avg(A) = " << A->average() << std::endl;
	*/

	/*
	A->transpose();
	Matrix<double> B = A->getBlock(0,0,n/2+1,n/2+1);
	B.removeCol(n/4);
	B.display();

	std::cout << "det(A) = " << A->determinant() << std::endl;

	Matrix<double> *C = new Matrix<double>(n,n,0);
	*C += *A;
	A->transpose();
	*C = *C + *A;

	C->display();

	*C *= -1.0;
	*C = *C * 2.0;
	C->display();
	delete(C);

	*A = *A * *C;
	A->display();
	*/

	// Matrix<double> *I;
	// MatrixGenerator::eye(n);
	// I->display();

	/*
	Matrix S = MatrixGenerator::swapTwoRows(n,1,3);
	S.display();

	Matrix B = S*A;
	B.display();
	*/

	/*Matrix T = MatrixGenerator::addOneToAnother(n,2,1,-2);
	A = T*A;
	A.display();*/

	// A.display();
	// Matrix B = A;
	// B = A + B;
	// B.display();

	/*
	A.transpose();
	Matrix B = A;
	A.transpose();
	
	Matrix C = A*B;

	A.display();
	B.display();
	C.display();

	Matrix b(3,1,1);

	Matrix x = B*b;
	x.display();*/

	/*
	Color red = Color(1.0,0.1,0.0);
	Color green = Color(0.0,1.0,0.1);
	Color blue = Color(0.1,0.0,1.0);
	
	Color rainbow = red + green;
	std::cout << rainbow << std::endl;
	rainbow += blue;
	std::cout << rainbow << std::endl;
	rainbow *= 0.5;
	std::cout << rainbow << std::endl;
	rainbow *= red;
	std::cout << rainbow << std::endl;
	rainbow = red * green;
	std::cout << rainbow << std::endl;
	*/

	return 0;
}