#include <iostream>
#include "Matrix.h"
#include "Matrix.cpp" // to avoid linking errors
#include "Color.h"
#include "CImg.h" // CImg Library

#include "IOcolorMatrix.h"

Matrix<Color> IOcolorMatrix::readFromFile(const char *const filename)
{
	cimg_library::CImg<unsigned char> img(filename);

	int imgWidth = img.width();
	int imgHeight = img.height();

	Matrix<Color> C(imgHeight,imgWidth,Color());

	for (int x = 0; x < imgWidth; x++){
    for (int y = 0; y < imgHeight; y++){
    	unsigned char *r = img.data(x, y, 0, 0);
    	unsigned char *g = img.data(x, y, 1, 0);
    	unsigned char *b = img.data(x, y, 2, 0);
    	Color pxColor = Color(*r,*g,*b);
    	pxColor *= 1.0/256.0;
    	C.setElement(y,x,pxColor);
    }
	}

	return C;
}

void IOcolorMatrix::writeToFile(Matrix<Color> C, const char *const filename)
{
  int imgHeight = C.getDimension()[0];
	int imgWidth = C.getDimension()[1];

	cimg_library::CImg<unsigned char> img(imgWidth, imgHeight, 1, 3);

	for (int x = 0; x < imgWidth; x++){
		for (int y = 0; y < imgHeight; y++){
			Color colorPx = C.getElement(y,x);
			double *rgbPx = colorPx.getRGB();
			double r = rgbPx[0]*255.0;
			double g = rgbPx[1]*255.0;
			double b = rgbPx[2]*255.0;
			img(x,y,0,0) = r;
			img(x,y,1,0) = g;
			img(x,y,2,0) = b;
		}
	}

	img.save(filename);
}