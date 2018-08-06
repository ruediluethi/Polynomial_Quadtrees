#include "Matrix.h"
#include "Color.h"

class IOcolorMatrix{
public:
	static Matrix<Color> readFromFile(const char *const);
	static void writeToFile(const Matrix<Color>, const char *const);
};