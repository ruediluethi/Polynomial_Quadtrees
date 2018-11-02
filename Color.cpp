#include <iostream>
#include <cassert>
#include <cmath>

#include "Color.h"

Color::Color(): r(0), g(0), b(0) {}
Color::Color(double red, double green, double blue): r(red), g(green), b(blue) {}

// Copy constructor
Color::Color(const Color& c): r(c.getRGB()[0]), g(c.getRGB()[1]), b(c.getRGB()[2]) {}

// Assignment operator
Color& Color::operator=(const Color& c)
{
	r = c.getRGB()[0];
	g = c.getRGB()[1];
	b = c.getRGB()[2];
	return *this;
}

//destructor
Color::~Color(){}

double* Color::getRGB() const
{
	double* rgb = new double[3];
	rgb[0] = r;
	rgb[1] = g;
	rgb[2] = b;
	return rgb;
}

void Color::setRGB(double red, double green, double blue)
{
	r = red;
	g = green;
	b = blue;
}

double Color::distance(const Color& other)
{
	double *rgbThis = getRGB();
	double *rgbOther = other.getRGB();

	double dr = rgbThis[0] - rgbOther[0];
	double dg = rgbThis[1] - rgbOther[1];
	double db = rgbThis[2] - rgbOther[2];

	return sqrt(2*dr*dr + 4*dg*dg + 3*db*db); // source: https://en.wikipedia.org/wiki/Color_difference
	// return 2*dr*dr + 4*dg*dg + 3*db*db;
}

std::ostream& operator<< (std::ostream& out, Color one)
{
	double *rgbOne = one.getRGB();
	double r = round(rgbOne[0]*100)/100.0;
	double g = round(rgbOne[1]*100)/100.0;
	double b = round(rgbOne[2]*100)/100.0;
    return out << r << "|" << g << "|" << b;
}

// addition
Color& operator+= (Color& one, Color other)
{
	double *rgbOne = one.getRGB();
	double *rgbOther = other.getRGB();

	double r = rgbOne[0] + rgbOther[0];
	double g = rgbOne[1] + rgbOther[1];
	double b = rgbOne[2] + rgbOther[2];

	one.setRGB(r,g,b);
	return one;
}

Color operator+ (Color one, Color other)
{
	Color sum = one;
	sum += other;
	return sum;
}

// difference
Color& operator-= (Color& one, Color other)
{
	double *rgbOne = one.getRGB();
	double *rgbOther = other.getRGB();

	double r = fabs(rgbOne[0] - rgbOther[0]);
	double g = fabs(rgbOne[1] - rgbOther[1]);
	double b = fabs(rgbOne[2] - rgbOther[2]);

	one.setRGB(r,g,b);
	return one;
}

Color operator- (Color one, Color other)
{
	Color diff = one;
	diff -= other;
	return diff;
}


// scalar multipilication
Color& operator*= (Color& one, double scalar)
{
	double *rgbOne = one.getRGB();

	double r = rgbOne[0] * scalar;
	double g = rgbOne[1] * scalar;
	double b = rgbOne[2] * scalar;

	one.setRGB(r,g,b);
	return one;
}

Color operator* (Color one, double scalar)
{
	Color prod = one;
	prod *= scalar;
	return prod;
}

// Color multiplication
Color& operator*= (Color& one, Color other)
{
	double *rgbOne = one.getRGB();
	double *rgbOther = other.getRGB();

	double r = rgbOne[0] * rgbOther[0];
	double g = rgbOne[1] * rgbOther[1];
	double b = rgbOne[2] * rgbOther[2];

	one.setRGB(r,g,b);
	return one;
}

Color operator* (Color one, Color other)
{
	Color prod = one;
	prod *= other;
	return prod;
}