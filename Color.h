// avoid double includes
#ifndef _COLOR_H_
#define _COLOR_H_

class Color{
public:
	Color();
	Color(double, double, double);

	// rule of three
	Color(const Color&); // copy
	Color& operator=(const Color&); // assignment
	~Color(); // destructor

	double* getRGB() const;
	void setRGB(double, double, double);
	double distance(const Color&);

private:
	double r;
	double g;
	double b;
};

// output operator
std::ostream& operator<< (std::ostream&, Color);

// addition
Color& operator+= (Color&, Color);
Color operator+ (Color, Color);

// difference
Color& operator-= (Color&, Color);
Color operator- (Color, Color);
// distance
//double operator- (Color, Color);

// multiplication with scalar
Color& operator*= (Color&, double);
Color operator* (Color, double);

// color multiplication
Color& operator*= (Color&, Color);
Color operator* (Color, Color);

#endif