#pragma once
#include <cmath>
#include <iostream>

class Point {
private:
	double xval, yval;
public:  // Constructor uses default arguments to allow calling with zero, one,
		 // or two values.
	Point(double x = 0.0, double y = 0.0) {
		xval = x;
		yval = y;
	}

	// Extractors.
	double x() { return xval; }
	double y() { return yval; }

	// Distance to another point.  Pythagorean thm.
	double dist(Point other) {
		double xd = xval - other.xval;
		double yd = yval - other.yval;
		return sqrt(xd*xd + yd*yd);
	}

	// Add or subtract two points.
	Point add(Point b)
	{
		return Point(xval + b.xval, yval + b.yval);
	}
	Point sub(Point b)
	{
		return Point(xval - b.xval, yval - b.yval);
	}
	Point mult(double scalar)
	{
		return Point(xval * scalar, yval * scalar);
	}
	// dot product
	double dot(Point b)
	{
		return xval * b.xval + yval * b.yval;
	}
	// apply Periodic Boundary Conditions - To print configurations in a sensible way, print the configurations with .pbc(). To calculate pair separations, use .sub().pbc()
	Point pbc(const double &boxl, const double &boxinv) {
		return Point(xval - static_cast<int>(xval * boxinv + ((xval >= 0.0) ? 0.5 : -0.5)) * boxl, yval - static_cast<int>(yval * boxinv + ((yval >= 0.0) ? 0.5 : -0.5)) * boxl);
	}
	// Print the point on the stream.  The class ostream is a base class
	// for output streams of various types.
	void print(std::ostream &strm)
	{
		strm << xval << ' ' << yval << '\n';
	}
};
