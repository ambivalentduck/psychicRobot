#ifndef MAT2_H
#define MAT2_H

#include "point.h"

class mat2
{
public:
	mat2();
	mat2(double a, double b, double c, double d);
	point operator /(point p);
	point operator *(point p);
	mat2 operator *(double d);
	mat2 operator +(mat2 m);
	mat2 transpose();
	double m[4];
};

#endif

