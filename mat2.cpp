#include "mat2.h"

mat2::mat2(double a, double b, double c, double d)
{
	m[0]=a;
	m[1]=b;
	m[2]=c;
	m[3]=d;
}

mat2::mat2()
{
	m[0]=0;
	m[1]=0;
	m[2]=0;
	m[3]=0;
}

point mat2::operator /(point p)
{
	double i=m[0]*m[3]-m[1]*m[2];
	if(i==0) return point();
	else return point((m[3]*p.p[0]-m[1]*p.p[1])/i, (m[0]*p.p[1]-m[2]*p.p[0])/i);
}

point mat2::operator *(point p)
{
	return point(m[0]*p.p[0]+m[1]*p.p[1], m[2]*p.p[0]+m[3]*p.p[1]);
}

mat2 mat2::transpose()
{
	return mat2(m[0],m[2],m[1],m[3]);
}
