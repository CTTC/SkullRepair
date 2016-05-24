#include "iGeometry.h"

iPoint const operator+(const iPoint &p1, const iPoint &p2)
{
	return iPoint(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
}

iPoint const operator-(const iPoint &p1, const iPoint &p2)
{
	return iPoint(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
}

iPoint const operator*(const double d, const iPoint &p)
{
	return iPoint(d * p.x, d * p.y, d * p.z);
}

iPoint const operator*(const iPoint &p, const double d)
{
	return iPoint(d * p.x, d * p.y, d * p.z);
}

iPoint const operator/(const iPoint &p, const double d)
{
	return iPoint(p.x / d, p.y / d, p.z / d);
}

double const operator*(const iPoint &p1, const iPoint &p2)
{
	return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z;
}

iPoint const operator^(const iPoint &p1, const iPoint &p2)
{
	return iPoint(p1.y*p2.z - p1.z*p2.y, p1.z*p2.x - p1.x*p2.z, p1.x*p2.y - p1.y*p2.x);
}

double GetAngle(iPoint p1, iPoint p2)
{
	return acos((p1*p2)/(p1.module()*p2.module()));
}

std::ostream& operator<<(std::ostream& cout, const iPoint& p)
{
	cout << "x:" << p.x << " y:" << p.y << " z:" << p.z;
	return cout;
}

double GetTheta(iPoint p1, iPoint p2, iPoint n)
{
	iPoint nn = p1^p2;
	double judge = nn*n;
	if (judge<0)
		return acos(p1*p2 / (p1.module()*p2.module()));
	else
		return 2 * PI - acos(p1*p2 / (p1.module()*p2.module()));
}