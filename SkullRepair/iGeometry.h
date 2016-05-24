#include <math.h>
#include <iostream>
#define sqr(x) ((x)*(x))
#include "vtkCommon.h"
#pragma once

#define PI 3.1415926535
#define SIGMA 1.0

struct iPoint
{
	double x;
	double y;
	double z;
	iPoint* normal;

	iPoint() :x(0), y(0), z(0) {}
	iPoint(double xx, double yy, double zz)
		:x(xx), y(yy), z(zz) {}
	iPoint(const iPoint& p)
		:x(p.x), y(p.y), z(p.z) {}
	iPoint(double *m)
		:x(m[0]), y(m[1]), z(m[2]) {}
	iPoint(vtkSmartPointer<vtkPolyData> &polydata, vtkIdType p)
	{
		double m[3];
		polydata->GetPoint(p, m);
		x = m[0];
		y = m[1];
		z = m[2];
	}

	iPoint(vtkSmartPointer<vtkPolyData> &polydata, vtkIdType p1, vtkIdType p2)
	{
		double m1[3], m2[3];
		polydata->GetPoint(p1, m1);
		polydata->GetPoint(p2, m2);
		x = m2[0] - m1[0];
		y = m2[1] - m1[1];
		z = m2[2] - m1[2];
	}


	friend std::ostream& operator<<(std::ostream& cout, const iPoint& p);//使用友元函数重载<<输出运算符
	const void operator=(const iPoint& p){ x = p.x; y = p.y; z = p.z; }//重载=号
	const double distance_square(const iPoint &p) { return sqr(x - p.x) + sqr(y - p.y) + sqr(z - p.z); }

	iPoint* getnormal() { return normal; }
	double module() { return sqrt(x*x + y*y + z*z); }
};

iPoint const operator+(const iPoint &p1, const iPoint &p2);
iPoint const operator-(const iPoint &p1, const iPoint &p2);
iPoint const operator*(const double d, const iPoint &p);
iPoint const operator*(const iPoint &p, const double d);
iPoint const operator/(const iPoint &p, const double d);
double const operator*(const iPoint &p1, const iPoint &p2);
iPoint const operator^(const iPoint &p1, const iPoint &p2);
double GetAngle(iPoint p1, iPoint p2);
double GetTheta(iPoint p1, iPoint p2, iPoint n);
