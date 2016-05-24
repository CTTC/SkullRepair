#include<iostream>
#include<cmath>
#include "iGeometry.h"
#pragma once
using namespace std;
class Matrix
{
public:                           //外部接口
	Matrix(int dims=2);    //构造函数
	~Matrix();      //析构函数
	void setMatrix(double *rmatr);
	void printM();

protected:     //保护数据成员
	int index;    //矩阵维数
	static double *MatrixA;  //矩阵存放数组首地址

};
class Linequ:public Matrix    //公有派生类Linequ定义
{
public:       //外部接口
	Linequ(int dims=2);   //构造函数
	~Linequ();     //析构函数
	void setLinequ(double *a,double *b);   //方程赋值
	void printL();    //显示方程
	int Solve();    //全主元高斯消去法求解方程
	void showX();    //显示方程的解
	iPoint getR(iPoint ri,vector<iPoint>& results);
	iPoint getS(iPoint ri,vector<iPoint>& results);
	double getSx(iPoint ri,vector<iPoint>& results);
private:
	static double *sums;    //方程右端项
	static double *solu;    //方程的解
};
