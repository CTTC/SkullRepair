#include<iostream>
#include<cmath>
#include "iGeometry.h"
#pragma once
using namespace std;
class Matrix
{
public:                           //�ⲿ�ӿ�
	Matrix(int dims=2);    //���캯��
	~Matrix();      //��������
	void setMatrix(double *rmatr);
	void printM();

protected:     //�������ݳ�Ա
	int index;    //����ά��
	static double *MatrixA;  //�����������׵�ַ

};
class Linequ:public Matrix    //����������Linequ����
{
public:       //�ⲿ�ӿ�
	Linequ(int dims=2);   //���캯��
	~Linequ();     //��������
	void setLinequ(double *a,double *b);   //���̸�ֵ
	void printL();    //��ʾ����
	int Solve();    //ȫ��Ԫ��˹��ȥ����ⷽ��
	void showX();    //��ʾ���̵Ľ�
	iPoint getR(iPoint ri,vector<iPoint>& results);
	iPoint getS(iPoint ri,vector<iPoint>& results);
	double getSx(iPoint ri,vector<iPoint>& results);
private:
	static double *sums;    //�����Ҷ���
	static double *solu;    //���̵Ľ�
};
