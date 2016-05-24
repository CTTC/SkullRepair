#include "linequ.h"
#include <fstream>

double* Matrix::MatrixA = 0;
double* Linequ::solu = 0;
double* Linequ::sums = 0;
void Matrix::setMatrix(double *rmatr)//设置矩阵
{
	for (int i=0;i<index*index;i++)     
	{
		*(MatrixA+i)=rmatr[i];//矩阵成员赋初值
	}
}
Matrix::Matrix(int dims)
{
	index=dims;                    //保护数据赋值
	MatrixA=new double [index*index];//动态内存分配
	cout<<"construct Matrix class "<<endl;//标记Matrix类的构造函数

}

Matrix::~Matrix()
{
	delete []MatrixA;    //内存释放
	cout<<"destruct Matrix class  "<<endl;//标记Matrix类的析构函数
}

void Matrix::printM()    //显示矩阵元素
{
	cout<<"The Matrix is:"<<endl;
	for(int i=0;i<index;i++)
	{
		for(int j=0;j<index;j++)
			cout<<*(MatrixA+i*index+j)<<" ";
		cout<<endl;
	}

}

Linequ::Linequ(int dims):Matrix(dims)   //派生类Linequ的构造函数
{
	//使用参数调用基类构造函数
	sums=new double[dims];   //动态分配内存
	solu=new double[dims];
	cout<<"construct Linequ class "<<endl;//标记Linequ类构造函数
}

Linequ::~Linequ()    //派生类Linequ的析构函数
{
	//系统默认调用基类析构函数
	delete []sums;
	delete []solu;//释放内存
	cout<<"destruct Linequ class"<<endl;  //标记Linequ类析构函数
}

void Linequ::setLinequ(double *a,double *b)//设置方程组
{
	setMatrix(a);//调用基类函数
	for (int i=0;i<index;i++)
		sums[i]=b[i];
}
void Linequ::printL()//显示方程组
{
	cout<<"The Line equation is:  "<<endl;
	for(int i=0;i<index;i++)
	{
		for(int j=0;j<index;j++)
			cout<<*(MatrixA+i*index+j)<<" ";
		cout<<"  "<<sums[i]<<endl;
	}
}

void Linequ::showX()//输出方程解
{
	cout<<"The Result is: "<<endl;
	ofstream fout("result.txt");
	for(int i=0;i<index;i++)
	{
		fout<<solu[i]<<endl;
	}
}

double Linequ::getSx(iPoint ri,vector<iPoint>& results)
{
	double ret = 0;
	int size = results.size();
	for (int i=0;i<size;i++)
	{
		iPoint p = results[i];
		ret += solu[i]*1./(sqrt(2.*PI)*SIGMA)*exp(-1./(2.*sqr(SIGMA))*(pow(ri.x-p.x,2)+pow(ri.y-p.y,2)+pow(ri.z-p.z,2)));
	}
	size = index;
	ret += solu[size-4]+solu[size-3]*ri.x+solu[size-2]*ri.y+solu[size-1]*ri.z;
	return ret;
}
iPoint Linequ::getS(iPoint ri,vector<iPoint>& results)
{
	double dx=0,dy=0,dz=0;
	int size = results.size();
	for (int i=0;i<size;++i)
	{
		iPoint p = results[i];
		dx += solu[i]*-1./(sqrt(2.*PI)*pow(SIGMA,3))*(ri.x-p.x)*exp(-1./(2.*sqr(SIGMA))*(pow(ri.x-p.x,2)+pow(ri.y-p.y,2)+pow(ri.z-p.z,2)));
		dy += solu[i]*-1./(sqrt(2.*PI)*pow(SIGMA,3))*(ri.y-p.y)*exp(-1./(2.*sqr(SIGMA))*(pow(ri.x-p.x,2)+pow(ri.y-p.y,2)+pow(ri.z-p.z,2)));
		dz += solu[i]*-1./(sqrt(2.*PI)*pow(SIGMA,3))*(ri.z-p.z)*exp(-1./(2.*sqr(SIGMA))*(pow(ri.x-p.x,2)+pow(ri.y-p.y,2)+pow(ri.z-p.z,2)));
	}
	size = index;
	dx += solu[size-3];
	dy += solu[size-2];
	dz += solu[size-1];
	iPoint retp(dx,dy,dz);
	return retp;
}
iPoint Linequ::getR(iPoint ri,vector<iPoint>& results)
{
	iPoint ret = ri,p,s;
	do 
	{
		ri = ret;
		s = getS(ri,results);
		ret = ri - (getSx(ri,results)/pow(s.module(),2))*s;
		p = ret - ri;
	} while (p.module()>=0.001);
	return ret;
}

int Linequ::Solve()       //全主元高斯消去法求解方程
{
	ifstream fin("result.txt");
	double data;
	if(fin)
	{
		for (int i=0;i<index;i++)
		{
			fin>>data;
			solu[i] = data;
		}
		return 1;
	}
	int *js,l,k,i,j,is,p,q;
	double d,t;
	js=new int[index];
	l=1;
	for(k=0;k<=index-2;k++) //消去过程
	{
		d = (double)k/(double)(index-2)*100.;
		cout<<"%"<<d<<endl;
		d=0.0;
		for(i=k;i<index-1;i++)
			for(j=k;j<index-1;j++)
			{
				t=fabs(MatrixA[i*index+j]);
				if(t>d)
				{
					d=t;
					js[k]=j;
					is=i;
				}
			}

			if(d+1.0==1.0)
				l=0;
			else
			{
				if(js[k]!=k)
					for(i=0;i<=index-1;i++)
					{
						p=i*index+k;
						q=i*index+js[k];
						t=MatrixA[p];
						MatrixA[p]=MatrixA[q];
						MatrixA[q]=t;
					}

					if(is!=k)
					{
						for(j=k;j<=index-1;j++)
						{
							p=k*index+j;
							q=is*index+j;
							t=MatrixA[p];
							MatrixA[p]=MatrixA[q];
							MatrixA[q]=t;
						}
						t=sums[k];
						sums[k]=sums[is];
						sums[is]=t;
					}
			}
			if(l==0)
			{
				delete [] js;
				cout<<"fail "<<endl;
				return(0);
			}

			d=MatrixA[k*index+k];
			for(j=k+1;j<=index-1;j++)
			{
				p=k*index+j;
				MatrixA[p]=MatrixA[p]/d;
			}
			sums[k]=sums[k]/d;
			for(i=k+1;i<=index-1;i++)
			{
				for (j=k+1;j<=index-1;j++)
				{
					p=i*index+j;
					MatrixA[p]=MatrixA[p]-MatrixA[i*index+k]*MatrixA[k*index+j];
				}
				sums[i]=sums[i]-MatrixA[i*index+k]*sums[k];
			}
	}

	d=MatrixA[(index-1)*index+index-1];
	if(fabs(d)+1.0==1.0)
	{
		delete [] js;
		cout<<"fail "<<endl;
		return (0);
	}

	solu[index-1]=sums[index-1]/d;    //回代过程
	for(i=index-2;i>=0;i--)
	{
		t=0.0;
		for(j=i+1;j<=index-1;j++)
			t=t+MatrixA[i*index+j]*solu[j];
		solu[i]=sums[i]-t;
	}
	js[index-1]=index-1;
	for(k=index-1;k>=0;k--)
	{
		if(js[k]!=k)
		{
			t=solu[k];
			solu[k]=solu[js[k]];
			solu[js[k]]=t;
		}
	}
	delete [] js;
	return(1);
}