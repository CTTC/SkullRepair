#pragma once

#include "vtkCommon.h"
#include "iGeometry.h"
#include <map>
#include <vector>
#include <ctime>

const iPoint VIEWPOINT0(0, 200, 100);
const iPoint VIEWPOINT1(150, 100, 100);
const int DEBUG = 1;

struct DisplayPack
{
	vtkSmartPointer<vtkPolyData> data;
	double r, g, b;
	float width, size;

	DisplayPack() {}
	DisplayPack(vtkSmartPointer<vtkPolyData> _data, double _r, double _g, double _b, float _width, float _size)
		: data(_data), r(_r), g(_g), b(_b), width(_width), size(_size) {}
	const void operator=(const DisplayPack& p)
	{
		data = p.data;
		r = p.r;
		g = p.g;
		b = p.b;
		width = p.width;
		size = p.size;
	}
};

vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
std::vector<iPoint> cellNormal;
std::vector< std::pair<vtkIdType, vtkIdType> > linedata;
std::vector< std::pair<vtkIdType, vtkIdType> > nextPoly;
std::map< std::pair<vtkIdType, vtkIdType>, vtkIdType> lineMap;
std::vector<vtkIdType> uBound;
std::vector<vtkIdType> lBound;
std::vector<iPoint> uMesh;
std::vector <iPoint> lMesh;
iPoint viewDir;

//===============================================================================================================================

void InputData(vtkSmartPointer<vtkPolyData> &polydata);

void GetCellNormal(vtkSmartPointer<vtkPolyData> polydata, std::vector<iPoint> &cellNormal);

void LableLine(vtkSmartPointer<vtkPolyData> polydata, std::vector< std::pair<vtkIdType, vtkIdType> > &linedata, std::vector< std::pair<vtkIdType, vtkIdType> > &nextPoly,
	std::map< std::pair<vtkIdType, vtkIdType>, vtkIdType> &lineMap);

void GetViewVect(iPoint &viewDir);

void Display(int num, ...);

void PushLine(vtkIdType p1, vtkIdType p2, vtkIdType cellId, std::map< std::pair<vtkIdType, vtkIdType>, vtkIdType> &lineMap,
	std::vector< std::pair<vtkIdType, vtkIdType> > &linedata, std::vector< std::pair<vtkIdType, vtkIdType> > &nextPoly);

void GetViewData(vtkSmartPointer<vtkPolyData> &viewdata);

void GetBeltData(vtkSmartPointer<vtkPolyData> polydata, std::vector<vtkIdType> belt, vtkSmartPointer<vtkPolyData> &beltdata);

void GetMeshData(std::vector<iPoint> mesh, vtkSmartPointer<vtkPolyData> &meshdata);