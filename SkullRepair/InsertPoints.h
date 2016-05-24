#include "iGeometry.h"
#include <vector>
#include "vtkCommon.h"
#include <set>
#pragma once

struct InsertPoints
{
	vtkSmartPointer<vtkPolyData> meshdata;
	vtkSmartPointer<vtkPoints> meshPoint;
	vtkSmartPointer<vtkCellArray> meshLine;

	struct ComAngle
	{
		bool operator () (const std::pair<vtkIdType, double> &p1, const std::pair<vtkIdType, double> &p2) const
		{
			return p1.second < p2.second;
		}
	};

	vtkIdType nPt;
	std::vector<vtkIdType> prev;
	std::vector<vtkIdType> succ;
	std::vector<double> angle;
	std::vector<bool> inSet;
	std::set<std::pair<vtkIdType, double>, ComAngle> angleSet;

	InsertPoints() {}

	vtkSmartPointer<vtkPolyData> Insert(std::vector<iPoint>& points, iPoint viewDir);
	double PointAngle(vtkIdType ptId, iPoint viewDir);
	void InsertAngle(vtkIdType ptId, iPoint viewDir);
	void AddLine(vtkIdType p1Id, vtkIdType p2Id);
};
