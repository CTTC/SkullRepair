#pragma once

#include "vtkCommon.h"
#include "iGeometry.h"
#include <map>
#include <vector>

struct EdgePicking
{
	EdgePicking() {}
	static void GetResult(vtkSmartPointer<vtkPolyData> polydata, std::vector<iPoint> cellNormal, std::vector< std::pair<vtkIdType, vtkIdType> > linedata,
		std::vector< std::pair<vtkIdType, vtkIdType> > nextPoly, std::map< std::pair<vtkIdType, vtkIdType>, vtkIdType> &lineMap, iPoint viewDir,
		std::vector<vtkIdType> &uBound, std::vector<vtkIdType> &lBound);

	static void FindEdge(vtkSmartPointer<vtkPolyData> polydata, std::vector<iPoint> cellNormal, std::vector< std::pair<vtkIdType, vtkIdType> > linedata, std::vector< std::pair<vtkIdType, vtkIdType> > nextPoly,
		iPoint viewDir, std::vector< std::pair<vtkIdType, vtkIdType> > &edgeLineId);

	static void GetTargetPoint(vtkSmartPointer<vtkPolyData> polydata, std::vector< std::pair<vtkIdType, vtkIdType> > edgeLineId, vtkIdType &targetPt);

	static void FindHole(vtkSmartPointer<vtkPolyData> polydata, std::vector< std::pair<vtkIdType, vtkIdType> > edgeLineId, vtkIdType targetPt,
		std::vector<vtkIdType> &holePtIndex, std::vector<vtkIdType> &holeSucc);

	static void GetConnect(std::vector< std::pair<vtkIdType, vtkIdType> > edgeLineId, std::vector<std::vector<vtkIdType>> &connect);
	
	static vtkIdType dfs(vtkIdType nowId, std::vector<std::vector<vtkIdType>> &connect, std::vector<int> &depth, std::vector<vtkIdType> &succ);

	static void GetUpperBound(vtkSmartPointer<vtkPolyData> polydata, std::vector< std::pair<vtkIdType, vtkIdType> > &linedata, std::map< std::pair<vtkIdType, vtkIdType>, vtkIdType> &lineMap,
		std::vector< std::pair<vtkIdType, vtkIdType> > &nextPoly, iPoint viewDir, std::vector<vtkIdType> &holePtIndex, std::vector<vtkIdType> &uBound);


	static vtkIdType ThirdPt(vtkSmartPointer<vtkPolyData> &polydata, vtkIdType cellId, vtkIdType p1, vtkIdType p2);

	static vtkIdType GetLineId(std::map< std::pair<vtkIdType, vtkIdType>, vtkIdType> &lineMap, vtkIdType p1, vtkIdType p2);
};