#include "EdgePicking.h"

void EdgePicking::GetResult(vtkSmartPointer<vtkPolyData> polydata, std::vector<iPoint> cellNormal, std::vector< std::pair<vtkIdType, vtkIdType> > linedata,
	std::vector< std::pair<vtkIdType, vtkIdType> > nextPoly, std::map< std::pair<vtkIdType, vtkIdType>, vtkIdType> &lineMap, iPoint viewDir,
	std::vector<vtkIdType> &uBound, std::vector<vtkIdType> &lBound)
{
	std::vector< std::pair<vtkIdType, vtkIdType> > edgeLineId;
	vtkIdType targetPt;

	std::vector<vtkIdType> holePtIndex;
	std::vector<vtkIdType> holeSucc;

	FindEdge(polydata, cellNormal, linedata, nextPoly, viewDir, edgeLineId);
	GetTargetPoint(polydata, edgeLineId, targetPt);
	FindHole(polydata, edgeLineId, targetPt, holePtIndex, holeSucc);

	GetUpperBound(polydata, linedata, lineMap, nextPoly, viewDir, holePtIndex, uBound);

	//std::cout << endl << targetPt << endl;
}

void EdgePicking::FindEdge(vtkSmartPointer<vtkPolyData> polydata, std::vector<iPoint> cellNormal, std::vector< std::pair<vtkIdType, vtkIdType> > linedata, std::vector< std::pair<vtkIdType, vtkIdType> > nextPoly,
	iPoint viewDir, std::vector< std::pair<vtkIdType, vtkIdType> > &edgeLineId)
{
	std::cout << "FindEdge started...";
	edgeLineId.clear();

	std::cout << viewDir << endl;

	std::vector<vtkIdType> ptUse;
	for (vtkIdType ptId = 0; ptId < polydata->GetNumberOfPoints(); ptId++) ptUse.push_back(-1);

	for (vtkIdType lineId = 0; lineId < linedata.size(); lineId++)
	{
		vtkIdType poly1 = nextPoly[lineId].first, poly2 = nextPoly[lineId].second;
		double dot1 = 0.0, dot2 = 0.0;
		dot1 = cellNormal[poly1] * viewDir;
		dot2 = cellNormal[poly2] * viewDir;

		if (dot1 * dot2 < -1e-7)
		{
			vtkIdType p1 = linedata[lineId].first, p2 = linedata[lineId].second;
			edgeLineId.push_back(std::make_pair(p1, p2));
		}
	}
	std::cout << "FindEdge ended" << endl;
}

void EdgePicking::GetTargetPoint(vtkSmartPointer<vtkPolyData> polydata, std::vector< std::pair<vtkIdType, vtkIdType> > edgeLineId, vtkIdType &targetPt)
{
	std::cout << "GetTargetPoint started...";
	targetPt = -1;
	for (std::vector<std::pair<vtkIdType, vtkIdType>>::iterator iter = edgeLineId.begin(); iter != edgeLineId.end() && targetPt == -1; iter++)
	{
		double m[3];
		polydata->GetPoint(iter->first, m);
		if (m[0] > 120 && m[1] < 100 && m[2] > 80 && targetPt == -1)
		{
			targetPt = iter->first;
		}
		polydata->GetPoint(iter->second, m);
		if (m[0] > 120 && m[1] < 100 && m[2] > 80 && targetPt == -1)
		{
			targetPt = iter->first;
		}
	}

	double m[3];
	polydata->GetPoint(targetPt, m);
	//std::cout << m[0] << ' ' << m[1] << ' ' << m[2] << endl;

	std::cout << "GetTargetPoint ended" << endl;
}

void EdgePicking::FindHole(vtkSmartPointer<vtkPolyData> polydata, std::vector< std::pair<vtkIdType, vtkIdType> > edgeLineId, vtkIdType targetPt,
	std::vector<vtkIdType> &holePtIndex, std::vector<vtkIdType> &holeSucc)
{
	std::cout << "FindHole started...";
	std::vector<std::vector<vtkIdType>> connect;
	connect.clear();
	connect.resize(polydata->GetNumberOfPoints());
	for (std::vector<std::vector<vtkIdType>>::iterator iter = connect.begin(); iter != connect.end(); iter++) iter->clear();
	GetConnect(edgeLineId, connect);

	std::vector<int> depth;
	holeSucc.clear();
	depth.clear();
	for (vtkIdType ptId = 0; ptId < polydata->GetNumberOfPoints(); ptId++)
	{
		depth.push_back(-1);
		holeSucc.push_back(-1);
	}
	depth[targetPt] = 0;
	vtkIdType depPt1 = dfs(targetPt, connect, depth, holeSucc);
	for (vtkIdType ptId = 0; ptId < polydata->GetNumberOfPoints(); ptId++)
	{
		depth[ptId] = -1;
		holeSucc[ptId] = -1;
	}
	depth[depPt1] = 0;
	vtkIdType depPt2 = dfs(depPt1, connect, depth, holeSucc);
	//succ[depPt2] = depPt1;

	holePtIndex.clear();
	for (vtkIdType ptId = depPt1; ptId != depPt2; ptId = holeSucc[ptId])
	{
		holePtIndex.push_back(ptId);
	}
	holePtIndex.push_back(depPt2);
	std::cout << "FindHole ended" << endl;
}

void EdgePicking::GetUpperBound(vtkSmartPointer<vtkPolyData> polydata, std::vector< std::pair<vtkIdType, vtkIdType> > &linedata, std::map< std::pair<vtkIdType, vtkIdType>, vtkIdType> &lineMap,
	std::vector< std::pair<vtkIdType, vtkIdType> > &nextPoly, iPoint viewDir, std::vector<vtkIdType> &holePtIndex, std::vector<vtkIdType> &uBound)
{
	std::cout << "GetUpperBound started...";
	std::vector<vtkIdType> prev, succ;
	prev.clear(); succ.clear();
	for (vtkIdType ptId = 0; ptId < polydata->GetNumberOfPoints(); ptId++)
	{
		succ.push_back(-1); prev.push_back(-1);
	}
	vtkIdType head = holePtIndex[0];
	vtkIdType uBdLen = holePtIndex.size();

	for (int i = 0; i < uBdLen; i++)
	{
		succ[holePtIndex[i]] = holePtIndex[(i + 1) % uBdLen];
		prev[holePtIndex[i]] = holePtIndex[(i + uBdLen - 1) % uBdLen];
	}

	std::vector<std::vector<vtkIdType>> connect;
	connect.clear();
	connect.resize(polydata->GetNumberOfPoints());
	GetConnect(linedata, connect);

	vtkIdType ptId = head;
	for (int nTime = 1000 * uBdLen; nTime > 0; nTime--)
	{
		vtkIdType p1 = prev[ptId], p2 = ptId, p3 = succ[ptId];
		vtkIdType l1 = GetLineId(lineMap, p1, p2), l2 = GetLineId(lineMap, p2, p3);

		int flag;

		//Type 1
		flag = 0;
		if (GetLineId(lineMap, p1, p3) != -1)
		{
			if (iPoint(viewDir)*iPoint(polydata, p2, p1) > 0 || iPoint(viewDir)*iPoint(polydata, p2, p3) > 0)
			{
				flag = 1;
				prev[p2] = succ[p2] = -1;
				succ[p1] = p3; prev[p3] = p1;
				uBdLen--;
				ptId = p1;
				continue;
			}
		}

		vtkIdType cellId, pNew;
		//Type 2
		flag = 0;
		vtkIdType pNew1 = ThirdPt(polydata, nextPoly[l1].first, p1, p2), pNew2 = ThirdPt(polydata, nextPoly[l1].second, p1, p2);
		pNew = (GetAngle(iPoint(viewDir), iPoint(polydata, p2, pNew1)) < GetAngle(iPoint(viewDir), iPoint(polydata, p2, pNew2))) ? pNew1 : pNew2;
		if (GetLineId(lineMap, pNew, p3) != -1 && succ[pNew] == -1)
		{
			if (iPoint(viewDir)*iPoint(polydata, p2, pNew) > 0)
			{
				flag = 1;
				prev[p2] = succ[p2] = -1;
				prev[pNew] = p1; succ[p1] = pNew;
				succ[pNew] = p3; prev[p3] = pNew;
				ptId = p1;
				continue;
			}
		}

		//Type 3
		if ((iPoint(viewDir) * iPoint(polydata, p2, pNew) > 0 && iPoint(viewDir)*iPoint(polydata, p1, pNew) > 0) && succ[pNew] == -1)
		{
			flag = 2;
			prev[pNew] = p1; succ[p1] = pNew;
			succ[pNew] = p2; prev[p2] = pNew;
			uBdLen++;
			ptId = pNew;
			continue;
		}

		ptId = succ[ptId];
	}
	//std::cout << uBdLen << endl;

	head = -1;
	for (vtkIdType ptId = 0; ptId < polydata->GetNumberOfPoints(); ptId++)
	{
		if (succ[ptId] != -1) head = ptId;
	}
	uBound.clear();
	for (vtkIdType ptId = 0; ptId < uBdLen; ptId++)
	{
		uBound.push_back(head);
		head = succ[head];
	}
	std::cout << "GetUpperBound ended" << endl;
}

void EdgePicking::GetConnect(std::vector< std::pair<vtkIdType, vtkIdType> > edgeLineId, std::vector<std::vector<vtkIdType>> &connect)
{
	for (std::vector< std::pair<vtkIdType, vtkIdType> >::iterator iter = edgeLineId.begin(); iter != edgeLineId.end(); iter++)
	{
		vtkIdType p1 = iter->first, p2 = iter->second;
		connect[p1].push_back(p2);
		connect[p2].push_back(p1);
	}
}

vtkIdType EdgePicking::dfs(vtkIdType nowId, std::vector<std::vector<vtkIdType>> &connect, std::vector<int> &depth, std::vector<vtkIdType> &succ)
{
	vtkIdType maxId = nowId;
	for (std::vector<vtkIdType>::iterator iter = connect[nowId].begin(); iter != connect[nowId].end(); iter++)
	{
		if (depth[*iter] == -1)
		{
			depth[*iter] = depth[nowId] + 1;
			vtkIdType maxChild = dfs(*iter, connect, depth, succ);
			if (depth[maxChild] > depth[maxId])
			{
				maxId = maxChild;
				succ[nowId] = *iter;
			}
		}
	}
	return maxId;
}

vtkIdType EdgePicking::GetLineId(std::map< std::pair<vtkIdType, vtkIdType>, vtkIdType> &lineMap, vtkIdType p1, vtkIdType p2)
{
	if (p1 < p2)
	{
		std::map <std::pair<vtkIdType, vtkIdType>, vtkIdType>::iterator iter = lineMap.find(std::make_pair(p1, p2));
		if (iter == lineMap.end()) return -1;
		else return lineMap[std::make_pair(p1, p2)];
	}
	else
	{
		std::map <std::pair<vtkIdType, vtkIdType>, vtkIdType>::iterator iter = lineMap.find(std::make_pair(p2, p1));
		if (iter == lineMap.end()) return -1;
		else return lineMap[std::make_pair(p2, p1)];

	}
}

vtkIdType EdgePicking::ThirdPt(vtkSmartPointer<vtkPolyData> &polydata, vtkIdType cellId, vtkIdType p1, vtkIdType p2)
{
	vtkIdType ptsId[3];
	vtkSmartPointer<vtkIdList> ptIdList = vtkSmartPointer<vtkIdList>::New();
	ptIdList = polydata->GetCell(cellId)->GetPointIds();
	for (vtkIdType i = 0; i < ptIdList->GetNumberOfIds(); i++)
	{
		ptsId[i] = ptIdList->GetId(i);
	}
	vtkIdType res = ptsId[0];
	if (res == p1 || res == p2) res = ptsId[1];
	if (res == p1 || res == p2) res = ptsId[2];
	return res;
}