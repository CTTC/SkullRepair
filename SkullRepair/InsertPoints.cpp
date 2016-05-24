#include "insertPoints.h"

vtkSmartPointer<vtkPolyData> InsertPoints::Insert(std::vector<iPoint>& points, iPoint viewDir)
{
	std::cout << "Insert started...";
	meshdata = vtkSmartPointer<vtkPolyData>::New();
	meshPoint = vtkSmartPointer<vtkPoints>::New();
	meshLine = vtkSmartPointer<vtkCellArray>::New();

	nPt = points.size();
	//std::cout << nPt << endl;

	prev.clear(); succ.clear();

	for (vtkIdType ptId = 0; ptId < nPt; ptId++)
	{
		meshPoint->InsertPoint(ptId, points[ptId].x, points[ptId].y, points[ptId].z);
		prev.push_back( (ptId == 0) ? nPt - 1 : ptId - 1 );
		succ.push_back( (ptId == nPt - 1) ? 0 : ptId + 1 );
	}
	for (vtkIdType ptId = 0; ptId < nPt; ptId++)
	{
		vtkSmartPointer<vtkLine> singleLine = vtkSmartPointer<vtkLine>::New();
		singleLine->GetPointIds()->SetId(0, ptId);
		singleLine->GetPointIds()->SetId(1, (ptId + 1) % nPt);
		meshLine->InsertNextCell(singleLine);
	}

	angle.clear();
	angleSet.clear();
	inSet.clear();
	for (vtkIdType ptId = 0; ptId < nPt; ptId++)
	{
		angle.push_back(2 * PI);
		inSet.push_back(false);
	}
	for (vtkIdType ptId = 0; ptId < nPt; ptId++)
	{
		InsertAngle(ptId, viewDir);
	}

	ofstream fout("E:/debug.txt");

	for (int COUNT = 0; COUNT < 55; COUNT++)
	{
		std::set<std::pair<vtkIdType, double>>::iterator now = angleSet.begin();

		fout << now->second / PI * 180.0;

		//std::cout << now->second / PI * 180.0 << ' ';

		if (now->second < PI / 180.0 * 85.0)
		{
			vtkIdType p0Id = now->first, p2Id = succ[now->first], p1Id = prev[now->first];

			AddLine(p1Id, p2Id);

			prev[now->first] = succ[now->first] = -1;
			prev[p2Id] = p1Id; succ[p1Id] = p2Id;

			angleSet.erase(now); angle[p0Id] = 2 * PI; inSet[p0Id] = false;
			InsertAngle(p1Id, viewDir);
			InsertAngle(p2Id, viewDir);

			fout << endl;
		}
		else
		{
			int ith = rand() % angleSet.size();
			//int ith = 0;
			now = angleSet.begin();
			while (ith > 0)
			{
				now++; ith--;
			}

			fout << " --> " << now->second / PI * 180.0 << endl;

			if (now->second < PI / 180.0 * 135.0)
			{
				vtkIdType p0Id = now->first, p2Id = succ[now->first], p1Id = prev[now->first];
				iPoint p0 = iPoint(meshPoint->GetPoint(p0Id)), p2 = iPoint(meshPoint->GetPoint(p2Id)), p1 = iPoint(meshPoint->GetPoint(p1Id));
				double theta = now->second / 2.0;

				iPoint v1 = p0 - p1, v2 = p2 - p0;
				iPoint pn = v1 ^ v2;
				iPoint px = v2;
				iPoint py = pn ^ px;
				px = px / px.module(); py = py / py.module();
				iPoint pnew = px * cos(theta) + py * sin(theta);

				double len = (v1.module() + v2.module()) / 2.0;
				pnew = pnew / pnew.module() * len;
				pnew = p0 + pnew;

				vtkIdType pnewId = nPt;
				nPt++;

				meshPoint->InsertPoint(pnewId, pnew.x, pnew.y, pnew.z);
				prev.push_back(-1); succ.push_back(-1);
				angle.push_back(2 * PI); inSet.push_back(false);

				AddLine(p1Id, pnewId);
				AddLine(p0Id, pnewId);
				AddLine(pnewId, p2Id);

				prev[p0Id] = succ[p0Id] = -1;
				prev[pnewId] = p1Id; succ[p1Id] = pnewId;
				succ[pnewId] = p2Id; prev[p2Id] = pnewId;

				angleSet.erase(now); angle[p0Id] = 2 * PI; inSet[p0Id] = false;
				InsertAngle(p1Id, viewDir);
				InsertAngle(pnewId, viewDir);
				InsertAngle(p2Id, viewDir);
			}
			else
			{
				vtkIdType p0Id = now->first, p2Id = succ[now->first], p1Id = prev[now->first];
				iPoint p0 = iPoint(meshPoint->GetPoint(p0Id)), p2 = iPoint(meshPoint->GetPoint(p2Id)), p1 = iPoint(meshPoint->GetPoint(p1Id));
				double theta = now->second / 3.0;

				iPoint v1 = p0 - p1, v2 = p2 - p0;
				iPoint pn = v1 ^ v2;
				iPoint px, py;

				px = v2;
				py = pn ^ px;
				px = px / px.module(); py = py / py.module();
				iPoint pnew1 = px * cos(theta) + py * sin(theta);

				px = pnew1;
				py = pn ^ px;
				px = px / px.module(); py = py / py.module();
				iPoint pnew2 = px * cos(theta) + py * sin(theta);

				double len = (v1.module() + v2.module()) / 2.0;
				pnew1 = pnew1 / pnew1.module() * len;
				pnew1 = p0 + pnew1;
				pnew2 = pnew2 / pnew2.module() * len;
				pnew2 = p0 + pnew2;

				vtkIdType pnew1Id = nPt;
				nPt++;
				vtkIdType pnew2Id = nPt;
				nPt++;

				meshPoint->InsertPoint(pnew1Id, pnew1.x, pnew1.y, pnew1.z);
				prev.push_back(-1); succ.push_back(-1);
				angle.push_back(2 * PI); inSet.push_back(false);
				meshPoint->InsertPoint(pnew2Id, pnew2.x, pnew2.y, pnew2.z);
				prev.push_back(-1); succ.push_back(-1);
				angle.push_back(2 * PI); inSet.push_back(false);

				AddLine(p1Id, pnew2Id);
				AddLine(pnew2Id, pnew1Id);
				AddLine(pnew1Id, p2Id);
				AddLine(p0Id, pnew1Id);
				AddLine(p0Id, pnew2Id);

				prev[p0Id] = succ[p0Id] = -1;
				prev[pnew2Id] = p1Id; succ[p1Id] = pnew2Id;
				succ[pnew2Id] = pnew1Id; prev[pnew1Id] = pnew2Id;
				succ[pnew1Id] = p2Id; prev[p2Id] = pnew1Id;

				angleSet.erase(now); angle[p0Id] = 2 * PI; inSet[p0Id] = false;
				InsertAngle(p1Id, viewDir);
				InsertAngle(pnew2Id, viewDir);
				InsertAngle(pnew1Id, viewDir);
				InsertAngle(p2Id, viewDir);
			}
		}
	}

	for (std::set<std::pair<vtkIdType, double>>::iterator iter = angleSet.begin(); iter != angleSet.end(); iter++)
	{
		fout << "(" << iter->first << ", " << iter->second / PI * 180.0 << ") ";
	}
	fout << endl;

	fout.close();

	meshdata->SetPoints(meshPoint);
	meshdata->SetLines(meshLine);

	std::cout << "Insert ended" << endl;

	return meshdata;
}

double InsertPoints::PointAngle(vtkIdType ptId, iPoint viewDir)
{
	//std::cout << ptId << endl;

	//std::cout << ptId << ' ' << prev.size() << ' ' << succ.size() << endl;

	iPoint p0 = iPoint(meshPoint->GetPoint(ptId));
	iPoint p1 = iPoint(meshPoint->GetPoint(prev[ptId]));
	iPoint p2 = iPoint(meshPoint->GetPoint(succ[ptId]));
	return GetTheta(p1 - p0, p2 - p0, viewDir);
}

void InsertPoints::InsertAngle(vtkIdType ptId, iPoint viewDir)
{
	if (inSet[ptId])
	{
		std::pair<vtkIdType, double> pr = std::make_pair(ptId, angle[ptId]);
		std::set<std::pair<vtkIdType, double>>::iterator iter = angleSet.find(pr);
		if (iter != angleSet.end()) angleSet.erase(iter);
		inSet[ptId] = false;
	}
	angle[ptId] = PointAngle(ptId, viewDir);
	if (angle[ptId] < PI)
	{
		std::pair<vtkIdType, double> pr = std::make_pair(ptId, angle[ptId]);
		angleSet.insert(pr);
		inSet[ptId] = true;
	}
}

void InsertPoints::AddLine(vtkIdType p1Id, vtkIdType p2Id)
{
	vtkSmartPointer<vtkLine> singleLine = vtkSmartPointer<vtkLine>::New();
	singleLine->GetPointIds()->SetId(0, p1Id);
	singleLine->GetPointIds()->SetId(1, p2Id);
	meshLine->InsertNextCell(singleLine);
}