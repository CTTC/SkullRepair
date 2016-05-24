#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#include "vtkSmartPointer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkSTLReader.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkCell.h"
#include "vtkLineSource.h"
#include "vtkLine.h"
#include "vtkCellArray.h"
#include <vtkVersion.h>
#include <vtkPoints.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetMapper.h>
#include <vtkTriangleFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include "insertPoints.h"
#include "vtkSTLReader.h"
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vector>
#include <fstream>
#include <map>
#include <set>

#include "EdgePicking.h"
#include "InsertPoints.h"
#include "iGeometry.h"
#include "main.h"
#include "linequ.h"
#include "stdarg.h"

void GetPointNormal(vtkSmartPointer<vtkPolyData> polydata, std::vector<iPoint> &ptNormal)
{
	std::cout << "GetPointNormal started...\n";
	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
	normalGenerator->SetInputData(polydata);
	normalGenerator->ComputePointNormalsOn();
	normalGenerator->ComputeCellNormalsOff();
	normalGenerator->Update();
	polydata = normalGenerator->GetOutput();

	vtkDataArray* normalsGeneric = polydata->GetPointData()->GetNormals(); //works

	//std::cout << normalsGeneric->GetNumberOfTuples();

	for (vtkIdType ptId = 0; ptId < normalsGeneric->GetNumberOfTuples(); ptId++)
	{
		double m[3] = {0,0,0};
		normalsGeneric->GetTuple(ptId, m);
		iPoint p(m);
		ptNormal.push_back(p);
	}
	std::cout << "GetPointNormal ended\n";
}



int main ()
{
	InputData(polydata);

	GetCellNormal(polydata, cellNormal);
	LableLine(polydata, linedata, nextPoly, lineMap);
	GetViewVect(viewDir);

	vtkSmartPointer<vtkPolyData> viewdata = vtkSmartPointer<vtkPolyData>::New();
	GetViewData(viewdata);

	EdgePicking::GetResult(polydata, cellNormal, linedata, nextPoly, lineMap, viewDir,
		uBound, lBound);

	vtkSmartPointer<vtkPolyData> udata = vtkSmartPointer<vtkPolyData>::New();
	GetBeltData(polydata, uBound, udata);
	/*
	std::vector<iPoint> uPoints;
	uPoints.clear();
	for (std::vector<vtkIdType>::iterator iter = uBound.begin(); iter != uBound.end(); iter++)
	{
		uPoints.push_back(iPoint(polydata->GetPoint(*iter)));
	}


	InsertPoints *hInsPt = new(InsertPoints);
	vtkSmartPointer<vtkPolyData> uMeshdata = hInsPt->Insert(uPoints, viewDir);
	delete(hInsPt);
	*/

	//reader
	std::string inputFilename = "e:/cranio.stl";

	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(inputFilename.c_str());
	reader->SetOutput(reader->GetOutput());
	reader->Update();

	vtkSmartPointer<vtkPolyDataMapper> readermapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	readermapper->SetInputConnection(reader->GetOutputPort());

	vtkSmartPointer<vtkActor> readeractor =
		vtkSmartPointer<vtkActor>::New();
	readeractor->SetMapper(readermapper);
	readeractor->GetProperty()->SetRepresentationToWireframe();
	readeractor->GetProperty()->SetLineWidth(1);

	//show edge points;
	std::set<vtkIdType> edge_point_id;
	std::vector<iPoint> ptnormals;
	std::vector<iPoint> testpt;
	std::vector<iPoint> edgept;
	GetPointNormal(polydata,ptnormals);

	int point_num = uBound.size();
	double x,y,z,mx=0,my=0,mz=0;
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> vertices =
		vtkSmartPointer<vtkCellArray>::New();

	vtkIdType pid[1];
	std::cout<<"size: "<<uBound.size()<<endl;
	for (int i = 0; i < uBound.size(); i++)
	{
		double pp[3];
		polydata->GetPoint(uBound[i],pp);
		x = pp[0];
		y = pp[1];
		z = pp[2];
		iPoint p(x,y,z);
		edgept.push_back(p);
		mx += x;
		my += y;
		mz += z;
		points->InsertPoint(i, x, y, z);
		pid[0] = i;
		vertices->InsertNextCell(1, pid);
		edge_point_id.insert(polydata->FindPoint(x,y,z));

	}
	mx /= point_num;
	my /= point_num;
	mz /= point_num;
	vtkSmartPointer<vtkPolyData> point =
		vtkSmartPointer<vtkPolyData>::New();
	point->SetPoints(points);
	point->SetVerts(vertices);

	vtkSmartPointer<vtkPolyDataMapper> pointmapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();

	pointmapper->SetInputData(point);

	vtkSmartPointer<vtkActor> pointactor =
		vtkSmartPointer<vtkActor>::New();
	pointactor->SetMapper(pointmapper);
	pointactor->GetProperty()->SetPointSize(10);
	pointactor->GetProperty()->SetColor(0, 0, 1);

	//middle point
	iPoint mp(mx,my,mz);
	for (int i=0;i<edgept.size();++i)
	{
		iPoint v = mp -  edgept[i];
		iPoint p = edgept[i] + 5*(v/v.module());
		testpt.push_back(p);
		p = edgept[i] + 10*(v/v.module());
		testpt.push_back(p);
		p = edgept[i] + 15*(v/v.module());
		testpt.push_back(p);
		p = edgept[i] + 20*(v/v.module());
		testpt.push_back(p);
		p = edgept[i] + 25*(v/v.module());
		testpt.push_back(p);
		p = edgept[i] + 30*(v/v.module());
		testpt.push_back(p);
	}
	vtkSmartPointer<vtkPolyData> middle_point =
		vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> middle_points =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> middle_vertices =
		vtkSmartPointer<vtkCellArray>::New();
	pid[0] = middle_points->InsertNextPoint(mx, my, mz);
	middle_vertices->InsertNextCell(1,pid);
	middle_point->SetPoints(middle_points);
	middle_point->SetVerts(middle_vertices);
	vtkSmartPointer<vtkPolyDataMapper> middle_pointmapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();

	middle_pointmapper->SetInputData(middle_point);

	vtkSmartPointer<vtkActor> middle_pointactor =
		vtkSmartPointer<vtkActor>::New();
	middle_pointactor->SetMapper(middle_pointmapper);
	middle_pointactor->GetProperty()->SetPointSize(5);
	middle_pointactor->GetProperty()->SetColor(1, 0, 0);
	//show cell
	std::set<vtkIdType> neighbors;
	std::set<vtkIdType> resultPoints;
	std::set<vtkIdType> tmpPoints;
	std::set<vtkIdType> reducedPoints;
	std::vector<iPoint> inPoints,outPoints;
	for(std::set<vtkIdType>::iterator it1 = edge_point_id.begin(); it1 != edge_point_id.end(); it1++)
	{
		iPoint p(polydata->GetPoint(*it1));
		double distance = p.distance_square(mp);
		iPoint v = p-mp;
		vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		polydata->GetPointCells(*it1,cellIds);
		for(vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
		{
			if(neighbors.find(cellIds->GetId(i))==neighbors.end())
			{
				polydata->GetCellPoints(cellIds->GetId(i),ptIds);
				bool isright = true;
				iPoint p1(polydata->GetPoint(ptIds->GetId(0)));
				iPoint p2(polydata->GetPoint(ptIds->GetId(1)));
				iPoint normal = p1^p2;
				for(vtkIdType j = 0; j < ptIds->GetNumberOfIds(); j++)
				{
					iPoint tmp(polydata->GetPoint(ptIds->GetId(j)));
					if(tmp.distance_square(mp)<distance)
						isright = false;
				}
				double angel = GetAngle(v,normal);
				if(isright && angel>PI/4 && angel<3*PI/4)
				{
					neighbors.insert(cellIds->GetId(i));
					for(vtkIdType j = 0; j < ptIds->GetNumberOfIds(); j++)
					{
						reducedPoints.insert(ptIds->GetId(j));
						if(edge_point_id.find(ptIds->GetId(j))==edge_point_id.end())
							tmpPoints.insert(ptIds->GetId(j));
					}
				}
			}
		}
	}
	resultPoints = tmpPoints;
	for(int i=0;i<4;++i)
	{
		std::set<vtkIdType> tmpPtIds = tmpPoints;
		tmpPoints.clear();
		for(std::set<vtkIdType>::iterator it1 = tmpPtIds.begin(); it1 != tmpPtIds.end(); it1++)
		{
			iPoint p(polydata->GetPoint(*it1));
			double distance = p.distance_square(mp);
			iPoint v = p-mp;
			vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
			vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
			polydata->GetPointCells(*it1,cellIds);
			for(vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
			{
				if(neighbors.find(cellIds->GetId(i))==neighbors.end())
				{
					polydata->GetCellPoints(cellIds->GetId(i),ptIds);
					bool isright = true;
					for(vtkIdType j = 0; j < ptIds->GetNumberOfIds(); j++)
					{
						iPoint tmp(polydata->GetPoint(ptIds->GetId(j)));
						if(tmp.distance_square(mp)<distance)
							isright = false;
					}
					if(isright)
					{
						neighbors.insert(cellIds->GetId(i));
						for(vtkIdType j = 0; j < ptIds->GetNumberOfIds(); j++)
						{
							if(reducedPoints.find(ptIds->GetId(j))==reducedPoints.end())
								tmpPoints.insert(ptIds->GetId(j));
							resultPoints.insert(ptIds->GetId(j));
							reducedPoints.insert(ptIds->GetId(j));
						}
					}
				}
			}
		}
	}
	for (std::set<vtkIdType>::iterator it1 = resultPoints.begin(); it1 != resultPoints.end(); it1++)
	{
		iPoint p(polydata->GetPoint(*it1));
		iPoint n = ptnormals[*it1];
		n = n/n.module();
		double l = 1.0;
		iPoint inpt,outpt;
		inpt = p - l*n;
		outpt = p + l*n;
		inPoints.push_back(inpt);
		outPoints.push_back(outpt);
	}
	vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
	ids->SetNumberOfComponents(1);
	for(std::set<vtkIdType>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
	{
		ids->InsertNextValue(*it1);
	}

	vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
	selectionNode->SetFieldType(vtkSelectionNode::CELL);
	selectionNode->SetContentType(vtkSelectionNode::INDICES);
	selectionNode->SetSelectionList(ids);

	vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
	selection->AddNode(selectionNode);

	vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
	extractSelection->SetInputConnection(0, reader->GetOutputPort());
	extractSelection->SetInputData(1, selection);
	extractSelection->Update();

	vtkSmartPointer<vtkDataSetMapper> neighborCellsMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	neighborCellsMapper->SetInputConnection(extractSelection->GetOutputPort());

	vtkSmartPointer<vtkActor> neighborCellsActor = vtkSmartPointer<vtkActor>::New();
	neighborCellsActor->SetMapper(neighborCellsMapper);
	neighborCellsActor->GetProperty()->SetColor(0,1,0);


	//for(std::set<vtkIdType>::iterator it1 = resultPoints.begin(); it1 != resultPoints.end(); it1++)
	//{
	//	pid[0] = middle_points->InsertNextPoint(polydata->GetPoint(*it1));
	//	middle_vertices->InsertNextCell(1,pid);
	//}
	/*for(int i=0;i<outPoints.size();++i)
	{
		pid[0] = middle_points->InsertNextPoint(outPoints[i].x,outPoints[i].y,outPoints[i].z);
		middle_vertices->InsertNextCell(1,pid);
	}*/
	//calculate Phi
	std::cout<<"calculate phi\n";
	std::vector<iPoint> results;
	for (std::set<vtkIdType>::iterator it1 = resultPoints.begin(); it1 != resultPoints.end(); it1++)
	{
		iPoint p(polydata->GetPoint(*it1));
		results.push_back(p);
	}
	for (int i=0;i<outPoints.size();++i)
		results.push_back(outPoints[i]);
	for (int i=0;i<inPoints.size();++i)
		results.push_back(inPoints[i]);
	int size = results.size();
	double** phis = new double*[size];
	for (int i=0;i<size;++i)
	{
		phis[i] = new double[size];
	}
	for (int i=0;i<size;++i)
	{
		for (int j=0;j<size;++j)
		{
			iPoint pi = results[i];
			iPoint pj = results[j];
			double phi = 1./(sqrt(2.*PI)*SIGMA)*exp(-1./(2.*sqr(SIGMA))*(pow(pi.x-pj.x,2)+pow(pi.y-pj.y,2)+pow(pi.z-pj.z,2)));
			phis[i][j] = phi;
		}
	}
	//ax=b
	std::cout<<"ax=b\n";
	size += 4; 
	static double** a = new double*[size];//size = 3N+4
	int i;
	int j;
	for (i=0;i<size;++i)
	{
		a[i] = new double[size];
	}
	static double* b = new double[size];
	static double* aa = new double[size*size];
	std::cout<<"clear to 0\n";
	//clear a b to 0
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
			a[i][j] = 0;
	}
	for(i=0;i<size;i++)
	{
		b[i] = 0;
	}
	std::cout<<"input a\n";
	//input a
	for(i=0;i<size-4;i++)
		for(j=0;j<size-4;j++)
			a[i][j] = phis[i][j];
	for(i=0;i<size-4;i++)
	{
		a[i][size-4] = 1;
		a[size-4][i] = 1;
	}
	for(j=0;j<size-4;j++)
	{
		a[size-3][j] = results[j].x;
		a[size-2][j] = results[j].y;
		a[size-1][j] = results[j].z;
		a[j][size-3] = results[j].x;
		a[j][size-2] = results[j].y;
		a[j][size-1] = results[j].z;
	}
	//input aa
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
			aa[i*size+j] = a[i][j];
	}
	std::cout<<"size:"<<size<<"\n";
	delete[] a;
	//input b
	std::cout<<"input b\n";
	size = results.size();
	std::cout<<"size:"<<size<<"\n";
	for( i=size/3;i<2*size/3;i++)
		b[i] = 1;
	for( i=2*size/3;i<size;i++)
		b[i] = -1;
	//solve
	Linequ equl(size+4);
	equl.setLinequ(aa,b);
	//equl.printL();
	if(equl.Solve())
		equl.showX();
	else
		std::cout<<"no solve Fail "<<std::endl;



	for (int i=0;i<testpt.size();++i)
	{
		iPoint p = equl.getR(testpt[i],results);
		pid[0] = middle_points->InsertNextPoint(p.x,p.y,p.z);
		middle_vertices->InsertNextCell(1,pid);
	}

	vtkSmartPointer<vtkDelaunay3D> del =  vtkSmartPointer<vtkDelaunay3D>::New();
	vtkSmartPointer<vtkPolyData> polydata2 = 
		vtkSmartPointer<vtkPolyData>::New();
	polydata2->SetPoints(middle_points);
	polydata2->SetVerts(middle_vertices);
	del->SetInputData(polydata2);
	del->BoundingTriangulationOn();
	del->SetTolerance(1);
	del->SetAlpha(100);
	del->BoundingTriangulationOff();
	// Generate a tetrahedral mesh from the input points. By
	// default, the generated volume is the convex hull of the points.
	//vtkDelaunay3D Mapper
	vtkDataSetMapper *delMapper = vtkDataSetMapper::New();
	//vtkSmartPointer<vtkDataSetMapper>delMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	delMapper->SetInputConnection(del->GetOutputPort());

	//vtkDelaunay3D Actor
	vtkSmartPointer<vtkActor>delActor = vtkSmartPointer<vtkActor>::New();
	delActor->SetMapper(delMapper);
	delActor->GetProperty()->SetColor(1, 0, 0);

	vtkSmartPointer<vtkRenderer> originalRenderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderer> delaunayRenderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderer> delaunayAlphaRenderer =
		vtkSmartPointer<vtkRenderer>::New();


	//display 
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

//
	renderer->AddActor(pointactor);
	renderer->AddActor(neighborCellsActor);
	renderer->AddActor(middle_pointactor);
	renderer->AddActor(readeractor);
	renderer->SetBackground(0,0,0); // Background color green
	renderer->ResetCamera();
	renderer->AddActor(delActor);

//
	renderWindow->Render();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}





void InputData(vtkSmartPointer<vtkPolyData> &polydata)
{
	std::cout << "InputData started...";
	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetOutput(reader->GetOutput());
	reader->SetFileName("e:/cranio.stl");
	reader->Update();
	polydata->DeepCopy(reader->GetOutput());
	std::cout << "InputData ended" << endl;
}

void GetCellNormal(vtkSmartPointer<vtkPolyData> polydata, std::vector<iPoint> &cellNormal)
{
	std::cout << "GetCellNormal started...";
	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
	normalGenerator->SetInputData(polydata);
	normalGenerator->ComputePointNormalsOff();
	normalGenerator->ComputeCellNormalsOn();
	normalGenerator->Update();
	polydata = normalGenerator->GetOutput();

	vtkDataArray* normalsGeneric = polydata->GetCellData()->GetNormals(); //works

	//std::cout << normalsGeneric->GetNumberOfTuples();

	for (vtkIdType cellId = 0; cellId < normalsGeneric->GetNumberOfTuples(); cellId++)
	{
		double m[3];
		normalsGeneric->GetTuple(cellId, m);
		cellNormal.push_back(iPoint(m));
	}
	std::cout << "GetCellNormal ended" << endl;
}

void LableLine(vtkSmartPointer<vtkPolyData> polydata, std::vector< std::pair<vtkIdType, vtkIdType> > &linedata, std::vector< std::pair<vtkIdType, vtkIdType> > &nextPoly,
	std::map< std::pair<vtkIdType, vtkIdType>, vtkIdType> &lineMap)
{
	std::cout << "LableLine started...";
	linedata.clear();
	lineMap.clear();

	for (vtkIdType cellId = 0; cellId < polydata->GetNumberOfCells(); cellId++)
	{
		vtkIdType ptsId[3];
		vtkSmartPointer<vtkIdList> ptIdList = vtkSmartPointer<vtkIdList>::New();
		ptIdList = polydata->GetCell(cellId)->GetPointIds();
		for (vtkIdType i = 0; i < ptIdList->GetNumberOfIds(); i++)
		{
			ptsId[i] = ptIdList->GetId(i);
		}
		PushLine(ptsId[0], ptsId[1], cellId, lineMap, linedata, nextPoly);
		PushLine(ptsId[1], ptsId[2], cellId, lineMap, linedata, nextPoly);
		PushLine(ptsId[2], ptsId[0], cellId, lineMap, linedata, nextPoly);
	}
	std::cout << "LableLine ended" << endl;
}

void GetViewVect(iPoint &viewDir)
{
	std::cout << "GetViewVect started...";
	viewDir = VIEWPOINT1 - VIEWPOINT0;
	viewDir = viewDir / viewDir.module();
	std::cout << "GetViewVect ended" << endl;
}


void Display(int num, ...)
{
	va_list argu;
	va_start(argu, num);

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(0.2, 0.3, 0.4);

	for (int i = 0; i < num; i++)
	{
		DisplayPack pack = va_arg(argu, DisplayPack);

		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData(pack.data);
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(pack.r, pack.g, pack.b);
		actor->GetProperty()->SetLineWidth(pack.width);
		//std::cout << pack.width << ' ';
		actor->GetProperty()->SetPointSize(pack.size);
		renderer->AddActor(actor);
	}
	va_end(argu);

	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(renderer);
	renWin->SetSize(500, 500);

	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);

	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	iren->SetInteractorStyle(style);

	iren->Initialize();
	iren->Start();
}

void PushLine(vtkIdType p1, vtkIdType p2, vtkIdType cellId, std::map< std::pair<vtkIdType, vtkIdType>, vtkIdType> &lineMap,
	std::vector< std::pair<vtkIdType, vtkIdType> > &linedata, std::vector< std::pair<vtkIdType, vtkIdType> > &nextPoly)
{
	std::map < std::pair<vtkIdType, vtkIdType>, vtkIdType > ::iterator l_it;
	std::pair<vtkIdType, vtkIdType> lineNow;

	if (p1 < p2)
		lineNow = std::make_pair(p1, p2);
	else
		lineNow = std::make_pair(p2, p1);

	l_it = lineMap.find(lineNow);
	if (l_it == lineMap.end())
	{
		lineMap[lineNow] = linedata.size();
		linedata.push_back(lineNow);
		nextPoly.push_back(std::make_pair(cellId, -1));
	}
	else
	{
		nextPoly[lineMap[lineNow]].second = cellId;
	}
}

void GetViewData(vtkSmartPointer<vtkPolyData> &viewdata)
{
	vtkSmartPointer<vtkPoints> viewPoint = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> viewLine = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkLine> singleLine = vtkSmartPointer<vtkLine>::New();
	viewPoint->InsertPoint(0, VIEWPOINT0.x, VIEWPOINT0.y, VIEWPOINT0.z);
	viewPoint->InsertPoint(1, VIEWPOINT1.x, VIEWPOINT1.y, VIEWPOINT1.z);
	singleLine->GetPointIds()->SetId(0, 0);
	singleLine->GetPointIds()->SetId(1, 1);
	viewLine->InsertNextCell(singleLine);
	viewdata->SetPoints(viewPoint);
	viewdata->SetLines(viewLine);
}

void GetBeltData(vtkSmartPointer<vtkPolyData> polydata, std::vector<vtkIdType> belt, vtkSmartPointer<vtkPolyData> &beltdata)
{
	vtkSmartPointer<vtkPoints> beltPoint = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> beltLine = vtkSmartPointer<vtkCellArray>::New();

	for (vtkIdType ptId = 0; ptId < belt.size(); ptId++)
	{
		double m[3];
		polydata->GetPoint(belt[ptId], m);
		beltPoint->InsertPoint(ptId, m);
	}
	for (int i = 0; i < belt.size(); i++)
	{
		vtkSmartPointer<vtkLine> singleLine = vtkSmartPointer<vtkLine>::New();
		singleLine->GetPointIds()->SetId(0, i);
		singleLine->GetPointIds()->SetId(1, (i+1)%belt.size());
		beltLine->InsertNextCell(singleLine);
	}
	beltdata->SetPoints(beltPoint);
	beltdata->SetLines(beltLine);
}

void GetMeshData(std::vector<iPoint> mesh, vtkSmartPointer<vtkPolyData> &meshdata)
{
	vtkSmartPointer<vtkPoints> meshPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> meshVertices = vtkSmartPointer<vtkCellArray>::New();
	for (int i = 0; i < mesh.size(); i++)
	{
		vtkIdType pid[1];

		meshPoints->InsertPoint(i, mesh[i].x, mesh[i].y, mesh[i].z);
		pid[0] = i;
		meshVertices->InsertNextCell(1, pid);
	}
	meshdata->SetPoints(meshPoints);
	meshdata->SetVerts(meshVertices);
}