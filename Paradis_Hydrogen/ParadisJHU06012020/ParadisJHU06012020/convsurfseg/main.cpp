#include "stdio.h"
#include "string"
#include "Tools.h"
#include "Point.h"
#include "list"

using namespace std;
using namespace SupportSystem;
using namespace EZ;

void ProcessGroup(const string& sFileName,const unsigned int& iMaxStep,list<Point>& loPoints)
{
	FILE* fpInput = fopen(sFileName.c_str(),"r");
	string sRead = "";
	unsigned int iStep = 0;
	double dX1 = 0.0;
	double dY1 = 0.0;
	double dZ1 = 0.0;
	double dX2 = 0.0;
	double dY2 = 0.0;
	double dZ2 = 0.0;
	double dTemp = 0.0;
	double dTolerance = 1.0E0;
	Point oPoint1;
	Point oPoint2;
	while(!feof(fpInput))
	{
		sRead = GetRealString(512,fpInput);
		sscanf(sRead.c_str(),"%d,%*d : (%*d,%*d) -> (%*d,%*d) : (%*f,%*f,%*f) / (%*f,%*f,%*f)\n",&iStep);
		if(iStep > iMaxStep)
		{
			break;
		}
		sRead = GetRealString(512,fpInput);
		sscanf(sRead.c_str(),"(%lf,%lf,%lf) -> (%lf,%lf,%lf) : (%*f,%*f,%*f) -> (%*f,%*f,%*f)\n",&dX1,&dY1,&dZ1,&dX2,&dY2,&dZ2);
		dTemp = (dX2 - dX1)*(dX2 - dX1) + (dY2 - dY1)*(dY2 - dY1) + (dZ2 - dZ1)*(dZ2 - dZ1);
		oPoint1.Set(dX1,dY1,dZ1);
		oPoint2.Set(dX2,dY2,dZ2);
		if(oPoint1.Distance(oPoint2) < dTolerance)
		{
			continue;
		}
		loPoints.push_back(oPoint1);
		loPoints.push_back(oPoint2);
	}
	fclose(fpInput);
	printf("%s : done\n",sFileName.c_str());
}

void WriteParaviewFile(list<Point>* ploPoints)
{
	unsigned int iPointsCount = (unsigned int)ploPoints->size();
	if(iPointsCount < 2)
	{
		printf("no points to write!!\n");
		return;
	}
	FILE* fpFile = fopen("output.vtk","w");
	fprintf(fpFile,"# vtk DataFile Version 1.0\n");
	fprintf(fpFile,"dislocation microstructure\n");
	fprintf(fpFile,"ASCII\n");
	fprintf(fpFile,"DATASET POLYDATA\n");
	fprintf(fpFile,"POINTS %d float\n",iPointsCount);
	list<Point>::iterator liPoints;
	for(liPoints = ploPoints->begin() ; liPoints != ploPoints->end() ; liPoints++)
	{
		fprintf(fpFile,"%e\t%e\t%e\n",(*liPoints).GetX(),(*liPoints).GetY(),(*liPoints).GetZ());
	}
	unsigned int iLinesCount = iPointsCount/2;
	fprintf(fpFile,"LINES %d %d\n",iLinesCount,3*iLinesCount);
	unsigned int i = 0;
	for(i = 0 ; i < iLinesCount ; i++)
	{
		fprintf(fpFile,"2\t%d\t%d\n",2*i,2*i + 1);
	}
	fprintf(fpFile,"POINT_DATA %d\n",iPointsCount);
	fprintf(fpFile,"scalars NodeID integer\n");
	fprintf(fpFile,"LOOKUP_TABLE default\n");
	for(i = 0 ; i < iPointsCount  ; i++)
	{
		fprintf(fpFile,"%d\n",1);
	}
	fprintf(fpFile,"CELL_DATA %d\n",iLinesCount);
	fprintf(fpFile,"SCALARS SegmentID integer\n");
	fprintf(fpFile,"LOOKUP_TABLE default\n");
	for(i = 0 ; i < iLinesCount ; i++)
	{
		fprintf(fpFile,"%d\n",2);
	}
	fclose(fpFile);
}

void CleanPointsList(list<Point>* ploPoints)
{
	list<Point>::iterator liPoints;
	liPoints = ploPoints->begin();
	Point oPoint1;
	Point oPoint2;
	double dTolerance = 1.0E-3;
	while(liPoints != ploPoints->end())
	{
		oPoint1 = (*liPoints);
		liPoints++;
		if(liPoints == ploPoints->end())
		{
			break;
		}
		oPoint2 = (*liPoints);
		liPoints--;
		if(oPoint1.Distance(oPoint2) < dTolerance)
		{
			// erase the two points
			liPoints = ploPoints->erase(liPoints);
			liPoints = ploPoints->erase(liPoints);
		}
		else
		{
			liPoints++;
		}
	}
}

int main(int argc,char** argv)
{
	if(argc != 5)
	{
		printf("error: few arguments\n");
		printf("usage: conv groups_count start_group input_filename_base max_step\n");
		return 1;
	}

	char cFileName[512];
	string sFileName;
	unsigned int iGroupOrder = atoi(argv[2]);
	unsigned int iMaxStep = atoi(argv[4]);
	list<Point> loPoints;
	while(iGroupOrder < atoi(argv[1]))
	{
		sprintf(cFileName,"%s_%d.txt",argv[3],iGroupOrder);
		sFileName = cFileName;
		ProcessGroup(sFileName,iMaxStep,loPoints);
		iGroupOrder = iGroupOrder + 1;
	}
	CleanPointsList(&loPoints);
	WriteParaviewFile(&loPoints);
	loPoints.clear();

	return 0;
}

