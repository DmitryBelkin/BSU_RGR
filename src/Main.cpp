#include "CocgComplex.h"

const char * pathToResidual1 = "../resources/OUT_data/nev1.txt";
const char * pathToResidual2 = "../resources/OUT_data/nev2.txt";
const char * pathToSlaeDir   = "../resources/IN_data";

void InputSlae(
	  vector < double > &di
	, vector < double > &ggl
	, vector < int    > &ig
	, vector < int    > &jg
	, vector < int    > &idi
	, vector < int    > &ijg
	, vector < double > &rightPart
	, int               &slaeDimension
	, int               &blockSize
	, double            &epsilon
	, int               &maxiter
	, vector < double > &check
	)
{
	char filename[100];
	sprintf_s(filename, "%s/kuslau", pathToSlaeDir);
	FILE *fp;
	fopen_s(&fp, filename, "r");
	fscanf_s(fp, "%d%le%d", &slaeDimension, &epsilon, &maxiter);
	fclose(fp);

	blockSize = slaeDimension / 2;

	ig.resize(blockSize + 1);
	sprintf_s(filename, "%s/ig", pathToSlaeDir);
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < blockSize + 1; ++i)
	{
		fread(&ig[i], sizeof(int), 1, fp);
		ig[i]--;
	}
	fclose(fp);

	sprintf_s(filename, "%s/idi", pathToSlaeDir);
	idi.resize(blockSize + 1);
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < blockSize + 1; ++i)
	{
		fread(&idi[i], sizeof(int), 1, fp);
		idi[i]--;
	}
	fclose(fp);

	sprintf_s(filename, "%s/jg", pathToSlaeDir);
	jg.resize(ig[blockSize]);
	fopen_s(&fp, filename, "rb");
	for (int i = 0, size = jg.size(); i < size; ++i)
	{
		fread(&jg[i], sizeof(int), 1, fp);
		jg[i]--;
	}
	fclose(fp);

	sprintf_s(filename, "%s/ijg", pathToSlaeDir);
	ijg.resize(ig[blockSize]+1);
	fopen_s(&fp, filename, "rb");
	for (int i = 0, size = ijg.size(); i < size; ++i)
	{
		fread(&ijg[i], sizeof(int), 1, fp);
		ijg[i]--;
	}
	fclose(fp);

	sprintf_s(filename, "%s/gg", pathToSlaeDir);
	ggl.resize(ijg[ig[blockSize]]);
	fopen_s(&fp, filename, "rb");
	for (int i = 0, size = ggl.size(); i < size; ++i)
	{
		fread(&ggl[i], sizeof(double), 1, fp);	
	}
	fclose(fp);

	sprintf_s(filename, "%s/di", pathToSlaeDir);
	di.resize(idi[blockSize]);
	fopen_s(&fp, filename, "rb");
	for (int i = 0, size = di.size(); i < size; ++i)
	{
		fread(&di[i], sizeof(double), 1, fp);
	}
	fclose(fp);

	sprintf_s(filename, "%s/pr", pathToSlaeDir);
	rightPart.resize(slaeDimension);
	fopen_s(&fp, filename, "rb");
	for (int i = 0, size = rightPart.size(); i < size; ++i)
	{
		fread(&rightPart[i], sizeof(double), 1, fp);
	}
	fclose(fp);

	sprintf_s(filename, "%s/y.txt", pathToSlaeDir);
	check.resize(slaeDimension);
	fopen_s(&fp, filename, "rb");
	for (int i = 0, size = check.size(); i < size; ++i)
	{
		fread(&check[i], sizeof(double), 1, fp);
	}
	fclose(fp);
}

void main()
{
	vector < double > di, ggl, rightPart, result, check;
	vector < int    > idi, ig, ijg, jg;
	int               maxiter, slaeDimension, blockSize;
	double            epsilon;

	InputSlae(di, ggl, ig, jg, idi, ijg, rightPart, slaeDimension, blockSize, epsilon, maxiter, check);

	FILE *fp;
	fopen_s(&fp, pathToResidual1, "w"); fclose(fp);
	fopen_s(&fp, pathToResidual2, "w"); fclose(fp);

	result.resize(slaeDimension);
	CocgComplex(ig, jg, ggl, di, ijg, idi, rightPart, blockSize, result, epsilon, maxiter);
}
