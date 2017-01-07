#include "CocgComplex.h"

const char * pathToResidual1 = "../resources/output/nev1.txt";
const char * pathToResidual2 = "../resources/output/nev2.txt";
const char * pathToSlaeDir   = "../resources/input";

void InputSlae(
	  vector < double > &di
	, vector < double > &gg
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
	gg.resize(ijg[ig[blockSize]]);
	fopen_s(&fp, filename, "rb");
	for (int i = 0, size = gg.size(); i < size; ++i)
	{
		fread(&gg[i], sizeof(double), 1, fp);	
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
	vector < double > di, gg, rightPart, result, check;
	vector < int    > idi, ig, ijg, jg;
	int               maxiter, slaeDimension, blockSize;
	double            epsilon;

	InputSlae(di, gg, ig, jg, idi, ijg, rightPart, slaeDimension, blockSize, epsilon, maxiter, check);

	printf("slaeDimension\t=\t%d\n", slaeDimension);
	printf("maxiter\t=\t%d\n", maxiter);
	printf("epsilon\t=\t%le\n", epsilon);
	//epsilon = 1e-16;
	FILE *fp;
	fopen_s(&fp, pathToResidual1, "w"); fclose(fp);
	fopen_s(&fp, pathToResidual2, "w"); fclose(fp);

	result.resize(slaeDimension);
	CocgComplex(ig, jg, gg, di, ijg, idi, rightPart, blockSize, result, epsilon, maxiter);

	// output
	fopen_s(&fp, "../resources/output/true_soulution.txt", "w");
	for (int i = 0, size = check.size(); i < size; ++i)
	{
		fprintf(fp, "%lf\n", check[i]);
	}
	fclose(fp);

	fopen_s(&fp, "../resources/output/result_soulution.txt", "w");
	for (int i = 0, size = result.size(); i < size; ++i)
	{
		fprintf(fp, "%lf\n", result[i]);
	}
	fclose(fp);

	fopen_s(&fp, "../resources/output/diff_soulutions.txt", "w");
	for (int i = 0, size = result.size(); i < size; ++i)
	{
		fprintf(fp, "%lf\n", check[i] - result[i]);
	}
	fclose(fp);
}
