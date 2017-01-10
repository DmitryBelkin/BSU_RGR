#include "CocgComplex.h"
#include <mkl.h>
#include <omp.h>

const char * pathToResidual1 = "../resources/output/nev1.txt";
const char * pathToResidual2 = "../resources/output/nev2.txt";
const char * pathToSlaeDir   = "../resources/input";

void InputSlae(
	  double *& di
	, double *& gg
	, int    *& ig
	, int    *& jg
	, int    *& idi
	, int    *& ijg
	, double *& rightPart
	, int     & slaeDimension
	, int     & blockSize
	, double  & epsilon
	, int     & maxiter
	, double *& check
	)
{
	char filename[100];
	sprintf_s(filename, "%s/kuslau", pathToSlaeDir);
	FILE *fp;
	fopen_s(&fp, filename, "r");
	fscanf_s(fp, "%d%le%d", &slaeDimension, &epsilon, &maxiter);
	fclose(fp);

	blockSize = slaeDimension / 2;

	ig = new int[blockSize+1];
	sprintf_s(filename, "%s/ig", pathToSlaeDir);
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < blockSize + 1; ++i)
	{
		fread(&ig[i], sizeof(int), 1, fp);
		ig[i]--;
	}
	fclose(fp);

	sprintf_s(filename, "%s/idi", pathToSlaeDir);
	idi = new int[blockSize+1];
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < blockSize + 1; ++i)
	{
		fread(&idi[i], sizeof(int), 1, fp);
		idi[i]--;
	}
	fclose(fp);

	sprintf_s(filename, "%s/jg", pathToSlaeDir);
	jg = new int[ig[blockSize]];
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < ig[blockSize]; ++i)
	{
		fread(&jg[i], sizeof(int), 1, fp);
		jg[i]--;
	}
	fclose(fp);

	sprintf_s(filename, "%s/ijg", pathToSlaeDir);
	ijg = new int[ig[blockSize]+1];
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < ig[blockSize]+1; ++i)
	{
		fread(&ijg[i], sizeof(int), 1, fp);
		ijg[i]--;
	}
	fclose(fp);

	sprintf_s(filename, "%s/gg", pathToSlaeDir);
	gg = new double[ijg[ig[blockSize]]];
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < ijg[ig[blockSize]]; ++i)
	{
		fread(&gg[i], sizeof(double), 1, fp);	
	}
	fclose(fp);

	sprintf_s(filename, "%s/di", pathToSlaeDir);
	di = new double[idi[blockSize]];
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < idi[blockSize]; ++i)
	{
		fread(&di[i], sizeof(double), 1, fp);
	}
	fclose(fp);

	sprintf_s(filename, "%s/pr", pathToSlaeDir);
	rightPart = new double[slaeDimension];
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < slaeDimension; ++i)
	{
		fread(&rightPart[i], sizeof(double), 1, fp);
	}
	fclose(fp);

	sprintf_s(filename, "%s/y.txt", pathToSlaeDir);
	check = new double[slaeDimension];
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < slaeDimension; ++i)
	{
		fread(&check[i], sizeof(double), 1, fp);
	}
	fclose(fp);
}

void main()
{
	double * di        = NULL;
	double * gg        = NULL;
	double * rightPart = NULL;
	double * result    = NULL;
	double * check     = NULL;
	int    * idi       = NULL;
	int    * ig        = NULL;
	int    * ijg       = NULL;
	int    * jg        = NULL;
	int      maxiter = 0, slaeDimension = 0, blockSize = 0;
	double   epsilon = 0;

	InputSlae(di, gg, ig, jg, idi, ijg, rightPart, slaeDimension, blockSize, epsilon, maxiter, check);

	printf("slaeDimension\t=\t%d\n", slaeDimension);
	printf("maxiter\t\t=\t%d\n", maxiter);
	printf("epsilon\t\t=\t%le\n", epsilon);

	FILE *fp;
	fopen_s(&fp, pathToResidual1, "w"); fclose(fp);
	fopen_s(&fp, pathToResidual2, "w"); fclose(fp);

	result = new double[slaeDimension];
	mkl_set_num_threads(2);
	omp_set_num_threads(2);
	CocgComplex(ig, jg, gg, di, ijg, idi, rightPart, blockSize, result, epsilon, maxiter);

	// output
	fopen_s(&fp, "../resources/output/true_soulution.txt", "w");
	for (int i = 0; i < slaeDimension; ++i)
	{
		fprintf(fp, "%lf\n", check[i]);
	}
	fclose(fp);

	fopen_s(&fp, "../resources/output/result_soulution.txt", "w");
	for (int i = 0; i < slaeDimension; ++i)
	{
		fprintf(fp, "%lf\n", result[i]);
	}
	fclose(fp);

	fopen_s(&fp, "../resources/output/diff_soulutions.txt", "w");
	for (int i = 0; i < slaeDimension; ++i)
	{
		fprintf(fp, "%lf\n", check[i] - result[i]);
	}
	fclose(fp);

	delete [] di;
	delete [] gg;
	delete [] rightPart;
	delete [] result;
	delete [] check;
	delete [] idi;
	delete [] ig;
	delete [] ijg;
	delete [] jg;
}
