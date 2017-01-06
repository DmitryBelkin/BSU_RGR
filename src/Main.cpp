#include "CocgComplex.h"

void InputSparseComplexBlockSLAE(
	char* dirname,
	vector < double > &di,
	vector < double > &ggl,
	vector < int    > &ig,
	vector < int    > &jg,
	vector < int    > &idi,
	vector < int    > &ijg,
	vector < double > &right_part,
	int               &FullTaskSize,
	int               &BlockSize,
	double            &eps,
	int               &MaxIter,
	vector < double > &VectorForCheckMatrixVectorSLAEMultiplication
	)
{
	char filename[500];
	//kuslau
	sprintf_s(filename, "%s/kuslau", dirname);
	FILE *fp;
	fopen_s(&fp, filename, "r");
	fscanf_s(fp, "%d%le%d", &FullTaskSize, &eps, &MaxIter);
	fclose(fp);

	BlockSize = FullTaskSize / 2;

	//ig
	ig.resize(BlockSize + 1);
	sprintf_s(filename, "%s/ig", dirname);
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < BlockSize + 1; ++i)
	{
		fread(&ig[i], sizeof(int), 1, fp);
		ig[i]--;
	}
	fclose(fp);

	//idi
	sprintf_s(filename, "%s/idi", dirname);
	idi.resize(BlockSize + 1);
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < BlockSize + 1; ++i)
	{
		fread(&idi[i], sizeof(int), 1, fp);
		idi[i]--;
	}
	fclose(fp);

	//jg
	sprintf_s(filename, "%s/jg", dirname);
	jg.resize(ig[BlockSize]);
	fopen_s(&fp, filename, "rb");
	for (int i = 0, size = jg.size(); i < size; ++i)
	{
		fread(&jg[i], sizeof(int), 1, fp);
		jg[i]--;
	}
	fclose(fp);
	//ijg
	sprintf_s(filename, "%s/ijg", dirname);
	ijg.resize(ig[BlockSize]+1);
	fopen_s(&fp, filename, "rb");
	for (int i = 0, size = ijg.size(); i < size; ++i)
	{
		fread(&ijg[i], sizeof(int), 1, fp);
		ijg[i]--;
	}
	fclose(fp);
	//ggl
	sprintf_s(filename, "%s/gg", dirname);
	ggl.resize(ijg[ig[BlockSize]]);
	fopen_s(&fp, filename, "rb");
	for (int i = 0, size = ggl.size(); i < size; ++i)
	{
		fread(&ggl[i], sizeof(double), 1, fp);	
	}
	fclose(fp);
	//di
	sprintf_s(filename, "%s/di", dirname);
	di.resize(idi[BlockSize]);
	fopen_s(&fp, filename, "rb");
	for (int i = 0, size = di.size(); i < size; ++i)
	{
		fread(&di[i], sizeof(double), 1, fp);
	}
	fclose(fp);
	//right part
	sprintf_s(filename, "%s/pr", dirname);
	right_part.resize(FullTaskSize);
	fopen_s(&fp, filename, "rb");
	for (int i = 0, size = right_part.size(); i < size; ++i)
	{
		fread(&right_part[i], sizeof(double), 1, fp);
	}
	fclose(fp);
	//VectorForCheckMatrixVectorSLAEMultiplication
	sprintf_s(filename, "%s/y.txt", dirname);
	VectorForCheckMatrixVectorSLAEMultiplication.resize(FullTaskSize);
	fopen_s(&fp, filename, "rb");
	for (int i = 0, size = VectorForCheckMatrixVectorSLAEMultiplication.size(); i < size; ++i)
	{
		fread(&VectorForCheckMatrixVectorSLAEMultiplication[i], sizeof(double), 1, fp);
	}
	fclose(fp);
}

void main()
{
	vector < double > Cdi, Cggl, Cright_part, Cresult;
	vector < int    > Cidi, Cig, Cijg, Cjg;
	vector < double > VectorForCheckMatrixVectorSLAEMultiplication;
	int    MaxIter, FullTaskSize, BlockSize;
	double eps;

	InputSparseComplexBlockSLAE("../resources/IN_data", Cdi, Cggl, Cig, Cjg,
		Cidi, Cijg, Cright_part, FullTaskSize, BlockSize, eps,
		MaxIter, VectorForCheckMatrixVectorSLAEMultiplication);

	Cresult.resize(FullTaskSize);

	FILE *fp;
	fopen_s(&fp, "../resources/OUT_data/nev1.txt", "w");
	fclose(fp);
	fopen_s(&fp, "../resources/OUT_data/nev2.txt", "w");
	fclose(fp);

	COCG(Cig, Cjg, Cggl, Cdi, 
		Cijg, Cidi, Cright_part, 
		BlockSize, Cresult,eps, MaxIter);
}