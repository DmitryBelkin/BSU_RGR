#include "IO.h"

void InputSparseComplexBlockSLAE(
	char* dirname,
	vector < double > &di,
	vector < double > &ggl,
	vector < int > &ig,
	vector < int > &jg,
	vector < int > &idi,
	vector < int > &ijg,
	vector < double > &righr_part,
	int &FullTaskSize,
	int &BlockSize,
	double &eps,
	int &MaxIter,
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
	for (int i = 0; i < BlockSize + 1; i++)
	{
		fread(&ig[i], sizeof(int), 1, fp);
		ig[i]--;
	}
	fclose(fp);

	//idi
	sprintf_s(filename, "%s/idi", dirname);
	idi.resize(BlockSize + 1);
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < BlockSize + 1; i++)
	{
		fread(&idi[i], sizeof(int), 1, fp);
		idi[i]--;
	}
	fclose(fp);

	//jg
	sprintf_s(filename, "%s/jg", dirname);
	jg.resize(ig[BlockSize]);
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < jg.size(); i++)
	{
		fread(&jg[i], sizeof(int), 1, fp);
		jg[i]--;
	}
	fclose(fp);
	//ijg
	sprintf_s(filename, "%s/ijg", dirname);
	ijg.resize(ig[BlockSize]+1);
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < ijg.size(); i++)
	{
		fread(&ijg[i], sizeof(int), 1, fp);
		ijg[i]--;
	}
	fclose(fp);
	//ggl
	sprintf_s(filename, "%s/gg", dirname);
	ggl.resize(ijg[ig[BlockSize]]);
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < ggl.size(); i++)
	{
		fread(&ggl[i], sizeof(double), 1, fp);	
	}
	fclose(fp);
	//di
	sprintf_s(filename, "%s/di", dirname);
	di.resize(idi[BlockSize]);
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < di.size(); i++)
	{
		fread(&di[i], sizeof(double), 1, fp);
	}
	fclose(fp);
	//right part
	sprintf_s(filename, "%s/pr", dirname);
	righr_part.resize(FullTaskSize);
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < righr_part.size(); i++)
	{
		fread(&righr_part[i], sizeof(double), 1, fp);
	}
	fclose(fp);
	//VectorForCheckMatrixVectorSLAEMultiplication
	sprintf_s(filename, "%s/y.txt", dirname);
	VectorForCheckMatrixVectorSLAEMultiplication.resize(FullTaskSize);
	fopen_s(&fp, filename, "rb");
	for (int i = 0; i < VectorForCheckMatrixVectorSLAEMultiplication.size(); i++)
	{
		fread(&VectorForCheckMatrixVectorSLAEMultiplication[i], sizeof(double), 1, fp);
	}
	fclose(fp);
}

	//void inputSparceBlockMatrix_from_txt(
	//	char* filename,
	//	vector < double > &di,
	//	vector < double > &ggl,
	//	vector < int > &ig,
	//	vector < int > &jg,
	//	vector < int > &idi,
	//	vector < int > &ijg,
	//	vector < double > &righr_part,
	//	int &FullTaskSize,
	//	int &BlockSize,
	//	double &eps,
	//	int &MaxIter
	//	)
	//{
	//	
	//	
	//	FILE *fp;
	//	fp = fopen_s(filename, "r");
	//	fscanf_s(fp, "%d%le%d", &FullTaskSize, &eps, &MaxIter);
	//	BlockSize = FullTaskSize / 2;
	//	//ig
	//	ig.resize(BlockSize + 1);	
	//	for (int i = 0; i < BlockSize + 1; i++)
	//	{
	//		fscanf_s(fp, "%d", &ig[i]);
	//		ig[i]--;
	//	}
	//	//idi		
	//	idi.resize(BlockSize + 1);	
	//	for (int i = 0; i < BlockSize + 1; i++)
	//	{
	//		fscanf_s(fp, "%d", &idi[i]);
	//		idi[i]--;
	//	}
	//	//jg	
	//	jg.resize(ig[BlockSize]);	
	//	for (int i = 0; i < jg.size(); i++)
	//	{
	//		fscanf_s(fp, "%d", &jg[i]);
	//		jg[i]--;
	//	}	
	//	//ijg		
	//	ijg.resize(ig[BlockSize] + 1);		
	//	for (int i = 0; i < ijg.size(); i++)
	//	{
	//		fscanf_s(fp, "%d", &ijg[i]);
	//		ijg[i]--;
	//	}	
	//	//ggl		
	//	ggl.resize(ijg[ig[BlockSize]]);		
	//	for (int i = 0; i < ggl.size(); i++)
	//	{
	//		fscanf_s(fp, "%le", &ggl[i]);
	//	}		
	//	//di		
	//	di.resize(idi[BlockSize]);		
	//	for (int i = 0; i < di.size(); i++)
	//	{
	//		fscanf_s(fp, "%le",&di[i]);
	//	}		
	//	//right part		
	//	righr_part.resize(FullTaskSize);		
	//	for (int i = 0; i < righr_part.size(); i++)
	//	{
	//		fscanf_s(fp, "%le", &righr_part[i]);
	//	}
	//	fclose(fp);
	//}

////////////////////////////////////////////////////////////////////////////////////////////////////
	void OutputItersResidual(
		double Residual, 
		int IterNumber,
		char *filename
		)
	{
		FILE *fp;
		if (IterNumber == 0)
		{
			fopen_s(&fp, filename, "w");
		}
		else
		{
			fopen_s(&fp, filename, "a");
		}
		fprintf(fp, "%d\t%.15le\n", IterNumber, Residual);
		fclose(fp);
	}