#include "LLT.h"
#include "MVFunctions.h"

void ResizeL(
	        int    *& LLT_ig
	,       int    *& LLT_jg
	,       int    *& LLT_ijg
	,       int    *& LLT_idi
	,       double *& LLT_gg
	,       double *& LLT_di
	,       int    *& jg
	,       int    *& ig
	, const int       blockSize
	)
{
	LLT_ig = new int[blockSize + 1];
	LLT_jg = new int[ig[blockSize]];

	CopyIntArray(ig, LLT_ig, blockSize + 1);
	CopyIntArray(jg, LLT_jg, ig[blockSize]);
	
	LLT_ijg = new int[ ig[blockSize] + 1 ];
	LLT_ijg[0] = 0;
	//��� ����� � ���������� ����� ������ 2
	for (int i = 1; i < ig[blockSize] + 1; ++i)
	{
		LLT_ijg[i] = LLT_ijg[i - 1] + 2;
	}
	LLT_gg  = new double[ LLT_ijg[ LLT_ig[blockSize] ] ];
	LLT_idi = new int   [ blockSize + 1 ];
	LLT_idi[0] = 0;
	for (int i = 1; i < blockSize + 1; ++i)
	{
		LLT_idi[i] = LLT_idi[i - 1] + 2;
	}

	LLT_di = new double[ LLT_idi[blockSize] ];
	for (int i = 0; i < LLT_idi[blockSize]; ++i)
	{
		LLT_di[i] = 0;
	}
	
}

void SqrtComplex(double *ab, double*xy)
{
	xy[0] = sqrt((ab[0] + sqrt(ab[0] * ab[0] + ab[1] * ab[1])) / 2);
	xy[1] = ab[1] / (2 * xy[0]);
}

void LLT_Factorization(
	        int    *& ig
	,       int    *& jg
	,       int    *& ijg
	,       int    *& idi
	,       double *& gg
	,       double *& di
	,       int    *& LLT_ig
	,       int    *& LLT_jg
	,       int    *& LLT_ijg
	,       int    *& LLT_idi
	,       double *& LLT_gg
	,       double *& LLT_di
	, const int       blockSize
	)
{
	ResizeL(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_gg, LLT_di, jg, ig, blockSize);
	double sum[2];
	int size;
	for (int i = 0; i < blockSize; ++i)
	{
		int ib0 = ig[i];
		int ib1 = ig[i + 1];
		//���� �� ������ ��� ����������
		//gg
		for (int m = ib0; m < ib1; ++m)
		{
			int j = jg[m];
			int jb0 = ig[j];
			int jb1 = ig[j + 1];
			sum[0] = 0;
			sum[1] = 0;
			int i_cur = ib0;
			//��������� ��������� ������ � ������-�������
			for (int k = jb0; k < jb1; ++k)
			{
				while (jg[i_cur] < jg[k])
				{
					i_cur++;
				}
				if (jg[i_cur] == jg[k])
				{
					size = LLT_ijg[j+1] - LLT_ijg[j];
					MultiplyBlock(&LLT_gg[LLT_ijg[k]], sum, &LLT_gg[LLT_ijg[i_cur]], size);
				}
			}
			//new elem 
			//SubtractComplexNumbers(&gg[ijg[m]], sum, sum);
			size = ijg[m + 1] - ijg[m];
			if (size == 1)
			{
				sum[0] = gg[ijg[m]] - sum[0];
			}
			else
			{
				sum[0] = gg[ijg[m]  ] - sum[0];
				sum[1] = gg[ijg[m]+1] - sum[1];
			}

			DivideComplexNumbers(sum, &LLT_di[LLT_idi[j]], sum);
			LLT_gg[LLT_ijg[m]  ] = sum[0];
			LLT_gg[LLT_ijg[m]+1] = sum[1];
		}
		//di
		sum[0] = 0;
		sum[1] = 0;
		for (int k = ib0; k < ib1; ++k)
		{
			size = LLT_ijg[k + 1] - LLT_ijg[k];
			MultiplyBlock(&LLT_gg[LLT_ijg[k]], sum, &LLT_gg[LLT_ijg[k]], size);
		}
		double tmp[2];
		size = idi[i + 1] - idi[i];
		if (size == 1)
		{
			sum[0] = di[idi[i]]-sum[0];
		}
		else
		{
			sum[0] = di[idi[i]  ] - sum[0];
			sum[1] = di[idi[i]+1] - sum[1];
		}
		SqrtComplex(sum, tmp);
		LLT_di[2 * i    ] = tmp[0];
		LLT_di[2 * i + 1] = tmp[1];
	}
}

void SLAE_Forward_Complex(
	        int    *& LLT_ig
	,       int    *& LLT_jg
	,       int    *& LLT_ijg
	,       int    *& LLT_idi
	,       double *& LLT_gg
	,       double *& LLT_di
	,       double *& rightPart
	,       double *& result
	, const int       blockSize
	)
{
	CopyDoubleArray(rightPart, result, 2 * blockSize);
	double tmp[2];
	for (int i = 0; i < blockSize; ++i)
	{
		int ib0 = LLT_ig[i];
		int ib1 = LLT_ig[i + 1];
		for (int m = ib0; m < ib1; ++m)
		{
			int j = LLT_jg[m];
			MultiplyComplexNumbers(&LLT_gg[LLT_ijg[m]], &result[2 * j], tmp);
			SubtractComplexNumbers(&result[2 * i], tmp, &result[2 * i]);
		}
		DivideComplexNumbers(&result[2 * i], &LLT_di[LLT_idi[i]], &result[2 * i]);
	}
}

void SLAE_Backward_Complex(
	        int    *& LLT_ig
	,       int    *& LLT_jg
	,       int    *& LLT_ijg
	,       int    *& LLT_idi
	,       double *& LLT_gg
	,       double *& LLT_di
	,       double *& rightPart
	,       double *& result
	, const int       blockSize
	)
{
	CopyDoubleArray(rightPart, result, 2 * blockSize);
	double tmp[2];
	for (int j = blockSize - 1; j >= 0; j--)
	{
		DivideComplexNumbers(&result[2 * j], &LLT_di[LLT_idi[j]], &result[2 * j]);
		int jb0 = LLT_ig[j];
		int jb1 = LLT_ig[j+1];
		for (int i = jb0; i < jb1; ++i)
		{
			int ib = LLT_jg[i];
			tmp[0] = 0;
			tmp[1] = 0;
			MultiplyBlock(&LLT_gg[LLT_ijg[i]], tmp, &result[2 * j], 2);
			SubtractComplexNumbers(&result[2 * ib], tmp, &result[2 * ib]);
		}
	}
}
