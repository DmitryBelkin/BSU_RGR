#include "MVFunctions.h"

double RealScalarProduct(vector < double > &vec1,vector < double > &vec2)
{
	double result = 0;
	for (int i = 0, size = vec1.size(); i < size; ++i)
	{
		result += vec1[i] * vec2[i];
	}
	return result;
}

void RealMultiplyVectorScalar(vector < double > &vec, double scalar, vector < double > &result)
{
	for (int i = 0, size = vec.size(); i < size; ++i)
	{
		result[i] = vec[i] * scalar;
	}
}

void MultiplyBlock(double *x, double *y, double *a, int size)
{
	if (size == 1)
	{
		y[0] += a[0] * x[0];
		y[1] += a[0] * x[1];
	}
	if (size == 2)
	{
		y[0] += a[0] * x[0] - a[1] * x[1];
		y[1] += a[0] * x[1] + a[1] * x[0];
	}
}

void MultiplyRarefiedMatrixOnVector(
	vector < int    > &ig ,
	vector < int    > &jg ,
	vector < double > &ggl,
	vector < double > &di ,
	vector < int    > &ijg,
	vector < int    > &idi,
	vector < double > &x  ,
	vector < double > &y  ,
	int nb
	)
{
	for (int i = 0; i < nb * 2; ++i) { y[i] = 0; }
	for (int i = 0; i < nb    ; ++i)
	{
		int size = idi[i + 1] - idi[i];
		MultiplyBlock(&x[i * 2], &y[i * 2], &di[idi[i]], size);
	}
	for (int i = 0; i < nb; ++i)
	{
		for (int j = ig[i]; j < ig[i + 1]; ++j)
		{
			int size = ijg[j + 1] - ijg[j];
			MultiplyBlock(&x[jg[j] * 2], &y[i     * 2], &ggl[ijg[j]], size);
			MultiplyBlock(&x[i     * 2], &y[jg[j] * 2], &ggl[ijg[j]], size);
		}
	}
}

void MultiplyComplexNumbers(double *x, double *y, double *result)
{
	result[0] = x[0] * y[0] - x[1] * y[1];
	result[1] = x[0] * y[1] + x[1] * y[0];
}

void DivideComplexNumbers(double *x, double *y, double *result)
{
	double tmp[2];
	tmp[0] = x[0];
	tmp[1] = x[1];

	result[0] = tmp[0] * y[0] + tmp[1] * y[1];
	result[1] = tmp[1] * y[0] - tmp[0] * y[1];

	double norm = y[0] * y[0] + y[1] * y[1];
	result[0] /= norm;
	result[1] /= norm;
}
void DiagonalPreconditioning(double * x, double *result)
{
	double complex_one[2];
	complex_one[0] = 1.0;
	complex_one[1] = 0.0;
	DivideComplexNumbers(complex_one, x, result);
}

void ComplexScalarConjugateProduct(vector < double > &x, vector < double > &y, double  *result, int nb)
{
	result[0] = 0.0;
	result[1] = 0.0;
	double temporaryResult[2];
	for (int i = 0; i < nb; ++i)
	{
		MultiplyComplexNumbers(&x[2 * i], &y[2 * i], temporaryResult);
		result[0] += temporaryResult[0];
		result[1] += temporaryResult[1];
	}
}

double Norm(vector < double > x, int nb)
{
	double result = 0;
	for (int i = 0; i < nb; ++i)
	{
		result += x[2 * i] * x[2 * i] + x[2 * i + 1] * x[2 * i + 1];
	}
	return sqrt(result);
}

void SubtractVectors(vector < double > &x, vector < double > &y, vector < double > &result)
{
	for (int i = 0, size = x.size(); i < size; ++i)
	{
		result[i] = x[i] - y[i];
	}
}

void SummVectors(vector < double > &x, vector < double > &y, vector < double > &result)
{
	for (int i = 0, size = x.size(); i < size; ++i)
	{
		result[i] = x[i] + y[i];
	}
}

void CopyVector(double *x, double *y, int size) // копирование комплексного вектора
{
	for (int i = 0; i < 2 * size; i++)
	{
		y[i] = x[i];
	}
}

void ComplexMultiplyVectorScalar(vector < double > &vec, double *skalar, vector < double > &result, int nb)
{
	for (int i = 0; i < nb; ++i)
	{
		MultiplyComplexNumbers(&vec[2 * i], skalar, &result[2 * i]);
	}
}

void SubtractComplexNumbers(double *x, double *y, double *result)
{
	result[0] = x[0] - y[0];
	result[1] = x[1] - y[1];
}

void MultDiOnVect(vector < double > &di, vector < double > &vec, vector < double > &res, int Nb)
{
	double tmp[2];
	for (int i = 0; i < Nb; ++i)
	{
		tmp[0] = 0;
		tmp[1] = 0;
		MultiplyBlock(&di[2 * i], tmp, &vec[2 * i], 2);
		res[2 * i    ] = tmp[0];
		res[2 * i + 1] = tmp[1];
	}
}

void MultVMatrixOnVector(vector < vector < double > > &V, vector < double > vec, vector < double > &res, int Nb, int m)
{
	for (int i = 0; i < Nb; ++i)
	{
		res[2 * i   ] = 0;
		res[2 * i +1] = 0;
		for (int j = 0; j < m; ++j)
		{
			MultiplyBlock(&V[j][2 * i], &res[2 * i], &vec[2 * j], 2);
		}
	}	
}
