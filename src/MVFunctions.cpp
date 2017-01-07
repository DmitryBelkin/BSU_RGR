#include "MVFunctions.h"

double RealScalarProduct(const vector < double > &vec1, const vector < double > &vec2)
{
	double result = 0;
	for (int i = 0, size = vec1.size(); i < size; ++i)
	{
		result += vec1[i] * vec2[i];
	}
	return result;
}

void ComplexScalarConjugateProduct(const vector < double > &x, const vector < double > &y, double  *result, const int blockSize)
{
	result[0] = 0.0;
	result[1] = 0.0;
	double temporaryResult[2];
	for (int i = 0; i < blockSize; ++i)
	{
		MultiplyComplexNumbers(&x[2 * i], &y[2 * i], temporaryResult);
		result[0] += temporaryResult[0];
		result[1] += temporaryResult[1];
	}
}

void RealMultiplyVectorScalar(const vector < double > &vec, const double scalar, vector < double > &result)
{
	for (int i = 0, size = vec.size(); i < size; ++i)
	{
		result[i] = vec[i] * scalar;
	}
}

void ComplexMultiplyVectorScalar(const vector < double > &vec, const double *skalar, vector < double > &result, const int blockSize)
{
	for (int i = 0; i < blockSize; ++i)
	{
		MultiplyComplexNumbers(&vec[2 * i], skalar, &result[2 * i]);
	}
}

void MultiplyBlock(const double *x, double *dest, const double *a, const int size)
{
	if (size == 1)
	{
		dest[0] += a[0] * x[0];
		dest[1] += a[0] * x[1];
	}
	if (size == 2)
	{
		dest[0] += a[0] * x[0] - a[1] * x[1];
		dest[1] += a[0] * x[1] + a[1] * x[0];
	}
}

void MultiplyComplexNumbers(const double *x, const double *y, double *result)
{
	result[0] = x[0] * y[0] - x[1] * y[1];
	result[1] = x[0] * y[1] + x[1] * y[0];
}

void DivideComplexNumbers(const double *x, const double *y, double *result)
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

void SubtractComplexNumbers(const double *x, const double *y, double *result)
{
	result[0] = x[0] - y[0];
	result[1] = x[1] - y[1];
}

void SubtractVectors(const vector < double > &x, const vector < double > &y, vector < double > &result)
{
	for (int i = 0, size = x.size(); i < size; ++i)
	{
		result[i] = x[i] - y[i];
	}
}

void SummVectors(const vector < double > &x, const vector < double > &y, vector < double > &result)
{
	for (int i = 0, size = x.size(); i < size; ++i)
	{
		result[i] = x[i] + y[i];
	}
}

void DiagonalPreconditioning(const double * x, double *result)
{
	double complex_one[2];
	complex_one[0] = 1.0;
	complex_one[1] = 0.0;
	DivideComplexNumbers(complex_one, x, result);
}

double Norm(const vector < double > x, const int blockSize)
{
	double result = 0;
	for (int i = 0; i < blockSize; ++i)
	{
		result += x[2 * i] * x[2 * i] + x[2 * i + 1] * x[2 * i + 1];
	}
	return sqrt(result);
}

void MultDiOnVect(const vector < double > &di, const vector < double > &vec, vector < double > &result, const int blockSize)
{
	double tmp[2];
	for (int i = 0; i < blockSize; ++i)
	{
		tmp[0] = 0;
		tmp[1] = 0;
		MultiplyBlock(&di[2 * i], tmp, &vec[2 * i], 2);
		result[2 * i    ] = tmp[0];
		result[2 * i + 1] = tmp[1];
	}
}

void MultVMatrixOnVector(const vector < vector < double > > &V, const vector < double > vec, vector < double > &result, const int blockSize, const int m)
{
	for (int i = 0; i < blockSize; ++i)
	{
		result[2 * i   ] = 0;
		result[2 * i +1] = 0;
		for (int j = 0; j < m; ++j)
		{
			MultiplyBlock(&V[j][2 * i], &result[2 * i], &vec[2 * j], 2);
		}
	}	
}

void MultiplyRarefiedMatrixOnVector(
	  const vector < int    > &ig
	, const vector < int    > &jg
	,       vector < double > &gg
	,       vector < double > &di
	, const vector < int    > &ijg
	, const vector < int    > &idi
	,       vector < double > &x
	,       vector < double > &y
	, const int blockSize
	)
{
	for (int i = 0; i < blockSize * 2; ++i) { y[i] = 0; }
	for (int i = 0; i < blockSize    ; ++i)
	{
		int size = idi[i + 1] - idi[i];
		MultiplyBlock(&x[i * 2], &y[i * 2], &di[idi[i]], size);
	}
	for (int i = 0; i < blockSize; ++i)
	{
		for (int j = ig[i]; j < ig[i + 1]; ++j)
		{
			int size = ijg[j + 1] - ijg[j];
			MultiplyBlock(&x[jg[j] * 2], &y[i     * 2], &gg[ijg[j]], size);
			MultiplyBlock(&x[i     * 2], &y[jg[j] * 2], &gg[ijg[j]], size);
		}
	}
}
