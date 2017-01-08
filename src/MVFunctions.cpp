#include "MVFunctions.h"

double RealScalarProduct(double *& vec1, double *& vec2, const int size)
{
	double result = 0;
	for (int i = 0; i < size; ++i)
	{
		result += vec1[i] * vec2[i];
	}
	return result;
}

void ComplexScalarConjugateProduct(double *& x, double *& y, double * result, const int blockSize)
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

void RealMultiplyArrayScalar(double *& vec, const double scalar, double *& result, const int size)
{
	for (int i = 0; i < size; ++i)
	{
		result[i] = vec[i] * scalar;
	}
}

void ComplexMultiplyArrayScalar(double *& vec, const double * skalar, double *& result, const int blockSize)
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

void MultiplyComplexNumbers(const double * x, const double * y, double * result)
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

void SubtractArrays(double *& x, double *& y, double *& result, const int size)
{
	for (int i = 0; i < size; ++i)
	{
		result[i] = x[i] - y[i];
	}
}

void SummArrays(double *& x, double *& y, double *& result, const int size)
{
	for (int i = 0; i < size; ++i)
	{
		result[i] = x[i] + y[i];
	}
}

void CopyDoubleArray(double *& source, double *& dest, const int size)
{
	for (int i = 0; i < size; ++i)
	{
		dest[i] = source[i];
	}
}

void CopyIntArray(int *& source, int *& dest, const int size)
{
	for (int i = 0; i < size; ++i)
	{
		dest[i] = source[i];
	}
}

void DiagonalPreconditioning(const double * x, double *result)
{
	double complex_one[2];
	complex_one[0] = 1.0;
	complex_one[1] = 0.0;
	DivideComplexNumbers(complex_one, x, result);
}

double Norm(double * x, const int blockSize)
{
	double result = 0;
	for (int i = 0; i < blockSize; ++i)
	{
		result += x[2 * i] * x[2 * i] + x[2 * i + 1] * x[2 * i + 1];
	}
	return sqrt(result);
}

void MultDiOnVect(double *& di, double *& vec, double *& result, const int blockSize)
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

void MultiplyRarefiedMatrixOnVector(
	        int    *& ig
	,       int    *& jg
	,       double *& gg
	,       double *& di
	,       int    *& ijg
	,       int    *& idi
	,       double *& x
	,       double *& y
	, const int       blockSize
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
