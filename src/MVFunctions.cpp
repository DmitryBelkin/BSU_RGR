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

void RealMultVector_by_Scalar(vector < double > &vec, double scalar, vector < double > &result)
{
	for (int i = 0, size = vec.size(); i < size; ++i)
	{
		result[i] = vec[i] * scalar;
	}
}

void mult_block(double *x, double *y, double *a, int size)
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

void mult_MV(
	vector < int    > &ig,
	vector < int    > &jg,
	vector < double > &ggl,
	vector < double > &di,
	vector < int    > &ijg,
	vector < int    > &idi,
	vector < double > &x,
	vector < double > &y,
	int nb
	)
{
	for (int i = 0; i < nb * 2; ++i) { y[i] = 0; }
	for (int i = 0; i < nb    ; ++i)
	{
		int size = idi[i + 1] - idi[i];
		mult_block(&x[i * 2], &y[i * 2], &di[idi[i]], size);
	}
	for (int i = 0; i < nb; ++i)
	{
		for (int j = ig[i]; j < ig[i + 1]; ++j)
		{
			int size = ijg[j + 1] - ijg[j];
			mult_block(&x[jg[j] * 2], &y[i     * 2], &ggl[ijg[j]], size);
			mult_block(&x[i     * 2], &y[jg[j] * 2], &ggl[ijg[j]], size);
		}
	}
}

void mult_complex_numbers(double *x, double *y, double *res)
{
	res[0] = x[0] * y[0] - x[1] * y[1];
	res[1] = x[0] * y[1] + x[1] * y[0];
}

void div_complex_numbers(double *x, double *y, double *res)
{
	double tmp[2];
	tmp[0] = x[0];
	tmp[1] = x[1];

	res[0] = tmp[0] * y[0] + tmp[1] * y[1];
	res[1] = tmp[1] * y[0] - tmp[0] * y[1];

	double norm = y[0] * y[0] + y[1] * y[1];
	res[0] /= norm;
	res[1] /= norm;
}
void di_preconditioning(double * x, double *res)
{
	double complex_one[2];
	complex_one[0] = 1.0;
	complex_one[1] = 0.0;
	div_complex_numbers(complex_one, x, res);
}

void conjugate_SK_mult(vector < double > &x, vector < double > &y, double  *res, int nb)
{
	res[0] = 0.0;
	res[1] = 0.0;
	double temp_res[2];
	for (int i = 0; i < nb; ++i)
	{
		mult_complex_numbers(&x[2 * i], &y[2 * i], temp_res);
		res[0] += temp_res[0];
		res[1] += temp_res[1];
	}
}

void norm(vector < double > x, double &res, int nb)
{
	res = 0;
	for (int i = 0; i < nb; ++i)
	{
		res += x[2 * i] * x[2 * i] + x[2 * i + 1] * x[2 * i + 1];
	}
	res = sqrt(res);
}

void sub_vec(vector < double > &x, vector < double > &y, vector < double > &res)
{
	for (int i = 0, size = x.size(); i < size; ++i)
	{
		res[i] = x[i] - y[i];
	}
}

void sum_vec(vector < double > &x, vector < double > &y, vector < double > &res)
{
	for (int i = 0, size = x.size(); i < size; ++i)
	{
		res[i] = x[i] + y[i];
	}
}

void copy_vec(double *x, double *y, int size) // копирование комплексного вектора
{
	for (int i = 0; i < 2 * size; i++)
	{
		y[i] = x[i];
	}
}

void mult_vector_by_scalar(vector < double > &x, double *skalar, vector < double > &res, int nb)
{
	for (int i = 0; i < nb; ++i)
	{
		mult_complex_numbers(&x[2 * i], skalar, &res[2 * i]);
	}
}

void mult_DI_matrix_vector(vector < double > &inverse_di, vector < double > &r, vector < double > &res, int nb)
{
	for (int i = 0; i < nb; ++i)
	{
		mult_complex_numbers(&inverse_di[2 * i], &r[2 * i], &res[2 * i]);
	}
}

void sub_complex_numbers(double *a, double *b, double *res)
{
	res[0] = a[0] - b[0];
	res[1] = a[1] - b[1];
}

void Complex_SK_mult(vector < double > &x, vector < double > &y, double *res, int nb)
{
	res[0] = 0;
	res[1] = 0;
	for (int i = 0; i < nb; ++i)
	{
		res[0] += x[2 * i] * y[2 * i    ] + x[2 * i + 1] * y[2 * i + 1];
		res[1] += x[2 * i] * y[2 * i + 1] - x[2 * i + 1] * y[2 * i    ];
	}
}

void MultDiOnVect(vector < double > &di, vector < double > &vec, vector < double > &res, int Nb)
{
	double tmp[2];
	for (int i = 0; i < Nb; ++i)
	{
		tmp[0] = 0;
		tmp[1] = 0;
		mult_block(&di[2 * i], tmp, &vec[2 * i], 2);
		res[2 * i    ] = tmp[0];
		res[2 * i + 1] = tmp[1];
	}
}

void MultLTMatrixOnVector(
	vector < int    > &LLT_ig,
	vector < int    > &LLT_jg,
	vector < int    > &LLT_ijg,
	vector < double > &LLT_ggl,
	vector < double > &LLT_di,
	vector < double > &vec,
	vector < double > &result,
	int Nb
	)
{
	double tmp[2];
	MultDiOnVect(LLT_di, vec, result, Nb);
	for (int i = 0; i < Nb; ++i)
	{
		int jb0 = LLT_ig[i    ];
		int jb1 = LLT_ig[i + 1];
		for (int m = jb0; m < jb1; ++m)
		{
			int i = LLT_jg[m];
			tmp[0] = 0;
			tmp[1] = 0;
			mult_block(&LLT_ggl[LLT_ijg[m]], tmp, &vec[2 * i], 2);
			result[2 * i    ] += tmp[0];
			result[2 * i + 1] += tmp[0];
		}
	}
}

void MultVMatrixOnVector(vector < vector < double > > &V,vector < double > vec,vector < double > &res,int Nb, int m)
{
	for (int i = 0; i < Nb; ++i)
	{
		res[2 * i   ] = 0;
		res[2 * i +1] = 0;
		for (int j = 0; j < m; ++j)
		{
			mult_block(&V[j][2 * i], &res[2 * i], &vec[2 * j], 2);
		}
	}	
}

void NormalizeVec(vector < double > &r, double norma, int Nb)
{
	for (int i = 0; i < 2 * Nb; ++i)
	{
		r[i] /= norma;
	}
}
