#include "COCG_for_Complex.h"
#include "LLT.h"
#include <omp.h>

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

void COCG(
	vector < int    > &ig        ,
	vector < int    > &jg        ,
	vector < double > &ggl       ,
	vector < double > &di        ,
	vector < int    > &ijg       ,
	vector < int    > &idi       ,
	vector < double > &right_part,
	int                nb        ,
	vector < double > &result    ,
	double             eps       ,
	int                MaxIter
	)
{
	//LLT sparce matrix factorization
	vector < double > LLT_ggl, LLT_di;
	vector < int > LLT_ig, LLT_jg, LLT_ijg, LLT_idi;
	LLT_Factorization(ig, jg, ijg, idi, ggl, di, LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, nb);

	vector < double > inverse_di;
	inverse_di.resize(2 * nb);
	vector < double > temp;
	vector < double > p, z, r;
	vector < double > Ap;
	double alfa[2], beta[2];
	double complex_number1[2], complex_number2[2];
	vector < double > y, s, temp_spline;
	double etta;
	p   .resize(2 * nb);
	z   .resize(2 * nb);
	r   .resize(2 * nb);
	Ap  .resize(2 * nb);
	temp.resize(2 * nb);
	temp_spline.resize(2 * nb);
	int i = 0;
	int FactorizationType = 2;//1 - di, 2 - LLT
	for (i = 0; i < 2 * nb; i++)
	{
		result[i] = 0.0;
	}
	for (i = 0; i < nb; i++)
	{
		int size = idi[i + 1] - idi[i];
		if (size == 2)
		{
			di_preconditioning(&di[idi[i]], &inverse_di[2 * i]);
		}
		else
		{
			inverse_di[2 * i] = 1.0 / di[idi[i]];
			inverse_di[2 * i + 1] = 0;
		}
	}

	double nev;
	double nevSecond;
	mult_MV(ig, jg, ggl, di, ijg, idi, result, temp, nb);
	s.resize(2 * nb);
	y.resize(2 * nb);
	sub_vec(right_part, temp, r);
	norm(r, nev, nb);
	int flag = 0;
	double norm0;
	norm(right_part, norm0, nb);
	nev /= norm0;
	nevSecond = nev;
	OutputItersResidual(nev, 0, "../resources/OUT_data/nev1.txt");
	OutputItersResidual(nevSecond, 0, "../resources/OUT_data/nev2.txt");
	if (FactorizationType == 1)
	{
		MultDiOnVect(inverse_di, r, z, nb);
	}
	else
	{
		SLAE_Forward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, r, z, nb);
		SLAE_Backward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, z, z, nb);
	}

/*	p = z;
	s = r;
	y = result;*/
	copy_vec(&z[0], &p[0], nb);
	copy_vec(&r[0], &s[0], nb);
	copy_vec(&result[0], &y[0], nb);
	int iter = 0;
	double t_start;
	double t_end;
	t_start = omp_get_wtime();
	while (nevSecond > eps && iter < MaxIter)
	{
		conjugate_SK_mult(r, z, complex_number1, nb);
		mult_MV(ig, jg, ggl, di, ijg, idi, p, Ap, nb);
		conjugate_SK_mult(Ap, p, complex_number2, nb);
		div_complex_numbers(complex_number1, complex_number2, alfa);

		mult_vector_by_scalar(p, alfa, temp, nb);
		sum_vec(result, temp, result);//xj+1

		mult_vector_by_scalar(Ap, alfa, temp, nb);
		sub_vec(r, temp, r);//rj+1
		//etta
		sub_vec(r, s, temp_spline);
		etta = -RealScalarProduct(temp_spline, s) / RealScalarProduct(temp_spline, temp_spline);
		if (etta > 1.0)
		{
			flag = 1;
		/*	y = result;
			s = r;*/
			copy_vec(&result[0], &y[0], nb);
			copy_vec(&r[0], &s[0], nb);
		}
		else if (etta > 0)
		{
			flag = 1;
			sub_vec(result, s, temp_spline);
			RealMultVector_by_Scalar(temp_spline, etta, temp_spline);
			sum_vec(y, temp_spline, y);

			sub_vec(r, s, temp_spline);
			RealMultVector_by_Scalar(temp_spline, etta, temp_spline);
			sum_vec(s, temp_spline, s);
		}

		
		if (FactorizationType == 1)
		{
			MultDiOnVect(inverse_di, r, z, nb);//zj+1*/
		}
		else
		{
			SLAE_Forward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, r, z, nb);
			SLAE_Backward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, z, z, nb);
		}
		/*//*/
		conjugate_SK_mult(r, z, complex_number2, nb);//beta
		div_complex_numbers(complex_number2, complex_number1, beta);

		mult_vector_by_scalar(p, beta, temp, nb);
		sum_vec(z, temp, p);

		norm(r, nev, nb);
		nev /= norm0;
		norm(s, nevSecond, nb);
		nevSecond /= norm0;		
		OutputItersResidual(nev, iter+1, "../resources/OUT_data/nev1.txt");
		OutputItersResidual(nevSecond, iter+1, "../resources/OUT_data/nev2.txt");
		printf("nev = %le\r", nev);
		iter++;
	}
	if (flag)
	{
		/*y = result;*/
		copy_vec(&result[0], &y[0], nb);
	}
	mult_MV(ig, jg, ggl, di, ijg, idi, result, temp, nb);
	sub_vec(temp, right_part, temp);
	norm(temp, nev, nb);
	nev /= norm0;	
	OutputItersResidual(nev, iter + 1, "../resources/OUT_data/nev1.txt");
	OutputItersResidual(nevSecond, iter + 1, "../resources/OUT_data/nev2.txt");
	printf("\nEnd nev = %le\n", nev);
	printf("\nIters\t=\t%d\n", iter);
	t_end = omp_get_wtime();
	FILE *fp;
	fopen_s(&fp, "Times.txt", "a");
	fprintf(fp, "%.15le\n", t_end - t_start);
	fclose(fp);
}
