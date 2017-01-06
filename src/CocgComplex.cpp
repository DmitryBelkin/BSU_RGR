#include "CocgComplex.h"
#include "LLT.h"
#include <omp.h>

const char * timeMeasurements = "../resources/timeMeasurements.txt";

void OutputIterationsAndResidual(const double residual, const int iteration, const char *filename)
{
	FILE *fp;
	fopen_s(&fp, filename, "a");
	fprintf(fp, "%d\t%.15le\n", iteration, residual);
	fclose(fp);
}

void COCG(
	vector < int    > &ig        ,
	vector < int    > &jg        ,
	vector < double > &ggl       ,
	vector < double > &di        ,
	vector < int    > &ijg       ,
	vector < int    > &idi       ,
	vector < double > &rightPart,
	int                nb        ,
	vector < double > &result    ,
	double             eps       ,
	int                maxiter
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
	for (i = 0; i < 2 * nb; ++i)
	{
		result[i] = 0.0;
	}
	for (i = 0; i < nb; i++)
	{
		int size = idi[i + 1] - idi[i];
		if (size == 2)
		{
			DiagonalPreconditioning(&di[idi[i]], &inverse_di[2 * i]);
		}
		else
		{
			inverse_di[2 * i] = 1.0 / di[idi[i]];
			inverse_di[2 * i + 1] = 0;
		}
	}

	double nev;
	double nevSecond;
	MultiplyRarefiedMatrixOnVector(ig, jg, ggl, di, ijg, idi, result, temp, nb);
	s.resize(2 * nb);
	y.resize(2 * nb);
	SubtractVectors(rightPart, temp, r);
	nev = Norm(r, nb);
	int flag = 0;
	double norm0;
	norm0 = Norm(rightPart, nb);
	nev /= norm0;
	nevSecond = nev;
	OutputIterationsAndResidual(nev, 0, "../resources/OUT_data/nev1.txt");
	OutputIterationsAndResidual(nevSecond, 0, "../resources/OUT_data/nev2.txt");
	if (FactorizationType == 1)
	{
		MultDiOnVect(inverse_di, r, z, nb);
	}
	else
	{
		SLAE_Forward_Complex (LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, r, z, nb);
		SLAE_Backward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, z, z, nb);
	}

/*	p = z;
	s = r;
	y = result;*/
	CopyVector(&z[0], &p[0], nb);
	CopyVector(&r[0], &s[0], nb);
	CopyVector(&result[0], &y[0], nb);
	int iter = 0;
	double timeStart = omp_get_wtime();
	while (nevSecond > eps && iter < maxiter)
	{
		ComplexScalarConjugateProduct(r, z, complex_number1, nb);
		MultiplyRarefiedMatrixOnVector(ig, jg, ggl, di, ijg, idi, p, Ap, nb);
		ComplexScalarConjugateProduct(Ap, p, complex_number2, nb);
		DivideComplexNumbers(complex_number1, complex_number2, alfa);

		ComplexMultiplyVectorScalar(p, alfa, temp, nb);
		SummVectors(result, temp, result);//xj+1

		ComplexMultiplyVectorScalar(Ap, alfa, temp, nb);
		SubtractVectors(r, temp, r);//rj+1
		//etta
		SubtractVectors(r, s, temp_spline);
		etta = -RealScalarProduct(temp_spline, s) / RealScalarProduct(temp_spline, temp_spline);
		if (etta > 1.0)
		{
			flag = 1;
		/*	y = result;
			s = r;*/
			CopyVector(&result[0], &y[0], nb);
			CopyVector(&r[0], &s[0], nb);
		}
		else if (etta > 0)
		{
			flag = 1;
			SubtractVectors(result, s, temp_spline);
			RealMultiplyVectorScalar(temp_spline, etta, temp_spline);
			SummVectors(y, temp_spline, y);

			SubtractVectors(r, s, temp_spline);
			RealMultiplyVectorScalar(temp_spline, etta, temp_spline);
			SummVectors(s, temp_spline, s);
		}

		
		if (FactorizationType == 1)
		{
			MultDiOnVect(inverse_di, r, z, nb);//zj+1*/
		}
		else
		{
			SLAE_Forward_Complex (LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, r, z, nb);
			SLAE_Backward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, z, z, nb);
		}

		ComplexScalarConjugateProduct(r, z, complex_number2, nb);//beta
		DivideComplexNumbers(complex_number2, complex_number1, beta);

		ComplexMultiplyVectorScalar(p, beta, temp, nb);
		SummVectors(z, temp, p);

		nev = Norm(r, nb);
		nev /= norm0;
		nevSecond = Norm(s, nb);
		nevSecond /= norm0;
		OutputIterationsAndResidual(nev, iter+1, "../resources/OUT_data/nev1.txt");
		OutputIterationsAndResidual(nevSecond, iter+1, "../resources/OUT_data/nev2.txt");
		printf("nev = %le\r", nev);
		iter++;
	}
	if (flag)
	{
		/*y = result;*/
		CopyVector(&result[0], &y[0], nb);
	}
	MultiplyRarefiedMatrixOnVector(ig, jg, ggl, di, ijg, idi, result, temp, nb);
	SubtractVectors(temp, rightPart, temp);
	nev = Norm(temp, nb);
	nev /= norm0;
	OutputIterationsAndResidual(nev      , iter + 1, "../resources/OUT_data/nev1.txt");
	OutputIterationsAndResidual(nevSecond, iter + 1, "../resources/OUT_data/nev2.txt");
	printf("\nEnd nev = %le\n", nev);
	printf("\nIters\t=\t%d\n", iter);
	double timeEnd = omp_get_wtime();
	FILE *fp;
	fopen_s(&fp, timeMeasurements, "a");
	fprintf(fp, "%.15le\n", timeEnd - timeStart);
	fclose(fp);
}
