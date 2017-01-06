#include "CocgComplex.h"
#include "LLT.h"
#include <omp.h>

// иначе LLT
#define DIAGONAL_FACTORIZATION 0
#define LLT_FACTORIZATION !DIAGONAL_FACTORIZATION

const char * timeMeasurements = "../resources/output/timeMeasurements.txt";

void OutputIterationsAndResidual(const double residual, const int iteration, const char *filename)
{
	FILE *fp;
	fopen_s(&fp, filename, "a");
	fprintf(fp, "%d\t%.15le\n", iteration, residual);
	fclose(fp);
}

void CocgComplex(
	  vector < int    > &ig
	, vector < int    > &jg
	, vector < double > &ggl
	, vector < double > &di
	, vector < int    > &ijg
	, vector < int    > &idi
	, vector < double > &rightPart
	, int                blockSize
	, vector < double > &result
	, double             epsilon
	, int                maxiter
	)
{
#if LLT_FACTORIZATION
	vector < double > LLT_ggl, LLT_di;
	vector < int    > LLT_ig, LLT_jg, LLT_ijg, LLT_idi;
	LLT_Factorization(ig, jg, ijg, idi, ggl, di, LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, blockSize);
#endif

	vector < double > inverseDi;
	inverseDi.resize(2 * blockSize);
	vector < double > temp;
	vector < double > p, z, r;
	vector < double > Ap;
	double alfa[2], beta[2];
	double complexNumber1[2], complexNumber2[2];
	vector < double > y, s, temp_spline;
	double etta;
	p   .resize(2 * blockSize);
	z   .resize(2 * blockSize);
	r   .resize(2 * blockSize);
	Ap  .resize(2 * blockSize);
	temp.resize(2 * blockSize);
	temp_spline.resize(2 * blockSize);
	int i = 0;
	for (i = 0; i < 2 * blockSize; ++i)
	{
		result[i] = 0.0;
	}
	for (i = 0; i < blockSize; ++i)
	{
		int size = idi[i + 1] - idi[i];
		if (size == 2)
		{
			DiagonalPreconditioning(&di[idi[i]], &inverseDi[2 * i]);
		}
		else
		{
			inverseDi[2 * i] = 1.0 / di[idi[i]];
			inverseDi[2 * i + 1] = 0;
		}
	}

	MultiplyRarefiedMatrixOnVector(ig, jg, ggl, di, ijg, idi, result, temp, blockSize);
	s.resize(2 * blockSize);
	y.resize(2 * blockSize);
	SubtractVectors(rightPart, temp, r);
	double residual = Norm(r, blockSize);
	double norm0 = Norm(rightPart, blockSize);
	residual /= norm0;
	double residualSecond = residual;
	OutputIterationsAndResidual(residual      , 0, "../resources/output/residual1.txt");
	OutputIterationsAndResidual(residualSecond, 0, "../resources/output/residual2.txt");

#if DIAGONAL_FACTORIZATION
	MultDiOnVect(inverseDi, r, z, blockSize);
#else
	SLAE_Forward_Complex (LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, r, z, blockSize);
	SLAE_Backward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, z, z, blockSize);
#endif

/*	p = z;
	s = r;
	y = result;*/
	CopyVector(&z[0], &p[0], blockSize);
	CopyVector(&r[0], &s[0], blockSize);
	CopyVector(&result[0], &y[0], blockSize);
	int flag = 0;
	int iteration = 0;
	const double timeStart = omp_get_wtime();
	while (residualSecond > epsilon && iteration < maxiter)
	{
		ComplexScalarConjugateProduct(r, z, complexNumber1, blockSize);
		MultiplyRarefiedMatrixOnVector(ig, jg, ggl, di, ijg, idi, p, Ap, blockSize);
		ComplexScalarConjugateProduct(Ap, p, complexNumber2, blockSize);
		DivideComplexNumbers(complexNumber1, complexNumber2, alfa);

		ComplexMultiplyVectorScalar(p, alfa, temp, blockSize);
		SummVectors(result, temp, result);//xj+1

		ComplexMultiplyVectorScalar(Ap, alfa, temp, blockSize);
		SubtractVectors(r, temp, r);//rj+1
		//etta
		SubtractVectors(r, s, temp_spline);
		etta = -RealScalarProduct(temp_spline, s) / RealScalarProduct(temp_spline, temp_spline);
		if (etta > 1.0)
		{
			flag = 1;
		/*	y = result;
			s = r;*/
			CopyVector(&result[0], &y[0], blockSize);
			CopyVector(&r[0], &s[0], blockSize);
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

#if DIAGONAL_FACTORIZATION
		MultDiOnVect(inverseDi, r, z, blockSize);//zj+1*/
#else
		SLAE_Forward_Complex (LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, r, z, blockSize);
		SLAE_Backward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, z, z, blockSize);
#endif

		ComplexScalarConjugateProduct(r, z, complexNumber2, blockSize);//beta
		DivideComplexNumbers(complexNumber2, complexNumber1, beta);

		ComplexMultiplyVectorScalar(p, beta, temp, blockSize);
		SummVectors(z, temp, p);

		residual        = Norm(r, blockSize);
		residual       /= norm0;
		residualSecond  = Norm(s, blockSize);
		residualSecond /= norm0;
		OutputIterationsAndResidual(residual      , iteration + 1, "../resources/output/residual1.txt");
		OutputIterationsAndResidual(residualSecond, iteration + 1, "../resources/output/residual2.txt");
		printf("residual = %le\r", residual);
		iteration++;
	}
	if (flag)
	{
		/*y = result;*/
		CopyVector(&result[0], &y[0], blockSize);
	}
	MultiplyRarefiedMatrixOnVector(ig, jg, ggl, di, ijg, idi, result, temp, blockSize);
	SubtractVectors(temp, rightPart, temp);
	residual  = Norm(temp, blockSize);
	residual /= norm0;
	OutputIterationsAndResidual(residual      , iteration + 1, "../resources/output/residual1.txt");
	OutputIterationsAndResidual(residualSecond, iteration + 1, "../resources/output/residual2.txt");
	printf("end residual\t=\t%le\n", residual);
	printf("iterations\t=\t%d\n", iteration);
	const double timeEnd = omp_get_wtime();

	printf("time spent\t=\t%.6lf sec\n", timeEnd - timeStart);
	FILE *fp;
	fopen_s(&fp, timeMeasurements, "a");
	fprintf(fp, "%.6lf\n", timeEnd - timeStart);
	fclose(fp);
}
