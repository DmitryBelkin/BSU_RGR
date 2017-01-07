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
	,       vector < int    > &jg
	,       vector < double > &gg
	,       vector < double > &di
	,       vector < int    > &ijg
	,       vector < int    > &idi
	,       vector < double > &rightPart
	, const int                blockSize
	,       vector < double > &result
	, const double             epsilon
	, const int                maxiter
	)
{
#if DIAGONAL_FACTORIZATION
	vector < double > di_1;
	di_1.resize(2 * blockSize);
	for (int i = 0; i < blockSize; ++i)
	{
		int size = idi[i + 1] - idi[i];
		if (size == 2) // если комплексное число
		{
			DiagonalPreconditioning(&di[idi[i]], &di_1[2 * i]);
		}
		else
		{
			di_1[2 * i    ] = 1.0 / di[idi[i]];
			di_1[2 * i + 1] = 0.0;
		}
	}
#else
	vector < double > LLT_gg, LLT_di;
	vector < int    > LLT_ig, LLT_jg, LLT_ijg, LLT_idi;
	LLT_Factorization(ig, jg, ijg, idi, gg, di, LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_gg, LLT_di, blockSize);
#endif

	vector < double > temp;
	vector < double > p, z, r;
	vector < double > Ap;
	double alfa[2], beta[2];
	double complexNumber1[2], complexNumber2[2];
	vector < double > y, s, test;
	p   .resize(2 * blockSize);
	z   .resize(2 * blockSize);
	r   .resize(2 * blockSize);
	Ap  .resize(2 * blockSize);
	temp.resize(2 * blockSize);
	test.resize(2 * blockSize);
	for (int i = 0; i < 2 * blockSize; ++i) { result[i] = 0.0; }

	MultiplyRarefiedMatrixOnVector(ig, jg, gg, di, ijg, idi, result, temp, blockSize);
	s.resize(2 * blockSize);
	y.resize(2 * blockSize);
	SubtractVectors(rightPart, temp, r);
	double residual = Norm(r, blockSize);
	double norm0    = Norm(rightPart, blockSize);
	residual       /= norm0;
	double residualSecond = residual;
	OutputIterationsAndResidual(residual      , 0, "../resources/output/residual1.txt");
	OutputIterationsAndResidual(residualSecond, 0, "../resources/output/residual2.txt");

#if DIAGONAL_FACTORIZATION
	MultDiOnVect(di_1, r, z, blockSize);
#else
	SLAE_Forward_Complex (LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_gg, LLT_di, r, z, blockSize);
	SLAE_Backward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_gg, LLT_di, z, z, blockSize);
#endif

	p = z;
	s = r;
	y = result;

	double etta;
	int    flag = 0;
	int    iteration = 0;

	const double timeStart = omp_get_wtime();

	while (residualSecond > epsilon && iteration < maxiter)
	{
		ComplexScalarConjugateProduct(r, z, complexNumber1, blockSize);
		MultiplyRarefiedMatrixOnVector(ig, jg, gg, di, ijg, idi, p, Ap, blockSize);
		ComplexScalarConjugateProduct(Ap, p, complexNumber2, blockSize);
		DivideComplexNumbers(complexNumber1, complexNumber2, alfa); // alpha_j

		ComplexMultiplyVectorScalar(p, alfa, temp, blockSize);      // alpha_j * p_j
		SummVectors(result, temp, result);                          // x_j+1

		ComplexMultiplyVectorScalar(Ap, alfa, temp, blockSize);     // alpha_j * A * p_j
		SubtractVectors(r, temp, r);                                // r_j+1

		SubtractVectors(r, s, test);
		etta = -RealScalarProduct(test, s) / RealScalarProduct(test, test);
		if (etta > 1.0)
		{
			flag = 1;
			y = result;
			s = r;
		}
		else if (etta > 0)
		{
			flag = 1;
			SubtractVectors(result, s, test);
			RealMultiplyVectorScalar(test, etta, test);
			SummVectors(y, test, y);

			SubtractVectors(r, s, test);
			RealMultiplyVectorScalar(test, etta, test);
			SummVectors(s, test, s);
		}

#if DIAGONAL_FACTORIZATION
		MultDiOnVect(di_1, r, z, blockSize);//zj+1*/
#else
		SLAE_Forward_Complex (LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_gg, LLT_di, r, z, blockSize);
		SLAE_Backward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_gg, LLT_di, z, z, blockSize);
#endif

		ComplexScalarConjugateProduct(r, z, complexNumber2, blockSize);
		DivideComplexNumbers(complexNumber2, complexNumber1, beta);     // beta_j

		ComplexMultiplyVectorScalar(p, beta, temp, blockSize);          // beta_j * p_j
		SummVectors(z, temp, p);                                        // p_j+1

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
		y = result;
	}
	MultiplyRarefiedMatrixOnVector(ig, jg, gg, di, ijg, idi, result, temp, blockSize);
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
