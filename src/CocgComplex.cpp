#include "CocgComplex.h"
#include "LLT.h"
#include <mkl.h>
#include <omp.h>

#define DIAGONAL_FACTORIZATION 1
#define LLT_FACTORIZATION !DIAGONAL_FACTORIZATION

const char * timeMeasurements = "../resources/output/timeMeasurements.txt";

void OutputIterationsAndResidual(const double residual1, const int iteration, const char *filename)
{
	FILE *fp;
	fopen_s(&fp, filename, "a");
	fprintf(fp, "%d\t%.15le\n", iteration, residual1);
	fclose(fp);
}

void CocgComplex(
	        int    *& ig
	,       int    *& jg
	,       double *& gg
	,       double *& di
	,       int    *& ijg
	,       int    *& idi
	,       double *& rightPart
	, const int       blockSize
	,       double *& result
	, const double    epsilon
	, const int       maxiter
	)
{
#if DIAGONAL_FACTORIZATION
	double * di_1 = NULL;
	di_1 = new double[2 * blockSize];
	for (int i = 0; i < blockSize; ++i)
	{
		int size = idi[i + 1] - idi[i];
		if (size == 2)
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
	double * LLT_gg  = NULL;
	double * LLT_di  = NULL;
	int    * LLT_ig  = NULL;
	int    * LLT_jg  = NULL;
	int    * LLT_ijg = NULL;
	int    * LLT_idi = NULL;
	LltFactorization(ig, jg, ijg, idi, gg, di, LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_gg, LLT_di, blockSize);
#endif

	double * temp = NULL;
	double * p    = NULL;
	double * z    = NULL;
	double * r    = NULL;
	double * Ap   = NULL;
	double * y    = NULL;
	double * s    = NULL;
	double * test = NULL;
	double alfa[2], beta[2];
	double complexNumber1[2], complexNumber2[2];
	p    = new double[2 * blockSize];
	z    = new double[2 * blockSize];
	r    = new double[2 * blockSize];
	Ap   = new double[2 * blockSize];
	temp = new double[2 * blockSize];
	test = new double[2 * blockSize];
	for (int i = 0; i < 2 * blockSize; ++i) { result[i] = 0.0; }

	MultiplyRarefiedMatrixOnVector(ig, jg, gg, di, ijg, idi, result, temp, blockSize);
	s = new double[2 * blockSize];
	y = new double[2 * blockSize];
	SubtractArrays(rightPart, temp, r, 2 * blockSize);
	const double normRightPart = Norm(rightPart, blockSize);
	double residual1 = Norm(r, blockSize) / normRightPart;
	double residual2 = residual1;
	//OutputIterationsAndResidual(residual1, 0, "../resources/output/residual1.txt");
	//OutputIterationsAndResidual(residual2, 0, "../resources/output/residual2.txt");

#if DIAGONAL_FACTORIZATION
	MultDiOnVect(di_1, r, z, blockSize);
#else
	ForwardSlae (LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_gg, LLT_di, r, z, blockSize);
	BackwardSlae(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_gg, LLT_di, z, z, blockSize);
#endif

	//CopyDoubleArray(z     , p, 2 * blockSize);
	//CopyDoubleArray(r     , s, 2 * blockSize);
	//CopyDoubleArray(result, y, 2 * blockSize);
	cblas_dcopy(2 * blockSize, z     , 1, p, 1);
	cblas_dcopy(2 * blockSize, r     , 1, s, 1);
	cblas_dcopy(2 * blockSize, result, 1, y, 1);

	double etta;
	int    flag      = 0;
	int    iteration = 0;

	const double timeStart = omp_get_wtime();

	while (residual2 > epsilon && iteration < maxiter)
	{
		ComplexScalarConjugateProduct (r, z, complexNumber1, blockSize);            // ( _(r_j), z_j )
		MultiplyRarefiedMatrixOnVector(ig, jg, gg, di, ijg, idi, p, Ap, blockSize); //     A * p_j
		ComplexScalarConjugateProduct(Ap, p, complexNumber2, blockSize);            // ( _(A * p_j), p_j )
		DivideComplexNumbers(complexNumber1, complexNumber2, alfa);                 // alpha_j

		ComplexMultiplyArrayScalar(p, alfa, temp, blockSize);                       // alpha_j * p_j
		SummArrays(result, temp, result, 2 * blockSize);                            // x_j+1

		ComplexMultiplyArrayScalar(Ap, alfa, temp, blockSize);                      // alpha_j * A * p_j
		SubtractArrays(r, temp, r, 2 * blockSize);                                  // r_j+1

		SubtractArrays(r, s, test, 2 * blockSize);
		//etta = -RealScalarProduct(test, s, 2 * blockSize) / RealScalarProduct(test, test, 2 * blockSize);
		etta = -cblas_ddot(2 * blockSize, test, 1, s, 1) / cblas_ddot(2 * blockSize, test, 1, test, 1);
		if (etta > 1.0)
		{
			flag = 1;
			//CopyDoubleArray(result, y, 2 * blockSize);
			//CopyDoubleArray(r     , s, 2 * blockSize);
			cblas_dcopy(2 * blockSize, result, 1, y, 1);
			cblas_dcopy(2 * blockSize, r     , 1, s, 1);
		}
		else if (etta > 0)
		{
			flag = 1;
			SubtractArrays(result, s, test, 2 * blockSize);
			//RealMultiplyArrayScalar(test, etta, test, 2 * blockSize);
			cblas_dscal (2 * blockSize, etta, test, 1);
			SummArrays(y, test, y, 2 * blockSize);

			SubtractArrays(r, s, test, 2 * blockSize);
			//RealMultiplyArrayScalar(test, etta, test, 2 * blockSize);
			cblas_dscal (2 * blockSize, etta, test, 1);
			SummArrays(s, test, s, 2 * blockSize);
		}

#if DIAGONAL_FACTORIZATION
		MultDiOnVect(di_1, r, z, blockSize);                                             // z_j+1
#else
		ForwardSlae (LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_gg, LLT_di, r, z, blockSize);
		BackwardSlae(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_gg, LLT_di, z, z, blockSize); // z_j+1
#endif

		ComplexScalarConjugateProduct(r, z, complexNumber2, blockSize);
		DivideComplexNumbers(complexNumber2, complexNumber1, beta);     // beta_j

		ComplexMultiplyArrayScalar(p, beta, temp, blockSize);           // beta_j * p_j
		SummArrays(z, temp, p, 2 * blockSize);                          // p_j+1

		residual1  = Norm(r, blockSize);
		residual1 /= normRightPart;
		residual2  = Norm(s, blockSize);
		residual2 /= normRightPart;
		//OutputIterationsAndResidual(residual1, iteration + 1, "../resources/output/residual1.txt");
		//OutputIterationsAndResidual(residual2, iteration + 1, "../resources/output/residual2.txt");
		//printf("residual\t=\t%le\r", residual1);
		iteration++;
	}

	if (flag)
	{
		//CopyDoubleArray(result, y, 2 * blockSize);
		cblas_dcopy(2 * blockSize, result, 1, y, 1);
	}
	MultiplyRarefiedMatrixOnVector(ig, jg, gg, di, ijg, idi, result, temp, blockSize);
	SubtractArrays(temp, rightPart, temp, 2 * blockSize);
	residual1  = Norm(temp, blockSize);
	residual1 /= normRightPart;
	//OutputIterationsAndResidual(residual1, iteration + 1, "../resources/output/residual1.txt");
	//OutputIterationsAndResidual(residual2, iteration + 1, "../resources/output/residual2.txt");
	//printf("final residual\t=\t%le\n", residual1);
	//printf("iterations\t=\t%d\n", iteration);
	const double timeEnd = omp_get_wtime();

	printf("time spent\t=\t%.6lf sec\n", timeEnd - timeStart);
	FILE *fp;
	fopen_s(&fp, timeMeasurements, "a");
	fprintf(fp, "%.6lf\n", timeEnd - timeStart);
	fclose(fp);

	delete [] temp;
	delete [] p;
	delete [] z;
	delete [] r;
	delete [] Ap;
	delete [] y;
	delete [] s;
	delete [] test;

#if DIAGONAL_FACTORIZATION
	delete [] di_1;
#else
	delete [] LLT_gg;
	delete [] LLT_di;
	delete [] LLT_ig;
	delete [] LLT_jg;
	delete [] LLT_ijg;
	delete [] LLT_idi;
#endif
}
