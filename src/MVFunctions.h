#pragma once
#include <stdio.h>

double RealScalarProduct            (double *& vec1, double *& vec2, const int size);
void   ComplexScalarConjugateProduct(double *& x, double *& y, double * result, const int blockSize);
void   RealMultiplyArrayScalar      (double *& vec, const double   scalar, double *& result, const int size);
void   ComplexMultiplyArrayScalar   (double *& vec, const double * skalar, double *& result, const int blockSize);
void   MultiplyBlock                (const double *x, double *dest, const double *a, const int size);

void   MultiplyComplexNumbers(const double * x, const double * y, double * result);
void   DivideComplexNumbers  (const double * x, const double * y, double * result);
void   SubtractComplexNumbers(const double * x, const double * y, double * result);

void   SubtractArrays (double *& x, double *& y, double *& result, const int size);
void   SummArrays     (double *& x, double *& y, double *& result, const int size);
void   CopyDoubleArray(double *& source, double *& dest, const int size);
void   CopyIntArray   (int    *& source, int    *& dest, const int size);

void   DiagonalPreconditioning(const double * x, double *result);

double Norm(double * x, const int blockSize);

void   MultDiOnVect(double *& di, double *& vec, double *& result, const int blockSize);

void   MultiplyRarefiedMatrixOnVector(
	        int    *& ig
	,       int    *& jg
	,       double *& gg
	,       double *& di
	,       int    *& ijg
	,       int    *& idi
	,       double *& x
	,       double *& y
	, const int       blockSize
	);
