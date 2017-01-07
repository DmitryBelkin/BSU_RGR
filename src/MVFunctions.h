#include <vector>
using namespace std;

double RealScalarProduct            (const vector < double > &vec1, const vector < double > &vec2);
void   ComplexScalarConjugateProduct(const vector < double > &x, const vector < double > &y, double  *result, const int blockSize);
void   RealMultiplyVectorScalar     (const vector < double > &vec, const double  scalar, vector < double > &result);
void   ComplexMultiplyVectorScalar  (const vector < double > &vec, const double *scalar, vector < double > &result, const int blockSize);
void   MultiplyBlock(const double *x, double *dest, const double *a, const int size);

void   MultiplyComplexNumbers(const double *x, const double *y, double *result);
void   DivideComplexNumbers  (const double *x, const double *y, double *result);
void   SubtractComplexNumbers(const double *x, const double *y, double *result);

void   SubtractVectors(const vector < double > &x, const vector < double > &y, vector < double > &result);
void   SummVectors    (const vector < double > &x, const vector < double > &y, vector < double > &result);

void   DiagonalPreconditioning(const double * x, double *result);

double Norm(const vector < double > x, const int blockSize);

void   MultDiOnVect(const vector < double > &di, const vector < double > &vec, vector < double > &result, const int blockSize);

void   MultiplyRarefiedMatrixOnVector(
	  const vector < int    > &ig
	, const vector < int    > &jg
	,       vector < double > &gg
	,       vector < double > &di
	, const vector < int    > &ijg
	, const vector < int    > &idi
	,       vector < double > &x
	,       vector < double > &y
	, const int blockSize
	);
