#include <vector>
using namespace std;

double RealScalarProduct            (vector < double > &vec1, vector < double > &vec2);
void   ComplexScalarConjugateProduct(vector < double > &x, vector < double > &y,double  *result,int nb);
void   RealMultiplyVectorScalar     (vector < double > &vec, double  scalar, vector < double > &result);
void   ComplexMultiplyVectorScalar  (vector < double > &vec, double *scalar, vector < double > &result, int nb);
void   MultiplyBlock                (double *x, double *y, double *a, int size);

void   MultiplyComplexNumbers(double *x, double *y, double *result);
void   DivideComplexNumbers  (double *x, double *y, double *result);
void   SubtractComplexNumbers(double *x, double *y, double *result);
void   CopyVector            (double *x, double *y, int size);

void   SubtractVectors(vector < double > &x, vector < double > &y, vector < double > &result);
void   SummVectors    (vector < double > &x, vector < double > &y, vector < double > &result);

void   DiagonalPreconditioning(double * x, double *result);

double Norm(vector < double > x, int nb);

void   MultDiOnVect(vector < double > &di, vector < double > &vec, vector < double > &res, int Nb);
void   MultVMatrixOnVector(vector < vector < double > > &V, vector < double > vec, vector < double > &res, int Nb, int m);

void   MultiplyRarefiedMatrixOnVector(
	vector < int    > &ig ,
	vector < int    > &jg ,
	vector < double > &ggl,
	vector < double > &di ,
	vector < int    > &ijg,
	vector < int    > &idi,
	vector < double > &x  ,
	vector < double > &y  ,
	int nb
	);
