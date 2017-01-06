#include "IO.h"
#include "COCG_for_Complex.h"

void main()
{
	vector < double > Cdi, Cggl, Cright_part, Cresult;
	vector < int    > Cidi, Cig, Cijg, Cjg;
	vector < double > VectorForCheckMatrixVectorSLAEMultiplication;
	int    MaxIter, FullTaskSize, BlockSize;
	double eps;
	
	InputSparseComplexBlockSLAE("../resources/IN_data", Cdi, Cggl, Cig, Cjg,
		Cidi, Cijg, Cright_part, FullTaskSize, BlockSize, eps,
		MaxIter, VectorForCheckMatrixVectorSLAEMultiplication);

	Cresult.resize(FullTaskSize);
	COCG(Cig, Cjg, Cggl, Cdi, 
		Cijg, Cidi, Cright_part, 
		BlockSize, Cresult,eps, MaxIter);

	getchar();
}