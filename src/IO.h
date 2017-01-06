#include <vector>
using namespace std;

//INPUT//
void InputSparseComplexBlockSLAE(
	char* dirname,
	vector < double > &di,
	vector < double > &ggl,
	vector < int    > &ig,
	vector < int    > &jg,
	vector < int    > &idi,
	vector < int    > &ijg,
	vector < double > &righr_part,
	int               &FullTaskSize,
	int               &BlockSize,
	double            &eps,
	int               &MaxIter,
	vector < double > &VectorForCheckMatrixVectorSLAEMultiplication
	);

void inputSparceBlockMatrix_from_txt(
	char* filename,
	vector < double > &di,
	vector < double > &ggl,
	vector < int    > &ig,
	vector < int    > &jg,
	vector < int    > &idi,
	vector < int    > &ijg,
	vector < double > &righr_part,
	int               &FullTaskSize,
	int               &BlockSize,
	double            &eps,
	int               &MaxIter
	);

///////////////////////////////////////////////
//OUTPUT//
void OutputItersResidual(double Residual, int IterNumber, char *filename);
