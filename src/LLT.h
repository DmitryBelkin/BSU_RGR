#include <vector>
using namespace std;

void LLT_Factorization(
	        vector < int    > &ig
	,       vector < int    > &jg
	,       vector < int    > &ijg
	,       vector < int    > &idi
	,       vector < double > &gg
	,       vector < double > &di
	,       vector < int    > &LLT_ig
	,       vector < int    > &LLT_jg
	,       vector < int    > &LLT_ijg
	,       vector < int    > &LLT_idi
	,       vector < double > &LLT_gg
	,       vector < double > &LLT_di
	, const int               blockSize
	);

void SLAE_Forward_Complex(
	        vector < int    > &LLT_ig
	,       vector < int    > &LLT_jg
	,       vector < int    > &LLT_ijg
	,       vector < int    > &LLT_idi
	,       vector < double > &LLT_gg
	,       vector < double > &LLT_di
	,       vector < double > &rightPart
	,       vector < double > &result
	, const int                blockSize
	);

void SLAE_Backward_Complex(
	        vector < int    > &LLT_ig
	,       vector < int    > &LLT_jg
	,       vector < int    > &LLT_ijg
	,       vector < int    > &LLT_idi
	,       vector < double > &LLT_gg
	,       vector < double > &LLT_di
	,       vector < double > &rightPart
	,       vector < double > &result
	, const int                blockSize
	);

void SqrtComplex(double *ab, double*xy);
