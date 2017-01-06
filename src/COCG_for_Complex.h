#include "MVFunctions.h"
#include "IO.h"

void COCG(
	vector < int > &ig,
	vector < int > &jg,
	vector < double > &ggl,
	vector < double > &di,
	vector < int > &ijg,
	vector < int > &idi,
	vector < double > &right_part,
	int nb,
	vector < double > &result,
	double eps,
	int MaxIter
	);