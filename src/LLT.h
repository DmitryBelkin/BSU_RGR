#pragma once
#include <stdio.h>

void LltFactorization(
	        int    *& ig
	,       int    *& jg
	,       int    *& ijg
	,       int    *& idi
	,       double *& gg
	,       double *& di
	,       int    *& LLT_ig
	,       int    *& LLT_jg
	,       int    *& LLT_ijg
	,       int    *& LLT_idi
	,       double *& LLT_gg
	,       double *& LLT_di
	, const int       blockSize
	);

void ForwardSlae(
	        int    *& LLT_ig
	,       int    *& LLT_jg
	,       int    *& LLT_ijg
	,       int    *& LLT_idi
	,       double *& LLT_gg
	,       double *& LLT_di
	,       double *& rightPart
	,       double *& result
	, const int       blockSize
	);

void BackwardSlae(
	        int    *& LLT_ig
	,       int    *& LLT_jg
	,       int    *& LLT_ijg
	,       int    *& LLT_idi
	,       double *& LLT_gg
	,       double *& LLT_di
	,       double *& rightPart
	,       double *& result
	, const int       blockSize
	);

void SqrtComplex(double *ab, double*xy);
