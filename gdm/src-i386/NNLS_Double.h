//
// NNLS_Double.h
//
#ifndef __NNLS_DOUBLE_H__
#define __NNLS_DOUBLE_H__

#if defined _M_X64


double *WeightedNNLSRegression( char *lpTmpFile, 
							    double *pEnvDataMatrix, long long nRows, long long nCols, 
								double *pRespVector, double *pDeviance, double *pWeights);

double CalcGDMDevianceDouble( double *pY, double *pU, double *pW, long long nLen );

double GetWeightedNULLDeviance( long long nRows, double *pRespVector, double *pWeights );

double *CopyEnvMatrixDouble( double *pMatrix, long long nRows, long long nCols );

double *nnlsFITDouble( double *pEnvDataMatrix, long long nRows, long long nCols, double *pRespVector, double *pWeights );


#elif defined _WIN32


double *WeightedNNLSRegression( char *lpTmpFile, 
							    double *pEnvDataMatrix, int nRows, int nCols, 
								double *pRespVector, double *pDeviance, double *pWeights);


double CalcGDMDevianceDouble( double *pY, double *pU, double *pW, int nLen );

double GetWeightedNULLDeviance( int nRows, double *pRespVector, double *pWeights );

double GetTestWeightedNULLDeviance( double *pMatrix, int nCols, int nRows, double *pRespVector, double *pWeights );

double *CopyEnvMatrixDouble( double *pMatrix, int nRows, int nCols );

double *nnlsFITDouble( double *pEnvDataMatrix, int nRows, int nCols, double *pRespVector, double *pWeights );

#else // Linux

double *WeightedNNLSRegression( char *lpTmpFile, 
							    double *pEnvDataMatrix, int nRows, int nCols, 
								double *pRespVector, double *pDeviance, double *pWeights);


double CalcGDMDevianceDouble( double *pY, double *pU, double *pW, int nLen );

double GetWeightedNULLDeviance( int nRows, double *pRespVector, double *pWeights );

double GetTestWeightedNULLDeviance( double *pMatrix, int nCols, int nRows, double *pRespVector, double *pWeights );

double *CopyEnvMatrixDouble( double *pMatrix, int nRows, int nCols );

double *nnlsFITDouble( double *pEnvDataMatrix, int nRows, int nCols, double *pRespVector, double *pWeights );


#endif

#endif  // __NNLS_DOUBLE_H__