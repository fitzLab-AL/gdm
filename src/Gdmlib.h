//
// Gdm4Rlib.h
//
#ifndef	__GDMLIB_H__
#define	__GDMLIB_H__

#include "stdafx.h"

extern "C" 
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These are specifically called from the R-Stats package
//
void GDM_FitFromTable( char **wspath, double *pData, 
                       int *pDoGeo, int *pPreds, 
 				       int *pRows, int *pCols, 
					   int *pSplines, double *pQuantiles,
					   double *pGDMDev, double *pNullDev, double *pExpDev, 
					   double *pIntercept, double *pCoeffs,
					   double *pY, double *pX, double *pE );



void GDM_PredictFromTable(double *pData, 
                          int *pDoGeo, int *pPreds, int *pRows, 
				          double *pQuantiles,  int *pSplines, double *pCoeffs,
				          double *pX);



void GetPredictorPlotData( double *pPredData, int *pLength,
					       double *pCoefficients,
					       double *pQuantiles,
					       int *pSplines );


void GDM_TransformFromTable( int *pRows, int *pCols, 
                             int *pDoGeo, int *pPreds, 
                             int *pSplines, double *pQuants, double *pCoeffs,
	                         double *pInData,
	                         double *pOutData );


void GetWordSize(int *p1);


//
// Contruct a regression matrix for the predictors
//
double *ConstructMatrix(int nDoGeo, double *pData, double *pQuants, int nPreds, int *pSplines, int nRows);


}

//
// define the constant data column indices
//
#define COL_RESPONSE 0
#define COL_WEIGHTS  1
#define COL_SITE1_X0 2
#define COL_SITE1_Y0 3
#define COL_SITE2_X1 4
#define COL_SITE2_Y1 5
#define LEADING_COLS 6

//
// Local support functions
//
#if defined _M_X64


//
// Note that pData represents a column major matrix
//
void GDM_PredictFromTable(double *pData, 
		                  int *pDoGeo, int *pPreds, int *pRows, 
					      double *pQuantiles, int *pSplines, double *pCoeffs,
					      double *pX);

//
// Populate pPredData as a transformed GDM predictor to plot in R-Stats
//
void GetPredictorPlotData( double *pPredData, int *pLength,
						   double *pCoefficients,
						   double *pQuantiles,
						   int *pSplines );

//
// Predict dissimilarity from spline data and coefficients
//
double CalcDissimilarity( double *pData, double *pCoeffs, long long nRows, int nCols, long long nIndex );


//
// calculate total number of splines
//
int GetTotalSplineCount(int *pSplines, int nPreds);


//
// Contruct a regression matrix for predictors
//
double *ConstructMatrix(int nDoGeo, double *pData, double *pQuants, int nPreds, int *pSplines, long long nRows);


//
// Calculate the GDM transform
//
double CalculateGDMTransform(double dValue, int nSplines, double *pQuants, double *pCoeffs);


//
// Calculate the I-Spline value for dVal given quantiles q1, q2, q3
//
double DoSplineCalc( double dVal, double q1, double q2, double q3 );


//
// Display the contents of the quantile vector
//
void ShowQuantiles(double *pQuants, int nPreds, int *pSplines);


//
// Write the Predictor Matrix to a comma delimited file
//
void DebugPredMatrix(char *pPath, double *pPredData, long long nRows, int nPreds, int *pSplines, int nCols);

#elif defined _WIN32


//
// Predict dissimilarity from spline data and coefficients
//
double CalcDissimilarity( double *pData, double *pCoeffs, int nRows, int nCols, int nIndex );


//
// calculate total number of splines
//
int GetTotalSplineCount(int *pSplines, int nPreds);


//
// Contruct a regression matrix for the predictors
//
double *ConstructMatrix(int nDoGeo, double *pData, double *pQuants, int nPreds, int *pSplines, int nRows);


//
// Calculate the GDM transform
//
double CalculateGDMTransform(double dValue, int nSplines, double *pQuants, double *pCoeffs);


//
// Calculate the I-Spline value for dVal given quantiles q1, q2, q3
//
double DoSplineCalc( double dVal, double q1, double q2, double q3 );



//
// Display the contents of the quantile vector
//
void ShowQuantiles(double *pQuants, int nPreds, int *pSplines);


//
// Write the Predictor Matrix to a comma delimited file
//
void DebugPredMatrix(char *pPath, double *pPredData, int nRows, int nPreds, int *pSplines, int nCols);

#else // Linux

//
// Predict dissimilarity from spline data and coefficients
//
double CalcDissimilarity( double *pData, double *pCoeffs, int nRows, int nCols, int nIndex );


//
// calculate total number of splines
//
int GetTotalSplineCount(int *pSplines, int nPreds);


//
// Contruct a regression matrix for the predictors
//
double *ConstructMatrix(int nDoGeo, double *pData, double *pQuants, int nPreds, int *pSplines, int nRows);


//
// Calculate the GDM transform
//
double CalculateGDMTransform(double dValue, int nSplines, double *pQuants, double *pCoeffs);


//
// Calculate the I-Spline value for dVal given quantiles q1, q2, q3
//
double DoSplineCalc( double dVal, double q1, double q2, double q3 );



//
// Display the contents of the quantile vector
//
void ShowQuantiles(double *pQuants, int nPreds, int *pSplines);


//
// Write the Predictor Matrix to a comma delimited file
//
void DebugPredMatrix(char *pPath, double *pPredData, int nRows, int nPreds, int *pSplines, int nCols);


#endif


#endif // __GDMBINLIB_H__
