//
// GdmTabLib
//
#include <string>
#include "stdafx.h"
#include "Gdmlib.h"
#include "NNLS_Double.h"

#include "Message.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include <unistd.h>
#include <Rcpp.h>
using namespace Rcpp;

//#ifdef _WIN32
//	#include <io.h>
//	#include <fcntl.h>
//#else
//	#include <unistd.h>
//#endif
//
//#include <sys/stat.h>


//
// return the number of bytes in a pointer to an int (32bit = 4 and 64bit = 8)
//
void GetWordSize(int *p1)
{
	*p1 = (int)(sizeof(int *));
}

#ifndef HAVE_STRLCAT
/*
 * '_cups_strlcat()' - Safely concatenate two strings.
 */

size_t                  /* O - Length of string */
strlcat(char       *dst,        /* O - Destination string */
              const char *src,      /* I - Source string */
          size_t     size)      /* I - Size of destination string buffer */
{
  size_t    srclen;         /* Length of source string */
  size_t    dstlen;         /* Length of destination string */


 /*
  * Figure out how much room is left...
  */

  dstlen = strlen(dst);
  size   -= dstlen + 1;

  if (!size)
    return (dstlen);        /* No room, return immediately... */

 /*
  * Figure out how much room is needed...
  */

  srclen = strlen(src);

 /*
  * Copy the appropriate amount...
  */

  if (srclen > size)
    srclen = size;

  memcpy(dst + dstlen, src, srclen);
  dst[dstlen + srclen] = '\0';

  return (dstlen + srclen);
}
#endif /* !HAVE_STRLCAT */

#ifndef HAVE_STRLCPY
/*
 * '_cups_strlcpy()' - Safely copy two strings.
 */

size_t                  /* O - Length of string */
strlcpy(char       *dst,        /* O - Destination string */
              const char *src,      /* I - Source string */
          size_t      size)     /* I - Size of destination string buffer */
{
  size_t    srclen;         /* Length of source string */


 /*
  * Figure out how much room is needed...
  */

  size --;

  srclen = strlen(src);

 /*
  * Copy the appropriate amount...
  */

  if (srclen > size)
    srclen = size;

  memcpy(dst, src, srclen);
  dst[srclen] = '\0';

  return (srclen);
}
#endif /* !HAVE_STRLCPY */

#if defined _M_X64

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Functions called from R -Stats ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Main GDM fitting functions called from R
//
//
// Note that pData represents a column major matrix
//

void GDM_FitFromTable(char **wspath,
				      double *pData,
	                  int *pDoGeo, int *pPreds,
				      int *pRows, int *pCols,
				      int *pSplines, double *pQuantiles,
				      double *pGDMDev, double *pNullDev, double *pExpDev,
				      double *pIntercept, double *pCoeffs,
				      double *pY, // observed
				      double *pX, // predicted
				      double *pE) // ecological dist
{
	//
	// Everything passed from R needs to be dereferenced...
	//
	int nDoGeo = *pDoGeo;
	int nPreds = *pPreds;
	long nRows = *pRows;
	//int nCols = *pCols;


	//
	// Get the response column
	//
	double *pResponse = &pData[COL_RESPONSE * nRows];
	if ( NULL == pResponse )
	{
		//Message("pResponse is NULL", "ERROR in GDM_FitFromTable");
		return;
	}


	//
	// Get the weights column
	//
	double *pWeights = &pData[COL_WEIGHTS * nRows];
	if ( NULL == pWeights )
	{
		//Message("pWeights is NULL", "ERROR in GDM_FitFromTable");
		return;
	}


	//
	// Get the predictor matrix
	//
	double *pPredData = ConstructMatrix(nDoGeo, pData, pQuantiles, nPreds, pSplines, (long long)nRows);
	if ( NULL == pPredData )
	{
		//Message("pPredData is NULL", "ERROR in GDM_FitFromTable");
		return;
	}
	// Sum the splines vector for the total number of splines
	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	// Message(nTotalSplines, "nTotalSplines");
	// pQuantiles will have a length of nPreds * nTotalSplines
	// ShowQuantiles(pQuantiles, nPreds, pSplines);
	//DebugPredMatrix("PredMatrix.csv", pPredData, nRows, nPreds, pSplines, nTotalSplines+1);


	//
	// Write a binary file image of the predictor matrix
	//file extension
	std::string bin = ".bin";
	char *cbin = new char[bin.length() + 1];
	strcpy(cbin, bin.c_str());

	char lpTmpFile[1024];
	strlcpy(lpTmpFile, *wspath, sizeof(lpTmpFile));
	strlcat(lpTmpFile, cbin, sizeof(lpTmpFile));
	//Rcpp::Rcout << lpTmpFile << std::endl;

	int h = _open( lpTmpFile, _O_BINARY | _O_CREAT | _O_TRUNC | _O_RDWR, S_IREAD | S_IWRITE );
	if ( h < 0 )
	{
		//Message("Cannot create binary image file", "ERROR in GDM_FitFromTable");
		if (pPredData) delete[] pPredData;
		return;
	}
	long nThis = 0;
	for ( int i=0; i<(nTotalSplines+1); i++ )
	{
		_write( h, &pPredData[nThis], nRows * sizeof( double ) );
		nThis += nRows;
	}
	_close( h );



	//
	// Do the matrix regression
	//
	double dGDMDeviance;
	double *pCoefficients = WeightedNNLSRegression( lpTmpFile,
		                                            pPredData,
											        nRows,
											        nTotalSplines+1,
											        pResponse,
											        &dGDMDeviance,
											        pWeights);

	//
	// remove the temporary matrix file if it exists (and it should!!) now that we don't need it
	//
	if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );


	//
	// create a NULL model and return the deviance
	//
	double dNullDeviance = GetWeightedNULLDeviance( nRows, pResponse, pWeights );


	//
	// calculate deviance explained as a percentage
	//
	double dDevianceExplained = ( 1.0 - ( dGDMDeviance / dNullDeviance ) ) * 100;


	//
	// display results
	//
	pGDMDev[0] = dGDMDeviance;
	pNullDev[0] = dNullDeviance;
	pExpDev[0] = dDevianceExplained;
	pIntercept[0] = pCoefficients[0];
	for ( int i=1; i<nTotalSplines+1; i++ )
	{
		pCoeffs[i-1] = pCoefficients[i];
	}


	//
	// copy to response data to a vector
	//
	for ( long long i=0; i<nRows; i++ )
	{
		pY[i] = pResponse[i];
		pE[i] = CalcDissimilarity( pPredData, pCoefficients, nRows, nTotalSplines+1, i );
		pX[i] = 1.0 - exp(-pE[i]);
	}


	//
	// Clean up
	//
	if (pPredData) delete[] pPredData;
	if ( pCoefficients ) delete[] pCoefficients;
}





//
// Main GDM predict function called from R
//
//
// Note that pData represents a column major matrix
//
void GDM_PredictFromTable(double *pData,
		                  int *pDoGeo, int *pPreds, int *pRows,
					      double *pQuantiles, int *pSplines, double *pCoeffs,
					      double *pX) // predicted
{
	//
	// Everything passed from R needs to be dereferenced...
	//
	int nDoGeo = *pDoGeo;
	int nPreds = *pPreds;
	long nRows  = *pRows;


	//
	// Get the predictor matrix
	//
	double *pPredData = ConstructMatrix(nDoGeo, pData, pQuantiles, nPreds, pSplines, (long long)nRows);
	if ( NULL == pPredData )
	{
		//Message("pPredData is NULL", "ERROR");
		return;
	}

	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	//Message(nTotalSplines, "nTotalSplines");
	// pQuantiles will have a length of nPreds * nTotalSplines
	//ShowQuantiles(pQuantiles, nPreds, pSplines);
	//DebugPredMatrix("PredMatrix.csv", pPredData, nRows, nPreds, pSplines, nTotalSplines+1);


	//
	// copy predicted data to pX
	//
	for ( long i=0; i<nRows; i++ )
	{
		pX[i] = 1.0 - exp(-CalcDissimilarity( pPredData, pCoeffs, nRows, nTotalSplines+1, i ));
	}


	//
	// Clean up
	//
	if (pPredData) delete[] pPredData;
}




//
// Populate pPredData as a transformed GDM predictor to plot in R-Stats
//
void GetPredictorPlotData( double *pPredData, int *pLength,
						   double *pCoefficients,
						   double *pQuantiles,
						   int *pSplines )
{
	int nLen  = *pLength;
	int nSplines = *pSplines;

	double dInc = fabs(pQuantiles[nSplines-1]-pQuantiles[0]) / nLen;
	double dVal = pQuantiles[0];
	for ( int i=0; i<nLen; i++ )
	{
		pPredData[i] = 0;
		for ( int j=0; j<nSplines; j++ )
		{
			if ( j==0 )
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j], pQuantiles[j], pQuantiles[j+1] );
			}

			else if ( j == nSplines-1 )
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j-1], pQuantiles[j], pQuantiles[j] );
			}

			else
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j-1], pQuantiles[j], pQuantiles[j+1] );
			}
		}
		dVal += dInc;
	}
}


//
//
// Transform a dataset using the coefficients and spline positions from a GDM model
//
//
void GDM_TransformFromTable(int *pRows, int *pCols,
	                        int *pDoGeo, int *pPreds,
							int *pSplines, double *pQuants, double *pCoeffs,
							double *pInData,
							double *pOutData)
{
	//
	// Everything passed from R needs to be dereferenced...
	//
	int nRows = *pRows;
	// Modified by DNL:
	//	int nCols = *pCols;
	int nDoGeo = *pDoGeo;
	int nPreds = *pPreds;

	int thisItem = 0;
	int thisQuant = 0;
	int thisCoeff = 0;

	if (nDoGeo)
	{
		int numSplines = pSplines[0];
		double *pQuantiles = &pQuants[0];
		double *pCoefficients = &pCoeffs[0];

		//
		// Setup a linear scalar to approximate the general slope of
		// the model curve to apply to the distance from the left edge.
		// Distances exceeding this are simply extrapolated by the linear scalar.
		//
		double dTranX = 0.0;
		double dTranY = 0.0;
		for ( int i=0; i<numSplines; i++ )
		{
			dTranY += pCoefficients[i];

			// get the maximum quantile value here
			if ( i == numSplines-1 )
			{
				dTranX = pQuantiles[i];
			}
		}
		double dScalar = dTranY / dTranX;

		//
		// now find the euclidean distance from the left edge and top edge and transform
		//
		// get min easting
		double *pXCol = &pInData[thisItem];
		double dMinX = pXCol[0];
		for ( int i=1; i<nRows; i++ )
		{
			if (pXCol[i] < dMinX) dMinX = pXCol[i];
		}

		// get min northing
		double *pYCol = &pInData[thisItem+nRows];
		double dMinY = pYCol[0];
		for ( int i=1; i<nRows; i++ )
		{
			if (pYCol[i] < dMinY) dMinY = pYCol[i];
		}


		// Calculate the transform for X
		for ( int i=0; i<nRows; i++ )
		{
			// get the data value
			double dDataVal = pInData[thisItem];

			// calculate the transform
			double dOutVal = dScalar * fabs(dDataVal-dMinX);

			// update output table
			pOutData[thisItem] = dOutVal;

			// update the data index
			++thisItem;
		}


		// Calculate the transform for Y
		for ( int i=0; i<nRows; i++ )
		{
			// get the data value
			double dDataVal = pInData[thisItem];

			// calculate the transform
			double dOutVal = dScalar * fabs(dDataVal-dMinY);

			// update output table
			pOutData[thisItem] = dOutVal;

			// update the data index
			++thisItem;
		}


		//
		// update indices and do environmental preds
		//
		thisQuant += numSplines;
		thisCoeff += numSplines;
		for (int pred=1; pred<nPreds; pred++)
		{
			int numSplines = pSplines[pred];
			double *pQuantiles = &pQuants[thisQuant];
			double *pCoefficients = &pCoeffs[thisCoeff];

			for ( int i=0; i<nRows; i++ )
			{
				// get the data value
				double dDataVal = pInData[thisItem];

				// calculate the transform
				double dOutVal = CalculateGDMTransform(dDataVal, numSplines, pQuantiles, pCoefficients);

				// update output table
				pOutData[thisItem] = dOutVal;

				// update the data index
				++thisItem;
			}

			//
			// update indices
			//
			thisQuant += numSplines;
			thisCoeff += numSplines;
		}

	}
	else
	{
		for (int pred=0; pred<nPreds; pred++)
		{
			int numSplines = pSplines[pred];
			double *pQuantiles = &pQuants[thisQuant];
			double *pCoefficients = &pCoeffs[thisCoeff];

			for ( int i=0; i<nRows; i++ )
			{
				// get the data value
				double dDataVal = pInData[thisItem];

				// calculate the transform
				double dOutVal = CalculateGDMTransform(dDataVal, numSplines, pQuantiles, pCoefficients);

				// update output table
				pOutData[thisItem] = dOutVal;

				// update the data index
				++thisItem;
			}

			//
			// update indices
			//
			thisQuant += numSplines;
			thisCoeff += numSplines;
		}
	}
}


//
// Calculate the GDM transform
//
double CalculateGDMTransform(double dValue, int nSplines, double *pQuants, double *pCoeffs)
{
	double dVal = 0.0;
	for ( int i=0; i<nSplines; i++ )
	{
		if (i == 0) // first spline
		{
			dVal += DoSplineCalc( dValue, pQuants[i], pQuants[i], pQuants[i+1] ) * pCoeffs[i];
		}

		else if (i == nSplines-1) // last spline
		{
			dVal += DoSplineCalc( dValue, pQuants[i-1], pQuants[i], pQuants[i] ) * pCoeffs[i];
		}

		else  // a middle spline
		{
			dVal += DoSplineCalc( dValue, pQuants[i-1], pQuants[i], pQuants[i+1] ) * pCoeffs[i];
		}
	}
	return(dVal);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Local Support Functions ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Predict dissimilarity from spline data and coefficients
//
double CalcDissimilarity( double *pData, double *pCoeffs, long long nRows, int nCols, long long nIndex )
{
	double dVal = 0.0;

	double *pTmp = &pData[nIndex];
	for ( int i=0; i<nCols; i++ )
	{
		dVal += *pTmp * pCoeffs[i];
		pTmp += nRows;
	}
	return(dVal);
}


//
// calculate total number of splines
//
int GetTotalSplineCount(int *pSplines, int nPreds)
{
	int nCount = 0;
	for ( int i=0; i<nPreds; i++ )
	{
		nCount += pSplines[i];
	}
	return(nCount);
}


//
// Construct a regression matrix for predictors
//
double *ConstructMatrix(int nDoGeo, double *pData, double *pQuants, int nPreds, int *pSplines, long long nRows)
{
	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	//Message(nTotalSplines, "nTotalSplines");

	//
	// Construct the predictor matrix with an extra column for the intercept
	//
	long long nSize = nRows * (nTotalSplines + 1);
	double *pPredData = new double [nSize];
	if ( NULL == pPredData )
	{
		//Message("Cannot allocate Predictor Data", "ERROR in ConstructMatrix");
		return(NULL);
	}
	for ( long long i=0; i<nSize; i++ ) pPredData[i] = 0.0;


	//
	// Initialise the intercept column
	//
	double *pLoc = &pPredData[0];
	for ( long long i=0; i<nRows; i++ ) pLoc[i] = 1.0;
	pLoc += nRows;   // start at the first column after the intercept


	//
	// Initialise for the geographic distance if we are using it
	//
	if (nDoGeo)
	{
		//
		// Initialise the matrix for the geographic followed by the environmental predictors
		//
		double dMin,dMid,dMax,dVal1,dVal2;
		int CurrentSpline = 0;

		for ( int pred=0; pred<nPreds; pred++ )
		{
			int d1Offset = ((LEADING_COLS+pred-1)*nRows);
			int d2Offset = ((LEADING_COLS+pred+nPreds-2)*nRows);

			int nThis = 0;
			for ( long long i=0; i<nRows; i++ )
			{
				if ( pred == 0 )
				{
					double distX = fabs( pData[(COL_SITE1_X0 * nRows)+i] - pData[(COL_SITE2_X1 * nRows)+i] );
					double distY = fabs( pData[(COL_SITE1_Y0 * nRows)+i] - pData[(COL_SITE2_Y1 * nRows)+i] );

					dVal1 = 0.0;
					dVal2 = sqrt( ( distX * distX ) + ( distY * distY ) );
				}
				else
				{
					dVal1 = pData[d1Offset+i];
					dVal2 = pData[d2Offset+i];
				}

				for ( int spl = 0; spl<pSplines[pred]; spl++ )
				{
					if ( spl == 0 )                              // first spline
					{
						dMin = pQuants[CurrentSpline+spl+0];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];
					}

					else if ( spl == (pSplines[pred]-1) )        // last spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+0];
					}

					else	                                     // a middle spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];
					}

					// calculate the spline values for each site
					double d1a = DoSplineCalc( dVal1, dMin, dMid, dMax );
					double d2a = DoSplineCalc( dVal2, dMin, dMid, dMax );

					// set the distance between sites for this spline
					pLoc[ (spl * nRows) + nThis ] = fabs( d2a - d1a );

				} // for ( int spl = 0; spl<pSplines[pred]; spl++ )

				// increment for the next row
				++nThis;

			} // for ( long long i=0; i<nRows; i++ )

			// increment pointer to the start of the next predictor's spline columns
			pLoc += pSplines[pred] * nRows;

			// increment the current quantile index for the next pred
			CurrentSpline += pSplines[pred];

		} // for ( int pred=0; pred<nPreds; pred++ )
	}


	else
	{
		//
		// Initialise the matrix for the environmental predictors only
		//
		double dMin,dMid,dMax;
		int CurrentSpline = 0;

		for ( int pred=0; pred<nPreds; pred++ )
		{
			int d1Offset = ((LEADING_COLS+pred)*nRows);
			int d2Offset = ((LEADING_COLS+nPreds+pred)*nRows);

			int nThis = 0;
			for ( long long i=0; i<nRows; i++ )
			{
				double dVal1 = pData[d1Offset+i];
				double dVal2 = pData[d2Offset+i];

				for ( int spl = 0; spl<pSplines[pred]; spl++ )
				{
					if ( spl == 0 )                              // first spline
					{
						dMin = pQuants[CurrentSpline+spl+0];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];

					}

					else if ( spl == (pSplines[pred]-1) )        // last spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+0];
					}

					else	                                     // a middle spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];
					}


					// calculate the spline values for each site
					double d1a = DoSplineCalc( dVal1, dMin, dMid, dMax );
					double d2a = DoSplineCalc( dVal2, dMin, dMid, dMax );

					// set the distance between sites for this spline
					pLoc[ (spl * nRows) + nThis ] = fabs( d2a - d1a );

				} // for ( int spl = 0; spl<pSplines[pred]; spl++ )

				// increment for the next row
				++nThis;

			} // for ( int i=0; i<nRows; i++ )

			// increment pointer to the start of the next predictor's spline columns
			pLoc += pSplines[pred] * nRows;

			// increment the current quantile index for the next pred
			CurrentSpline += pSplines[pred];

		} // for ( int pred=0; pred<nPreds; pred++ )

	} // environmental predictors only

	return(pPredData);
}




//
// Calculate the I-Spline value for dVal given quantiles q1, q2, q3
//
double DoSplineCalc( double dVal, double q1, double q2, double q3 )
{
	if ( dVal <= q1 ) return(0.0);

	else if ( dVal >= q3 ) return( 1.0 );

	else if ( ( q1 < dVal ) && ( dVal < q2 ) )
		return( ( ( ( dVal - q1 ) * ( dVal - q1 ) ) / ( ( q2 - q1 ) * ( q3 - q1 ) ) ) );

	else
		return( ( 1.0 - ( ( ( q3 - dVal ) * ( q3 - dVal ) ) / ( ( q3 - q2 ) * ( q3 - q1) ) ) ) );
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Local Debugging Functions /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Display the contents of the quantile vector
//
void ShowQuantiles(double *pQuants, int nPreds, int *pSplines)
{
	char buff[1024];
	double *pTmp = &pQuants[0];
	for ( int i=0; i<nPreds; i++ )
	{
		//sprintf( buff, "Quant %d: ", i+1 );
		//snprintf( buff, sizeof(buff), "Quant %d: ", i+1 );
		snprintf( buff, 10000, "Quant %d: ", i+1 );
		for ( int j=0; j<pSplines[i]; j++ )
		{

			// Modified by DNL:
			//sprintf( buff, "%s %lf ", buff, *pTmp);
			//sprintf( buff, "%s %f ", buff, *pTmp);
			//snprintf( buff, sizeof(buff), "%s %f ", buff, *pTmp);
			snprintf( buff, 10000, "%s %f ", buff, *pTmp);
			++pTmp;
		}
		//Message(buff);
	}
}




//
// Write the Predictor Matrix to a comma delimited file
//
void DebugPredMatrix(char *pPath, double *pPredData, long long nRows, int nPreds, int *pSplines, int nCols)
{
	FILE *fp = fopen(pPath, "w+t");

	// write the header
	fprintf( fp, "Intercept,");
	for ( int p=0; p<nPreds; p++ )
	{
		for ( int s=0; s<pSplines[p]; s++ )
		{
			fprintf( fp, "Pred%dSpline%d", p+1,s+1);
			if ( s < pSplines[p]-1 )
				fprintf( fp, "," );
		}
		if ( p < nPreds-1 )
			fprintf( fp, "," );
		else
			fprintf( fp, "\n" );

	}

	for ( long long i=0; i<nRows; i++ )
	{
		for ( int j=0; j<nCols; j++ )
		{
			// Modified by DNL:
			//fprintf( fp, "%lf", pPredData[(j*nRows)+i]);
			fprintf( fp, "%f", pPredData[(j*nRows)+i]);

			if ( j < nCols-1 )
				fprintf(fp, "," );
			else
				fprintf(fp, "\n" );
		}
	}

	if (fp) fclose(fp);
}




#elif defined _WIN32


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Functions called from R -Stats ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Main GDM fitting functions called from R
//
//
// Note that pData represents a column major matrix
//
void GDM_FitFromTable(char **wspath,
				      double *pData,
	                  int *pDoGeo, int *pPreds,
				      int *pRows, int *pCols,
				      int *pSplines, double *pQuantiles,
				      double *pGDMDev, double *pNullDev, double *pExpDev,
				      double *pIntercept, double *pCoeffs,
				      double *pY, // observed
				      double *pX, // predicted
				      double *pE) // ecological dist
{
	//
	// Everything passed from R needs to be dereferenced...
	//
	int nDoGeo = *pDoGeo;
	int nPreds = *pPreds;
	int nRows = *pRows;
	//int nCols = *pCols;


	//
	// Get the response column
	//
	double *pResponse = &pData[COL_RESPONSE * nRows];
	if ( NULL == pResponse )
	{
		//Message("pResponse is NULL", "ERROR in GDM_FitFromTable");
		return;
	}


	//
	// Get the weights column
	//
	double *pWeights = &pData[COL_WEIGHTS * nRows];
	if ( NULL == pWeights )
	{
		//Message("pWeights is NULL", "ERROR in GDM_FitFromTable");
		return;
	}


	//
	// Get the predictor matrix
	//
	double *pPredData = ConstructMatrix(nDoGeo, pData, pQuantiles, nPreds, pSplines, nRows);
	if ( NULL == pPredData )
	{
		//Message("pPredData is NULL", "ERROR in GDM_FitFromTable");
		return;
	}
	// Sum the splines vector for the total number of splines
	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	// Message(nTotalSplines, "nTotalSplines");
	// pQuantiles will have a length of nPreds * nTotalSplines
	// ShowQuantiles(pQuantiles, nPreds, pSplines);
	//DebugPredMatrix("PredMatrix.csv", pPredData, nRows, nPreds, pSplines, nTotalSplines+1);


	//
	// Write a binary file image of the predictor matrix
	//file extension
	std::string bin = ".bin";
	char *cbin = new char[bin.length() + 1];
	strcpy(cbin, bin.c_str());

	char lpTmpFile[1024];
	strlcpy(lpTmpFile, *wspath, sizeof(lpTmpFile));
	strlcat(lpTmpFile, cbin, sizeof(lpTmpFile));
	//Rcpp::Rcout << lpTmpFile << std::endl;

	int h = _open( lpTmpFile, _O_BINARY | _O_CREAT | _O_TRUNC | _O_RDWR, S_IREAD | S_IWRITE );
	if ( h < 0 )
	{
		//Message("Cannot create binary image file", "ERROR in GDM_FitFromTable");
		if (pPredData) delete[] pPredData;
		return;
	}
	int nSize = nRows * (nTotalSplines+1);
	_write( h, pPredData, nSize * sizeof( double ) );
	_close( h );



	//
	// Do the matrix regression
	//
	double dGDMDeviance;
	double *pCoefficients = WeightedNNLSRegression( lpTmpFile,
		                                            pPredData,
											        nRows,
											        nTotalSplines+1,
											        pResponse,
											        &dGDMDeviance,
											        pWeights);

	//
	// remove the temporary matrix file if it exists (and it should!!) now that we don't need it
	//
	if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );


	//
	// create a NULL model and return the deviance
	//
	double dNullDeviance = GetWeightedNULLDeviance( nRows, pResponse, pWeights );


	//
	// calculate deviance explained as a percentage
	//
	double dDevianceExplained = ( 1.0 - ( dGDMDeviance / dNullDeviance ) ) * 100;


	//
	// display results
	//
	pGDMDev[0] = dGDMDeviance;
	pNullDev[0] = dNullDeviance;
	pExpDev[0] = dDevianceExplained;
	pIntercept[0] = pCoefficients[0];
	for ( int i=1; i<nTotalSplines+1; i++ )
	{
		pCoeffs[i-1] = pCoefficients[i];
	}


	//
	// copy to response data to a vector
	//
	for ( int i=0; i<nRows; i++ )
	{
		pY[i] = pResponse[i];
		pE[i] = CalcDissimilarity( pPredData, pCoefficients, nRows, nTotalSplines+1, i );
		pX[i] = 1.0 - exp(-pE[i]);
	}


	//
	// Clean up
	//
	if (pPredData) delete[] pPredData;
	if ( pCoefficients ) delete[] pCoefficients;
}





//
// Main GDM predict function called from R
//
//
// Note that pData represents a column major matrix
//
void GDM_PredictFromTable(double *pData,
		                  int *pDoGeo, int *pPreds, int *pRows,
					      double *pQuantiles, int *pSplines, double *pCoeffs,
					      double *pX) // predicted
{
	//
	// Everything passed from R needs to be dereferenced...
	//
	int nDoGeo = *pDoGeo;
	int nPreds = *pPreds;
	int nRows  = *pRows;


	//
	// Get the predictor matrix
	//
	double *pPredData = ConstructMatrix(nDoGeo, pData, pQuantiles, nPreds, pSplines, nRows);
	if ( NULL == pPredData )
	{
		//Message("pPredData is NULL", "ERROR");
		return;
	}

	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	//Message(nTotalSplines, "nTotalSplines");
	// pQuantiles will have a length of nPreds * nTotalSplines
	//ShowQuantiles(pQuantiles, nPreds, pSplines);
	//DebugPredMatrix("PredMatrix.csv", pPredData, nRows, nPreds, pSplines, nTotalSplines+1);


	//
	// copy predicted data to pX
	//
	for ( int i=0; i<nRows; i++ )
	{
		pX[i] = 1.0 - exp(-CalcDissimilarity( pPredData, pCoeffs, nRows, nTotalSplines+1, i ));
	}


	//
	// Clean up
	//
	if (pPredData) delete[] pPredData;
}




//
// Populate pPredData as a transformed GDM predictor to plot in R-Stats
//
void GetPredictorPlotData( double *pPredData, int *pLength,
						   double *pCoefficients,
						   double *pQuantiles,
						   int *pSplines )
{
	int nLen  = *pLength;
	int nSplines = *pSplines;

	double dInc = fabs(pQuantiles[nSplines-1]-pQuantiles[0]) / nLen;
	double dVal = pQuantiles[0];
	for ( int i=0; i<nLen; i++ )
	{
		pPredData[i] = 0;
		for ( int j=0; j<nSplines; j++ )
		{
			if ( j==0 )
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j], pQuantiles[j], pQuantiles[j+1] );
			}

			else if ( j == nSplines-1 )
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j-1], pQuantiles[j], pQuantiles[j] );
			}

			else
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j-1], pQuantiles[j], pQuantiles[j+1] );
			}
		}
		dVal += dInc;
	}
}


//
//
// Transform a dataset using the coefficients and spline positions from a GDM model
//
//
void GDM_TransformFromTable(int *pRows, int *pCols,
	                        int *pDoGeo, int *pPreds,
							int *pSplines, double *pQuants, double *pCoeffs,
							double *pInData,
							double *pOutData)
{
	//
	// Everything passed from R needs to be dereferenced...
	//
	int nRows = *pRows;
	// Modified by DNL:
	//int nCols = *pCols;
	int nDoGeo = *pDoGeo;
	int nPreds = *pPreds;

	int thisItem = 0;
	int thisQuant = 0;
	int thisCoeff = 0;

	if (nDoGeo)
	{
		int numSplines = pSplines[0];
		double *pQuantiles = &pQuants[0];
		double *pCoefficients = &pCoeffs[0];

		//
		// Setup a linear scalar to approximate the general slope of
		// the model curve to apply to the distance from the left edge.
		// Distances exceeding this are simply extrapolated by the linear scalar.
		//
		double dTranX = 0.0;
		double dTranY = 0.0;
		for ( int i=0; i<numSplines; i++ )
		{
			dTranY += pCoefficients[i];

			// get the maximum quantile value here
			if ( i == numSplines-1 )
			{
				dTranX = pQuantiles[i];
			}
		}
		double dScalar = dTranY / dTranX;

		//
		// now find the euclidean distance from the left edge and top edge and transform
		//
		// get min easting
		double *pXCol = &pInData[thisItem];
		double dMinX = pXCol[0];
		for ( int i=1; i<nRows; i++ )
		{
			if (pXCol[i] < dMinX) dMinX = pXCol[i];
		}

		// get min northing
		double *pYCol = &pInData[thisItem+nRows];
		double dMinY = pYCol[0];
		for ( int i=1; i<nRows; i++ )
		{
			if (pYCol[i] < dMinY) dMinY = pYCol[i];
		}


		// Calculate the transform for X
		for ( int i=0; i<nRows; i++ )
		{
			// get the data value
			double dDataVal = pInData[thisItem];

			// calculate the transform
			double dOutVal = dScalar * fabs(dDataVal-dMinX);

			// update output table
			pOutData[thisItem] = dOutVal;

			// update the data index
			++thisItem;
		}


		// Calculate the transform for Y
		for ( int i=0; i<nRows; i++ )
		{
			// get the data value
			double dDataVal = pInData[thisItem];

			// calculate the transform
			double dOutVal = dScalar * fabs(dDataVal-dMinY);

			// update output table
			pOutData[thisItem] = dOutVal;

			// update the data index
			++thisItem;
		}


		//
		// update indices and do environmental preds
		//
		thisQuant += numSplines;
		thisCoeff += numSplines;
		for (int pred=1; pred<nPreds; pred++)
		{
			int numSplines = pSplines[pred];
			double *pQuantiles = &pQuants[thisQuant];
			double *pCoefficients = &pCoeffs[thisCoeff];

			for ( int i=0; i<nRows; i++ )
			{
				// get the data value
				double dDataVal = pInData[thisItem];

				// calculate the transform
				double dOutVal = CalculateGDMTransform(dDataVal, numSplines, pQuantiles, pCoefficients);

				// update output table
				pOutData[thisItem] = dOutVal;

				// update the data index
				++thisItem;
			}

			//
			// update indices
			//
			thisQuant += numSplines;
			thisCoeff += numSplines;
		}

	}
	else
	{
		for (int pred=0; pred<nPreds; pred++)
		{
			int numSplines = pSplines[pred];
			double *pQuantiles = &pQuants[thisQuant];
			double *pCoefficients = &pCoeffs[thisCoeff];

			for ( int i=0; i<nRows; i++ )
			{
				// get the data value
				double dDataVal = pInData[thisItem];

				// calculate the transform
				double dOutVal = CalculateGDMTransform(dDataVal, numSplines, pQuantiles, pCoefficients);

				// update output table
				pOutData[thisItem] = dOutVal;

				// update the data index
				++thisItem;
			}

			//
			// update indices
			//
			thisQuant += numSplines;
			thisCoeff += numSplines;
		}
	}
}


//
// Calculate the GDM transform
//
double CalculateGDMTransform(double dValue, int nSplines, double *pQuants, double *pCoeffs)
{
	double dVal = 0.0;
	for ( int i=0; i<nSplines; i++ )
	{
		if (i == 0) // first spline
		{
			dVal += DoSplineCalc( dValue, pQuants[i], pQuants[i], pQuants[i+1] ) * pCoeffs[i];
		}

		else if (i == nSplines-1) // last spline
		{
			dVal += DoSplineCalc( dValue, pQuants[i-1], pQuants[i], pQuants[i] ) * pCoeffs[i];
		}

		else  // a middle spline
		{
			dVal += DoSplineCalc( dValue, pQuants[i-1], pQuants[i], pQuants[i+1] ) * pCoeffs[i];
		}
	}
	return(dVal);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Local Support Functions ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Predict dissimilarity from spline data and coefficients
//
double CalcDissimilarity( double *pData, double *pCoeffs, int nRows, int nCols, int nIndex )
{
	double dVal = 0.0;

	double *pTmp = &pData[nIndex];
	for ( int i=0; i<nCols; i++ )
	{
		dVal += *pTmp * pCoeffs[i];
		pTmp += nRows;
	}
	return(dVal);
}


//
// calculate total number of splines
//
int GetTotalSplineCount(int *pSplines, int nPreds)
{
	int nCount = 0;
	for ( int i=0; i<nPreds; i++ )
	{
		nCount += pSplines[i];
	}
	return(nCount);
}



//
// Construct a regression matrix for predictors
//
double *ConstructMatrix(int nDoGeo, double *pData, double *pQuants, int nPreds, int *pSplines, int nRows)
{
	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	//Message(nTotalSplines, "nTotalSplines");

	//
	// Construct the predictor matrix with an extra column for the intercept
	//
	int nSize = nRows * (nTotalSplines + 1);
	double *pPredData = new double [nSize];
	if ( NULL == pPredData )
	{
		//Message("Cannot allocate Predictor Data", "ERROR in ConstructMatrix");
		return(NULL);
	}
	for ( int i=0; i<nSize; i++ ) pPredData[i] = 0.0;


	//
	// Initialise the intercept column
	//
	double *pLoc = &pPredData[0];
	for ( int i=0; i<nRows; i++ ) pLoc[i] = 1.0;
	pLoc += nRows;   // start at the first column after the intercept


	//
	// Initialise for the geographic distance if we are using it
	//
	if (nDoGeo)
	{
		//
		// Initialise the matrix for the geographic followed by the environmental predictors
		//
		double dMin,dMid,dMax,dVal1,dVal2;
		int CurrentSpline = 0;

		for ( int pred=0; pred<nPreds; pred++ )
		{
			int d1Offset = ((LEADING_COLS+pred-1)*nRows);
			int d2Offset = ((LEADING_COLS+pred+nPreds-2)*nRows);

			int nThis = 0;
			for ( int i=0; i<nRows; i++ )
			{
				if ( pred == 0 )
				{
					double distX = fabs( pData[(COL_SITE1_X0 * nRows)+i] - pData[(COL_SITE2_X1 * nRows)+i] );
					double distY = fabs( pData[(COL_SITE1_Y0 * nRows)+i] - pData[(COL_SITE2_Y1 * nRows)+i] );

					dVal1 = 0.0;
					dVal2 = sqrt( ( distX * distX ) + ( distY * distY ) );
				}
				else
				{
					dVal1 = pData[d1Offset+i];
					dVal2 = pData[d2Offset+i];
				}

				for ( int spl = 0; spl<pSplines[pred]; spl++ )
				{
					if ( spl == 0 )                              // first spline
					{
						dMin = pQuants[CurrentSpline+spl+0];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];
					}

					else if ( spl == (pSplines[pred]-1) )        // last spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+0];
					}

					else	                                     // a middle spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];
					}

					// calculate the spline values for each site
					double d1a = DoSplineCalc( dVal1, dMin, dMid, dMax );
					double d2a = DoSplineCalc( dVal2, dMin, dMid, dMax );

					// set the distance between sites for this spline
					pLoc[ (spl * nRows) + nThis ] = fabs( d2a - d1a );

				} // for ( int spl = 0; spl<pSplines[pred]; spl++ )

				// increment for the next row
				++nThis;

			} // for ( int i=0; i<nRows; i++ )

			// increment pointer to the start of the next predictor's spline columns
			pLoc += pSplines[pred] * nRows;

			// increment the current quantile index for the next pred
			CurrentSpline += pSplines[pred];

		} // for ( int pred=0; pred<nPreds; pred++ )
	}


	else
	{
		//
		// Initialise the matrix for the environmental predictors only
		//
		double dMin,dMid,dMax;
		int CurrentSpline = 0;

		for ( int pred=0; pred<nPreds; pred++ )
		{
			int d1Offset = ((LEADING_COLS+pred)*nRows);
			int d2Offset = ((LEADING_COLS+nPreds+pred)*nRows);

			int nThis = 0;
			for ( int i=0; i<nRows; i++ )
			{
				double dVal1 = pData[d1Offset+i];
				double dVal2 = pData[d2Offset+i];

				for ( int spl = 0; spl<pSplines[pred]; spl++ )
				{
					if ( spl == 0 )                              // first spline
					{
						dMin = pQuants[CurrentSpline+spl+0];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];

					}

					else if ( spl == (pSplines[pred]-1) )        // last spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+0];
					}

					else	                                     // a middle spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];
					}


					// calculate the spline values for each site
					double d1a = DoSplineCalc( dVal1, dMin, dMid, dMax );
					double d2a = DoSplineCalc( dVal2, dMin, dMid, dMax );

					// set the distance between sites for this spline
					pLoc[ (spl * nRows) + nThis ] = fabs( d2a - d1a );

				} // for ( int spl = 0; spl<pSplines[pred]; spl++ )

				// increment for the next row
				++nThis;

			} // for ( int i=0; i<nRows; i++ )

			// increment pointer to the start of the next predictor's spline columns
			pLoc += pSplines[pred] * nRows;

			// increment the current quantile index for the next pred
			CurrentSpline += pSplines[pred];

		} // for ( int pred=0; pred<nPreds; pred++ )

	} // environmental predictors only

	return(pPredData);
}



//
// Calculate the I-Spline value for dVal given quantiles q1, q2, q3
//
double DoSplineCalc( double dVal, double q1, double q2, double q3 )
{
	if ( dVal <= q1 ) return(0.0);

	else if ( dVal >= q3 ) return( 1.0 );

	else if ( ( q1 < dVal ) && ( dVal < q2 ) )
		return( ( ( ( dVal - q1 ) * ( dVal - q1 ) ) / ( ( q2 - q1 ) * ( q3 - q1 ) ) ) );

	else
		return( ( 1.0 - ( ( ( q3 - dVal ) * ( q3 - dVal ) ) / ( ( q3 - q2 ) * ( q3 - q1) ) ) ) );
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Local Debugging Functions /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Display the contents of the quantile vector
//
void ShowQuantiles(double *pQuants, int nPreds, int *pSplines)
{
	char buff[1024];
	double *pTmp = &pQuants[0];
	for ( int i=0; i<nPreds; i++ )
	{
		//sprintf( buff, "Quant %d: ", i+1 );
		//snprintf( buff, sizeof(buff), "Quant %d: ", i+1 );
		snprintf( buff, 10000, "Quant %d: ", i+1 );
		for ( int j=0; j<pSplines[i]; j++ )
		{
			// Modified by DNL:
			//sprintf( buff, "%s %lf ", buff, *pTmp);
			//sprintf( buff, "%s %f ", buff, *pTmp);
			//snprintf( buff, sizeof(buff), "%s %f ", buff, *pTmp);
			snprintf( buff, 10000, "%s %f ", buff, *pTmp);
			++pTmp;
		}
		//Message(buff);
	}
}




//
// Write the Predictor Matrix to a comma delimited file
//
void DebugPredMatrix(char *pPath, double *pPredData, int nRows, int nPreds, int *pSplines, int nCols)
{
	FILE *fp = fopen(pPath, "w+t");

	// write the header
	fprintf( fp, "Intercept,");
	for ( int p=0; p<nPreds; p++ )
	{
		for ( int s=0; s<pSplines[p]; s++ )
		{
			fprintf( fp, "Pred%dSpline%d", p+1,s+1);
			if ( s < pSplines[p]-1 )
				fprintf( fp, "," );
		}
		if ( p < nPreds-1 )
			fprintf( fp, "," );
		else
			fprintf( fp, "\n" );

	}

	for ( int i=0; i<nRows; i++ )
	{
		for ( int j=0; j<nCols; j++ )
		{
			// Modified by DNL:
			//fprintf( fp, "%lf", pPredData[(j*nRows)+i]);
			fprintf( fp, "%f", pPredData[(j*nRows)+i]);

			if ( j < nCols-1 )
				fprintf(fp, "," );
			else
				fprintf(fp, "\n" );
		}
	}

	if (fp) fclose(fp);
}


#else // Linux


#define PERMS 0644 // RW for owner R for others

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Functions called from R -Stats ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Main GDM fitting functions called from R
//
//
// Note that pData represents a column major matrix
//
void GDM_FitFromTable(char **wspath,
				      double *pData,
	                  int *pDoGeo, int *pPreds,
				      int *pRows, int *pCols,
				      int *pSplines, double *pQuantiles,
				      double *pGDMDev, double *pNullDev, double *pExpDev,
				      double *pIntercept, double *pCoeffs,
				      double *pY, // observed
				      double *pX, // predicted
				      double *pE) // ecological dist
{
	//
	// Everything passed from R needs to be dereferenced...
	//
	int nDoGeo = *pDoGeo;
	int nPreds = *pPreds;
	int nRows = *pRows;
	//int nCols = *pCols;


	//
	// Get the response column
	//
	double *pResponse = &pData[COL_RESPONSE * nRows];
	if ( NULL == pResponse )
	{
		//Message("pResponse is NULL", "ERROR in GDM_FitFromTable");
		return;
	}


	//
	// Get the weights column
	//
	double *pWeights = &pData[COL_WEIGHTS * nRows];
	if ( NULL == pWeights )
	{
		//Message("pWeights is NULL", "ERROR in GDM_FitFromTable");
		return;
	}


	//
	// Get the predictor matrix
	//
	double *pPredData = ConstructMatrix(nDoGeo, pData, pQuantiles, nPreds, pSplines, nRows);
	if ( NULL == pPredData )
	{
		//Message("pPredData is NULL", "ERROR in GDM_FitFromTable");
		return;
	}
	// Sum the splines vector for the total number of splines
	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	// Message(nTotalSplines, "nTotalSplines");
	// pQuantiles will have a length of nPreds * nTotalSplines
	// ShowQuantiles(pQuantiles, nPreds, pSplines);
	//DebugPredMatrix("PredMatrix.csv", pPredData, nRows, nPreds, pSplines, nTotalSplines+1);


	//
	// Write a binary file image of the predictor matrix
	//file extension
	std::string bin = ".bin";
	char *cbin = new char[bin.length() + 1];
	strcpy(cbin, bin.c_str());

	char lpTmpFile[1024];
	strlcpy(lpTmpFile, *wspath, sizeof(lpTmpFile));
	strlcat(lpTmpFile, cbin, sizeof(lpTmpFile));
	//Rcpp::Rcout << lpTmpFile << std::endl;

	int h = creat( lpTmpFile, PERMS );
	if ( h < 0 )
	{
		//Message("Cannot create binary image file", "ERROR in GDM_FitFromTable");
		if (pPredData) delete[] pPredData;
		return;
	}
	int nThis = 0;
	// for ( int i=0; i<(nTotalSplines+1); i++ )
	// {
	// 	write( h, &pPredData[nThis], nRows * sizeof( double ) );
	// 	nThis += nRows;
	// }
	// close( h );

	for ( int i=0;i<(nTotalSplines+1); i++ )
    {
        ssize_t szTotalToWrite = nRows * sizeof( double );
        ssize_t szLeftToWrite = szTotalToWrite;
        ssize_t szBytesWritten = 0;
        while( szLeftToWrite > 0 )
            {
                szBytesWritten = write( h, &pPredData[nThis] + szTotalToWrite - szLeftToWrite,  szLeftToWrite);
                if ( szBytesWritten == -1 )
                    {
                        close(h);
                        return; // an error occurred, maybe keep trying if it was -EAGAIN?
                    }
                else
                    {
                        szLeftToWrite -= szBytesWritten;
                    }
            }
        nThis += nRows;
    }
	close(h);



	//
	// Do the matrix regression
	//
	double dGDMDeviance;
	double *pCoefficients = WeightedNNLSRegression( lpTmpFile,
		                                            pPredData,
											        nRows,
											        nTotalSplines+1,
											        pResponse,
											        &dGDMDeviance,
											        pWeights);

	//
	// remove the temporary matrix file if it exists (and it should!!) now that we don't need it
	//
	//if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
	if ( ( access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );


	//
	// create a NULL model and return the deviance
	//
	double dNullDeviance = GetWeightedNULLDeviance( nRows, pResponse, pWeights );


	//
	// calculate deviance explained as a percentage
	//
	double dDevianceExplained = ( 1.0 - ( dGDMDeviance / dNullDeviance ) ) * 100;


	//
	// display results
	//
	pGDMDev[0] = dGDMDeviance;
	pNullDev[0] = dNullDeviance;
	pExpDev[0] = dDevianceExplained;
	pIntercept[0] = pCoefficients[0];
	for ( int i=1; i<nTotalSplines+1; i++ )
	{
		pCoeffs[i-1] = pCoefficients[i];
	}


	//
	// copy to response data to a vector
	//
	for ( int i=0; i<nRows; i++ )
	{
		pY[i] = pResponse[i];
		pE[i] = CalcDissimilarity( pPredData, pCoefficients, nRows, nTotalSplines+1, i );
		pX[i] = 1.0 - exp(-pE[i]);
	}


	//
	// Clean up
	//
	if (pPredData) delete[] pPredData;
	if ( pCoefficients ) delete[] pCoefficients;
}





//
// Main GDM predict function called from R
//
//
// Note that pData represents a column major matrix
//
void GDM_PredictFromTable(double *pData,
		                  int *pDoGeo, int *pPreds, int *pRows,
					      double *pQuantiles, int *pSplines, double *pCoeffs,
					      double *pX) // predicted
{
	//
	// Everything passed from R needs to be dereferenced...
	//
	int nDoGeo = *pDoGeo;
	int nPreds = *pPreds;
	int nRows  = *pRows;


	//
	// Get the predictor matrix
	//
	double *pPredData = ConstructMatrix(nDoGeo, pData, pQuantiles, nPreds, pSplines, nRows);
	if ( NULL == pPredData )
	{
		//Message("pPredData is NULL", "ERROR");
		return;
	}

	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	//Message(nTotalSplines, "nTotalSplines");
	// pQuantiles will have a length of nPreds * nTotalSplines
	//ShowQuantiles(pQuantiles, nPreds, pSplines);
	//DebugPredMatrix("PredMatrix.csv", pPredData, nRows, nPreds, pSplines, nTotalSplines+1);


	//
	// copy predicted data to pX
	//
	for ( int i=0; i<nRows; i++ )
	{
		pX[i] = 1.0 - exp(-CalcDissimilarity( pPredData, pCoeffs, nRows, nTotalSplines+1, i ));
	}


	//
	// Clean up
	//
	if (pPredData) delete[] pPredData;
}




//
// Populate pPredData as a transformed GDM predictor to plot in R-Stats
//
void GetPredictorPlotData( double *pPredData, int *pLength,
						   double *pCoefficients,
						   double *pQuantiles,
						   int *pSplines )
{
	int nLen  = *pLength;
	int nSplines = *pSplines;

	double dInc = fabs(pQuantiles[nSplines-1]-pQuantiles[0]) / nLen;
	double dVal = pQuantiles[0];
	for ( int i=0; i<nLen; i++ )
	{
		pPredData[i] = 0;
		for ( int j=0; j<nSplines; j++ )
		{
			if ( j==0 )
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j], pQuantiles[j], pQuantiles[j+1] );
			}

			else if ( j == nSplines-1 )
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j-1], pQuantiles[j], pQuantiles[j] );
			}

			else
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j-1], pQuantiles[j], pQuantiles[j+1] );
			}
		}
		dVal += dInc;
	}
}


//
//
// Transform a dataset using the coefficients and spline positions from a GDM model
//
//
void GDM_TransformFromTable(int *pRows, int *pCols,
	                        int *pDoGeo, int *pPreds,
							int *pSplines, double *pQuants, double *pCoeffs,
							double *pInData,
							double *pOutData)
{
	//
	// Everything passed from R needs to be dereferenced...
	//
	int nRows = *pRows;
	//int nCols = *pCols;
	int nDoGeo = *pDoGeo;
	int nPreds = *pPreds;

	int thisItem = 0;
	int thisQuant = 0;
	int thisCoeff = 0;

	if (nDoGeo)
	{
		int numSplines = pSplines[0];
		double *pQuantiles = &pQuants[0];
		double *pCoefficients = &pCoeffs[0];

		//
		// Setup a linear scalar to approximate the general slope of
		// the model curve to apply to the distance from the left edge.
		// Distances exceeding this are simply extrapolated by the linear scalar.
		//
		double dTranX = 0.0;
		double dTranY = 0.0;
		for ( int i=0; i<numSplines; i++ )
		{
			dTranY += pCoefficients[i];

			// get the maximum quantile value here
			if ( i == numSplines-1 )
			{
				dTranX = pQuantiles[i];
			}
		}
		double dScalar = dTranY / dTranX;

		//
		// now find the euclidean distance from the left edge and top edge and transform
		//
		// get min easting
		double *pXCol = &pInData[thisItem];
		double dMinX = pXCol[0];
		for ( int i=1; i<nRows; i++ )
		{
			if (pXCol[i] < dMinX) dMinX = pXCol[i];
		}

		// get min northing
		double *pYCol = &pInData[thisItem+nRows];
		double dMinY = pYCol[0];
		for ( int i=1; i<nRows; i++ )
		{
			if (pYCol[i] < dMinY) dMinY = pYCol[i];
		}


		// Calculate the transform for X
		for ( int i=0; i<nRows; i++ )
		{
			// get the data value
			double dDataVal = pInData[thisItem];

			// calculate the transform
			double dOutVal = dScalar * fabs(dDataVal-dMinX);

			// update output table
			pOutData[thisItem] = dOutVal;

			// update the data index
			++thisItem;
		}


		// Calculate the transform for Y
		for ( int i=0; i<nRows; i++ )
		{
			// get the data value
			double dDataVal = pInData[thisItem];

			// calculate the transform
			double dOutVal = dScalar * fabs(dDataVal-dMinY);

			// update output table
			pOutData[thisItem] = dOutVal;

			// update the data index
			++thisItem;
		}


		//
		// update indices and do environmental preds
		//
		thisQuant += numSplines;
		thisCoeff += numSplines;
		for (int pred=1; pred<nPreds; pred++)
		{
			int numSplines = pSplines[pred];
			double *pQuantiles = &pQuants[thisQuant];
			double *pCoefficients = &pCoeffs[thisCoeff];

			for ( int i=0; i<nRows; i++ )
			{
				// get the data value
				double dDataVal = pInData[thisItem];

				// calculate the transform
				double dOutVal = CalculateGDMTransform(dDataVal, numSplines, pQuantiles, pCoefficients);

				// update output table
				pOutData[thisItem] = dOutVal;

				// update the data index
				++thisItem;
			}

			//
			// update indices
			//
			thisQuant += numSplines;
			thisCoeff += numSplines;
		}

	}
	else
	{
		for (int pred=0; pred<nPreds; pred++)
		{
			int numSplines = pSplines[pred];
			double *pQuantiles = &pQuants[thisQuant];
			double *pCoefficients = &pCoeffs[thisCoeff];

			for ( int i=0; i<nRows; i++ )
			{
				// get the data value
				double dDataVal = pInData[thisItem];

				// calculate the transform
				double dOutVal = CalculateGDMTransform(dDataVal, numSplines, pQuantiles, pCoefficients);

				// update output table
				pOutData[thisItem] = dOutVal;

				// update the data index
				++thisItem;
			}

			//
			// update indices
			//
			thisQuant += numSplines;
			thisCoeff += numSplines;
		}
	}
}


//
// Calculate the GDM transform
//
double CalculateGDMTransform(double dValue, int nSplines, double *pQuants, double *pCoeffs)
{
	double dVal = 0.0;
	for ( int i=0; i<nSplines; i++ )
	{
		if (i == 0) // first spline
		{
			dVal += DoSplineCalc( dValue, pQuants[i], pQuants[i], pQuants[i+1] ) * pCoeffs[i];
		}

		else if (i == nSplines-1) // last spline
		{
			dVal += DoSplineCalc( dValue, pQuants[i-1], pQuants[i], pQuants[i] ) * pCoeffs[i];
		}

		else  // a middle spline
		{
			dVal += DoSplineCalc( dValue, pQuants[i-1], pQuants[i], pQuants[i+1] ) * pCoeffs[i];
		}
	}
	return(dVal);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Local Support Functions ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Predict dissimilarity from spline data and coefficients
//
double CalcDissimilarity( double *pData, double *pCoeffs, int nRows, int nCols, int nIndex )
{
	double dVal = 0.0;

	double *pTmp = &pData[nIndex];
	for ( int i=0; i<nCols; i++ )
	{
		dVal += *pTmp * pCoeffs[i];
		pTmp += nRows;
	}
	return(dVal);
}


//
// calculate total number of splines
//
int GetTotalSplineCount(int *pSplines, int nPreds)
{
	int nCount = 0;
	for ( int i=0; i<nPreds; i++ )
	{
		nCount += pSplines[i];
	}
	return(nCount);
}


//
// Construct a regression matrix for predictors
//
double *ConstructMatrix(int nDoGeo, double *pData, double *pQuants, int nPreds, int *pSplines, int nRows)
{
	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	//Message(nTotalSplines, "nTotalSplines");

	//
	// Construct the predictor matrix with an extra column for the intercept
	//
	long long nSize = nRows * (nTotalSplines + 1);
	double *pPredData = new double [nSize];
	if ( NULL == pPredData )
	{
		//Message("Cannot allocate Predictor Data", "ERROR in ConstructMatrix");
		return(NULL);
	}
	for ( long long i=0; i<nSize; i++ ) pPredData[i] = 0.0;


	//
	// Initialise the intercept column
	//
	double *pLoc = &pPredData[0];
	for ( long long i=0; i<nRows; i++ ) pLoc[i] = 1.0;
	pLoc += nRows;   // start at the first column after the intercept


	//
	// Initialise for the geographic distance if we are using it
	//
	if (nDoGeo)
	{
		//
		// Initialise the matrix for the geographic followed by the environmental predictors
		//
		double dMin,dMid,dMax,dVal1,dVal2;
		int CurrentSpline = 0;

		for ( int pred=0; pred<nPreds; pred++ )
		{
			int d1Offset = ((LEADING_COLS+pred-1)*nRows);
			int d2Offset = ((LEADING_COLS+pred+nPreds-2)*nRows);

			int nThis = 0;
			for ( long long i=0; i<nRows; i++ )
			{
				if ( pred == 0 )
				{
					double distX = fabs( pData[(COL_SITE1_X0 * nRows)+i] - pData[(COL_SITE2_X1 * nRows)+i] );
					double distY = fabs( pData[(COL_SITE1_Y0 * nRows)+i] - pData[(COL_SITE2_Y1 * nRows)+i] );

					dVal1 = 0.0;
					dVal2 = sqrt( ( distX * distX ) + ( distY * distY ) );
				}
				else
				{
					dVal1 = pData[d1Offset+i];
					dVal2 = pData[d2Offset+i];
				}

				for ( int spl = 0; spl<pSplines[pred]; spl++ )
				{
					if ( spl == 0 )                              // first spline
					{
						dMin = pQuants[CurrentSpline+spl+0];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];
					}

					else if ( spl == (pSplines[pred]-1) )        // last spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+0];
					}

					else	                                     // a middle spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];
					}

					// calculate the spline values for each site
					double d1a = DoSplineCalc( dVal1, dMin, dMid, dMax );
					double d2a = DoSplineCalc( dVal2, dMin, dMid, dMax );

					// set the distance between sites for this spline
					pLoc[ (spl * nRows) + nThis ] = fabs( d2a - d1a );

				} // for ( int spl = 0; spl<pSplines[pred]; spl++ )

				// increment for the next row
				++nThis;

			} // for ( long long i=0; i<nRows; i++ )

			// increment pointer to the start of the next predictor's spline columns
			pLoc += pSplines[pred] * nRows;

			// increment the current quantile index for the next pred
			CurrentSpline += pSplines[pred];

		} // for ( int pred=0; pred<nPreds; pred++ )
	}


	else
	{
		//
		// Initialise the matrix for the environmental predictors only
		//
		double dMin,dMid,dMax;
		int CurrentSpline = 0;

		for ( int pred=0; pred<nPreds; pred++ )
		{
			int d1Offset = ((LEADING_COLS+pred)*nRows);
			int d2Offset = ((LEADING_COLS+nPreds+pred)*nRows);

			int nThis = 0;
			for ( long long i=0; i<nRows; i++ )
			{
				double dVal1 = pData[d1Offset+i];
				double dVal2 = pData[d2Offset+i];

				for ( int spl = 0; spl<pSplines[pred]; spl++ )
				{
					if ( spl == 0 )                              // first spline
					{
						dMin = pQuants[CurrentSpline+spl+0];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];

					}

					else if ( spl == (pSplines[pred]-1) )        // last spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+0];
					}

					else	                                     // a middle spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];
					}


					// calculate the spline values for each site
					double d1a = DoSplineCalc( dVal1, dMin, dMid, dMax );
					double d2a = DoSplineCalc( dVal2, dMin, dMid, dMax );

					// set the distance between sites for this spline
					pLoc[ (spl * nRows) + nThis ] = fabs( d2a - d1a );

				} // for ( int spl = 0; spl<pSplines[pred]; spl++ )

				// increment for the next row
				++nThis;

			} // for ( int i=0; i<nRows; i++ )

			// increment pointer to the start of the next predictor's spline columns
			pLoc += pSplines[pred] * nRows;

			// increment the current quantile index for the next pred
			CurrentSpline += pSplines[pred];

		} // for ( int pred=0; pred<nPreds; pred++ )

	} // environmental predictors only

	return(pPredData);
}




//
// Calculate the I-Spline value for dVal given quantiles q1, q2, q3
//
double DoSplineCalc( double dVal, double q1, double q2, double q3 )
{
	if ( dVal <= q1 ) return(0.0);

	else if ( dVal >= q3 ) return( 1.0 );

	else if ( ( q1 < dVal ) && ( dVal < q2 ) )
		return( ( ( ( dVal - q1 ) * ( dVal - q1 ) ) / ( ( q2 - q1 ) * ( q3 - q1 ) ) ) );

	else
		return( ( 1.0 - ( ( ( q3 - dVal ) * ( q3 - dVal ) ) / ( ( q3 - q2 ) * ( q3 - q1) ) ) ) );
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Local Debugging Functions /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Display the contents of the quantile vector
//
// void ShowQuantiles(double *pQuants, int nPreds, int *pSplines)
// {
// 	char buff[1024];
// 	double *pTmp = &pQuants[0];
// 	for ( int i=0; i<nPreds; i++ )
// 	{
// 		sprintf( buff, "Quant %d: ", i+1 );
//
// 		for ( int j=0; j<pSplines[i]; j++ )
// 		{
// 			// Modified by MCF (2/03/2020):
// 			//sprintf( buff, "%s %lf ", buff, *pTmp);
// 			sprintf( buff, "%s %f ", buff, *pTmp);
// 			++pTmp;
// 		}
// 		//Message(buff);
// 	}
// }




//
// Write the Predictor Matrix to a comma delimited file
//
void DebugPredMatrix(char *pPath, double *pPredData, int nRows, int nPreds, int *pSplines, int nCols)
{
	FILE *fp = fopen(pPath, "w+t");

	// write the header
	fprintf( fp, "Intercept,");
	for ( int p=0; p<nPreds; p++ )
	{
		for ( int s=0; s<pSplines[p]; s++ )
		{
			fprintf( fp, "Pred%dSpline%d", p+1,s+1);
			if ( s < pSplines[p]-1 )
				fprintf( fp, "," );
		}
		if ( p < nPreds-1 )
			fprintf( fp, "," );
		else
			fprintf( fp, "\n" );

	}

	for ( int i=0; i<nRows; i++ )
	{
		for ( int j=0; j<nCols; j++ )
		{
			fprintf( fp, "%lf", pPredData[(j*nRows)+i]);

			if ( j < nCols-1 )
				fprintf(fp, "," );
			else
				fprintf(fp, "\n" );
		}
	}

	if (fp) fclose(fp);
}
#endif
