/* ====================================================================
*           R-ScaLAPACK version 0.4.x:  ScaLAPACK interface to R
*              Oak Ridge National Laboratory, Oak Ridge TN.
*        Authors: David Bauer, Guruprasad Kora, Nagiza. F. Samatova, 
*                            Srikanth Yoginath.
*     Contact: Nagiza F. Samatova; (865) 241-4351; samatovan@ornl.gov
*     Contact: Guruprasad Kora; (865) 576-6210; koragh@ornl.gov
*                 Computer Science and Mathematics Division
*             Oak Ridge National Laboratory, Oak Ridge TN 37831 
*                   (C) 2004 All Rights Reserved
*
*                              NOTICE
*
* Permission to use, copy, modify, and distribute this software and
* its documentation for any purpose and without fee is hereby granted
* provided that the above copyright notice appear in all copies and
* that both the copyright notice and this permission notice appear in
* supporting documentation.
*
* Neither the Oak Ridge National Laboratory nor the Authors make any
* representations about the suitability of this software for any
* purpose.  This software is provided ``as is'' without express or
* implied warranty.
*
* RScaLAPACK (http://www.aspect-sdm.org/Parallel-R) was funded
* as part of the Scientific Data Management Center
* (http://sdm.lbl.gov/sdmcenter) under the Department of Energy's 
* Scientific Discovery through Advanced Computing (DOE SciDAC) program
* (http://www.scidac.org ). 
* ======================================================================
*
*       Based on:
*
*       -- ScaLAPACK example code --
*       University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*       and University of California, Berkeley.
*
*       Written by Antoine Petitet, August 1995 (petitet@cs.utk.edu)
*
* ==============================================================================
*/

#include "ParallelAgent.h"

void printContents (double *temp, int size)
{
	int i=0;
	for (i=0; i<size; i++)
		D_Rprintf(("vector [%d] - %f \n", i, *(temp+i)));
}

void PAdistData(double gMatrix[], int pgInfo[], int numRows, int numCols){

	int rowBlockSize, colBlockSize;
	int numProcessRows, numProcessCols;
	int currRow=0, currCol=0;
	int j=0;
	double *sendBuffer;
	int startIndex=0;
	int vectorSize;
	int one=1;
	int iCol=0, iRows=0;
	int sendToRank;
	int reqId=0;

	MPI_Request *requests;	
	MPI_Status *statuses;	
	rowBlockSize = pgInfo[4];
	colBlockSize = pgInfo[5];
	numProcessRows = pgInfo[6];
	numProcessCols = pgInfo[7];
/*
	used for non-blocking communication
	requests = (MPI_Request *) malloc (sizeof(MPI_Request) * rowBlockSize * iceil(numRows,rowBlockSize) * iceil(numCols,rowBlockSize));
	statuses = (MPI_Status *) malloc (sizeof(MPI_Status) * rowBlockSize * iceil(numRows,rowBlockSize) * iceil(numCols,rowBlockSize));
	
*/
	/* Loop over the column blocks */
	while (iCol < numCols)
	{	 
		j = 0;

		while (j < min( colBlockSize, numCols-iCol))
		{
			iRows =0;

			/* Loop over the block of rows */
			while  (iRows <  numRows )
			{
				D_Rprintf((" In third loop .. \n"));	
				vectorSize = min (rowBlockSize, numRows-iRows);

				sendToRank = currRow * numProcessCols + currCol;

				sendBuffer = gMatrix;

				startIndex = iRows + j * numRows + iCol * numRows ;

				/* 
				used for non-blocking communication

				PA_SendVectorToCR( &vectorSize, &one, sendBuffer+startIndex, &rowBlockSize, &sendToRank, requests+reqId);
				reqId++; 
				*/

				PA_SendVectorToCR( &vectorSize, &one, sendBuffer+startIndex, &rowBlockSize, &sendToRank);

				/*printContents (sendBuffer+startIndex, vectorSize);*/
				D_Rprintf(("Sent to rank - %d \n", sendToRank));

				currRow = (currRow+1) % numProcessRows;

				iRows = iRows + rowBlockSize;


			}

			j++;
			currRow = 0;	
		}

		currCol = (currCol+1) % numProcessCols;

		iCol = iCol + colBlockSize;	

	}

	D_Rprintf(("reqId : %d \n", reqId));	

	D_Rprintf(("Finished data distribution \n"));	

/* 
	used for non-blocking communication

	MPI_Waitall(reqId, requests, statuses); 
	free (requests);
	free (statuses);
*/	

}

void PAcollectData(double *retResult, int pgInfo[], int numRows, int numCols){

	int rowBlockSize, colBlockSize;
	int numProcessRows, numProcessCols;
	int currRow=0, currCol=0;
	int i=0,j=0;
	double *recvBuffer;
	int startIndex=0;
	int vectorSize;
	int one=1;
	int iColBlock=0, iRowBlock=0;
	int recvFromRank=0;
	int reqId=0;

/*	MPI_Request *requests;	
	MPI_Status *statuses;	
*/
	rowBlockSize = pgInfo[4];
	colBlockSize = pgInfo[5];
	numProcessRows = pgInfo[6];
	numProcessCols = pgInfo[7];
	recvBuffer = retResult;

	iColBlock = min( colBlockSize, numCols);


/* 
	used for non-blocking communication
	requests = (MPI_Request *) malloc (sizeof(MPI_Request) * rowBlockSize * iceil(numRows,rowBlockSize) * iceil(numCols,rowBlockSize));
	statuses = (MPI_Status *) malloc (sizeof(MPI_Status) * rowBlockSize * iceil(numRows,rowBlockSize) * iceil(numCols,rowBlockSize));
*/
	/* Handle the first block of column separately */
	while(i < iColBlock){
		currRow=0;

		iRowBlock = min(rowBlockSize, numRows);

		recvFromRank = currRow * numProcessCols + currCol;

		vectorSize = min (rowBlockSize, numRows);
		
		startIndex = i * numRows;

		/* 
		used for non-blocking communication

		PA_RecvVectorFromCR (&vectorSize, &one, recvBuffer+startIndex, &rowBlockSize, &recvFromRank, requests+reqId);
		reqId++;
		*/

		PA_RecvVectorFromCR (&vectorSize, &one, recvBuffer+startIndex, &rowBlockSize, &recvFromRank);
		currRow = ((currRow+1) % numProcessRows);

		/* Loop over the remainig block of rows */

		j =  iRowBlock;
		
		while ( j < numRows ){

			vectorSize = min (rowBlockSize, numRows-j);
			
			recvFromRank = currRow * numProcessCols + currCol;
			
			startIndex = i*numRows+j ;

			/* 
			used for non-blocking communication

			PA_RecvVectorFromCR (&vectorSize, &one, recvBuffer+startIndex, &rowBlockSize, &recvFromRank, requests+reqId);
			reqId++; 
			*/

			PA_RecvVectorFromCR (&vectorSize, &one, recvBuffer+startIndex, &rowBlockSize, &recvFromRank);

			j = j + rowBlockSize;

			currRow = ((currRow+1) % numProcessRows);

		}

		i++;
	}
			
	currCol = ((currCol+1) % numProcessCols);

	/* Handle the first block of column separately */

	while (iColBlock < numCols){

		iRowBlock = min(rowBlockSize,numRows);

		for (i=0;i < min (colBlockSize, numCols-iColBlock);i++) {
	
			currRow=0;

			vectorSize = min (rowBlockSize, numRows);

			recvFromRank = currRow * numProcessCols + currCol;
			
			startIndex = numRows*iColBlock + i*numRows ;
		
			/* 
			used for non-blocking communication

			PA_RecvVectorFromCR (&vectorSize, &one, recvBuffer+startIndex, &rowBlockSize, &recvFromRank, requests+reqId);
			reqId++; */

			PA_RecvVectorFromCR (&vectorSize, &one, recvBuffer+startIndex, &rowBlockSize, &recvFromRank);

			j = iRowBlock;

			currRow = ((currRow+1) % numProcessRows);

			/* Loop over the remainig block of rows */
			while ( j < numRows ){

				vectorSize = min (rowBlockSize, numRows-j);
				
				recvFromRank = currRow * numProcessCols + currCol;

				startIndex = numRows*iColBlock + i*numRows + j ;
				  
				/* 
				used for non-blocking communication

				PA_RecvVectorFromCR (&vectorSize, &one, recvBuffer+startIndex, &rowBlockSize, &recvFromRank, requests+reqId);
				reqId++;
				*/

				PA_RecvVectorFromCR (&vectorSize, &one, recvBuffer+startIndex, &rowBlockSize, &recvFromRank);

				j = j + rowBlockSize;

				currRow = ((currRow+1) % numProcessRows);
			}	


		}	
	
		currCol = ((currCol+1) % numProcessCols);
	
		iColBlock = iColBlock + colBlockSize;

	}
	D_Rprintf(("Finished data Collection  \n"));	

	/* 
	used for non-blocking communication

	MPI_Waitall(reqId, requests, statuses);
	
	free (requests);
	free (statuses);
	*/	
}
