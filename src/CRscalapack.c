/*===========================================================================
 *           R-ScaLAPACK version 0.2:  ScaLAPACK interface to R
 *              Oak Ridge National Laboratory, Oak Ridge TN.
 *      Authors: David Bauer, Nagiza. F. Samatova, Srikanth Yoginath
 *     Contact: Nagiza F. Samatova; (865) 241-4351; samatovan@ornl.gov
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
 * R-ScaLAPACK (http://www.aspect-sdm.org/R-ScaLAPACK) was funded
 * as part of the Scientific Data Management Center
 * (http://sdm.lbl.gov/sdmcenter) under the Department of Energy's 
 * Scientific Discovery through Advanced Computing (DOE SciDAC) program
 * (http://www.scidac.org ). 
=============================================================================*/
#include "CRscalapack.h"
#include <stdlib.h>

#ifdef DONT_SPAWN_R
#define Rprintf printf
#define AsInt(x) (SEXP)(x)
#endif

/* #define DEBUG_RSCALAPACK */

#ifndef DEBUG_RSCALAPACK
#define D_Rprintf(x)
#else
#define D_Rprintf(x) Rprintf x
#endif

/* **** CR_Exec ****
 * CR - stands for Child R-scalapack
 * CR_Exec is the main function that is invoked by all
 * the spawned processes 
 */
SEXP CR_Exec() {
	MPI_Comm mcParent;

	int iMyRank;
	int iNumProcs;
	int rel_flag = 0;
	int exitTemp =0;

	int ipGridAndDims[10]={0,0,0,0,0,0,0,0,0,0};	/* See CR_GetInputParams
													 * function for contents.
													 */

	D_Rprintf ((" SPAWNED PROCESS --- starts here !! \n"));

	/* Performing MPI_Init through BLACS */
	if (CR_InitializeEnv(&iMyRank, &iNumProcs) != 0)	{
		Rprintf("ERROR[1]: Initializing MPI thru BLACS ... FAILED .. EXITING !!\n");
		return AsInt(1);
	}

	D_Rprintf (("Initializing MPI thru BLACS successful \n"));

	/* Get the Parent Communicator */
	if (MPI_Comm_get_parent(&mcParent) != MPI_SUCCESS) {
		Rprintf("ERROR[2]: Getting Parent Comm ... FAILED .. EXITING !!\n");
		return AsInt(2);
	}

	if (mcParent == MPI_COMM_NULL) {
		char *cpNoParent = "No Parent found.  This program is not designed"
		" to be run by itself.\nThis program is used by the R-ScaLAPACK"
		" package as the computational engine.\n";
		fprintf(stderr, cpNoParent);
		return AsInt (3);
	}

	while (rel_flag != 1){

		/* Get the broadcasted input parameters from the Parallel Agent */
		if (CR_GetInputParams(mcParent,ipGridAndDims) != 0) {
			Rprintf ("ERROR[4]: FAILED while getting Input Parameters \n");
			return AsInt(4);
		}

		D_Rprintf(("Success:  Got Input Parameters from ParallelAgent\n"));

		rel_flag = ipGridAndDims[9];

		D_Rprintf(("relflag : %d \n", rel_flag));


		/* Call the requested ScaLAPACK Function */
		if (CR_CallScalapackFn(ipGridAndDims, iMyRank)!=0) {
			Rprintf ("ERROR[4]: FAILED while executing scalapack function  - CR_CallScalapackFn() \n");
			return AsInt(6);
		}

		D_Rprintf(("%d: After calling Fortran function (%d)\n", iMyRank, 
					ipGridAndDims[8]));

		/* If the function ID is 0 and the release flag is 1, then
		 * the function is sla.gridExit, so exit.
		 */ 
		if (ipGridAndDims[8] == 0 && rel_flag == 1){
			F77_CALL(blacs_exit_)(&exitTemp);   
			return AsInt(1); 
		}

	}	/*  Endof  while (...) */

	return AsInt(0);
}

/* Uses BLACS call to perform MPI_Init */
int CR_InitializeEnv(int *ipMyRank, int *ipNumProcs) {
	F77_CALL(blacs_pinfo_)(ipMyRank, ipNumProcs);
	return 0;
}

/* Receive - I/P Data dimensions and Process Grid Specifications from the parent
* 1. No. of rows in matrix A
* 2. No. of cols in matrix A
* 3. No. of rows in matrix B
* 4. No. of cols in matrix B
* 5. MB - Row Block size for matrix A
* 5. NB - Col Block size for matrix A
* 7. NPROW - Number of Process rows in the Process Grid - Row Block Size
* 8. NPCOL - Number of Process cols in the Process Grid - Col Block Size
* 9. Function id
*/
int CR_GetInputParams(MPI_Comm mcParent, int *ipGridAndDims) {

     MPI_Comm intercomm;

     if(MPI_Intercomm_merge(mcParent, 1, &intercomm)!= MPI_SUCCESS)
          return -1;
        
     if(MPI_Bcast(ipGridAndDims,10, MPI_INT, 0, intercomm) != MPI_SUCCESS) 
          return -2;
     else
          return 0;
}
/* ****  CR_CallScalapackFun  ****
 * This function calls the C wrappers around the Fortran wrappers around the
 * ScaLAPACK functions.
 * It basically does nothing but switch on the Function ID and call the
 * corresponding function.
 */ 
int CR_CallScalapackFn (int *ipGridAndDims, int iMyRank) {
	int iFunction = ipGridAndDims[8];

	switch(iFunction) {
		case 0:
			/* sla.gridInit, sla.gridExit, etc. */
			break;
		case 1:
			/* Solve */
			CRSF_solve (ipGridAndDims, iMyRank);
			break;
		case 2:
			/* SVD */
			CRSF_svd (ipGridAndDims, iMyRank);
			break;
		case 3:
			/* QR Decomposition */
			CRSF_qr(ipGridAndDims, iMyRank);
			break;
		case 4:
			/* Choleski Factorization */
			CRSF_chol(ipGridAndDims, iMyRank);
			break;
		case 5:
			/* Invert a (symmetric, positive definite, square) matrix from
			 * its Choleski decomposition */
			CRSF_chol2inv(ipGridAndDims, iMyRank);
			break;
		case 6:
			/* Eigen Value */
			CRSF_eigen(ipGridAndDims, iMyRank);
			break;
		default:
			Rprintf("%d:%s:%d: ASSERT ERROR: Reached unreachable default in switch statement!\n", iMyRank, __FILE__, __LINE__);
			return -1;

	}	/* Endof  switch (iFunction) */

	return 0;
}


/* Note:  CR_RecvVectorFromPA is translated to "cr_recvvectorfrompa__" 
 * for Fortran compatibility */
void CR_RecvVectorFromPA( int *ib, int *ia, double *A, int *mb ){
	MPI_Datatype GEMAT;
	MPI_Comm parent;

	MPI_Comm_get_parent (&parent);

	MPI_Type_vector(*ia, *ib, *mb, MPI_DOUBLE,&GEMAT);
	MPI_Type_commit (&GEMAT);

	MPI_Recv (A, 1, GEMAT,0 ,5000 ,parent, MPI_STATUS_IGNORE);

	MPI_Type_free (&GEMAT);
}

/* Note: CR_SendVectorToPA is translated to "cr_sendvectortopa__"
 *  for Fortran compatibility */
void CR_SendVectorToPA (int *ib, int *ia, double *work, int *mb){
	MPI_Datatype GEMAT;
	MPI_Comm parent;

	MPI_Comm_get_parent (&parent);

	MPI_Type_vector(*ia, *ib, *mb, MPI_DOUBLE,&GEMAT);
	MPI_Type_commit (&GEMAT);

	MPI_Send ( work, 1, GEMAT, 0, 15000 , parent);

	MPI_Type_free (&GEMAT);
}

/* Note:  CR_SendIntToPA is translated to "cr_sendinttopa__" 
 * for Fortran compatibility */
void CR_SendIntToPA(int *info, int *infoDim, int *tag) {
	MPI_Comm parent;
	MPI_Comm_get_parent (&parent);
	MPI_Send ( info, *infoDim, MPI_INT, 0, *tag , parent);
}

/* Note:  CR_SendDoubleToPA is translated to "cr_senddouble2pa__"
 *  for Fortran compatibility */
void CR_SendDoubleToPA(double *data, int *infoDim, int *tag) {
	MPI_Comm parent;
	MPI_Comm_get_parent (&parent);
	MPI_Send ( data, *infoDim, MPI_DOUBLE, 0, *tag , parent);
}

/* ****  Start CRSFs  -- Child R-ScaLAPACK Functions ****  */

/* **** CRSF_qr ****
 * This function is a C interface to the fortran implemented 
 * scalapack driver function "callpdgeqrf" that performs
 * qr decomposition on the input data.  
 */
int CRSF_qr(int dim[], int iMyRank) {
	int qrRank, sizeA, iMemSize;
	double *dpRet_val = NULL;
	double *dpWork = NULL;

	qrRank = min(dim[0], dim[1]);
	sizeA = dim[0]*dim[1];

	/* Calculate proper size for memory allocation.
	 *  ( sizeof(TAU array) ).
	 */

	dpRet_val = (double *) malloc(sizeof(double) * (qrRank) );
	memset(dpRet_val, 0xcc, sizeof(double) * (qrRank) );

	/* TODO Calculate proper iMemSize; */
	iMemSize = sizeA * 2 + 1024;
	dpWork = (double *) malloc(sizeof(double) * iMemSize);
	D_Rprintf(("%d: Preparing to call Fortran QR function\n", iMyRank));
	F77_CALL(callpdgeqrf)(dim, &sizeA, dpWork, &iMemSize,
			dpRet_val);
	D_Rprintf(("%d: After calling Fortran QR function\n", iMyRank));

	free(dpWork);
	free(dpRet_val);
	return 0;
}

/* **** CRSF_chol **** 
 * This function is a C interface to the fortran implemented 
 * scalapack driver function "callpdpotrf" that performs
 * Choleski factorization on the input data.
 */   
int CRSF_chol(int dim[], int iMyRank) {
	int iMemSize;
	double *dpWork = NULL;

	iMemSize =  dim[0] * dim[1] + dim [5] + 1;
	dpWork = (double *) malloc(sizeof(double) * iMemSize);

	F77_CALL(callpdpotrf)(dim, dpWork,&iMemSize);

	free (dpWork);
	return 0;
}

/* **** CRSF_chol2inv **** 
 * This function is a C interface to the fortran implemented 
 * scalapack driver function "callpdpotri" that performs
 * inverting a matrix from its Choleski Factorization
 */ 
int CRSF_chol2inv(int dim[], int iMyRank) {
	int iMemSize;
	double *dpWork = NULL;

	iMemSize = dim[0]*dim[1] + dim[5] + 1;
	dpWork = (double *) malloc(sizeof(double) * iMemSize);
	F77_CALL(callpdpotri)(dim, dpWork, &iMemSize);

	free(dpWork);
	return 0;
}

/* **** CRSF_eigen **** 
 * This function is a C interface to the fortran implemented 
 * scalapack driver function "callpdsyevd" that performs
 * eigen value decomposition or Spectral decomposition 
 */
int CRSF_eigen(int dim[], int iMyRank) {
	int iMemSize;
	double *dpRet_val = NULL;

	iMemSize = dim[0];
	if (iMyRank == 0){
		dpRet_val = (double *) malloc(sizeof(double) * iMemSize);
	}

	F77_CALL(callpdsyevd)(dim, dpRet_val);

	if (iMyRank == 0) 
		free(dpRet_val);
	return 0;
}

/* **** CRSF_svd ****
 *  This function is a C interface to the fortran implemented 
 *  scalapack driver function "callpdgesvd" that performs
 *  singular value decomposition 
 */
int CRSF_svd(int dim[], int iMyRank) {
	int iMemSize;
	double *dpWork = NULL;
	double *dpRet_val = NULL;

	int ipZero[] = { 0, 1, 2, 3 };
	int NPRow = dim[6];
	int NPCol = dim[7];
	int MyRow = iMyRank / NPCol;
	int MyCol = iMyRank % NPCol;

	if (dim[2] == 0 && dim[3] != 4) {
		/* A number of left or right vectors to return was specified.
		 ** Store that information in dim[2] and dim[3], because those
		 ** are already getting passed to the SVD function and aren't
		 ** being used for anything else.
		 **/
		dim[2] = (dim[3] & 0x01);
		dim[3] = (dim[3] & 0x02);
	} else {
		dim[2] = 1;		dim[3] = 1;
	}

	/* Calculate the size of the workspace. */
	{
		int iRowA = F77_CALL(numroc)(dim, dim+4, &MyRow, ipZero, &NPRow);
		int iRows = F77_CALL(numroc)(dim+1, dim+4, &MyRow, ipZero, &NPRow);
		int iNCol = F77_CALL(numroc)(dim+1, dim+5, &MyCol, ipZero, &NPCol);
		int iNQ = F77_CALL(numroc)(dim+1, dim+5, &MyCol, &iRowA, &NPCol);

		/* printf("Mem sizes: %d, %d, %d, %d\n", iRowA, iRows, iNCol, iNQ); */

		iMemSize = (max(1, iRowA) * 2 + max(1, iRows)) * iNCol + dim[1];
		iMemSize += max(iRowA, max(iRows, iNCol));
		iMemSize *= 2;
		iMemSize += max(NPRow*dim[4]*dim[1], 2 + 6*max(dim[0], dim[1]) +
				(dim[5]+1)*(iRowA+iNCol+1));
		
		/* printf("SVD: Calculated mem1: %d\n", iMemSize); */
		
		iMemSize += max( (dim[4]*(dim[4]-1))/2,
						 (2*iRows)*dim[4])+dim[4]*dim[4];
						 
		/* printf("SVD: Calculated mem2: %d\n", iMemSize); */

		/* TODO  BUG:  Fudge factor for memory allocation */
		iMemSize += 2*(iRowA + iNQ + max(1, iRowA)*iNCol*3 +
			max(0, NPRow*dim[4]*(dim[0] - dim[1])));

		/* printf("SVD: Fudged mem: %d\n", iMemSize); */

	}

	dpWork = (double *) malloc(sizeof(double) * iMemSize);

	if (iMyRank == 0 ){
		dpRet_val = (double *) malloc(sizeof(double)*min(dim[0],dim[1]));
	}

	F77_CALL(callpdgesvd)(dim, dpRet_val, dpWork, &iMemSize);

	free(dpWork);

	if (iMyRank == 0 )
		free(dpRet_val);
	return 0;
}

/* **** CRSF_solve ****
 *  This function is a C interface to the fortran implemented 
 *  scalapack driver function "callpdgesv" that solves the 
 *  equations A * X = B for 'X'
 */
int CRSF_solve(int dim[], int iMyRank) {
	int iMemSize;
	double *dpWork = NULL;

	int ipZero[] = { 0, 1, 2, 3 };
	int NPRow = dim[6];
	int NPCol = dim[7];
	int MyRow = iMyRank / NPCol;
	int MyCol = iMyRank % NPCol;

	int iNP = F77_CALL(numroc)(dim, dim+5, &MyRow, ipZero, &NPRow);
	int iLIPIV = (int) ( (4.0 * (iNP + dim[5]) + 7.0) / 8.0 );
	int iLLD = max(1, iNP);

	iMemSize = (2*F77_CALL(numroc)(dim, dim+5,&MyCol,ipZero,&NPCol)*iLLD +
			2*F77_CALL(numroc)(dim+3,dim+5,&MyCol,ipZero,&NPCol)*iLLD +
			max(max(iNP, iLIPIV), dim[5]) + dim[5] );
	/* TODO  WARNING  BUG !!!
	 * The amount of memory that I calculate Fortran as needing is a few bytes
	 * shy of what it calculates it needs.  (38*8 bytes shy for a 9x9 matrix;
	 * 65*8 bytes shy for a 2k x 2k matrix.  In order to protect against this
	 * problem, I allocate a few extra bytes.  The smallest safe variable is
	 * the last term in the allocation.
	 */

	iMemSize += 2 + max(max(iNP, iLIPIV), dim[5]);

	dpWork = (double *) malloc(sizeof(double) * iMemSize);
	memset(dpWork, 0xcc, sizeof(double) * iMemSize);

	D_Rprintf(("%d:  Allocated a work vector of size %d bytes\n",
				iMyRank, sizeof(double) * iMemSize));

	/* Rprintf("%d:  iNP = %d, iLLD = %d, MyRow = %d, MyCol = %d\n",
	 * iMyRank, iNP, iLLD, MyRow, MyCol);
	 *
	 * Rprintf("%d:  %d + %d + %d\n", iMyRank, 2*F77_CALL(numroc)
	 * (dim, dim+5,&MyCol,ipZero,&NPCol)*iLLD, 2*F77_CALL(numroc)
	 * (dim+3,dim+5,&MyCol,ipZero,&NPCol)*iLLD, max(max(iNP, iLIPIV), dim[5]) );
	 */

	F77_CALL(callpdgesv)(dim, dpWork, &iMemSize);

	free(dpWork);
	return 0;
}

/* ****  End CRSFs  -- Child R-ScaLAPACK Functions ****  */
