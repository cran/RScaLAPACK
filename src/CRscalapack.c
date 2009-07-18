/*======================================================================
 *           R-ScaLAPACK version 0.4.x:  ScaLAPACK interface to R
 *              Oak Ridge National Laboratory, Oak Ridge TN.
 *        Authors: David Bauer, Guruprasad Kora, Nagiza. F. Samatova, 
 *                            Srikanth Yoginath.
 *     Contact: Nagiza F. Samatova; (865) 241-4351; samatovan@ornl.gov
 *     Contact: Guruprasad Kora; (865) 576-6210; koragh@ornl.gov
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
=========================================================================*/
#include "CRscalapack.h"
#include <stdlib.h>

#define DONT_SPAWN_R

#ifdef DONT_SPAWN_R
#define Rprintf printf
#endif


#define AsInt(x) (SEXP)(x)

#if 0
#define DEBUG_RSCALAPACK
#endif

#ifndef DEBUG_RSCALAPACK
#define D_Rprintf(x)
#else
#define D_Rprintf(x) Rprintf x
#endif

MPI_Comm intercomm;

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

#ifdef OMPI
	dlopen("libmpi.so.0", RTLD_GLOBAL | RTLD_LAZY);
#endif


	/* Performing MPI_Init through BLACS */
	if (CR_InitializeEnv(&iMyRank, &iNumProcs) != 0)	{
		Rprintf("ERROR[1]: Initializing MPI thru BLACS ... FAILED .. EXITING !!\n");
		return AsInt(1);
	}



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


		rel_flag = ipGridAndDims[9];


		/* Call the requested ScaLAPACK Function */
		if (CR_CallScalapackFn(ipGridAndDims, iMyRank)!=0) {
			Rprintf ("ERROR[4]: FAILED while executing scalapack function  - CR_CallScalapackFn() \n");
			return AsInt(6);
		}


		/* If the function ID is 0 and the release flag is 1, then
		 * the function is sla.gridExit, so exit.
		 */ 
		if (ipGridAndDims[8] == 0 && rel_flag == 1){

			MPI_Comm_free( &intercomm );
			MPI_Finalize();
			return AsInt(1); 
		}

	}	/*  Endof  while (...) */


	//MPI_Barrier(intercomm);
	MPI_Comm_free( &intercomm );
	MPI_Finalize();


	return AsInt(0);
}

/* Uses BLACS call to perform MPI_Init */
int CR_InitializeEnv(int *ipMyRank, int *ipNumProcs) {
#if MPICH
	F77_CALL(blacs_pinfo)(ipMyRank, ipNumProcs);
#else
	F77_CALL( blacs_pinfo )( ipMyRank, ipNumProcs);
#endif
	return 0;
}

/* Receive - I/P Data dimensions and Process Grid Specifications from the parent
* 1. No. of rows in matrix A
* 2. No. of cols in matrix A
* 3. No. of rows in matrix B
* 4. No. of cols in matrix B
* 5. MB - Row Block size for matrix A
* 6. NB - Col Block size for matrix A
* 7. NPROW - Number of Process rows in the Process Grid - Row Block Size
* 8. NPCOL - Number of Process cols in the Process Grid - Col Block Size
* 9. Function id
* 10. Relaease Flag 
*/
int CR_GetInputParams(MPI_Comm mcParent, int *ipGridAndDims) {

	MPI_Comm parent;


	if (MPI_Comm_get_parent(&parent) != MPI_SUCCESS) {
		Rprintf("ERROR[2]: Getting Parent Comm ... FAILED .. EXITING !!\n");
		return AsInt(2);
	}

	if(MPI_Intercomm_merge(parent, 1, &intercomm)!= MPI_SUCCESS)
		return -1;

	if(MPI_Bcast(ipGridAndDims,10, MPI_INT, 0, intercomm) != MPI_SUCCESS) 
	{
		D_Rprintf(("Child: Broadcast error\n"));
		return -2;
	}
	else
	{
		return 0;
	}
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
		case 7:
			/* Matrix Multiplication */
			/*D_Rprintf(("Got every thing, Need to run multiplication alg.\n"));*/
			CRSF_multiply(ipGridAndDims, iMyRank);
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

	MPI_Type_vector(*ia, *ib, *ib, MPI_DOUBLE,&GEMAT);
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

	MPI_Type_vector(*ia, *ib, *ib, MPI_DOUBLE,&GEMAT);
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


/* Note:  CR_CheckFailFlag is translated to "crcheckfailflag__"
 *  for Fortran compatibility */
void CR_CheckFailFlag ( int* failFlag){
	int myFailFlag = *failFlag;
	MPI_Allreduce(&myFailFlag, failFlag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
}

/* ****  Start CRSFs  -- Child R-ScaLAPACK Functions ****  */

/* **** CRSF_qr ****
 * This function is a C interface to the fortran implemented 
 * scalapack driver function "callpdgeqrf" that performs
 * qr decomposition on the input data.  
 */
int CRSF_qr(int dim[], int iMyRank) {

	int iMemSize = 0;
	double *dpWork = NULL;

	int ipZero[] = { 0, 1, 2, 3 };
	int NPRow = dim[6];
	int NPCol = dim[7];
	int MyRow = iMyRank / NPCol;
	int MyCol = iMyRank % NPCol;

	int rowOfA = dim[0];
	int colOfA = dim[1];
	int rowBlockSize = dim[4];
	int colBlockSize = dim[5];

	/* Calculate required memory size */
	int localRowSizeOfA = F77_CALL(numroc)(&rowOfA, &rowBlockSize, &MyRow, ipZero, &NPRow);
	int localColSizeOfA = F77_CALL(numroc)(&colOfA, &colBlockSize, &MyCol, ipZero, &NPCol);
	int localSizeOfA = localRowSizeOfA * localColSizeOfA;

	int tauParam = min (rowOfA, colOfA);
	int lTau = F77_CALL(numroc)(&tauParam, &colBlockSize, &MyCol, ipZero, &NPCol);
	
	int iZero = 0;
	int iA = 1;
	int jA = 1;
	int iROff = rowBlockSize;
	int iCOff = colBlockSize;
	int iARow = F77_CALL(indxg2p)(&iA, &rowBlockSize, &MyRow,&iZero, &NPRow);
	int iACol = F77_CALL(indxg2p)(&jA, &colBlockSize, &MyCol,&iZero, &NPCol);
	int mp0Param = rowOfA + iROff;
	int mP0 = F77_CALL(numroc)(&mp0Param, &rowBlockSize, &MyRow, &iARow, &NPRow);
	int nQ0Param = colOfA + iCOff;
	int nQ0 = F77_CALL(numroc)(&nQ0Param, &colBlockSize, &MyCol, &iACol, &NPCol);
	int lWork = colBlockSize * (mP0 + nQ0 + colBlockSize);
	
	iMemSize = localSizeOfA + lTau + lWork;

	dpWork = (double *) malloc(sizeof(double) * iMemSize);
	memset(dpWork, 0xcc, sizeof(double) * iMemSize);

	D_Rprintf (("After allocating memory .. \n "));
	
	F77_CALL(callpdgeqrf)(dim, dpWork, &iMemSize);

	D_Rprintf (("AFTER FORTRAN FUNCTION EXECUTION \n "));

	free (dpWork);

	return 0;
}

/* **** CRSF_chol **** 
 * This function is a C interface to the fortran implemented 
 * scalapack driver function "callpdpotrf" that performs
 * Choleski factorization on the input data.
 */   
int CRSF_chol(int dim[], int iMyRank) {

	int iMemSize = 0;
	double *dpWork = NULL;

	int ipZero[] = { 0, 1, 2, 3 };
	int NPRow = dim[6];
	int NPCol = dim[7];
	int MyRow = iMyRank / NPCol;
	int MyCol = iMyRank % NPCol;

	int rowOfA = dim[0];
	int colOfA = dim[1];
	int rowBlockSize = dim[4];
	int colBlockSize = dim[5];

	/* Calculate required memory size */
	int localRowSizeOfA = F77_CALL(numroc)(&rowOfA, &rowBlockSize, &MyRow, ipZero, &NPRow);
	int localColSizeOfA = F77_CALL(numroc)(&colOfA, &colBlockSize, &MyCol, ipZero, &NPCol);
	int localSizeOfA = localRowSizeOfA * localColSizeOfA;
	int workSpace = max (rowBlockSize, colBlockSize);

	iMemSize = localSizeOfA + workSpace;
	
	dpWork = (double *) malloc(sizeof(double) * iMemSize);
	memset(dpWork, 0xcc, sizeof(double) * iMemSize);

	D_Rprintf (("After allocating memory .. \n "));
	
	F77_CALL(callpdpotrf)(dim, dpWork, &iMemSize);

	D_Rprintf (("AFTER FORTRAN FUNCTION EXECUTION \n "));

	free (dpWork);

	return 0;
}

/* **** CRSF_chol2inv **** 
 * This function is a C interface to the fortran implemented 
 * scalapack driver function "callpdpotri" that performs
 * inverting a matrix from its Choleski Factorization
 */ 
int CRSF_chol2inv(int dim[], int iMyRank) {

	int iMemSize = 0;
	double *dpWork = NULL;

	int ipZero[] = { 0, 1, 2, 3 };
	int NPRow = dim[6];
	int NPCol = dim[7];
	int MyRow = iMyRank / NPCol;
	int MyCol = iMyRank % NPCol;

	int rowOfA = dim[0];
	int colOfA = dim[1];
	int rowBlockSize = dim[4];
	int colBlockSize = dim[5];

	/* Calculate required memory size */
	int localRowSizeOfA = F77_CALL(numroc)(&rowOfA, &rowBlockSize, &MyRow, ipZero, &NPRow);
	int localColSizeOfA = F77_CALL(numroc)(&colOfA, &colBlockSize, &MyCol, ipZero, &NPCol);
	
	int localSizeOfA = localRowSizeOfA * localColSizeOfA;
	int workSpace = max (rowBlockSize, colBlockSize);

	iMemSize = localSizeOfA + workSpace;
	
	dpWork = (double *) malloc(sizeof(double) * iMemSize);
	memset(dpWork, 0xcc, sizeof(double) * iMemSize);

	D_Rprintf (("After allocating memory .. \n "));
	
	F77_CALL(callpdpotri)(dim, dpWork, &iMemSize);

	D_Rprintf (("AFTER FORTRAN FUNCTION EXECUTION \n "));

	free (dpWork);

	return 0;
}

/* **** CRSF_eigen **** 
 * This function is a C interface to the fortran implemented 
 * scalapack driver function "callpdsyevd" that performs
 * eigen value decomposition or Spectral decomposition 
 */
int CRSF_eigen(int dim[], int iMyRank) {

	int iMemSize = 0;
	double *dpWork = NULL;

	int ipZero[] = { 0, 1, 2, 3 };
	int NPRow = dim[6];
	int NPCol = dim[7];
	int MyRow = iMyRank / NPCol;
	int MyCol = iMyRank % NPCol;

	int rowOfA = dim[0];
	int colOfA = dim[1];
	int rowBlockSize = dim[4];
	int colBlockSize = dim[5];

	/* Calculate required memory size */
	int sizeOfLWORK;
	
	int tRILWMIN;
	int np, nq;
	
	int localRowSizeOfA = F77_CALL(numroc)(&rowOfA, &rowBlockSize, &MyRow, ipZero, &NPRow);
	int localColSizeOfA = F77_CALL(numroc)(&colOfA, &colBlockSize, &MyCol, ipZero, &NPCol);

	int localSizeOfA = localRowSizeOfA * localColSizeOfA;
	int localSizeOfZ = localRowSizeOfA * localColSizeOfA;
	int localSizeOfW = colOfA;
	int sizeOfLIWORK = 7 * colOfA + 8 * NPCol +2;


	
	np =  F77_CALL(numroc)(&colOfA, &colBlockSize, &MyRow, ipZero, &NPRow);
	nq =  F77_CALL(numroc)(&colOfA, &colBlockSize, &MyCol, ipZero, &NPCol);
	tRILWMIN = 3 * colOfA + max ( colBlockSize * (np +1) , 3 * colBlockSize);
	sizeOfLWORK = max (1 + 6 * colOfA + 2*np*nq, tRILWMIN ) + 2 * colOfA;


#if 0
	fprintf(stdout, "%d: localSizeOfA = %d, localSizeOfZ=%d, localSizeOfW = %d, sizeOfLIWORK=%d\n", iMyRank, localSizeOfA, localSizeOfZ, localSizeOfW, sizeOfLIWORK );fflush(stdout);
#endif


	iMemSize = localSizeOfA + localSizeOfZ + localSizeOfW + sizeOfLIWORK + sizeOfLWORK ;

	dpWork = (double *) malloc(sizeof(double) * iMemSize);

	if (dpWork == NULL)
	{
		fprintf(stdout, "%s:%d - Error allocating memory.\n", __FILE__, __LINE__ );
		fflush(stdout);

	}
	memset(dpWork, 0xcc, sizeof(double) * iMemSize);

	D_Rprintf (("After allocating memory .. \n "));
	
	F77_CALL(callpdsyevd)(dim, dpWork, &iMemSize);

	D_Rprintf (("AFTER FORTRAN FUNCTION EXECUTION \n "));

	free (dpWork);

	return 0;
}

/* **** CRSF_svd ****
 *  This function is a C interface to the fortran implemented 
 *  scalapack driver function "callpdgesvd" that performs
 *  singular value decomposition 
 */
int CRSF_svd(int dim[], int iMyRank) {
	int iMemSize = 0;
	double *dpWork = NULL;

	int ipZero[] = { 0, 1, 2, 3 };
	int NPRow = dim[6];
	int NPCol = dim[7];
	int MyRow = iMyRank / NPCol;
	int MyCol = iMyRank % NPCol;

	int rowOfA = dim[0];
	int colOfA = dim[1];
	int rowBlockSize = dim[4];
	int colBlockSize = dim[5];

	/* Calculate the required mem size */
	int localRowSizeOfA = 0;
	int localColSizeOfA = 0;
	
	int localSizeOfA = 0;
	int singularVectorSize = 0;
	int localSizeOfU = 0;
	int localSizeOfVT = 0;
	int localWorkSpaceSize = 0;

	D_Rprintf ((" In CRSF_svd func .. 1 \n"));

	/* Get the values */
	localRowSizeOfA = F77_CALL(numroc)(&rowOfA, &rowBlockSize, &MyRow, ipZero, &NPRow);
	localColSizeOfA = F77_CALL(numroc)(&colOfA, &colBlockSize, &MyCol, ipZero, &NPCol);

	localSizeOfA = localRowSizeOfA * localColSizeOfA;
	singularVectorSize = colOfA;
	localSizeOfU = localRowSizeOfA * localRowSizeOfA ;
	localSizeOfVT = localColSizeOfA * localColSizeOfA ;


	if (iMyRank == 0)
	{
		/* Calculate workspace size and Broadcast */
		int wATOBD =0;
		int ictxt;
		int mp0 = 0, nq0 = 0;
		int temp1=-1, temp2=0;
		int lldA = max (1, localRowSizeOfA);
		
		int wBDTOSVD;
		int size = min (rowOfA, colOfA);
		int sizeP = F77_CALL(numroc)(&size, &rowBlockSize, &MyRow, ipZero, &NPRow);
		int sizeQ = F77_CALL(numroc)(&size, &colBlockSize, &MyCol, ipZero, &NPRow);
		int wANTU = 1;
		int NRU = localRowSizeOfA;
		int wANTVT = 1;
		int NCVT = localColSizeOfA;
		int wBDSQR;
		int wPDORMBRQLN;
		int wPDORMBRPRT;

		D_Rprintf ((" Process 0 ... in cntrl statement ..after init \n"));

#if MPICH
		F77_CALL(blacs_get)(&temp1, &temp2, &ictxt);
#else
		F77_CALL(blacs_get)(&temp1, &temp2, &ictxt);
#endif
		/* Calculate the value for WATOBD */

		D_Rprintf ((" Process 0 ... in cntrl statement ..1 \n"));
		mp0 = F77_CALL(numroc)(&rowOfA, &rowBlockSize, &MyRow, &ictxt , &NPRow);
		nq0 = F77_CALL(numroc)(&rowOfA, &colBlockSize, &MyCol, &lldA , &NPRow);

		wATOBD = rowBlockSize * (mp0 + nq0 + 1) + nq0 ;

		D_Rprintf ((" Process 0 ... in cntrl statement ..2 \n"));
		/* Calculate the value for WBDTOSVD */
		
		wBDSQR = max (1, 2*size + (2*size - 4) * max (wANTU, wANTVT)) ;
		wPDORMBRQLN = max (colBlockSize *(colBlockSize -1)/2 , (sizeQ + mp0)*colBlockSize) + colBlockSize * colBlockSize;
		wPDORMBRPRT= max (rowBlockSize *(rowBlockSize -1)/2 , (sizeP + nq0)*rowBlockSize) + rowBlockSize * rowBlockSize;

		wBDTOSVD = size * (wANTU * NRU + wANTVT * NCVT) + max (wBDSQR,max (wANTU * wPDORMBRQLN , wANTVT * wPDORMBRPRT));

		localWorkSpaceSize = 2 + 6* max( rowOfA, colOfA) + max (wATOBD, wBDTOSVD);

		D_Rprintf ((" Process 0 ... in cntrl statement ..3 \n"));

	}

	D_Rprintf (("  after cntrl statement ..rank = %d  \n", iMyRank));
	
	if (MPI_Bcast (&localWorkSpaceSize, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS){
		Rprintf ("SVD: Failed during MPI_Bcast call ... Rank: %d Exiting \n", iMyRank); 
		exit (-1);
	} 
	
	iMemSize = localSizeOfA + singularVectorSize + localSizeOfU + localSizeOfVT + localWorkSpaceSize;
	
	D_Rprintf ((" After Broadcasting iMemSize = %d .. localWspace = %d ... rank = %d \n", iMemSize,localWorkSpaceSize, iMyRank));
	
	dpWork = (double *) malloc(sizeof(double) * iMemSize);
	memset(dpWork, 0xcc, sizeof(double) * iMemSize);

	D_Rprintf (("After allocating memory .. \n "));
	
	F77_CALL(callpdgesvd)(dim, dpWork, &iMemSize);

	D_Rprintf (("AFTER FORTRAN FUNCTION EXECUTION \n "));

	free(dpWork);

	return 0;

}


/* **** CRSF_solve ****
 *  This function is a C interface to the fortran implemented 
 *  scalapack driver function "callpdgesv" that solves the 
 *  equations A * X = B for 'X'
 */
int CRSF_solve(int dim[], int iMyRank) {

	int iMemSize = 0;
	double *dpWork = NULL;

	int ipZero[] = { 0, 1, 2, 3 };
	int NPRow = dim[6];
	int NPCol = dim[7];
	int MyRow = iMyRank / NPCol;
	int MyCol = iMyRank % NPCol;

	int rowOfA = dim[0];
	int colOfA = dim[1];
	int rowOfB = dim[2];
	int colOfB = dim[3];
	int rowBlockSize = dim[4];
	int colBlockSize = dim[5];

	/* Calculate the MemSize */	
	int localRowSizeOfA = F77_CALL(numroc)(&rowOfA, &rowBlockSize, &MyRow, ipZero, &NPRow);
	int localColSizeOfA = F77_CALL(numroc)(&colOfA, &colBlockSize, &MyCol, ipZero, &NPCol);

	int localRowSizeOfB = F77_CALL(numroc)(&rowOfB, &rowBlockSize, &MyRow, ipZero, &NPRow);
	int localColSizeOfB = F77_CALL(numroc)(&colOfB, &colBlockSize, &MyCol, ipZero, &NPCol);


	int localSizeOfA = localRowSizeOfA * localColSizeOfA;
	int localSizeOfB = localRowSizeOfB * localColSizeOfB;
	int localSizeOfPiv = localRowSizeOfA + rowBlockSize;
	int localWorkSize = colBlockSize; 

	D_Rprintf(("Child: inside CRSF_solve \n "));
	

	iMemSize = localSizeOfA + localSizeOfB + localSizeOfPiv + localWorkSize;

	dpWork = (double *) malloc(sizeof(double) * iMemSize);
	memset(dpWork, 0xcc, sizeof(double) * iMemSize);

	D_Rprintf(("Child: After allocating memory .. \n "));
	
	F77_CALL(callpdgesv)(dim, dpWork, &iMemSize);

	D_Rprintf (("Child: AFTER FORTRAN FUNCTION EXECUTION \n "));

	free (dpWork);

	return 0;
}

/* **** CRSF_multiply ****
 *  This function is a C interface to the fortran implemented 
 *  scalapack driver function "callpdgemm" that performs matrix
 *  multiplication.
 */
int CRSF_multiply(int dim[], int iMyRank) {

	int iMemSize;
	/*int ii;*/
	int localRowSizeOfA;
	int localColSizeOfA;
	int localRowSizeOfB;
	int localColSizeOfB;
	int localRowSizeOfC;
	int localColSizeOfC;
	double *dpWork = NULL;

	int ipZero[] = { 0, 1, 2, 3 };
	int NPRow = dim[6];
	int NPCol = dim[7];
	int MyRow = iMyRank / NPCol;
	int MyCol = iMyRank % NPCol;

	int rowOfA = dim[0];
	int colOfA = dim[1];
	int rowOfB = dim[2];
	int colOfB = dim[3];
	int rowBlockSize = dim[4];
	int colBlockSize = dim[5];

	int iLLD_A = 0;
	int iLLD_B = 0;
	int iLLD_C = 0;

	localRowSizeOfA = F77_CALL(numroc)(&rowOfA, &rowBlockSize, &MyRow, ipZero, &NPRow);
	localColSizeOfA = F77_CALL(numroc)(&colOfA, &colBlockSize, &MyCol, ipZero, &NPCol);

	localRowSizeOfB = F77_CALL(numroc)(&rowOfB, &rowBlockSize, &MyRow, ipZero, &NPRow);
	localColSizeOfB = F77_CALL(numroc)(&colOfB, &colBlockSize, &MyCol, ipZero, &NPCol);

	localRowSizeOfC = localRowSizeOfA;
	localColSizeOfC = localColSizeOfB;

	if ( localColSizeOfA > 0 )
		iLLD_A = max( 1, localRowSizeOfA );
	else
		iLLD_A = 1;


	if ( localColSizeOfB > 0 )
		iLLD_B = max( 1, localRowSizeOfB );
	else
		iLLD_B = 1;

	if ( localColSizeOfC > 0 )
		iLLD_C = max( 1, localRowSizeOfC );
	else
		iLLD_C = 1;


	/* I am calculating memory for only three matrices A, B & C(output). */

	iMemSize = iLLD_A * localColSizeOfA + iLLD_B * localColSizeOfB + iLLD_C * localColSizeOfC + colBlockSize + 1;

	D_Rprintf (( "%d/%d: Mem size for me = %d\n", MyRow, MyCol, iMemSize ));



	dpWork = (double *) malloc(sizeof(double) * iMemSize);
	memset(dpWork, 0xcc, sizeof(double) * iMemSize);

	D_Rprintf (("%d:  Allocated a work vector of size %d bytes\n",
				iMyRank, sizeof(double) * iMemSize));

	F77_CALL(callpdgemm)(dim, dpWork, &iMemSize);

	free(dpWork);
	return 0;
}


/* ****  End CRSFs  -- Child R-ScaLAPACK Functions ****  */
