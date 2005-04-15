/*======================================================================
 *           R-ScaLAPACK version 0.4.x:  ScaLAPACK interface to R
 *              Oak Ridge National Laboratory, Oak Ridge TN.
 *        Authors: David Bauer, Guruprasad Kora, Nagiza. F. Samatova, 
 *                            Srikanth Yoginath.
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
 * RScaLAPACK (http://www.aspect-sdm.org/Parallel-R) was funded
 * as part of the Scientific Data Management Center
 * (http://sdm.lbl.gov/sdmcenter) under the Department of Energy's 
 * Scientific Discovery through Advanced Computing (DOE SciDAC) program
 * (http://www.scidac.org ). 
=========================================================================*/
/*  ParallelAgent.c   The C code for the ParallelAgent of the R-ScaLAPACK
 *  package.
 */
#include "ParallelAgent.h"

/*#define DEBUG_RSCALAPACK*/

#ifndef DEBUG_RSCALAPACK
#define D_Rprintf(x)
#else
#define D_Rprintf(x) Rprintf x
#endif

#define max(a,b) ((a) > (b) ? (a) : (b))

static MPI_Comm childComm;
static int iGlobalNumChildren = 0;

/* ****  PA_Init ****
 * Check to see if MPI is initialized, and if it is not, then initialize it.
 * Returns 0 for success or non-zero on error.
 */
int PA_Init() {
	int flag;

	if( MPI_Initialized(&flag) != MPI_SUCCESS){
		Rprintf("ERROR[1]: Failed in call MPI_Initialized \n");
		return 1;
	}

	if (flag){
		D_Rprintf((" MPI has already been Initialized \n"));
		return 0;
	} else {
		MPI_Init(NULL, NULL);
		MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
		return 0;
	}
}

/* **** PA_Exit ****
 *  Calls MPI_Finalize.  The only time this function should be called is
 *  when R is unloading the R-ScaLAPACK library and exiting.
 */
SEXP PA_Exit() {
	/*MPI_Comm_free(&childComm);*/
	MPI_Finalize();
	return R_NilValue;
}

/* **** PA_Exec ****
 * The ParallelAgent's equivelant of "main".  This function spawns the
 * children processes, sends them the data, and gets back the results.
 */
SEXP PA_Exec(SEXP scriptLocn, SEXP sxInputVector) {
	int iFunction;
	int iSpawnFlag = 1;
	int iNumProcs;
	int ipDims[10]= { 0,0,0,0,0,0,0,0,0,0 };
	double *dpA = NULL;
	double *dpB = NULL;
	int iStatus;
	int returnValue;
	SEXP sRet;
#ifndef DONT_SPAWN_R
	char *cpProgram = "R";   /* Program to call; Usually "R" */

	char *child_args[] = {
		"BATCH",
		"--no-save",
		CHAR(STRING_ELT((scriptLocn),0)),
		"abc.out",
		NULL };

#else
	/* There are four ways to call the child processes (all go through
	 * MPI_COMM_SPAWN):
	 * 1.  Spawn a copy of R and let it call the child process function through
	 *        an R script.  This approach simplified the MPI-BLACS interation,
	 *        at the cost of major overhead.
	 * 2.  Spawn a shell script (dfCRDriver) which runs the child process
	 *        driver.  Spawning the program directly is better, but there is a
	 *        problem with needing the shared library in the library path.
	 * 3.  Spawn a shell (sh) which first adds the path to the scalapack.so
	 *        library to LD_LIBRARY_PATH and then calls the driver program.
	 *        This eliminates the need for the extra script (dfCRDriver), while
	 *        also being several thousandths of a second faster.
	 * 4.  Spawn the driver program directly.  This requires building the
	 *        program differently so that it doesn't need to scalapack.so
	 *        shared library.  This is the preferred/current method.
	 */
	char *cpProgram;    /* =R_PACKAGE_DIR"/exec/dfCRDriver.sh"; */
	char *child_args[] = { NULL, NULL };
	int iLength;

	/* The scriptLocn (script location) variable contains the path to the
	 * executable directory (followed by the script name).  Extract the path,
	 * and use it for the executable's path.
	 */
	cpProgram = CHAR(STRING_ELT((scriptLocn), 0));
	iLength = strrchr(cpProgram, '/') - cpProgram;
	if (iLength < 0) {
		Rprintf("Path to script is not complete.  Unable to continue.\n");
		return R_NilValue;
	}
	cpProgram = (char *) malloc(sizeof(char) * (iLength + 12));
	if (cpProgram == NULL) {
		Rprintf("Memory allocation (%d bytes) failed!\n", sizeof(char) *
				(iLength + 12));
		return R_NilValue;
	}
	*(cpProgram) = '\0';
	strncat(cpProgram, CHAR(STRING_ELT((scriptLocn),0)), iLength);
	strncat(cpProgram, "/CRDriver", 10);

	D_Rprintf(("Child process: \"%s\" \"%s\" \"%s\"\n", cpProgram, child_args[0], child_args[1]));
#endif    /* Endof  If DONT_SPAWN_R is defined */

	/*  Begin by unpacking the input vector into all of the seperate variables
	 *  that get their values from it */
	if (PA_UnpackInput(sxInputVector, ipDims, &dpA, &dpB, &iNumProcs,
				&iFunction, &iSpawnFlag) != 0) {
		free(cpProgram);
		return R_NilValue;
	}

	/*  Initialize MPI (if it is already initialized, it won't be
	 *  initialized again).	*/
	if (PA_Init() != 0){
		Rprintf(" ERROR[1]: Failed while intializing MPI \n");
		free(cpProgram);
		return R_NilValue;
	}

	if (iSpawnFlag != 0  && iGlobalNumChildren != 0) {
		Rprintf(" Error:  Attempt to spawn a new grid without releasing the previous grid.\n");
		return R_NilValue;
	}

	if(iSpawnFlag == 0 && iGlobalNumChildren == 0){
		Rprintf(" Error: Process Grid not present and Spawn option is set FALSE \n");
		return R_NilValue;
	}

	/* Begin:  Spawn the child processes */
	if ( iSpawnFlag != 0 ) {
		/* Begin:  Spawn the child processes */
		D_Rprintf(("Preparing to spawn %d child processes.\n", iNumProcs));
		iStatus = MPI_Comm_spawn(cpProgram, child_args, iNumProcs, MPI_INFO_NULL,
				0, MPI_COMM_WORLD, &childComm, MPI_ERRCODES_IGNORE);
		free(cpProgram);
		if (iStatus != MPI_SUCCESS) {
			Rprintf(" ERROR:  Failed to spawn (%d) child processes.\n", iNumProcs);
			return R_NilValue;
		}

		D_Rprintf(("SPAWNING SUCCESSFUL\n"));
		/* End:  Spawn the child processes */
		iGlobalNumChildren = iNumProcs;
	}


	/* SPECIAL for SVD */
	/* If the function is SVD, the child process needs to know the nu,nv
	 * parameters.
	 */

	if (iFunction == 2) {
		ipDims[2] = (int) dpB[0];
		ipDims[3] = (int) dpB[1];
	}
	
	/* DATA DISTRIBUTION */
	/* The data is distributed by the PA to all of the child processes. */
	if ((returnValue = PA_SendData(ipDims, dpA, dpB)) == 0)	{
		D_Rprintf (("SUCCESS[1]: DATA SENT TO CHILD PROCESSES.\n"));
	} else {	/* The send data failed, */ 
		Rprintf("ERROR [1] : DATA COULD NOT BE SENT TO CHILD PROCESSES.\n");
		return R_NilValue;
	}

	/* If the release flag == 1, then the grid will be (or was) released. */
	if (ipDims[9] == 1) iGlobalNumChildren = 0;

	/* If the function is sla.gridInit or sla.gridExit, just return. */
	if (iFunction == 0) return R_NilValue;

	/* GET BACK THE RESULT */
	sRet = PA_RecvResult(ipDims);

	return sRet;
}

/* ****  PA_UnpackInput  ****
 * The input parameters/data to the ParallelAgent is given as an R vector.
 * This function unpacks the vector, putting the pieces into the appropriate
 * variables.
 */
int PA_UnpackInput(SEXP sxInputVector, int *ipDims, double **dppA,
		double **dppB, int *ipNumProcs, int *ipFunction, int *ipSpawnFlag) {
	SEXP s;
	int iMB;
	int ipReleaseFlag;

	/* First parameter is the first matrix. --populates ipDims[0] and ipDims[1] */
	s = VECTOR_PTR(sxInputVector)[0];
	if (TYPEOF(s) != REALSXP) {
		Rprintf("1st parameter (Matrix A) is not an array of doubles.\n");
		return -1;
	}

	if (PA_GetTwoDims(s, ipDims) > 2) {
		Rprintf("1st parameter (Matrix A) has too many dimensions.\n");
		return -2;
	}

	/* If the object is one dimensional, set the second dimension to length 1 */
	if (ipDims[1] == 0) ipDims[1] = 1;

	/* Save a pointer to the first matrix */
	*dppA = REAL(s);

	/* Second parameter is the second matrix. --populates ipDims[2] and ipDims[3] */
	s = VECTOR_PTR(sxInputVector)[1];
	if (TYPEOF(s) != REALSXP) {
		Rprintf("2nd parameter (Matrix B) is not an array of doubles.\n");
		return -3;
	}

	if (PA_GetTwoDims(s, ipDims + 2) > 2) {
		Rprintf("2nd parameter (Matrix B) has too many dimensions.\n");
		return -4;
	}

	/* If the object is one dimensional, set the second dimension to length 1,
	 * unless the length is zero (i.e, there is no second matrix).    */
	if (ipDims[3] == 0 && LENGTH(s) != 0) ipDims[3] = 1;

	/* Save a pointer to the second matrix */
	*dppB = REAL(s);

	/* Third Parameter is number of Process Rows  -- populates ipDims[6] = NPROWS*/

	s = VECTOR_PTR(sxInputVector)[2];
	if (TYPEOF(s) != INTSXP) {
		Rprintf("Third parameter (number of row processors) is not an integer.\n");
		return -5;
	}
	if (LENGTH(s) != 1) {
		Rprintf("First parameter (number of row processors) is not a single number.\n");
		return -6;
	}
	ipDims[6] = INTEGER(s)[0];

	/* Fourth Parameter is number of Process Cols -- populates ipDims[7] = NPCOLS */
	s = VECTOR_PTR(sxInputVector)[3];
	if (TYPEOF(s) != INTSXP) {
		Rprintf("Fourth parameter (number of col processors) is not an integer.\n");
		return -7;
	}
	if (LENGTH(s) != 1) {
		Rprintf("Fourth parameter (number of col processors) is not a single number.\n");
		return -8;
	}
	ipDims[7] = INTEGER(s)[0];

	*ipNumProcs = ipDims[6] * ipDims[7];  /* iNumProcs = iNPRows * iNPCols  */

	/* Fifth Parameter is Block Size -- populates ipDims[4] =ipDims [5]= MB */

	s = VECTOR_PTR(sxInputVector)[4];
	if (TYPEOF(s) != INTSXP) {
		Rprintf("Fifth parameter (row block size of LHS matrix) is not an integer.\n");
		return -9;
	}
	if (LENGTH(s) != 1) {
		Rprintf("Fifth parameter (row block size of LHS matrix) is not a single number.\n");
		return -10;
	}
	iMB = INTEGER(s)[0];
	/* !!! If the block size is larger than the input size, reduce the block
	 * size down to the size of the input.  Needed for the eigen function.
	 */
	if (iMB > ipDims[0] && iMB > ipDims[1] && iMB > ipDims[2] &&
			iMB > ipDims[3])
		iMB = max(ipDims[0], max(ipDims[1], max(ipDims[2], ipDims[3])));


	/* Once upon a time, the block dimensions were taken as two inputs.  Now,
	 * the blocksize is forced to be square.
	 */
	ipDims[4] = ipDims[5] = iMB;

	/* Sixth Parameter is function identifier  -- populates ipDims[8] = functionID */
	s = VECTOR_PTR(sxInputVector)[5];
	if (TYPEOF(s) != INTSXP) {
		Rprintf("Sixth parameter (function) is not an integer.\n");
		return -11;
	}
	if (LENGTH(s) != 1) {
		Rprintf("Sixth parameter (function) is not a single number.\n");
		return -12;
	}

	*ipFunction = INTEGER(s)[0];
	if (*ipFunction < 0 || *ipFunction > 7) {
		Rprintf("Error:  Unknown function ID (%d).\n", *ipFunction);
		return -13;
	}

	ipDims[8]= *ipFunction;

	/* Seventh parameter --- Instruction to Release Grid or not -- populates ipDims[9]=ReleaseFlag */
	/* Populating ipDims array is complete ... */
	s = VECTOR_PTR(sxInputVector)[6];
	if (TYPEOF(s) != INTSXP) {
		Rprintf("Seventh parameter (function) is not an integer.\n");
		return -11;
	}

	ipReleaseFlag = INTEGER(s)[0];
	if(ipReleaseFlag != 0 && ipReleaseFlag != 1) {
		Rprintf ("Warning: Proper value for Release Flag= %d not used \n \t Release Flag is set to 1 \n", ipReleaseFlag );
		ipReleaseFlag = 1;  
	}

	ipDims[9]= ipReleaseFlag;

	/* ipSpawnFlag is the information sent from R to ParallelAgent requesting
	 * it either to spawn the processes or not to spawn (use the existing ones)
	 */

	/* Eighth parameter --- Instruction to Spawn Grid or not*/
	s = VECTOR_PTR(sxInputVector)[7];
	if (TYPEOF(s) != INTSXP) {
		Rprintf("Sixth parameter (function) is not an integer.\n");
		return -11;
	}
	*ipSpawnFlag = INTEGER(s)[0];

	D_Rprintf(("spawn flag : %d \n", *ipSpawnFlag));

	return 0;
}

/* ****  PA_SendData  ****
 *  This function sends the parameters to child 0 and then distributes the
 *  data to all of the child processes.  The data is distributed in the
 *  2D block/cyclic pattern required by the ScaLAPACK library.
 */
int PA_SendData(int ipDims[], double dpA[], double dpB[]) {
	int iWorksize, iMatrixNum;
	int iFunction;
	MPI_Comm intercomm;

	iFunction = ipDims[8];

	/* Broadcast the Data Dimension to the child Processes */
	PA_ErrorHandler(MPI_Intercomm_merge(childComm, 0, &intercomm));
	PA_ErrorHandler(MPI_Bcast(ipDims,10,MPI_INT, 0, intercomm)); 

	/* If the function was sla.gridInit or sla.gridExit, then there is
	 * no data to distribute.
	 */
	if (iFunction == 0) {
		return 0;
	} else {	/* Otherwise, distribute data as usual. */

		iWorksize = ipDims[5];
		iMatrixNum = 1;
		F77_CALL (padistdata)(dpA,ipDims,&iWorksize, &iMatrixNum);

		if (ipDims[2] != 0 && ipDims[8] != 2){
			iMatrixNum = 2;
			F77_CALL (padistdata)(dpB,ipDims,&iWorksize, &iMatrixNum);
		}

		return 0;
	}
}
/* ****  PA_RecvResult  ****
 *  After the computations are done, receive the result back from the children
 *  processes.  The parameters of the return values are sent by child 0, and
 *  the data is collected (from the 2D block/cyclic distribution) from all
 *  child processes.
 */
SEXP PA_RecvResult(int ipDims[])
{
	int iStatus;
	int iNumMatrices, iWorksize;
	int ipOutDims[3];
	SEXP sxRet = R_NilValue;
	SEXP sxTmp;
	int i;

	iWorksize = ipDims[5];
	D_Rprintf(("Preparing to receive .... "));

	/* First, get the number of matrices in the return object. */
	iStatus = MPI_Recv(&iNumMatrices,1,MPI_INT,0,202,childComm,MPINULL);
	if (iStatus != MPI_SUCCESS) {
		PA_ErrorHandler(iStatus);
		return R_NilValue;
	}

	D_Rprintf(("%d matrices\n", iNumMatrices));

	/* Setup the R vector to hold the R matrix objects */
	PROTECT( sxRet = allocVector(VECSXP, iNumMatrices) );

	/* Main loop */
	for (i = 0; i < iNumMatrices; i++) {
		/* Get the dimensions of matrix i */
		iStatus = MPI_Recv(ipOutDims, 3, MPI_INT, 0, 300+i, childComm, MPINULL);
		if (iStatus != MPI_SUCCESS) {
			PA_ErrorHandler(iStatus);
			UNPROTECT( 1 );
			return R_NilValue;
		}

		D_Rprintf(("Dimensions of matrix %d: [%d, %d], tag=%d\n", i,
					ipOutDims[1], ipOutDims[2], ipOutDims[0]));

		/* If the matrix is empty (size zero), handle the case seperately.
		 * SVD can return empty matrices.
		 */
		if (ipOutDims[1] == 0 || ipOutDims[2] == 0) {
			sxTmp = coerceVector(R_NilValue, NILSXP);
			SET_VECTOR_ELT(sxRet, i, sxTmp);
			continue;
		}

		/* If the matrix is actually a vector, handle that case. */
		if (ipOutDims[1] == 0) ipOutDims[1] = 1;
		if (ipOutDims[2] == 0) ipOutDims[2] = 1;

		/* Make an R array to hold the matrix */
		PROTECT( sxTmp = allocVector(REALSXP, ipOutDims[1]*ipOutDims[2]) );

		D_Rprintf(("Preparing to receive..."));
		/* Get matrix i */
		/* ipOutDims[0] is a flag saying whether to receive the matrix via
		 * a normal MPI call or via the PACollectData function */
		if (ipOutDims[0] == 1) {
			iStatus = MPI_Recv(REAL(sxTmp), ipOutDims[1]*ipOutDims[2],
					MPI_DOUBLE,	0, 400+i, childComm, MPINULL);
			if (iStatus != MPI_SUCCESS) {
				PA_ErrorHandler(iStatus);
				return R_NilValue;
			}
		} else {
			F77_CALL(pacollectdata)(REAL(sxTmp), ipOutDims+1, ipOutDims+2,
					ipDims, &iWorksize);
		}

		D_Rprintf(("Received.\n"));
		/* Set the dimensions of the matrix */
		PA_SetDim(sxTmp, 2, ipOutDims+1);
		/* Save the matrix into the vector */
		SET_VECTOR_ELT(sxRet, i, sxTmp);
		UNPROTECT( 1 );
	}

	UNPROTECT( 1 );

	return sxRet;
}

/* ****  PA_GetTwoDims  ****
 * Get the dimensions of a matrix (SEXP array).
 * If there are more than 2 dimensions, stop.
 * Returns the number of dimensions, or a negative integer in case of
 * error.
 */
int PA_GetTwoDims(SEXP s, int *ipDims) {
	int i, j;
	SEXP sxTmp = getAttrib(s, R_DimSymbol);

	if (sxTmp != R_NilValue) {
		if (TYPEOF(sxTmp) != INTSXP) {
			Rprintf("Error: Dim tag did not contain an integer array.\n");
			return -1;
		}

		i = LENGTH(sxTmp);
		/* If there are more than 2 dimensions, then quit. */
		if (i > 2) return i;
		for (j = 0; j < i; j++) {
			ipDims[j] = INTEGER(sxTmp)[j];
		}

		return i;
	}

	/* No dimensions tag, so it is one dimensional: use the length. */
	ipDims[0] = LENGTH(s);

	return 1;
}

/*  Given a SEXP object, the number of dimensions, and the array of dimensions,
 *  this function sets the dimensions of the SEXP object to the given
 *  dimensions.
 *  Returns 0 upon success and nonzero upon failure.
 */

/* ****  PA_SetDim  ****
 * Given an R array and dimensions, set the dimensions of the R array to the
 * dimensions given.
 */
int PA_SetDim(SEXP s, int iNumDim, int *ipMyDims) {
	int j;
	int iSum;
	SEXP sxDim;

	if (s != R_NilValue) {
		if (TYPEOF(s) != INTSXP && TYPEOF(s) != REALSXP) {
			Rprintf("Error:  Cannot give dimensions to non-array object.\n");
			return -1;
		}

		iSum = 1;
		for (j = 0; j < iNumDim; j++) iSum *= ipMyDims[j];

		if (iSum != LENGTH(s)) {
			Rprintf("Error:  Dimensions do not fit length of object.\n");
			return -2;
		}

		PROTECT(sxDim = allocVector(INTSXP, iNumDim));

		for (j = 0; j < iNumDim; j++) {
			INTEGER(sxDim)[j] = ipMyDims[j];
		}

		setAttrib(s, R_DimSymbol, sxDim);
		UNPROTECT( 1 );
	}

	return 0;
}

/* ****  PA_SendVectorToCR ****
 * MPI/C clone of the BLACS function to send a matrix, except it uses the
 * global variable childComm as the childCommunicator.
 */
/* Note:  PA_SendVectorToCR is translated to "pasend_" for Fortran compatibility */
void PA_SendVectorToCR(int *ib, int *ia, double *work, int *mb, int *s2rank) {
	MPI_Datatype GEMAT;

	PA_ErrorHandler(MPI_Type_vector(*ia, *ib, *mb, MPI_DOUBLE, &GEMAT));
	PA_ErrorHandler(MPI_Type_commit (&GEMAT));

	PA_ErrorHandler(MPI_Send ( work, 1, GEMAT, *s2rank, 5000 , childComm));

	PA_ErrorHandler(MPI_Type_free (&GEMAT));
}

/* ****  PA_RecvVectorFromCR  ****
 * MPI/C clone of the BLACS function to receive a matrix, except it uses the
 * global variable childComm as the intercommunicator.
 */
/* Note:  PA_RecvVectorFromCR is translated to "parecv_" for Fortran compatibility */
void PA_RecvVectorFromCR (int *ib, int *ia, double *A, int *mb, int *fromRank) {
	MPI_Datatype GEMAT;

	PA_ErrorHandler(MPI_Type_vector(*ia, *ib, *mb, MPI_DOUBLE, &GEMAT));
	PA_ErrorHandler(MPI_Type_commit (&GEMAT));

	PA_ErrorHandler(MPI_Recv (A, 1, GEMAT, *fromRank ,15000 ,childComm, MPINULL));

	PA_ErrorHandler(MPI_Type_free (&GEMAT));
}


/* **** PA_ErrorHandler ****
 * PA_ErrorHandler sends back the MPI error messages on their occurrence  
 */ 
int PA_ErrorHandler(int errcode)
{
	int errmsglen;
	char errmsg[MPI_MAX_ERROR_STRING];

	if (errcode != MPI_SUCCESS) {
		MPI_Error_string(errcode, errmsg, &errmsglen);
		/* error(errmsg); */
		Rprintf("MPI Error: \"%s\"\n", errmsg);
	}

	return errcode;
}
