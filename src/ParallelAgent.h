/*======================================================================
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
=========================================================================*/
#ifndef _PARALLELAGENT_H_
#define _PARALLELAGENT_H_
#include "mpi.h"
#define USE_RINTERNALS
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include <R_ext/PrtUtil.h>

#include <dlfcn.h>


#define DONT_SPAWN_R

#ifdef DONT_SPAWN_R
#define Rprintf printf
#endif

#define CHECKFAULT_TAG 1202	/* Tag used to check fault before data send  */
#define NUMMATRICES_TAG 202	/* Tag used to number of Results to collect  */
#define RECVRESULT_TAG 300	/* Tag used to receive results */
#define RECVVECORMAT_TAG 400	/* Tag used to communicate, whether to receive vector or matrix*/


#define MPINULL MPI_STATUS_IGNORE

#if 0
#define DEBUG_RSCALAPACK
#endif

#ifndef DEBUG_RSCALAPACK
#define D_Rprintf(x)
#else
#define D_Rprintf(x) Rprintf x
#endif

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define iceil(a,b) ((a+b-1)/b)

int PA_ErrorHandler(int errcode);
int erreturn(int errcode);

void PAdistData(double *, int *, int , int);
void PAcollectData(double *, int *, int, int);

/*
Performance degraded with non-blocking send and receive
void PA_SendVectorToCR (int *ib, int *ia, double *work, int *mb, int *s2rank, MPI_Request *request);
void PA_RecvVectorFromCR (int *ib, int *ia, double *A, int *mb, int *fromRank, MPI_Request *request);
*/

void PA_SendVectorToCR (int *ib, int *ia, double *work, int *mb, int *s2rank);
void PA_RecvVectorFromCR (int *ib, int *ia, double *A, int *mb, int *fromRank);

int PA_Init();
SEXP PA_Exit();
SEXP PA_SpawnProcs(SEXP NumProcs, SEXP scriptLocn );
int PA_UnpackInput(SEXP sxInputVector, int *ipDims, double **dppA,
        double **dppB, int *ipNumProcs, int *ipFunction, int *ipRelFlag);
SEXP PA_Exec(SEXP scriptLocn, SEXP sxInputVector);
int PA_SendData(int [], double [], double []);
SEXP PA_RecvResult(int []);
SEXP PA_GridInfo();

int PA_GetTwoDims(SEXP , int *);
int PA_SetDim(SEXP ,int , int *);
int PA_CheckFaultPriorRun();

#endif
