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
#ifndef _PARALLELAGENT_H_
#define _PARALLELAGENT_H_
#include "mpi.h"
#define USE_RINTERNALS
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include <R_ext/PrtUtil.h>

#define DONT_SPAWN_R

#define PC_SEND_DATA_DIM 501

#define MPINULL MPI_STATUS_IGNORE
#define PA_RecvVectorFromCR parecvvectorfromcr_
#define PA_SendVectorToCR pasendvectortocr_

#define max(a,b) ((a) > (b) ? (a) : (b))

int PA_ErrorHandler(int errcode);
int erreturn(int errcode);

void F77_NAME(padistdata)(double *, int *, int *, int *);
void F77_NAME(pacollectdata)(double *, int *, int *, int *, int *);

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

int PA_GetTwoDims(SEXP , int *);
int PA_SetDim(SEXP ,int , int *);

#endif
