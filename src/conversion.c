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
#include "ParallelAgent.h"
#include <sys/types.h>
#include <unistd.h>

SEXP AsInt (int x)
{
	SEXP sexp_x;
	PROTECT (sexp_x = allocVector (INTSXP, 1));
	INTEGER (sexp_x)[0] = x;
	UNPROTECT (1);
	return sexp_x;
}

/*AsIntArray(int *x, int length)
 * {
 * 	int i;
 * 	SEXP sexp_x;
 * 	PROTECT (sexp_x = allocVector(INTSXP, length));
 *
 * 	for ( i=0; i<length;i++)
 * 		INTEGER (sexp_x)[0] = x[i];
 *
 * 	return sexp_x;
 * }
 *
 *AsRealArray(double *cdouble, int length)
 *{
 * 	int i;
 * 	SEXP sexp_cdouble;
 * 	PROTECT (sexp_cdouble = allocVector(INTSXP, length));
 *
 * 	for ( i=0; i<length;i++)
 * 		INTEGER (sexp_cdouble)[0] = cdouble[i];
 *
 * 	return sexp_cdouble;
 * }
 *
 * */

SEXP AsReal (double cdouble)
{
	SEXP sexp_cdouble;
	PROTECT (sexp_cdouble = allocVector(REALSXP, 1));
	REAL (sexp_cdouble)[0] = cdouble;
	UNPROTECT(1);
	return sexp_cdouble;
}

