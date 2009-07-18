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
/* Driver program for the Child-R-ScaLAPACK process.
 * By design, this program is spawned by MPI_Comm_Spawn.  The main function
 * simply calls the CR_Exec function - the child process's equivelant to the
 * main function.
 */
#include "CRscalapack.h"

int main() {

	return (int) CR_Exec();
	
	}

