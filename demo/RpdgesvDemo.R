# ====================================================================
#           R-ScaLAPACK version 0.4.x:  ScaLAPACK interface to R
#              Oak Ridge National Laboratory, Oak Ridge TN.
#        Authors: David Bauer, Guruprasad Kora, Nagiza. F. Samatova, 
#                            Srikanth Yoginath.
#     Contact: Nagiza F. Samatova; (865) 241-4351; samatovan@ornl.gov
#                 Computer Science and Mathematics Division
#             Oak Ridge National Laboratory, Oak Ridge TN 37831 
#                   (C) 2004 All Rights Reserved
#
#                              NOTICE
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby granted
# provided that the above copyright notice appear in all copies and
# that both the copyright notice and this permission notice appear in
# supporting documentation.
#
# Neither the Oak Ridge National Laboratory nor the Authors make any
# representations about the suitability of this software for any
# purpose.  This software is provided ``as is'' without express or
# implied warranty.
#
# RScaLAPACK (http://www.aspect-sdm.org/Parallel-R) was funded
# as part of the Scientific Data Management Center
# (http://sdm.lbl.gov/sdmcenter) under the Department of Energy's 
# Scientific Discovery through Advanced Computing (DOE SciDAC) program
# (http://www.scidac.org ). 
# ======================================================================
RpdgesvDemo <- function () {

	library(RScaLAPACK)

        print ("****************************************************")
        print ("            DEMO: SOLVE SYSTEM EQNS                 ")
        print ("****************************************************")

	x <- c(6 ,3 ,0 ,0 ,3 , 0, 0, -3, -1, 1, 1, 0, -1, 0, 11, 0, 0, 10, 0, 0, 0, -11, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 8, 0, -10)

	dim(x) <- c(6,6)

	y <- c( 72, 0, 160, 0, 0, 0)

	dim(y) <- c(6,1)

	print (x)

	print (y)

	print (" SOLVING using SCALAPACK FUNCTION ")

	print (" PROCESS GRID SPECIFICATIONS .. in order ")

	print (" NPROWS = 2  NPCOLS = 2 ")

	print (" BLOCK SIZE LHS:  MB = 2   NB = 2 ")

	print (" CALLING SCALAPACK FUNCTION .. sla.solve(x,y,2,2,2) ")

	a <- sla.solve(x,y,2,2,2)

	print(a)

	print (" SOLVING USING CORRESPONDING LAPACK FUNCTION .. solve(x,y) ")

	b <- solve(x,y)

	print (b)

	print (" COMPARE BOTH RESULTS !! ")

	print (" EXPERIMENT WITH DIFFERENT PROCESS GRID PARAMETERS ")

       print ("****************************************************")

}

