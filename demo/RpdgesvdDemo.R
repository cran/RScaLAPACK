#  =======================================================================
#           R-ScaLAPACK version 0.3.x:  R interface to ScaLAPACK
#              Oak Ridge National Laboratory, Oak Ridge TN.
#      Authors: David Bauer, Nagiza. F. Samatova, Srikanth Yoginath
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
# R-ScaLAPACK (http://www.aspect-sdm.org/R-ScaLAPACK) was funded
# as part of the Scientific Data Management Center
# (http://sdm.lbl.gov/sdmcenter) under the Department of Energy's 
# Scientific Discovery through Advanced Computing (DOE SciDAC) program
# (http://www.scidac.org ). 
# ========================================================================
RpdgesvdDemo <- function () {

	library(RScaLAPACK)

        print ("****************************************************")
        print ("           DEMO: SVD DECOMPOSITON                   ")
        print ("****************************************************")

	x <- c(6 ,3 ,0 ,0 ,3 , 0, 0, -3, -1, 1, 1, 0, -1, 0, 11, 0, 0, 10, 0, 0, 0, -11, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 8, 0, -10)

	dim(x) <- c(4,9)

	print (x)

	print (" SVD Decomposition using ScaLAPACK function")

	print (" Process grid specifications .. in order ")

	print (" NPROWS = 2  NPCOLS = 2 ")

	print (" Block size:  MB = 2   NB = 2 ")

	print (" Calling ScaLAPACK function .. sla.svd(x,NULL,NULL,2,2,2) ")

	a <- sla.svd(x,NULL,NULL,2,2,2)

	print(a)

	print (" SVD Decomposition using corrosponding LAPACK function .. La.svd(x) ")

	b <- La.svd(x)

	print (b)

	print (" Compare the results !! ")

	print (" Experiment with different process grid parameters and different block sizes. (A block size of about 64 by 64 is recommended.)")

        print ("****************************************************")
}

