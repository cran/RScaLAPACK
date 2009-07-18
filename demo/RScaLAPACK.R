# ====================================================================
#           R-ScaLAPACK version 0.4.x:  ScaLAPACK interface to R
#              Oak Ridge National Laboratory, Oak Ridge TN.
#        Authors: David Bauer, Guruprasad Kora, Nagiza. F. Samatova, 
#                            Srikanth Yoginath.
#     Contact: Nagiza F. Samatova; (865) 241-4351; samatovan@ornl.gov
#     Contact: Guruprasad Kora; (865) 576-6210; koragh@ornl.gov
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

slaSolveDemo <- function () {

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

slaSvdDemo <- function () {

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


slaQrDemo <- function () {

        print ("****************************************************")
        print ("           DEMO: QR DECOMPOSITON                    ")
        print ("****************************************************")   

	x <- c(6 ,3 ,0 ,0 ,3 , 0, 0, -3, -1, 1, 1, 0, -1, 0, 11, 0, 0, 10, 0, 0, 0, -11, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 8, 0, -10)

	dim(x) <- c(4,9)

	print (x)

	print (" QR Decomposition using ScaLAPACK function")

	print (" Process grid specifications .. in order ")

	print (" NPROWS = 2  NPCOLS = 2 ")

	print (" Block size:  MB = 2   NB = 2 ")

	print (" Calling ScaLAPACK function .. sla.qr(x,2,2,2) ")

	a <- sla.qr(x,2,2,2)

	print(a)

	print (" QR Decomposition using corrosponding LAPACK function .. qr(x) ")

	b <- qr(x)

	print (b)

	print (" Compare the results !! ")

	print (" Experiment with different process grid parameters and different block sizes. (A block size of about 64 by 64 is recommended.)")

        print ("****************************************************")   
}


slaEigenDemo <- function () {

        print ("****************************************************")
        print ("           DEMO: EIGEN DECOMPOSITION                ")
        print ("****************************************************")


        x <- c(4.16, -3.12,  0.56, -0.10, -3.12,  5.03, -0.83,  1.18,  0.56, -0.83,  0.76,  0.34, -0.10,  1.18,  0.34,  1.18)

	dim(x) <- c(4,4)

	print (x)

	print (" EIGEN DECOMPOSITION using ScaLAPACK function")

	print (" Process grid specifications .. in order ")

	print (" NPROWS = 2  NPCOLS = 2 ")

	print (" Block size:  MB = 2   NB = 2 ")

	print (" Calling ScaLAPACK function .. sla.eigen(x) ")

	a <- sla.eigen(x)

	print(a)

	print (" Eigen Decomoposition using corrosponding LAPACK function .. eigen(x) ")

	b <- eigen(x)

	print (b)

	print (" Compare the results !! ")

	print (" Experiment with different process grid parameters and different block sizes.")

       print ("****************************************************")

}


slaCholDemo <- function () {

        print ("****************************************************")
        print ("           DEMO: CHOLESKY FACTORIZATION             ")
        print ("****************************************************")

	x <- c(4.16, -3.12,  0.56, -0.10, -3.12,  5.03, -0.83,  1.18,  0.56, -0.83,  0.76,  0.34, -0.10,  1.18,  0.34,  1.18)

	dim(x) <- c(4,4)

	print (x)

	print (" SOLVING using SCALAPACK FUNCTION ")

	print (" PROCESS GRID SPECIFICATIONS .. in order ")

	print (" NPROWS = 2  NPCOLS = 2 ")

	print (" BLOCK SIZE LHS:  MB = 2   NB = 2 ")

	print (" CALLING SCALAPACK FUNCTION .. sla.chol(x,2,2,2) ")

	a <- sla.chol(x,2,2,2)

	print(a)

	print (" SOLVING USING CORRESPONDING LAPACK FUNCTION .. La.chol(x) ")

	b <- La.chol(x)

	print (b)

	print (" COMPARE BOTH RESULTS !! ")

	print (" EXPERIMENT WITH DIFFERENT PROCESS GRID PARAMETERS and INPUT DATA")

        print ("****************************************************")
    
}


slaChol2invDemo <- function () {

        print ("****************************************************")
        print (" DEMO:INVERT A MATRIX FROM CHOLESKY DECOMPOSITION   ")
        print ("****************************************************")


	x <- c(4.16, -3.12,  0.56, -0.10, -3.12,  5.03, -0.83,  1.18,  0.56, -0.83,  0.76,  0.34, -0.10,  1.18,  0.34,  1.18)

	dim(x) <- c(4,4)

	print (x)

	print (" SOLVING using SCALAPACK FUNCTION ")

	print (" PROCESS GRID SPECIFICATIONS .. in order ")

	print (" NPROWS = 2  NPCOLS = 2 ")

	print (" BLOCK SIZE LHS:  MB = 2   NB = 2 ")

	print (" CALLING SCALAPACK FUNCTION .. sla.chol(x,2,2,2) ")

	a <- sla.chol(x,2,2,2)

	print(a)
          
        print (" CALLING SCALAPACK Function .... sla.chol2inv(a,2,2,2) ")

        aC2inv <- sla.chol2inv(a)

        print(aC2inv)

        print (" VERIFYING THE RESULT:  x %*% aC2inv ")

        ver <- x %*% aC2inv

        print(ver)

	print (" USING CORRESPONDING LAPACK FUNCTION .. La.chol(x) ")

	b <- La.chol(x)

	print (b)

        print ("CALLING ... chol2inv(b)")

        bC2inv <- chol2inv(b)

        print (" VERIFYING RESULT : ")

        ver <- x %*% bC2inv

        print(ver)
 
	print (" COMPARE BOTH RESULTS !! ")

	print (" EXPERIMENT WITH DIFFERENT PROCESS GRID PARAMETERS and INPUT DATA")
}


slaMultiplyDemo <- function () {

        print ("****************************************************")
        print ("            DEMO: Matrix Multiplication             ")
        print ("****************************************************")

	x <- matrix(rnorm(36),6,6)

	y <- matrix(rnorm(36),6,6)

	print (x)

	print (y)

	print (" Matrix Multiplication using SCALAPACK FUNCTION ")

	print (" PROCESS GRID SPECIFICATIONS .. in order ")

	print (" NPROWS = 2  NPCOLS = 2 ")

	print (" BLOCK SIZE LHS:  MB = 2   NB = 2 ")

	print (" CALLING SCALAPACK FUNCTION .. sla.multiply(x,y,2,2,2) ")

	a <- sla.multiply(x,y,2,2,2)

	print(a)

	print (" SOLVING USING CORRESPONDING LAPACK FUNCTION .. x %*% y ")

	b <- x %*% y

	print (b)

	print (" COMPARE BOTH RESULTS !! ")

	print (" EXPERIMENT WITH DIFFERENT PROCESS GRID PARAMETERS ")

       print ("****************************************************")

}

