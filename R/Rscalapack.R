#  =======================================================================
#           R-ScaLAPACK version 0.2:  ScaLAPACK interface to R
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

sla.pdgesv<- function () {
	abc<-.Call("CR_Exec", PACKAGE="RScaLAPACK")
	if (abc == 0) {
		print("YOUR EXAMPLE WORKS FINE ;-)");
	}
	else {
		print("ERROR IN sla.pdgesv - .Call returned nonzero")
	}	
}


sla.solve <- function (A=0, B=NULL, NPROWS=0, NPCOLS=0, MB=64,RFLAG=1, SPAWN=1 ){

# Save the dimensions of the matrices
	dimA = dim(A)
	dimB = dim(B)

# If Matrix B wasn't given, assume the identity matrix as the RHS
	if (is.null(B)) {
		B = diag(1,dimA[1], dimA[1])
		dimB = c(dimA[1], dimA[1])
	}

# If the B is not a matrix, make it one
	if (length(dimB) == 0 || length(dimB) == 1) {
		dimB = c(length(B), 1)
	}

# Check the matrices
	if (dimA[1] != dimA[2]) {
		stop("Coefficient matrix must be square.")
	}
	if (B != 0 && dimA[1] != dimB[1]) {
		stop("Right-hand side matrix present, but dimensions do not match the coefficient matrix.")
	}

# Force the matrices to be REALSXPs
	A <- as.real(A)
	B <- as.real(B)

# Restore the dimensions
	dim(A)<-dimA
	dim(B)<-dimB

	t = sla.ProcessRowsColumns(NPROWS, NPCOLS)
	if (is.null(t)) return
	NPROWS=t[1]
	NPCOLS=t[2]

	inputVector <- list(A,B,as.integer(NPROWS),as.integer(NPCOLS),as.integer(MB), as.integer(1),as.integer(RFLAG),as.integer(SPAWN))

	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector)
    if (length(x) == 1) {
        attributes(x) = list(names="x")
        return (x$x)
    } else  return (x$x)
}

sla.svd <- function (A=0, nu=NULL, nv=NULL, NPROWS=0, NPCOLS=0, MB=64, RFLAG=1, SPAWN=1){

# Save the dimensions of the matrices
	dimA = dim(A)

# Force the matrices to be REALSXPs
	A <- as.real(A)

# Restore the dimensions
	dim(A)<-dimA

	if (is.null(nu) && is.null(nv)) {
		B = as.real(NULL)
	} else {
		if (is.null(nv)) {
			B = c(nu, min(dim(A)))
		} else if (is.null(nu)) {
			B = c(min(dim(A)), nv)
		} else {
			B = c(nu, nv)
		}
		if ((B[1] != 0 && B[1] != min(dim(A))) || (B[2] != 0 && B[2] !=
													min(dim(A)))) {
			stop("sla.svd only accepts 0 or min(n,p) for nu and nv")
		}
	}

	t = sla.ProcessRowsColumns(NPROWS, NPCOLS)
	if (is.null(t)) return
	NPROWS=t[1]
	NPCOLS=t[2]

	inputVector <- list(A,B,as.integer(NPROWS),as.integer(NPCOLS),as.integer(MB), as.integer(2), as.integer(RFLAG),as.integer(SPAWN))

	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector)

    if (length(x) == 1) {
        attributes(x) = list(names="x")
        }

	if (length(x) == 3) {
		attributes(x) = list(names=c("d", "u", "vt"))
		dim(x$d) = NULL
		}

#	if (length(x$u) == 0) {
#		x$u = NULL
#		}
	
#	if (length(x$vt) == 0) {
#		x$vt = NULL
#		}

	return (x)
}

sla.qr <- function (A=0, NPROWS=0, NPCOLS=0, MB=48, RFLAG = 1, SPAWN = 1 ){

# Save the dimensions of the matrix
	dimA = dim(A)

# Force the matrix to be REALSXPs
	A <- as.real(A)

# Restore the dimensions
	dim(A)<-dimA

	t = sla.ProcessRowsColumns(NPROWS, NPCOLS)
	if (is.null(t)) return
	NPROWS=t[1]
	NPCOLS=t[2]

	inputVector <- list(A,as.real(NULL),as.integer(NPROWS),as.integer(NPCOLS),as.integer(MB), as.integer(3), as.integer(RFLAG),as.integer(SPAWN))

	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector)

    if (length(x) == 1) {
        attributes(x) = list(names="qr")
        }

	if (length(x) == 2) {
		attributes(x) = list(names=c("qr", "tau"))
		}
	if (length(x) == 3) {
		attr(x, "names") = c("qr", "rank", "qraux")
		attr(x, "useLAPACK") = TRUE
		attr(x, "class") = "qr"
		dim(x$qraux) = NULL
		}
	
	return (x)
}


sla.chol <- function (A=0, NPROWS=0, NPCOLS=0, MB=64, RFLAG = 1, SPAWN = 1 ){

# Save the dimensions of the matrices
	dimA = dim(A)

# Force the matrices to be REALSXPs
	A <- as.real(A)

# Restore the dimensions
	dim(A)<-dimA

	t = sla.ProcessRowsColumns(NPROWS, NPCOLS)
	if (is.null(t)) return
	NPROWS=t[1]
	NPCOLS=t[2]

	inputVector <- list(A,as.real(NULL),as.integer(NPROWS),as.integer(NPCOLS),as.integer(MB), as.integer(4), as.integer(RFLAG),as.integer(SPAWN))

	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector)
     
    if (length(x) == 1) {
        attributes(x) = list(names="x")
        }

     x$x[lower.tri(x$x)]<-0

     return (x$x)
}

sla.chol2inv <- function (A=0, NPROWS=0, NPCOLS=0, MB=48, RFLAG = 1, SPAWN = 1 ){

# Save the dimensions of the matrix
	dimA = dim(A)

# Force the matrix to be REALSXPs
	A <- as.real(A)

# Restore the dimensions
	dim(A)<-dimA

	t = sla.ProcessRowsColumns(NPROWS, NPCOLS)
	if (is.null(t)) return
	NPROWS=t[1]
	NPCOLS=t[2]

	inputVector <- list(A,as.real(NULL),as.integer(NPROWS),as.integer(NPCOLS),as.integer(MB), as.integer(5), as.integer(RFLAG),as.integer(SPAWN))

	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector)

    if (length(x) == 1) {
        attributes(x) = list(names="x")
        }
	x$x[lower.tri(x$x)]=t(x$x)[lower.tri(x$x)]
	return (x$x)
}

sla.eigen <- function (A=0, NPROWS=0, NPCOLS=0, MB=48, RFLAG = 1, SPAWN = 1 ){

# Save the dimensions of the matrix
	dimA = dim(A)

# Force the matrix to be REALSXPs
	A <- as.real(A)

# Restore the dimensions
	dim(A)<-dimA

	t = sla.ProcessRowsColumns(NPROWS, NPCOLS)
	if (is.null(t)) return
	NPROWS=t[1]
	NPCOLS=t[2]

	inputVector <- list(A,as.real(NULL),as.integer(NPROWS),as.integer(NPCOLS),as.integer(MB), as.integer(6),as.integer(RFLAG),as.integer(SPAWN))

	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector)

    if (length(x) == 1) {
        attributes(x) = list(names="values")
        }

	if (length(x) == 2) {
		attributes(x) = list(names=c("values", "vectors"))
		dim(x$values) = NULL
		}
	
	return (x)
}

sla.ProcessRowsColumns <- function(NPROWS, NPCOLS) {
    NPROWS = as.integer(NPROWS)
    NPCOLS = as.integer(NPCOLS)
                                                                                
    if (NPROWS < 0 || NPCOLS < 0) {
        stop("Bad grid given.  Rows and columns must be >= zero")
    }
                                                                                
# If the user didn't specify a grid, but a grid was found in the root
# environment space, use that (but give a warning.)
    if (NPROWS == 0 && exists(".RscalaGrid") && is.numeric(.RscalaGrid)) {
        NPROWS = .RscalaGrid[1];
        NPCOLS = .RscalaGrid[2];
        if (is.na(NPCOLS)) NPCOLS = 0
        if (is.na(NPROWS) || NPROWS < 1 || NPCOLS < 0) {
            print("Bad grid found in .RscalaGrid, ignoring")
            NPROWS = NPCOLS = 0
        } else {
            print("Using default grid from .RscalaGrid")
        }
    }
                                                                                
# If the user didn't specify a number of rows, use a single process
# If the user didn't specify a number of columns, assume the number
# of rows given should actually be the number of columns.
    if (NPROWS == 0) {
        NPROWS = 1
        NPCOLS = 1
        } else if (NPCOLS == 0) {
        NPCOLS = floor(sqrt(NPROWS))
        while (NPROWS / NPCOLS != floor(NPROWS/NPCOLS)) {
            NPCOLS = NPCOLS + 1
        }
        NPROWS = NPROWS / NPCOLS
                                                                                
                                                                                
        # Scalapack is faster with these switched.
        t = NPROWS
        NPROWS = NPCOLS
        NPCOLS = t
    }
                                                                                
    return( c(NPROWS, NPCOLS) )
}

sla.gridInit <- function (NPROCS=0) {
	inputVector <- list(as.real(NULL), as.real(NULL), as.integer(NPROCS),
		as.integer(1), as.integer(0), as.integer(0), as.integer(0),
		as.integer(1))
	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector)
		print (" Process Grid Initialized " ) 
	}

sla.gridExit <- function() {
	inputVector <- list(as.real(NULL),as.real(NULL),as.integer(0),as.integer(0),as.integer(0), as.integer(0), as.integer(1),as.integer(0))
	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector)
        print (" Process Grid Released ... " )  
}
