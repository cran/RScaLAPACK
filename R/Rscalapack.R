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

sla.pdgesv<- function () {
	abc<-.Call("CR_Exec", PACKAGE="RScaLAPACK")
	if (abc == 0) {
		print("YOUR EXAMPLE WORKS FINE ;-)");
	}
	else {
		print("ERROR IN sla.pdgesv - .Call returned nonzero")
	}	
}

sla.solve <- function (A, B=NULL, NPROWS=0, NPCOLS=0, MB=16,RFLAG=1, SPAWN=1 ){

# Check for illigal values

	if ( NPROWS < 0 )
	{
		stop("NPROWS cannot be a negative value.");
	}

	if ( NPCOLS < 0 )
	{
		stop("NPCOLS cannot be a negative value.");
	}

	if ( MB < 0 )
	{
		stop("MB cannot be a negative value.");
	}

	if ( RFLAG < 0 )
	{
		stop("RFLAG cannot be a negative value.");
	}

	if ( SPAWN < 0 )
	{
		stop("SPAWN cannot be a negative value.");
	}


# Check for unsupported datatype

	if (is.complex(A) || is.complex(B) ) {
		stop("Complex datatype is not supported in this version.");
	}

# Save the dimensions of the matrices
	dimA = dim(A)
	dimB = dim(B)

	if ( dimA[1] < NPROWS || dimA[1] < NPCOLS || dimA[2] < NPROWS || dimA[2] < NPCOLS )
	{
		stop("Input matrix is smaller than process grid");
	}

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

# Check block size for proper distribution
	MB = sla.checkBlockSize(dimA, MB)

	inputVector <- list(A,B,as.integer(NPROWS),as.integer(NPCOLS),as.integer(MB), as.integer(1),as.integer(RFLAG),as.integer(SPAWN))

	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector)
	
        attributes(x) = list(names="x")
        return (x$x)
}

sla.svd <- function (A, nu=NULL, nv=NULL, NPROWS=0, NPCOLS=0, MB=16, RFLAG=1, SPAWN=1){

# Check for illigal values

	if ( NPROWS < 0 )
	{
		stop("NPROWS cannot be a negative value.");
	}

	if ( NPCOLS < 0 )
	{
		stop("NPCOLS cannot be a negative value.");
	}

	if ( MB < 0 )
	{
		stop("MB cannot be a negative value.");
	}

	if ( RFLAG < 0 )
	{
		stop("RFLAG cannot be a negative value.");
	}

	if ( SPAWN < 0 )
	{
		stop("SPAWN cannot be a negative value.");
	}



# Check for unsupported datatype

	if (is.complex(A) ) {
		print("Complex datatype is not supported in this version.");
	}

# Save the dimensions of the matrices
	dimA = dim(A)

	if ( dimA[1] < NPROWS || dimA[1] < NPCOLS || dimA[2] < NPROWS || dimA[2] < NPCOLS )
	{
		stop("Input matrix is smaller than process grid");
	}

# Force the matrices to be REALSXPs
	A <- as.real(A)

# Restore the dimensions
	dim(A)<-dimA

	if (is.null(nu) && is.null(nv)) {
#		B = as.real(NULL)
		temp = as.real(min(dim(A)))
		B = c(temp, temp)
		#print (B)
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

# Check block size for proper distribution
	MB = sla.checkBlockSize(dimA, MB)

	inputVector <- list(A,B,as.integer(NPROWS),as.integer(NPCOLS),as.integer(MB), as.integer(2), as.integer(RFLAG),as.integer(SPAWN))

	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector)

	if (length(x) == 3) {
		attributes(x) = list(names=c("d", "u", "vt"))
		dim(x$d) = NULL
		}

	return (x)
}

sla.qr <- function (A, NPROWS=0, NPCOLS=0, MB=16, RFLAG = 1, SPAWN = 1 ){


# Check for illigal values

	if ( NPROWS < 0 )
	{
		stop("NPROWS cannot be a negative value.");
	}

	if ( NPCOLS < 0 )
	{
		stop("NPCOLS cannot be a negative value.");
	}

	if ( MB < 0 )
	{
		stop("MB cannot be a negative value.");
	}

	if ( RFLAG < 0 )
	{
		stop("RFLAG cannot be a negative value.");
	}

	if ( SPAWN < 0 )
	{
		stop("SPAWN cannot be a negative value.");
	}


 
# Check for unsupported datatype

    if (is.complex(A) ) {
        print("Complex datatype is not supported in this version.");
    }


# Save the dimensions of the matrix
	dimA = dim(A)
	if ( dimA[1] < NPROWS || dimA[1] < NPCOLS || dimA[2] < NPROWS || dimA[2] < NPCOLS )
	{
		stop("Input matrix is smaller than process grid");
	}

# Force the matrix to be REALSXPs
	A <- as.real(A)

# Restore the dimensions
	dim(A)<-dimA

	t = sla.ProcessRowsColumns(NPROWS, NPCOLS)
	if (is.null(t)) return
	NPROWS=t[1]
	NPCOLS=t[2]

# Check block size for proper distribution
	MB = sla.checkBlockSize(dimA, MB)
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

sla.chol <- function (A, NPROWS=0, NPCOLS=0, MB=16, RFLAG = 1, SPAWN = 1 ){

# Check for illigal values

	if ( NPROWS < 0 )
	{
		stop("NPROWS cannot be a negative value.");
	}

	if ( NPCOLS < 0 )
	{
		stop("NPCOLS cannot be a negative value.");
	}

	if ( MB < 0 )
	{
		stop("MB cannot be a negative value.");
	}

	if ( RFLAG < 0 )
	{
		stop("RFLAG cannot be a negative value.");
	}

	if ( SPAWN < 0 )
	{
		stop("SPAWN cannot be a negative value.");
	}


# Check for unsupported datatype

    if (is.complex(A)) {
        print("Complex datatype is not supported in this version.");
    }




# Save the dimensions of the matrices
	dimA = dim(A)

	if ( dimA[1] < NPROWS || dimA[1] < NPCOLS || dimA[2] < NPROWS || dimA[2] < NPCOLS )
	{
		stop("Input matrix is smaller than process grid");
	}

# Check for square matrix

	if ( dimA[1] != dimA[2] )
	{
		stop("non-square matrix in 'sla.eigen'");
	}


# Force the matrices to be REALSXPs
	A <- as.real(A)

# Restore the dimensions
	dim(A)<-dimA

	t = sla.ProcessRowsColumns(NPROWS, NPCOLS)
	if (is.null(t)) return
	NPROWS=t[1]
	NPCOLS=t[2]

# Check block size for proper distribution
	MB = sla.checkBlockSize(dimA, MB)

	inputVector <- list(A,as.real(NULL),as.integer(NPROWS),as.integer(NPCOLS),as.integer(MB), as.integer(4), as.integer(RFLAG),as.integer(SPAWN))

	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector)
     
    if (length(x) == 1) {
        attributes(x) = list(names="x")
        }

     x$x[lower.tri(x$x)]<-0

     return (x$x)
}

sla.chol2inv <- function (A, NPROWS=0, NPCOLS=0, MB=16, RFLAG = 1, SPAWN = 1 ){


# Check for illigal values

	if ( NPROWS < 0 )
	{
		stop("NPROWS cannot be a negative value.");
	}

	if ( NPCOLS < 0 )
	{
		stop("NPCOLS cannot be a negative value.");
	}

	if ( MB < 0 )
	{
		stop("MB cannot be a negative value.");
	}

	if ( RFLAG < 0 )
	{
		stop("RFLAG cannot be a negative value.");
	}

	if ( SPAWN < 0 )
	{
		stop("SPAWN cannot be a negative value.");
	}



# Check for unsupported datatype

    if (is.complex(A) ) {
        print("Complex datatype is not supported in this version.");
    }


# Save the dimensions of the matrix
	dimA = dim(A)

	if ( dimA[1] < NPROWS || dimA[1] < NPCOLS || dimA[2] < NPROWS || dimA[2] < NPCOLS )
	{
		stop("Input matrix is smaller than process grid");
	}

# Force the matrix to be REALSXPs
	A <- as.real(A)

# Restore the dimensions
	dim(A)<-dimA

	t = sla.ProcessRowsColumns(NPROWS, NPCOLS)
	if (is.null(t)) return
	NPROWS=t[1]
	NPCOLS=t[2]

# Check block size for proper distribution
	MB = sla.checkBlockSize(dimA, MB)

	inputVector <- list(A,as.real(NULL),as.integer(NPROWS),as.integer(NPCOLS),as.integer(MB), as.integer(5), as.integer(RFLAG),as.integer(SPAWN))

	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector)

    if (length(x) == 1) {
        attributes(x) = list(names="x")
        }
	x$x[lower.tri(x$x)]=t(x$x)[lower.tri(x$x)]
	return (x$x)
}

sla.eigen <- function (A, NPROWS=0, NPCOLS=0, MB=16, RFLAG = 1, SPAWN = 1 ){

# Check for illigal values

	if ( NPROWS < 0 )
	{
		stop("NPROWS cannot be a negative value.");
	}

	if ( NPCOLS < 0 )
	{
		stop("NPCOLS cannot be a negative value.");
	}

	if ( MB < 0 )
	{
		stop("MB cannot be a negative value.");
	}

	if ( RFLAG < 0 )
	{
		stop("RFLAG cannot be a negative value.");
	}

	if ( SPAWN < 0 )
	{
		stop("SPAWN cannot be a negative value.");
	}


# Check for unsupported datatype

    if (is.complex(A) ) {
        print("Complex datatype is not supported in this version.");
    }


# Save the dimensions of the matrix
	dimA = dim(A)

	if ( dimA[1] < NPROWS || dimA[1] < NPCOLS || dimA[2] < NPROWS || dimA[2] < NPCOLS )
	{
		stop("Input matrix is smaller than process grid");
	}

# Force the matrix to be REALSXPs
	A <- as.real(A)

# Restore the dimensions
	dim(A)<-dimA

	t = sla.ProcessRowsColumns(NPROWS, NPCOLS)
	if (is.null(t)) return
	NPROWS=t[1]
	NPCOLS=t[2]

# Check for square matrix

	if ( dimA[1] != dimA[2] )
	{
		stop("non-square matrix in 'sla.eigen'");
	}

# Check block size for proper distribution
	MB = sla.checkBlockSize(dimA, MB)

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
        NPROWS = .RscalaGrid$NPROWS
        NPCOLS = .RscalaGrid$NPCOLS
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
#        t = NPROWS
#        NPROWS = NPCOLS
#        NPCOLS = t
    }
                                                                                
    return( c(NPROWS, NPCOLS) )
}

sla.gridInfo <- function () {
	x<-.External("PA_GridInfo", PACKAGE="RScaLAPACK");
	
	invisible(x);
}

sla.gridInit <- function (NPROWS=1, NPCOLS=1, BLOCKSIZE=16) {

# Check for illigal values

	if ( NPROWS < 0 )
	{
		stop("NPROWS cannot be a negative value.");
	}

	if ( NPCOLS < 0 )
	{
		stop("NPCOLS cannot be a negative value.");
	}

	if ( BLOCKSIZE < 0 )
	{
		stop("BLOCKSIZE cannot be a negative value.");
	}


	if (exists(".RscalaGrid", envir=.GlobalEnv)) {
		stop("Attempt to spawn a new grid without releasing the previous grid.");	
	}
	inputVector <- list(as.real(NULL), as.real(NULL), as.integer(NPROWS),
		as.integer(NPCOLS), as.integer(BLOCKSIZE), as.integer(0), as.integer(0),
		as.integer(1))

	.RscalaGrid = list(NPROWS = NPROWS, NPCOLS = NPCOLS, BLOCKSIZE = BLOCKSIZE)
	assign(".RscalaGrid", .RscalaGrid, env = .GlobalEnv);

	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector)
		print ("RScaLAPACK:Process Grid Initialized " ) 
	}

sla.gridExit <- function() {
	inputVector <- list(as.real(NULL),as.real(NULL),as.integer(0),as.integer(0),as.integer(0), as.integer(0), as.integer(1),as.integer(0));
	x <- PA.exec("RScaLAPACK", "pdgesv.R",inputVector);
	rm( ".RscalaGrid", envir=.GlobalEnv);
    print ("RScaLAPACK:Process Grid Released " )  
}

sla.checkBlockSize <-function(tdimA, tMB) {
	temp = (tdimA[1] * tdimA[2]) / tMB
	if (temp > 32768 ){
		temp1 = (tdimA[1] * tdimA[2]) / 32768
		return (temp1)
	} else {
		smallerDim = ifelse( tdimA[1] < tdimA[2], tdimA[1], tdimA[2]);
		temp1 = ifelse( (tMB > smallerDim), (smallerDim-1), tMB );

		if ( temp1 == tMB )
		{
			return(tMB);
		}
		else
		{
			print("Warning:Block size is higher than the lowest dimension, forcing it to lowest dim");
			return(temp1);
		}
	}
}
	
sla.multiply <- function( A, B=NULL, NPROWS=0, NPCOLS=0, MB=16, RFLAG=1, SPAWN=1 ) {


# Check for illigal values

	if ( NPROWS < 0 )
	{
		stop("NPROWS cannot be a negative value.");
	}

	if ( NPCOLS < 0 )
	{
		stop("NPCOLS cannot be a negative value.");
	}

	if ( MB < 0 )
	{
		stop("MB cannot be a negative value.");
	}

	if ( RFLAG < 0 )
	{
		stop("RFLAG cannot be a negative value.");
	}

	if ( SPAWN < 0 )
	{
		stop("SPAWN cannot be a negative value.");
	}



# Check for unsupported datatype

    if (is.complex(A) || is.complex(B) ) {
        print("Complex datatype is not supported in this version.");
    }


#Check for availability of second matrix
        if ( is.matrix(B) == FALSE )
        {
                print("Warning:Second matrix is missing, Assuming second matrix as A")
                B = A
        }

#Get matrix sizes
        dimA = dim(A)
        dimB = dim(B)

		if ( dimA[1] < NPROWS || dimA[1] < NPCOLS || dimA[2] < NPROWS || dimA[2] < NPCOLS )
		{
			stop("Input matrix is smaller than process grid");
		}
        
#Get metrix dimensions
        dimAlength = length(dimA)
        dimBlength = length(dimB)
        
#Check for 2D matrix
        if( dimAlength != 2 || dimBlength != 2 )
        {       
                stop("One (or both) of the matrices is not a 2D matrix.\n3D matrices are not yes supported.")
        }

#Check for crossproduct rule conformity, Cik = Aij*Bjk
        if( dimA[2] != dimB[1] ){ 
                stop( "Given matrix sizes do not conform to matrix multiplication rules" )
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

# Check block size for proper distribution
	MB = sla.checkBlockSize(dimA, MB)

	inputVector <- list(A,B,as.integer(NPROWS), as.integer(NPCOLS), as.integer(MB), as.integer(7),as.integer(RFLAG),as.integer(SPAWN))

	#call 
	x <- PA.exec("RScaLAPACK", "pdgesv.R", inputVector)

	#check for the returned results
	
	#return the answer
        return (x)

}

sla.quit <- function (){
	PA.exit()
	q()
}
