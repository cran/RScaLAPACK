\name{sla.chol}
\title{Choleski Factorization of a real, symmetric, positive-definite, square-matrix using ScaLAPACK}
\usage{
sla.chol(A, NPROWS, NPCOLS, MB, RFLAG, SPAWN)
}
\alias{sla.chol}
\arguments{
	\item{A}{A real symmetric, positive definite matrix}
	
	\item{NPROWS}{Number of Process Rows in the Process Grid. }

	\item{NPCOLS}{Number of Process Cols in the Process Grid. }

	\item{MB}{Block Size. }

    \item{RFLAG}{ Flag saying whether the Process Grid should be released after
computation.}
                                                                                
    \item{SPAWN}{ Flag saying whether a new Process Grid should be spawned.}

}

\description{
     Compute the Choleski factorization of a real symmetric
     positive-definite square matrix using ScaLAPACK.
}

\details{
	If the number of processor rows and columns are both zero, one processor
	is used.  If the number of processor rows is nonzero and the number of
	processor columns is zero, then the processor rows is taken to be a 
	number of processors, and a grid is made accordingly.

}

\references{
	\url{http://www.netlib.org/scalapack/}
}

\seealso{
        \code{\link{sla.solve}} the SCALAPACK version of the R function \code{\link{solve}} 

        \code{\link{sla.qr}} the SCALAPACK version of the R function \code{\link{qr}}

        \code{\link{sla.svd}} the SCALAPACK version of the R function \code{\link{La.svd}} 
  
	\code{\link{La.chol}} the LAPACK function computing Choleski Factiorization in R.
}

\examples{
library(RScaLAPACK)
c(4.16,0,0,0,-3.12,5.03,0,0,0.56,-0.83,0.76,0,-0.10,1.18,0.34,1.18)->x
dim(x)=c(4,4)
sla.chol(x)
}

\concept{parallel}
\keyword{algebra}
