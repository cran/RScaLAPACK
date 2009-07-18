\name{sla.svd}
\title{Singular Value Decomposition of a Matrix}
\usage{
sla.svd(A, nu=NULL, nv=NULL, NPROWS=0, NPCOLS=0, MB=16, RFLAG=1, SPAWN=1)
}
\alias{sla.svd}
\arguments{
	\item{A}{A matrix whose SVD decomposition is to be computed }
	
	\item{nu}{The number of left singular vectors to be computed.  Must be either \code{0} or \code{min(nrow(x), ncol(x))}.}

	\item{nv}{The number of right singular vectors to be computed.  Must be either \code{0} or \code{min(nrow(x), ncol(x))}.}
	
	\item{NPROWS}{Number of Process Rows in the Process Grid. }

	\item{NPCOLS}{Number of Process Cols in the Process Grid. }

	\item{MB}{Block Size. }

    \item{RFLAG}{ Flag saying whether the Process Grid should be released after
computation.}
                                                                                
    \item{SPAWN}{ Flag saying whether a new Process Grid should be spawned.}

}

\description{
     Compute the singular-value decomposition of a rectangular matrix.
}

\details{
	If the number of processor rows and columns are both zero, one processor
	is used.  If the number of processor rows is nonzero and the number of
	processor columns is zero, then the processor rows is taken to be a 
	number of processors, and a grid is made accordingly.

}

\references{
	\url{http://www.netlib.org/scalapack/}

	\url{http://mathworld.wolfram.com/SingularValueDecomposition.html}
}

\seealso{
	\code{\link{sla.solve}} the SCALAPACK version of the R function \code{\link{solve}}

	\code{\link{sla.qr}} the SCALAPACK version of the R function \code{\link{qr}}

	\code{\link{La.svd}} the LAPACK implementation of singular value decomposition.
}

\examples{
library(RScaLAPACK)
rnorm(64*64)->x
dim(x)=c(64,64)
sla.svd(x, NPROWS=2, NPCOLS=2, MB=16)
}
\author{Nagiza Samatova (\email{samatovan@ornl.gov}), Guruprasad Kora (\email{koragh@ornl.gov}), Srikanth Yoginath (\email{yoginathsb@ornl.gov}), David Bauer (\email{bauerda@ornl.gov})}
\concept{parallel}
\keyword{algebra}
