\name{sla.eigen}
\title{Spectral (Eigen) Decomposition of a Matrix}
\usage{
sla.eigen(A, NPROWS, NPCOLS, MB, RFLAG, SPAWN)
}
\alias{sla.eigen}
\arguments{
	\item{A}{A symmetric, square matrix whose spectral (eigen) decomposition is to be computed }
	
	\item{NPROWS}{Number of Process Rows in the Process Grid. }

	\item{NPCOLS}{Number of Process Cols in the Process Grid. }

	\item{MB}{Block Size. }

    \item{RFLAG}{ Flag saying whether the Process Grid should be released after
computation.}
                                                                                
    \item{SPAWN}{ Flag saying whether a new Process Grid should be spawned.}

}

\description{
     Computes the eigenvalues and eigenvectores of a symmetric square matrix.
}

\details{
	If the number of processor rows and columns are both zero, one processor
	is used.  If the number of processor rows is nonzero and the number of
	processor columns is zero, then the processor rows is taken to be a 
	number of processors, and a grid is made accordingly.

}

\references{
	\url{http://www.netlib.org/scalapack/}

	\url{http://mathworld.wolfram.com/EigenDecomposition.html}
}

\seealso{
	\code{\link{eigen}} the normal R version of the function
}

\examples{
library(RScaLAPACK)
rnorm(100)->x
dim(x)=c(10,10)
x[lower.tri(x)]=t(x)[lower.tri(x)]
sla.eigen(x)
}

\concept{parallel}
\keyword{algebra}
