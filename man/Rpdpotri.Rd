\name{sla.chol2inv}
\title{Invert a Matrix from the Cholesky Decomposition (Using ScaLAPACK)}
\usage{
sla.chol2inv(A, NPROWS, NPCOLS, MB, RFLAG, SPAWN)
}
\alias{sla.chol2inv}
\arguments{
	\item{A}{a matrix.  The upper triangle should hold the Choleski Decomposition to compute the inverse of.}
	
	\item{NPROWS}{Number of Process Rows in the Process Grid. }

	\item{NPCOLS}{Number of Process Cols in the Process Grid. }

	\item{MB}{Block Size. }

    \item{RFLAG}{ Flag saying whether the Process Grid should be released after
computation.}
                                                                                
    \item{SPAWN}{ Flag saying whether a new Process Grid should be spawned.}

}

\description{
	Invert a symmetric, positive definite, square matrix from its Choleski decomposition using ScaLAPACK.
}

\details{
	If the number of processor rows and columns are both zero, one processor
	is used.  If the number of processor rows is nonzero and the number of
	processor columns is zero, then the processor rows is taken to be a 
	number of processors, and a grid is made accordingly.
}

\value{
	The inverse of the decomposed matrix is returned.
	}

\references{
	\url{http://mathworld.wolfram.com/CholeskyDecomposition.html}

	\url{http://www.netlib.org/scalapack/}
}

\seealso{
	\code{chol2inv} for the normal R version.
}

\examples{
library(RScaLAPACK)
rnorm(100)->x
dim(x)=c(10,10)
# First, generate a symmetric, positive definite matrix
sla.chol2inv(x)->y
print(y)
# Invert
sla.solve(y)->y2
# Perform a Choleski decomposition
sla.chol(y2)->x

}

\concept{parallel}
\keyword{algebra}
