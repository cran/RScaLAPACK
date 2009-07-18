\name{sla.chol2inv}
\title{Invert a Matrix from the Cholesky Decomposition (Using ScaLAPACK)}
\usage{
sla.chol2inv(A, NPROWS=0, NPCOLS=0, MB=16, RFLAG = 1, SPAWN = 1 )
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
rnorm(10000)->x
dim(x)=c(100,100)
# Generate a symmetric, positive definite matrix
sla.chol2inv(x, NPROWS=2, NPCOLS=2)->y

}
\author{Nagiza Samatova (\email{samatovan@ornl.gov}), Guruprasad Kora (\email{koragh@ornl.gov}), Srikanth Yoginath (\email{yoginathsb@ornl.gov}), David Bauer (\email{bauerda@ornl.gov})}
\concept{parallel}
\keyword{algebra}
