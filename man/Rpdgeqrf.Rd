\name{sla.qr}
\title{The QR Decomposition of a Matrix (Using ScaLAPACK)}
\usage{
sla.qr(A, NPROWS=0, NPCOLS=0, MB=16, RFLAG = 1, SPAWN = 1 )
}
\alias{sla.qr}
\arguments{
	\item{A}{A numeric matrix whose QR decomposition is to be computed.}
	
	\item{NPROWS}{Number of Process Rows in the Process Grid. }

	\item{NPCOLS}{Number of Process Cols in the Process Grid. }

	\item{MB}{Block Size. }

	\item{RFLAG}{ Flag saying whether the Process Grid should be released after computation.}

	\item{SPAWN}{ Flag saying whether a new Process Grid should be spawned.}

}

\description{
	The function \code{sla.qr} computes the QR decomposition of a matrix using
	the ScaLAPACK library routine PDGEQRF.
}

\details{
	If the number of processor rows and columns are both zero, one processor
	is used.  If the number of processor rows is nonzero and the number of
	processor columns is zero, then the processor rows is taken to be a 
	number of processors, and a grid is made accordingly.
}

\value{
	The QR decomposition of the matrix as computed by PDGEQRF is returned.
	\item{qr}{a matrix with the same dimensions as \code{A}.
	The upper triangle contains the \eqn{\bold{R}} of the decomposition and
	the lower triangle contains a representation of the \eqn{\bold{Q}} of
	the decomposition (stored in the LAPACK compact form).}
	\item{qraux}{a vector of length \code{min(ncol(A),nrow(A)} which contains
	information needed to extract \eqn{\bold{Q}}.}
	\item{rank}{the length of vector \eqn{\bold{qraux}}.}
	}

\references{
	\url{http://mathworld.wolfram.com/QRDecomposition.html}

	\url{http://www.netlib.org/scalapack/}
}

\seealso{
	\code{\link{qr}} for the normal solving function.
}

\examples{
library(RScaLAPACK)
rnorm(10000)->x
dim(x)=c(100,100)
sla.qr(x, NPROWS=2,NPCOLS=2)
}
\author{Nagiza Samatova (\email{samatovan@ornl.gov}), Guruprasad Kora (\email{koragh@ornl.gov}), Srikanth Yoginath (\email{yoginathsb@ornl.gov}), David Bauer (\email{bauerda@ornl.gov})}
\concept{parallel}
\keyword{algebra}
