\name{sla.solve}
\title{Solve a System of Equations using Scalapack functions}
\usage{

sla.solve(A, B=NULL, NPROWS=0, NPCOLS=0, MB=16,RFLAG=1, SPAWN=1 )

}
\alias{sla.solve}
\arguments{
	\item{A}{A square numeric matrix containing the coefficients of the
	linear system.}
	
	\item{B}{A numeric vector or matrix giving the RHS of the linear system. }

	\item{NPROWS}{Number of Process Rows in the Process Grid. }

	\item{NPCOLS}{Number of Process Cols in the Process Grid. }

	\item{MB}{Block Size. }

    \item{RFLAG}{ Flag saying whether the Process Grid should be released after
computation.}
                                                                                
    \item{SPAWN}{ Flag saying whether a new Process Grid should be spawned.}

}

\description{
	This generic function solves the equation \code{a \%*\% x = b} 
	for \code{x} using the ScaLAPACK library, where \code{b} can be
	either a vector or a matrix.
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
	\code{\link{solve}} for the normal solving function. 

	\code{\link{sla.qr}} the SCALAPACK version of the R function \code{\link{qr}}

}

\examples{
library(RScaLAPACK)
rnorm(10000)->x
dim(x)=c(100,100)
c(1,0*(1:99))->y
dim(y)=c(100,1)
sla.solve(x,y, NPROWS=2, NPCOLS=2)
}
\author{Nagiza Samatova (\email{samatovan@ornl.gov}), Guruprasad Kora (\email{koragh@ornl.gov}), Srikanth Yoginath (\email{yoginathsb@ornl.gov}), David Bauer (\email{bauerda@ornl.gov})}
\concept{parallel}
\keyword{algebra}
