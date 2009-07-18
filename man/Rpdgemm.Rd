\name{sla.multiply}
\title{Matrix Multiplication using Scalapack function}
\usage{

sla.multiply(A, B=NULL, NPROWS=0, NPCOLS=0, MB=16, RFLAG=1, SPAWN=1 )

}
\alias{sla.multiply}
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
	Given the input matrices, the generic function returns the result of 
        multiplication operation performed using the ScaLAPACK library routine.
 
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
	The matrix-matrix multiplication function in R a\%*\%b. 

}

\examples{
rnorm(10000)->x
dim(x)=c(100,100)
y<-matrix(rnorm(100*100),100,100)
z<-sla.multiply(x,y, NPROWS=2, NPCOLS=2, MB=10);
}
\author{Nagiza Samatova (\email{samatovan@ornl.gov}), Guruprasad Kora (\email{koragh@ornl.gov}), Srikanth Yoginath (\email{yoginathsb@ornl.gov}), David Bauer (\email{bauerda@ornl.gov})}
\concept{parallel}
\keyword{algebra}
