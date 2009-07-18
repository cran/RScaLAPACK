\name{sla.gridExit}
\alias{sla.gridExit}
\title{Releases the Process Grid for ScaLAPACK Computations}
\usage{
  sla.gridExit()
}
\description{
	Releases the grid used by sla.* functions.
}
\details{
	This function sends out a request to all the processes, telling them to exit.  The user is then free to spawn a new grid (of a different size).  It does NOT call \code{MPI_Finalize} or anything like that.
}

\references{
  \url{http://www.netlib.org/scalapack/}
}

\seealso{
  \code{\link{sla.gridInit}} spawns the processes to form a Process Grid,
}

\examples{

x = matrix(rnorm(128*128), 128, 128)
sla.gridInit(NPROCS=4)
sla.solve(x, NPROWS=4)->a
sla.gridExit()

}
\author{Nagiza Samatova (\email{samatovan@ornl.gov}), Guruprasad Kora (\email{koragh@ornl.gov}), Srikanth Yoginath (\email{yoginathsb@ornl.gov}), David Bauer (\email{bauerda@ornl.gov})}
\concept{parallel}
\keyword{algebra}
