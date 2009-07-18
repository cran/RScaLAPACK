\name{sla.gridInit}
\alias{sla.gridInit}
\title{Initialize the Process Grid for ScaLAPACK Computations}
\usage{
sla.gridInit(NPROWS=1, NPCOLS=1, BLOCKSIZE=16)
}
\arguments{
	\item{NPROWS}{ Number of Process Rows in the process grid }

	\item{NPCOLS}{ Number of Process Columns in the process grid }

	\item{BLOCKSIZE}{Block Size. }
}
\description{
	Sets up grid for use by the sla.* functions by spawning the requested number of processes.
}
\details{
	The requested number of processes are spawned and left waiting for further instruction.  If a grid a already present, the function instead prints and error and returns.  The grid can then be used by sla.* functions by setting the \code{SPAWN} parameter of those function to zero (most of the time, \code{RFLAG} should also be set to zero).  Functions do not have to use the full grid, and the shape doesn't matter.  For example if user spawns 6 processes by "sla.gridInit(6)",
 it can be used as a "2 x 2", "1 x 6", "2 x 3", "3 x 2", etc process grid.
}

\note{
	This was implemented for two reasons.  One, it slightly reduces the overhead of the call.  Second, it allows many ScaLAPACK calls without needing to open lots of sockets (as LAM/MPI does for each process spawned).  We strongly recommend that you release the process grid before exiting R, although in our tests the MPI daemon usually closed the child processes when the R session was closed.
}

\references{
  \url{http://www.netlib.org/scalapack/}
}

\seealso{
  \code{\link{sla.gridExit}} for releasing the Process Grid,
}

\examples{

x = matrix(rnorm(128*128), 128, 128)
sla.gridInit(NPROWS=2, NPCOLS=2)
sla.solve(x, NPROWS=2, NPCOLS=2, RFLAG=FALSE, SPAWN=FALSE)->a
sla.gridExit()

}
\author{Nagiza Samatova (\email{samatovan@ornl.gov}), Guruprasad Kora (\email{koragh@ornl.gov}), Srikanth Yoginath (\email{yoginathsb@ornl.gov}), David Bauer (\email{bauerda@ornl.gov})}
\concept{parallel}
\keyword{algebra}
