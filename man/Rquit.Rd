\name{sla.quit}
\alias{sla.quit}
\title{Releases the Parallel Environment and quits R}
\usage{
  sla.quit()
}
\description{
	Always call this function to quit R after using RScLAPACK library. This function is used to release the parallel environment and quit R.
}
\details{
	Calls \code{MPI_Finalize} and quits R.
}

\references{
  \url{http://www.aspect-sdm.org/Parallel-R}
}
\author{Nagiza Samatova (\email{samatovan@ornl.gov}), Guruprasad Kora (\email{koragh@ornl.gov}), Srikanth Yoginath (\email{yoginathsb@ornl.gov}), David Bauer (\email{bauerda@ornl.gov})}
\concept{parallel}
\keyword{internal}
