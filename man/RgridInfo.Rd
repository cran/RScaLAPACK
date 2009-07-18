\name{sla.gridInfo}
\alias{sla.gridInfo}
\title{Retrieve the Existing Process Grid Information}
\usage{
sla.gridInfo()
}
\description{
	Retrieve the existing process grid information stored in the environment.
}
\details{
	Retrieve the existing process grid information stored in the environment.
}

\references{
  \url{http://www.netlib.org/scalapack/}
}

\seealso{
  \code{\link{sla.gridInit}} for initializing the Process Grid,
  \code{\link{sla.gridExit}} for releasing the Process Grid,
}

\examples{

x = matrix(rnorm(128*128), 128, 128)
sla.gridInit(NPROWS=2, NPCOLS=2)
sla.gridInfo()
sla.gridExit()

}
\author{Nagiza Samatova (\email{samatovan@ornl.gov}), Guruprasad Kora (\email{koragh@ornl.gov}), Srikanth Yoginath (\email{yoginathsb@ornl.gov}), David Bauer (\email{bauerda@ornl.gov})}
\concept{parallel}
\keyword{algebra}
