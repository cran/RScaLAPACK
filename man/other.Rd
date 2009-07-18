\name{sla}
\alias{PA.exec}
\alias{PA.exit}
\alias{sla.pdgesv}
\alias{sla.ProcessRowsColumns}
\alias{sla.checkBlockSize}
\alias{sla}
\alias{scalapack}
\alias{.RscalaGrid}
\title{R-ScaLAPACK}
\description{ The R-ScaLAPACK package provides a backwards compatible interface to the ScaLACPACK version of several R functions which presently use LIN/LAPACK.}
\details{
	The ScaLAPACK library performs many common linear algebra functions in parallel on distributed or shared memory architectures.  It is based on the LAPACK library.  The R-ScaLAPACK package provides an interface to the most common ScaLAPACK functions available in an R compatible form.

	All of the functions are named "sla.*", and the take the same (or as similar as possible) parameters as the normal R version.  All functions have the additional parameters of the process grid and block size.  A default process grid may be placed in the variable ".RscalaGrid", which will be used by the sla.* functions when no grid is given in their parameters.}
	
\seealso{
	 \code{\link{sla.chol}}
	 
	 \code{\link{sla.chol2inv}}
	 
	 \code{\link{sla.eigen}}
	 
	 \code{\link{sla.prcomp}}
	 
	 \code{\link{sla.qr}}
	 
	 \code{\link{sla.solve}}

	 \code{\link{sla.svd}}
	 }
\author{Nagiza Samatova (\email{samatovan@ornl.gov}), Guruprasad Kora (\email{koragh@ornl.gov}), Srikanth Yoginath (\email{yoginathsb@ornl.gov}), David Bauer (\email{bauerda@ornl.gov})}
\keyword{internal}
\concept{parallel}
\keyword{algebra}
