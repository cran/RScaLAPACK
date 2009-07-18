\name{sla.prcomp}
\title{Principal Components Analysis (Using ScaLAPACK)}
\usage{
sla.prcomp(x, retx=TRUE, center=TRUE, scale.=FALSE, tol=NULL, NPROWS, NPCOLS, MB)
}
\alias{sla.prcomp}
\arguments{
	\item{x}{A numeric matrix which provides the data for the principal components analysis.}

	\item{retx}{a logical value indicating whether the rotated variables should be returned.}

	\item{center}{a logical value indicating whether the variables should be shifted to be zero centered.  Alternately, a vector of length ncols(x).  The value is passed to \code{scale}.}

	\item{scale.}{a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place.  The default is \code{FALSE} for old compatibility, but in general scaling is advisable.  Alternately, a vector of length ncols(x).  The value is passed to \code{scale}.}

	\item{tol}{a value indicating the magnitude below which components should be omitted.  The default is to omit no components.}
	
	\item{NPROWS}{Number of Process Rows in the Process Grid. }

	\item{NPCOLS}{Number of Process Cols in the Process Grid. }

	\item{MB}{Block Size. }
}

\description{
	The function \code{sla.prcomp} performs a principal components analysis on the give data matrix (using \code{sla.svd}) and returns the results as a object of class \code{prcomp}.}

\details{
	If the number of processor rows and columns are both zero, one processor
	is used.  If the number of processor rows is nonzero and the number of
	processor columns is zero, then the processor rows is taken to be a 
	number of processors, and a grid is made accordingly.
}

\value{
	\item{sdev}{the standard deviations of the principal components.}
	\item{rotation}{the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors).}
	\item{x}{if \code{retx} is true, the value of the rotated data (the original data multiplied by the \code{rotation} matrix).}
	}

\references{
	\url{http://www.google.com/search?q=principal+components+analysis}

	\url{http://www.netlib.org/scalapack/}
}

\seealso{
	\code{\link{prcomp}} for the normal PCA function.
}

\examples{
data(USArrests)
sla.prcomp(USArrests, scale=TRUE, NPROWS=2, NPCOLS=2, MB=2)
}
\author{Nagiza Samatova (\email{samatovan@ornl.gov}), Guruprasad Kora (\email{koragh@ornl.gov}), Srikanth Yoginath (\email{yoginathsb@ornl.gov}), David Bauer (\email{bauerda@ornl.gov})}
\concept{parallel}
\keyword{algebra}
