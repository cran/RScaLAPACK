\name{sla.princomp}
\title{Principal Components Analysis (Using ScaLAPACK)}
\usage{
sla.princomp(x, \dots)
\method{sla.princomp}{formula}(formula, data = NULL, subset, na.action, NPROWS=0, NPCOLS=0, MB=16, \dots)
\method{sla.princomp}{default}(x, cor = FALSE, scores = TRUE, covmat = NULL, subset = rep(TRUE, nrow(as.matrix(x))), NPROWS=0, NPCOLS=0, MB=16, \dots)
}
\alias{sla.princomp}
\alias{sla.princomp.default}
\alias{sla.princomp.formula}
\arguments{
	\item{formula}{a formula with no response variable.}
	\item{data}{an optional data frame containing the variables in the
		formula \code{formula}. By default the variables are taken from
			\code{environment(formula)}.}
	\item{x}{a matrix or data frame which provides the data for the
		principal components analysis.}
	\item{subset}{an optional vector used to select rows (observations) of the
		data matrix \code{x}.}
	\item{na.action}{a function which indicates what should happen
		when the data contain \code{NA}s.  The default is set by
			the \code{na.action} setting of \code{\link{options}}, and is
			\code{\link{na.fail}} if that is unset. The \dQuote{factory-fresh}
		default is \code{\link{na.omit}}.}
	\item{cor}{a logical value indicating whether the calculation should
		use the correlation matrix or the covariance matrix.}
	\item{scores}{a logical value indicating whether the score on each
		principal component should be calculated.}
	\item{covmat}{a covariance matrix, or a covariance list as returned by
		\code{\link{cov.wt}}, \code{\link[MASS]{cov.mve}} or \code{\link[MASS]{cov.mcd}}.
			If supplied, this is used rather than the covariance matrix of
			\code{x}.}
	\item{\dots}{arguments passed to or from other methods. If \code{x} is
		a formula one might specify \code{cor} or \code{scores}.}

	\item{NPROWS}{Number of Process Rows in the Process Grid. }

	\item{NPCOLS}{Number of Process Cols in the Process Grid. }

	\item{MB}{Block Size. }
}

\description{
	The function \code{sla.prcomp} performs a principal components analysis on the give data matrix (using \code{sla.eigen}) and returns the results as a object of class \code{princomp}.}

\details{
	If the number of processor rows and columns are both zero, one processor
	is used.  If the number of processor rows is nonzero and the number of
	processor columns is zero, then the processor rows is taken to be a 
	number of processors, and a grid is made accordingly.
}

\value{
  \code{sla.princomp} returns a list with class \code{"princomp"}
  containing the following components:
  \item{sdev}{the standard deviations of the principal components.}
  \item{loadings}{the matrix of variable loadings (i.e., a matrix
    whose columns contain the eigenvectors).  This is of class
    \code{"loadings"}: see \code{\link{loadings}} for its \code{print}
    method.}
  \item{center}{the means that were subtracted.}
  \item{scale}{the scalings applied to each variable.}
  \item{n.obs}{the number of observations.}
  \item{scores}{if \code{scores = TRUE}, the scores of the supplied
    data on the principal components.}
  \item{call}{the matched call.}
  \item{na.action}{If relevant.}
}

\references{
  \code{\link{summary.princomp}}, \code{\link{screeplot}},
  \code{\link{biplot.princomp}},
  \code{\link{prcomp}}, \code{\link{cor}}, \code{\link{cov}},
  \code{\link{eigen}}.

	\url{http://www.google.com/search?q=principal+components+analysis}

	\url{http://www.netlib.org/scalapack/}
}

\seealso{
	\code{\link{prcomp}} for the normal PCA function.
}

\examples{
data(USArrests)
sla.princomp(USArrests, cor=TRUE, NPROWS=2, NPCOLS=2, MB=2)
}
\author{Nagiza Samatova (\email{samatovan@ornl.gov}), Guruprasad Kora (\email{koragh@ornl.gov}), Srikanth Yoginath (\email{yoginathsb@ornl.gov}), David Bauer (\email{bauerda@ornl.gov})}
\concept{parallel}
\keyword{multivariate}
