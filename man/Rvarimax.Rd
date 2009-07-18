\name{sla.varimax}
\alias{sla.promax}
\alias{sla.varimax}
\title{Rotation Methods for Factor Analysis}
\description{
  These functions \sQuote{rotate} loading matrices in factor analysis using ScaLAPACK.
}
\usage{
sla.varimax(x, normalize = TRUE, eps = 1e-5, NPROWS=0, NPCOLS=0, MB=16)
sla.promax(x, m = 4, NPROWS=0, NPCOLS=0, MB=16)
}
\arguments{
	\item{x}{A loadings matrix, with \eqn{p} rows and \eqn{k < p} columns}
	\item{m}{The power used the target for \code{promax}.  Values of 2 to
		4 are recommended.}
	\item{normalize}{logical. Should Kaiser normalization be performed?
		If so the rows of \code{x} are re-scaled to unit length before
			rotation, and scaled back afterwards.}
	\item{eps}{The tolerance for stopping: the relative change in the sum
		of singular values.}
	\item{NPROWS}{Number of Process Rows in the Process Grid. }
	\item{NPCOLS}{Number of Process Cols in the Process Grid. }
	\item{MB}{Block Size. }
}
\details{
  These seek a \sQuote{rotation} of the factors \code{x \%*\% T} that
  aims to clarify the structure of the loadings matrix.  The matrix
  \code{T} is a rotation (possibly with reflection) for \code{varimax},
  but a general linear transformation for \code{promax}, with the
  variance of the factors being preserved.
}
\value{
  A list with components
  \item{loadings}{The \sQuote{rotated} loadings matrix,
    \code{x \%*\% rotmat}.}
  \item{rotmat}{The \sQuote{rotation} matrix.}
}
\references{
  Hendrickson, A. E. and White, P. O. (1964) Promax: a quick method for
  rotation to orthogonal oblique structure. \emph{British Journal of
    Statistical Psychology}, \bold{17}, 65--70.

  Horst, P. (1965) \emph{Factor Analysis of Data Matrices.} Holt,
  Rinehart and Winston.  Chapter 10.

  Kaiser, H. F. (1958) The varimax criterion for analytic rotation in
  factor analysis. \emph{Psychometrika} \bold{23}, 187--200.

  Lawley, D. N. and Maxwell, A. E. (1971) \emph{Factor Analysis as a
    Statistical Method}. Second edition. Butterworths.
}
\seealso{
	\code{\link{varimax}} for the normal varmax,
	\code{\link{promax}} for the normal promax,
	\code{\link{sla.factanal}},
	\code{\link{factanal}}, \code{\link{Harman74.cor}}.}

\examples{
## varimax with normalize = T is the default

v1 <- c(1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,4,5,6)
v2 <- c(1,2,1,1,1,1,2,1,2,1,3,4,3,3,3,4,6,5)
v3 <- c(3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,5,4,6)
v4 <- c(3,3,4,3,3,1,1,2,1,1,1,1,2,1,1,5,6,4)
v5 <- c(1,1,1,1,1,3,3,3,3,3,1,1,1,1,1,6,4,5)
v6 <- c(1,1,1,2,1,3,3,3,4,3,1,1,1,2,1,6,5,4)
m1 <- cbind(v1,v2,v3,v4,v5,v6)
fa <- sla.factanal(m1, factors=3, rotation="sla.promax", NPROWS=2, NPCOLS=2, MB=2 )
sla.varimax(fa$loadings, normalize = FALSE, NPROWS=2, NPCOLS=2, MB=2)
sla.promax(fa$loadings, NPROWS=2, NPCOLS=2, MB=2)
}
\concept{parallel}
\keyword{multivariate}
\author{Nagiza Samatova (\email{samatovan@ornl.gov}), Guruprasad Kora (\email{koragh@ornl.gov}), Srikanth Yoginath (\email{yoginathsb@ornl.gov}), David Bauer (\email{bauerda@ornl.gov})}
