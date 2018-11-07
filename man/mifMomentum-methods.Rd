\name{mifMomentum-methods}
\docType{methods}
\alias{mifMomentum-methods}
\alias{logLik,mifMomentum-method}
\alias{logLik-mifMomentum}
\alias{conv.rec}
\alias{conv.rec,mifMomentum-method}
\alias{conv.rec-mifMomentum}
\alias{conv.rec,mifMomentumList-method}
\alias{conv.rec-mifMomentumList}
\alias{plot-mifMomentum}
\alias{plot,mifMomentum-method}
\alias{plot-mifMomentumList}
\alias{plot,mifMomentumList-method}
\alias{mifMomentumList-class}
\alias{c-mifMomentum}
\alias{c,mifMomentum-method}
\alias{c-mifMomentumList}
\alias{c,mifMomentumList-method}
\alias{[-mifMomentumList}
\alias{[,mifMomentumList-method}
\alias{compare.mifMomentum}
\title{Methods of the "mifMomentum" class}
\description{Methods of the \code{mifMomentum} class.}
\usage{
\S4method{logLik}{mifMomentum}(object, \dots)
\S4method{conv.rec}{mifMomentum}(object, pars, transform = FALSE, \dots)
\S4method{conv.rec}{mifMomentumList}(object, \dots)
\S4method{plot}{mifMomentum}(x, y, \dots)
\S4method{plot}{mifMomentumList}(x, y, \dots)
\S4method{c}{mifMomentum}(x, \dots, recursive = FALSE)
\S4method{c}{mifMomentumList}(x, \dots, recursive = FALSE)
compare.mifMomentum(z)
}
\arguments{
  \item{object}{The \code{mifMomentum} object.}
  \item{pars}{Names of parameters.}
  \item{x}{The \code{mifMomentum} object.}
  \item{y, recursive}{Ignored.}
  \item{z}{A \code{mifMomentum} object or list of \code{mifMomentum} objects.}
  \item{transform}{
    optional logical;
    should the parameter transformations be applied?
    See \code{\link[=coef-pomp]{coef}} for details.
  }
  \item{\dots}{
    Further arguments (either ignored or passed to underlying functions).
  }
}
\section{Methods}{
  \describe{
    \item{conv.rec}{
      \code{conv.rec(object, pars = NULL)} returns the columns of the convergence-record matrix corresponding to the names in \code{pars}.
      By default, all rows are returned.
    }
    \item{logLik}{
      Returns the value in the \code{loglik} slot.
    }
    \item{c}{
      Concatenates \code{mifMomentum} objects into an \code{mifMomentumList}.
    }
    \item{plot}{
      Plots a series of diagnostic plots.
    }
    \item{compare.mifMomentum}{
      Deprecated: use \code{plot} instead.
    }
  }
}
\author{{Dao Nguyen \email{dxnguyen at olemiss dot edu}}, {Xin Dang \email{xdang at olemiss dot edu}}, {Duc Anh Doan \email{ddoan at olemiss dot edu}}}
\seealso{\code{\link{mifMomentum}}, \code{\link{pfilter2}}}
\keyword{models}
\keyword{ts}
