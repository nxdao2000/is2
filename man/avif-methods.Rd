\name{avif-methods}
\docType{methods}
\alias{avif-methods}
\alias{logLik,avif-method}
\alias{logLik-avif}
\alias{conv.rec}
\alias{conv.rec,avif-method}
\alias{conv.rec-avif}
\alias{conv.rec,avifList-method}
\alias{conv.rec-avifList}
\alias{plot-avif}
\alias{plot,avif-method}
\alias{plot-avifList}
\alias{plot,avifList-method}
\alias{avifList-class}
\alias{c-avif}
\alias{c,avif-method}
\alias{c-avifList}
\alias{c,avifList-method}
\alias{[-avifList}
\alias{[,avifList-method}
\alias{compare.avif}
\title{Methods of the "avif" class}
\description{Methods of the \code{avif} class.}
\usage{
\S4method{logLik}{avif}(object, \dots)
\S4method{conv.rec}{avif}(object, pars, transform = FALSE, \dots)
\S4method{conv.rec}{avifList}(object, \dots)
\S4method{plot}{avif}(x, y, \dots)
\S4method{plot}{avifList}(x, y, \dots)
\S4method{c}{avif}(x, \dots, recursive = FALSE)
\S4method{c}{avifList}(x, \dots, recursive = FALSE)
compare.avif(z)
}
\arguments{
  \item{object}{The \code{avif} object.}
  \item{pars}{Names of parameters.}
  \item{x}{The \code{avif} object.}
  \item{y, recursive}{Ignored.}
  \item{z}{A \code{avif} object or list of \code{avif} objects.}
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
      Concatenates \code{avif} objects into an \code{avifList}.
    }
    \item{plot}{
      Plots a series of diagnostic plots.
    }
    \item{compare.avif}{
      Deprecated: use \code{plot} instead.
    }
  }
}
\author{{Dao Nguyen \email{dxnguyen at olemiss dot edu}}, {Xin Dang \email{xdang at olemiss dot edu}}, {Duc Anh Doan \email{ddoan at olemiss dot edu}}}
\seealso{\code{\link{avif}}, \code{\link{pfilter2}}}
\keyword{models}
\keyword{ts}
