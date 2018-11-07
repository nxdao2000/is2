\name{aif-methods}
\docType{methods}
\alias{aif-methods}
\alias{logLik,aif-method}
\alias{logLik-aif}
\alias{conv.rec}
\alias{conv.rec,aif-method}
\alias{conv.rec-aif}
\alias{conv.rec,aifList-method}
\alias{conv.rec-aifList}
\alias{plot-aif}
\alias{plot,aif-method}
\alias{plot-aifList}
\alias{plot,aifList-method}
\alias{aifList-class}
\alias{c-aif}
\alias{c,aif-method}
\alias{c-aifList}
\alias{c,aifList-method}
\alias{[-aifList}
\alias{[,aifList-method}
\alias{compare.aif}
\title{Methods of the "aif" class}
\description{Methods of the \code{aif} class.}
\usage{
\S4method{logLik}{aif}(object, \dots)
\S4method{conv.rec}{aif}(object, pars, transform = FALSE, \dots)
\S4method{conv.rec}{aifList}(object, \dots)
\S4method{plot}{aif}(x, y, \dots)
\S4method{plot}{aifList}(x, y, \dots)
\S4method{c}{aif}(x, \dots, recursive = FALSE)
\S4method{c}{aifList}(x, \dots, recursive = FALSE)
compare.aif(z)
}
\arguments{
  \item{object}{The \code{aif} object.}
  \item{pars}{Names of parameters.}
  \item{x}{The \code{aif} object.}
  \item{y, recursive}{Ignored.}
  \item{z}{A \code{aif} object or list of \code{aif} objects.}
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
      Concatenates \code{aif} objects into an \code{aifList}.
    }
    \item{plot}{
      Plots a series of diagnostic plots.
    }
    \item{compare.aif}{
      Deprecated: use \code{plot} instead.
    }
  }
}
\author{{Dao Nguyen \email{dxnguyen at olemiss dot edu}}, {Xin Dang \email{xdang at olemiss dot edu}}, {Duc Anh Doan \email{ddoan at olemiss dot edu}}}
\seealso{\code{\link{aif}}, \code{\link{pfilter2}}}
\keyword{models}
\keyword{ts}
