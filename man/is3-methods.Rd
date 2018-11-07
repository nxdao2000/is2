\name{is3-methods}
\docType{methods}
\alias{is3-methods}
\alias{logLik,is3-method}
\alias{logLik-is3}
\alias{conv.rec}
\alias{conv.rec,is3-method}
\alias{conv.rec-is3}
\alias{conv.rec,is3List-method}
\alias{conv.rec-is3List}
\alias{plot-is3}
\alias{plot,is3-method}
\alias{plot-is3List}
\alias{plot,is3List-method}
\alias{is3List-class}
\alias{c-is3}
\alias{c,is3-method}
\alias{c-is3List}
\alias{c,is3List-method}
\alias{[-is3List}
\alias{[,is3List-method}
\alias{compare.is3}
\title{Methods of the "is3" class}
\description{Methods of the \code{is3} class.}
\usage{
\S4method{logLik}{is3}(object, \dots)
\S4method{conv.rec}{is3}(object, pars, transform = FALSE, \dots)
\S4method{conv.rec}{is3List}(object, \dots)
\S4method{plot}{is3}(x, y, \dots)
\S4method{plot}{is3List}(x, y, \dots)
\S4method{c}{is3}(x, \dots, recursive = FALSE)
\S4method{c}{is3List}(x, \dots, recursive = FALSE)
compare.is3(z)
}
\arguments{
  \item{object}{The \code{is3} object.}
  \item{pars}{Names of parameters.}
  \item{x}{The \code{is3} object.}
  \item{y, recursive}{Ignored.}
  \item{z}{A \code{is3} object or list of \code{is3} objects.}
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
      Concatenates \code{is3} objects into an \code{is3List}.
    }
    \item{plot}{
      Plots a series of diagnostic plots.
    }
    \item{compare.is3}{
      Deprecated: use \code{plot} instead.
    }
  }
}
\author{{Dao Nguyen \email{dxnguyen at olemiss dot edu}}, {Xin Dang \email{xdang at olemiss dot edu}}, {Duc Anh Doan \email{ddoan at olemiss dot edu}}}
\seealso{\code{\link{is3}}, \code{\link{pfilter2}}}
\keyword{models}
\keyword{ts}
