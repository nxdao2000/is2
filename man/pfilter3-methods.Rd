\name{pfilter3-methods}
\alias{logLik,pfilter3d.pomp-method}
\alias{logLik-pfilter3d.pomp}
\alias{$,pfilter3d.pomp-method}
\alias{$-pfilter3d.pomp}
\alias{pred.mean}
\alias{pred.mean,pfilter3d.pomp-method}
\alias{pred.mean-pfilter3d.pomp}
\alias{pred.var}
\alias{pred.var,pfilter3d.pomp-method}
\alias{pred.var-pfilter3d.pomp}
\alias{filter.mean}
\alias{filter.mean,pfilter3d.pomp-method}
\alias{filter.mean-pfilter3d.pomp}
\alias{eff.sample.size}
\alias{eff.sample.size,pfilter3d.pomp-method}
\alias{eff.sample.size-pfilter3d.pomp}
\alias{cond.logLik}
\alias{cond.logLik,pfilter3d.pomp-method}
\alias{cond.logLik-pfilter3d.pomp}
\alias{as,pfilter3d.pomp-method}
\alias{as.data.frame.pfilter3d.pomp}
\alias{coerce,pfilter3d.pomp,data.frame-method}
\title{Methods of the "pfilter3d.pomp" class}
\description{Methods of the "pfilter3d.pomp" class.}
\usage{
\S4method{logLik}{pfilter3d.pomp}(object, \dots)
\S4method{pred.mean}{pfilter3d.pomp}(object, pars, \dots)
\S4method{pred.var}{pfilter3d.pomp}(object, pars, \dots)
\S4method{filter.mean}{pfilter3d.pomp}(object, pars, \dots)
\S4method{eff.sample.size}{pfilter3d.pomp}(object, \dots)
\S4method{cond.logLik}{pfilter3d.pomp}(object, \dots)
\S4method{as}{pfilter3d.pomp}(object, class)
\S4method{coerce}{pfilter3d.pomp,data.frame}(from, to = "data.frame", strict = TRUE)
\S3method{as.data.frame}{pfilter3d.pomp}(x, row.names, optional, \dots)
}
\arguments{
  \item{object, x}{
    An object of class \code{pfilter3d.pomp} or inheriting class \code{pfilter3d.pomp}.
  }
  \item{pars}{Names of parameters.}
  \item{class}{
    character;
    name of the class to which \code{object} should be coerced.
  }
  \item{from, to}{
    the classes between which coercion should be performed.
  }
  \item{strict}{
    ignored.
  }
  \item{row.names, optional}{
    ignored.
  }
  \item{\dots}{
    Additional arguments unused at present.
  }
}
\author{{Dao Nguyen \email{dxnguyen at olemiss dot edu}}, {Xin Dang \email{xdang at olemiss dot edu}}, {Duc Anh Doan \email{ddoan at olemiss dot edu}}}
\seealso{\code{\link{is2}}}
\seealso{\code{\link{pfilter3}}}
\keyword{models}
\keyword{ts}
