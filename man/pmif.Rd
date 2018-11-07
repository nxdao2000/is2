\name{Particle Iterated Filtering}
\title{The particle Iterated Filtering algorithm}
\docType{methods}
\alias{pmif}
\alias{pmif,missing-method}
\alias{pmif,ANY-method}
\alias{pmif,pmif-method}
\alias{pmif-pmif}
\alias{pmif,pfilterd.pomp-method}
\alias{pmif-pfilterd.pomp}
\alias{pmif,pomp-method}
\alias{pmif-pomp}
\alias{continue,pmif-method}
\alias{continue-pmif}
\alias{pmif-class}
\alias{pmif-methods}
\alias{logLik,pmif-method}
\alias{logLik-pmif}
\alias{logLik,pmifList-method}
\alias{logLik-pmifList}
\alias{conv.rec,pmif-method}
\alias{conv.rec-pmif}
\alias{conv.rec,pmifList-method}
\alias{conv.rec-pmifList}
\alias{filter.traj,pmif-method}
\alias{filter.traj,pmifList-method}
\alias{filter.traj-pmif}
\alias{filter.traj-pmifList}
\alias{plot-pmif}
\alias{plot,pmif-method}
\alias{plot-pmifList}
\alias{plot,pmifList-method}
\alias{pmifList-class}
\alias{covmat,pmif-method}
\alias{covmat-pmif}
\alias{covmat,pmifList-method}
\alias{covmat-pmifList}
\alias{c-pmif}
\alias{c,pmif-method}
\alias{c-pmifList}
\alias{c,pmifList-method}
\alias{[-pmifList}
\alias{[,pmifList-method}
\description{
  The Particle MCMC algorithm for estimating the parameters of a partially-observed Markov process.
  Running \code{pmif} causes a particle random-walk Metropolis-Hastings Markov chain algorithm to run for the specified number of proposals.
}
\usage{
\S4method{pmif}{pomp}(object, Nmif = 1, start, proposal, Np,
    tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"), \dots)
\S4method{pmif}{pfilterd.pomp}(object, Nmif = 1, Np, tol, \dots)
\S4method{pmif}{pmif}(object, Nmif, start, proposal, Np, tol,
    max.fail = Inf, verbose = getOption("verbose"), \dots)
\S4method{continue}{pmif}(object, Nmif = 1, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp}.
  }
  \item{Nmif}{
    The number of pmif iterations to perform.
  }
  \item{start}{
    named numeric vector;
    the starting guess of the parameters.
  }
  \item{proposal}{
    optional function that draws from the proposal distribution.
    Currently, the proposal distribution must be symmetric for proper inference:
    it is the user's responsibility to ensure that it is.
    Several functions that construct appropriate proposal function are provided:
    see \link{MCMC proposal functions} for more information.
  }
  \item{Np}{
    a positive integer;
    the number of particles to use in each filtering operation.
  }
  \item{tol}{
    numeric scalar; particles with log likelihood below \code{tol} are considered to be \dQuote{lost}.
    A filtering failure occurs when, at some time point, all particles are lost.
  }
  \item{max.fail}{
    integer; maximum number of filtering failures permitted.
    If the number of failures exceeds this number, execution will terminate with an error.
  }
  \item{verbose}{
    logical; if TRUE, print progress reports.
  }
  \item{\dots}{
    additional arguments that override the defaults.
  }
}
\value{
  An object of class \code{pmif}.
}
\section{Re-running pmif Iterations}{
  To re-run a sequence of pmif iterations, one can use the \code{pmif} method on a \code{pmif} object.
  By default, the same parameters used for the original pmif run are re-used (except for \code{tol}, \code{max.fail}, and \code{verbose}, the defaults of which are shown above).
  If one does specify additional arguments, these will override the defaults.
}
\section{Continuing pmif Iterations}{
  One can continue a series of pmif iterations from where one left off using the \code{continue} method.
  A call to \code{pmif} to perform \code{Nmif=m} iterations followed by a call to \code{continue} to perform \code{Nmif=n} iterations will produce precisely the same effect as a single call to \code{pmif} to perform \code{Nmif=m+n} iterations.
  By default, all the algorithmic parameters are the same as used in the original call to \code{pmif}.
  Additional arguments will override the defaults.
}
\section{Details}{
  \code{pmif} implements an MCMC algorithm in which the true likelihood of the data is replaced by an unbiased estimate computed by a particle filter.
  This gives an asymptotically correct Bayesian procedure for parameter estimation (Andrieu and Roberts, 2009).

  \strong{Note} that \code{pmif} does not make use of any parameter transformations supplied by the user.
}
\section{Methods}{
  \describe{
    \item{c}{
      Concatenates \code{pmif} objects into a \code{pmifList}.
    }
    \item{\code{conv.rec(object, pars)}}{
      returns the columns of the convergence-record matrix corresponding to the names in \code{pars} as an object of class \code{\link[coda]{mcmc}} or \code{\link[coda]{mcmc.list}}.
    }
    \item{\code{filter.traj(object, vars)}}{
      returns filter trajectories from a \code{pmif} or \code{pmifList} object.
    }
    \item{plot}{
      Diagnostic plots.
    }
    \item{logLik}{
      Returns the value in the \code{loglik} slot.
    }
    \item{coef}{
      Returns the last state of the MCMC chain.
      As such, it's not very useful for inference.
    }
    \item{\code{covmat(object, start, thin, expand)}}{
      computes the empirical covariance matrix of the MCMC samples beginning with iteration \code{start} and thinning by factor \code{thin}.
      It expands this by a factor \code{expand^2/n}, where \code{n} is the number of parameters estimated.
      By default, \code{expand=2.38}.
      The intention is that the resulting matrix is a suitable input to the proposal function \code{\link{mvn.rw}}.
    }
  }
}
\examples{
\dontrun{
library(pomp)

pompExample(ou2)

pmif(
      pomp(ou2,dprior=Csnippet("
   lik = dnorm(alpha_2,-0.5,1,1) + dnorm(alpha_3,0.3,1,1);
   lik = (give_log) ? lik : exp(lik);"),
           paramnames=c("alpha.2","alpha.3")),
      Nmif=2000,Np=500,verbose=TRUE,
      proposal=mvn.rw.adaptive(rw.sd=c(alpha.2=0.01,alpha.3=0.01),
        scale.start=200,shape.start=100)) -> chain
continue(chain,Nmif=2000,proposal=mvn.rw(covmat(chain))) -> chain
plot(chain)
chain <- pmif(chain)
plot(chain)

library(coda)
trace <- window(conv.rec(chain,c("alpha.2","alpha.3")),start=2000)
rejectionRate(trace)
effectiveSize(trace)
autocorr.diag(trace)

summary(trace)
plot(trace)

heidel.diag(trace)
geweke.diag(trace)
}
}
\seealso{
  \code{\link{pomp}}, \code{\link{pfilter}}, \link{MCMC proposal distributions},
  and the tutorials on the \href{https://kingaa.github.io/pomp}{package website}.
}
\references{
  C. Andrieu, A. Doucet and R. Holenstein,
  Particle Markov chain Monte Carlo methods,
  J. R. Stat. Soc. B, to appear, 2010.

  C. Andrieu and G.O. Roberts,
  The pseudo-marginal approach for efficient computation,
  Ann. Stat. 37:697-725, 2009.
}
\author{{Dao Nguyen \email{dxnguyen at olemiss dot edu}}, {Xin Dang \email{xdang at olemiss dot edu}}, {Duc Anh Doan \email{ddoan at olemiss dot edu}}}
\keyword{ts}
