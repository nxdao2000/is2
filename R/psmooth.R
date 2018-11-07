## particle filtering codes

setClass(
  "psmooth.pomp",
  contains="pomp",
  slots=c(
    pred.mean="array",
    pred.var="array",
    filter.mean="array",
    filter.traj="array",
    paramMatrix="array",
    indices="vector",
    eff.sample.size="numeric",
    cond.loglik="numeric",
    saved.states="list",
    saved.params="list",
    Np="integer",
    tol="numeric",
    nfail="integer",
    loglik="numeric",
    phats="numeric",
    covhats="array",
    pcovhats="array",
    lag = "numeric"
  ),
  prototype=prototype(
    pred.mean=array(data=numeric(0),dim=c(0,0)),
    pred.var=array(data=numeric(0),dim=c(0,0)),
    filter.mean=array(data=numeric(0),dim=c(0,0)),
    filter.traj=array(data=numeric(0),dim=c(0,0,0)),
    paramMatrix=array(data=numeric(0),dim=c(0,0)),
    indices=integer(0),
    eff.sample.size=numeric(0),
    cond.loglik=numeric(0),
    saved.states=list(),
    saved.params=list(),
    Np=as.integer(NA),
    tol=as.double(NA),
    nfail=as.integer(NA),
    loglik=as.double(NA),
    phats=numeric(0),
    covhats=array(data=numeric(0),dim=c(0,0)),
    pcovhats=array(data=numeric(0),dim=c(0,0,0,0)),
    lag = as.integer(NA)
  )
)
ancestor<-function(plist, t, lag, index){
    for (i in 0:(lag-1)){
        index=plist[[t-i]][index]
    }
    return(index)
}

smoothing<-function(aparticles, xparticles, pparticles, wparticles, nt, ntimes, lag, nvars, npars, rw, Np){
    at=rep(0,Np)        #current parent index
    bt=rep(0,Np)
    nlength<-length(rw)
    phat<-rep(0,npars)  #smoothed par
    pcovhat<-matrix(0,npars,lag) #covariance
    lcovhat<-array(0,dim=c(npars,npars,lag)) 
  
    if(lag>0){
        kk=nt+lag
        if(nt==ntimes){
        }
        else{
            if(kk<=ntimes){
                at=unlist(lapply(aparticles[[kk]],ancestor,  plist=aparticles, t=kk,lag=lag ))
                for ( jj in 1: npars ){
                    phat[jj]=sum(pparticles[[nt]][jj,at]*wparticles[[kk]])
                }
                for ( jj in 1: npars ){
                    for(nn in 1:(kk-nt)){
                        bt=unlist(lapply(aparticles[[kk]],ancestor,  plist=aparticles, t=kk,lag=lag+1-nn )) 
                        pcovhat[jj,nn] =sum(pparticles[[nt+nn]][jj,bt]*wparticles[[kk]])
                    }
                }
                for ( jj in 1: npars ){
                    for (ll in 1: npars){
                        for(nn in 1:lag){
                            bt=unlist(lapply(aparticles[[kk]],ancestor,  plist=aparticles, t=kk,lag=lag+1-nn ))
                            lcovhat[jj,ll,nn] = sum((pparticles[[nt+nn]][ll,bt]-pcovhat[ll,nn])*(pparticles[[nt]][jj,at]-phat[jj])*wparticles[[kk]])/(Np-1)
                        }
                    }
                }
            }  
            else{
                at=unlist(lapply(aparticles[[ntimes]],ancestor,  plist=aparticles, t=ntimes,lag=(ntimes-nt) ))
                kk=ntimes
                for ( jj in 1: npars ){
                    phat[jj]=sum(pparticles[[nt]][jj,at]*wparticles[[ntimes]])
                }
                for ( jj in 1: npars ){
                    for(nn in 1:(kk-nt)){
                        bt=unlist(lapply(aparticles[[kk]],ancestor,  plist=aparticles, t=kk,lag=nn )) 
                        pcovhat[jj,nn] =sum(pparticles[[nt]][jj,bt]*wparticles[[kk]])
                    }
                }
                for ( jj in 1: npars ){
                    for (ll in 1: npars){
                        for(nn in 1:(kk-nt)){
                            bt=unlist(lapply(aparticles[[kk]],ancestor,  plist=aparticles, t=kk,lag=nn ))
                            lcovhat[jj,ll,nn] = sum((pparticles[[nt]][ll,bt]-pcovhat[ll,nn])*(pparticles[[ntimes]][jj,at]-phat[jj])*wparticles[[kk]])/(Np-1)
                        }
                    }
                }
            }
        }
    }
    return(lcovhat)
}

psmooth.internal <- function (object, params, Np,
                              tol, max.fail,
                              pred.mean = FALSE,
                              pred.var = FALSE,
                              filter.mean = FALSE,
                              filter.traj = FALSE,
                              cooling, cooling.m, .wn=FALSE,.corr=FALSE, rd.sw, .transform,
                              verbose = FALSE,
                              save.states = FALSE,
                              save.params = FALSE, lag,
                              .getnativesymbolinfo = TRUE) {

  ep <- paste0("in ",sQuote("psmooth"),": ")
  if (missing(lag)) lag <- 0
  corr <- as.logical(.corr)
  wn <- as.logical(.wn)
  transform <- as.logical(.transform)
    
  object <- as(object,"pomp")
  pompLoad(object,verbose=verbose)

  gnsi.rproc <- gnsi.dmeas <- as.logical(.getnativesymbolinfo)
  pred.mean <- as.logical(pred.mean)
  pred.var <- as.logical(pred.var)
  filter.mean <- as.logical(filter.mean)
  filter.traj <- as.logical(filter.traj)
  verbose <- as.logical(verbose)
  save.states <- as.logical(save.states)
  save.params <- as.logical(save.params)
  
  
  
  if (length(params)==0)
    stop(ep,sQuote("params")," must be specified",call.=FALSE)

  if (missing(tol))
    stop(ep,sQuote("tol")," must be specified",call.=FALSE)

  one.par <- FALSE
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1

  if (missing(Np)) {
    if (is.matrix(params)) {
      Np <- ncol(params)
    } else {
      stop(ep,sQuote("Np")," must be specified",call.=FALSE)
    }
  }
  if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
      error = function (e) {
        stop(ep,"if ",sQuote("Np")," is a function, ",
             "it must return a single positive integer",call.=FALSE)
      }
    )
  }
  if (length(Np)==1)
    Np <- rep(Np,times=ntimes+1)
  else if (length(Np)!=(ntimes+1))
    stop(ep,sQuote("Np")," must have length 1 or length ",ntimes+1,call.=FALSE)
  if (any(Np<=0))
    stop(ep,"number of particles, ",sQuote("Np"),", must always be positive",call.=FALSE)
  if (!is.numeric(Np))
    stop(ep,sQuote("Np")," must be a number, a vector of numbers, or a function",call.=FALSE)
  Np <- as.integer(Np)
  if (is.matrix(params)) {
    if (!all(Np==ncol(params)))
      stop(ep,"when ",sQuote("params")," is provided as a matrix, do not specify ",
        sQuote("Np"),"!",call.=FALSE)
  }

  if (NCOL(params)==1) {        # there is only one parameter vector
    one.par <- TRUE
    coef(object) <- params     # set params slot to the parameters
    params <- as.matrix(params)
  }

  paramnames <- rownames(params)
  if (is.null(paramnames))
    stop(ep,sQuote("params")," must have rownames",call.=FALSE)

  init.x <- init.state(object,params=params,nsim=Np[1L])
  statenames <- rownames(init.x)
  nvars <- nrow(init.x)
  x <- init.x

  ## set up storage for saving samples from filtering distributions
  if (save.states | filter.traj) {
    xparticles <- setNames(vector(mode="list",length=ntimes),time(object))
  }
  if (save.params) {
    pparticles <- setNames(vector(mode="list",length=ntimes),time(object))
  } else {
    pparticles <- list()
  }
  if (filter.traj) {
    pedigree <- vector(mode="list",length=ntimes+1)
  }
  rw.names <- names(rw.sd)
  sigma <- rw.sd
  npars <- length(rw.names)
  nvars <- nrow(x)

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  nfail <- 0

  ## set up storage for prediction means, variances, etc.
  if (pred.mean) {
    pred.m <- matrix(
      data=0,
      nrow=nvars,
      ncol=ntimes,
      dimnames=list(
        variable=statenames,
        time=time(object))
    )
  } else {
    pred.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (pred.var) {
    pred.v <- matrix(
      data=0,
      nrow=nvars,
      ncol=ntimes,
      dimnames=list(
        variable=statenames,
        time=time(object))
    )
  } else {
    pred.v <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.mean) {
    filt.m <- matrix(
      data=0,
      nrow=nvars,
      ncol=ntimes,
      dimnames=list(
        variable=statenames,
        time=time(object))
    )
  } else {
    filt.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.traj) {
    filt.t <- array(
      data=0,
      dim=c(nvars,1,ntimes+1),
      dimnames=list(
        variable=statenames,
        rep=1,
        time=times)
    )
  } else {
    filt.t <- array(data=numeric(0),dim=c(0,0,0))
  }


    ##########################################
    # Fixed-lag Smoothing 
    ##########################################
    if ((lag<0)||(lag>ntimes))
        stop("Lag, ",sQuote("lag"),", must greater than 0 and less than ntimes",call.=FALSE)
    
    npars<- length(paramnames)
    phats<-rep(0,npars)
    names(phats)<-paramnames
  
    covhats <- array(
        0,
        dim=c(npars,npars)
    )
    pcovhats <- array(
        0,
        dim=c(0,0,0, 0)
    )
    if (lag>0 && !corr){
        asparticles <- vector(mode="list",length=(lag+1))
        xsparticles <- vector(mode="list",length=lag)
        psparticles <- vector(mode="list",length=lag)
    }
    
    if (lag>0 && corr){
        aparticles <- vector(mode="list",length=ntimes)
        wparticles <- vector(mode="list",length=ntimes)
    }
  
    ##########################################
  

  for (nt in seq_len(ntimes)) { ## main loop
    sigma1 <- sigma

    ## transform the parameters if necessary
    if (transform) tparams <- partrans(object,params,dir="fromEstimationScale",
                                       .getnativesymbolinfo=ptsi.for)
    ptsi.for <- FALSE


    ## advance the state variables according to the process model
    X <- tryCatch(
      rprocess(
        object,
        xstart=x,
        times=times[c(nt,nt+1)],
        params=params,
        offset=1,
        .getnativesymbolinfo=gnsi.rproc
      ),
      error = function (e) {
        stop(ep,"process simulation error: ",
             conditionMessage(e),call.=FALSE)
      }
    )
    gnsi.rproc <- FALSE

    if (pred.var) { ## check for nonfinite state variables and parameters
      problem.indices <- unique(which(!is.finite(X),arr.ind=TRUE)[,1L])
      if (length(problem.indices)>0) {  # state variables
        stop(
          ep,"non-finite state variable(s): ",
          paste(rownames(X)[problem.indices],collapse=', '),
          call.=FALSE
        )
      }
    }

    ## determine the weights
    weights <- tryCatch(
      dmeasure(
        object,
        y=object@data[,nt,drop=FALSE],
        x=X,
        times=times[nt+1],
        params=params,
        log=FALSE,
        .getnativesymbolinfo=gnsi.dmeas
      ),
      error = function (e) {
        stop(ep,"error in calculation of weights: ",
             conditionMessage(e),call.=FALSE)
      }
    )
    if (!all(is.finite(weights))) {
      first <- which(!is.finite(weights))[1L]
      datvals <- object@data[,nt]
      weight <- weights[first]
      states <- X[,first,1L]
      params <- if (one.par) params[,1L] else params[,first]
      msg <- nonfinite_dmeasure_error(time=times[nt+1],lik=weight,datvals,states,params)
      stop(ep,msg,call.=FALSE)
    }
    gnsi.dmeas <- FALSE

    ## compute prediction mean, prediction variance, filtering mean,
    ## effective sample size, log-likelihood
    ## also do resampling if filtering has not failed
    xx <- tryCatch(
      .Call(
        pfilter2_computations,
        x=X,
        params=params,
        Np=Np[nt+1],
        rw_sd=numeric(0),
        predmean=pred.mean,
        predvar=pred.var,
        filtmean=filter.mean,
        trackancestry=filter.traj,
        onepar=one.par,
        weights=weights,
        tol=tol
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE) # nocov
      }
    )
    all.fail <- xx$fail
    loglik[nt] <- xx$loglik
    eff.sample.size[nt] <- xx$ess

    x <- xx$states
        #random walk change to white noise for lag>0
    if(lag>0 && wn){
      params <- params
    } else{
      params <- xx$params
    }

    if (pred.mean)
      pred.m[,nt] <- xx$pm
    if (pred.var)
      pred.v[,nt] <- xx$pv
    if (filter.mean)
      filt.m[,nt] <- xx$fm
    if (filter.traj)
      pedigree[[nt]] <- xx$ancestry

    if (all.fail) { ## all particles are lost
      nfail <- nfail+1
      if (verbose)
        message("filtering failure at time t = ",times[nt+1])
      if (nfail>max.fail)
        stop(ep,"too many filtering failures",call.=FALSE)
    }

    if (save.states | filter.traj) {
      xparticles[[nt]] <- x
      dimnames(xparticles[[nt]]) <- setNames(dimnames(xparticles[[nt]]),c("variable","rep"))
    }

    if (save.params) {
      pparticles[[nt]] <- params
      dimnames(pparticles[[nt]]) <- setNames(dimnames(pparticles[[nt]]),c("variable","rep"))
    }

    if (verbose && (nt%%5==0))
      cat("psmooth timestep",nt,"of",ntimes,"finished\n")
      
     ##########################################
    if (lag>0 && !corr){
      if(nt<(lag+1)){
        xsparticles[[nt]] <- x
        psparticles[[nt]] <- params
        asparticles[[nt+1]] <- xx$pa
      }
      if(nt>lag && nt<=ntimes){
	    index<-unlist(ancestor(plist=asparticles,t=lag+1, lag=lag, 1:(Np[1])))
        psparticles[[1]][!is.finite(psparticles[[1]])] <- 0 
                
        C<-cov.wt(t(psparticles[[1]][,index]),wt=xx$weight)
        phats<-phats+C$center
        covhats<-covhats+C$cov/(Np[1]-1)
	    if (lag>1){                
		  for (i in 1:(lag-1)){
            psparticles[[i]]<-psparticles[[i+1]]
		    asparticles[[i+1]]<-asparticles[[i+2]]  
          }
	    }
        psparticles[[lag]]<-params
	    asparticles[[lag+1]]<-xx$pa
      }
      if(nt==ntimes){
	    index<-unlist(ancestor(plist=asparticles,t=lag+1, lag=lag, 1:(Np[1])))
        psparticles[[1]][!is.finite(psparticles[[1]])] <- 0 
        C<-cov.wt(t(psparticles[[1]][,index]),wt=xx$weight)
        phats<-phats+C$center
        covhats<-covhats+C$cov/(Np[1]-1)
	    if (lag>1){                
		  for (i in 1:(lag-1)){
            psparticles[[i]]<-psparticles[[i+1]]
		    asparticles[[i+1]]<-asparticles[[i+2]]  
          }
	    }
      }
    }
    if(lag>0 && corr){
      pparticles[[nt]] <- xx$params
      wparticles[[nt]] <-xx$weight
      if(nt==1)
        aparticles[[1]] <- 0   
      if (nt<ntimes)
        aparticles[[nt+1]] <- xx$pa  
    }
    if(lag>0 && corr){
      pparticles[[nt]] <- xx$params
      wparticles[[nt]] <-xx$weight
      if(nt==1)
        aparticles[[1]] <- 0   
      if (nt<ntimes)
        aparticles[[nt+1]] <- xx$pa  
    }
 

  }## end of main loop
    ###################################################################
    # fixed lag smoothing
    ###################################################################
  if(lag>0 && corr){
    pcovhats <- array(
          0,
          dim=c(npars,npars,lag, ntimes)
    )
      
    Np<-Np[1]
    #pcovhat<-matrix(0,npars,lag) #covariance
    results<-lapply(1:ntimes,smoothing,aparticles=aparticles,xparticles=xparticles,pparticles=pparticles,wparticles=wparticles,
                                ntimes=ntimes,lag=lag, nvars=nvars, npars=npars, sigma,Np=Np)
    for(i in 1:ntimes){
      pcovhats[,,,i]=results[[i]]
    }
    # Clean up
    gc()
  }
 
  if (filter.traj) { ## select a single trajectory
    if (max(weights)>0) {
      b <- sample.int(n=length(weights),size=1L,
                      prob=weights,replace=TRUE)
    } else {
      b <- sample.int(n=length(weights),size=1L,
                      replace=TRUE)
    }
    filt.t[,1L,ntimes+1] <- xparticles[[ntimes]][,b]
    for (nt in seq.int(from=ntimes-1,to=1L,by=-1L)) {
      b <- pedigree[[nt+1]][b]
      filt.t[,1L,nt+1] <- xparticles[[nt]][,b]
    }
    if (times[2L] > times[1L]) {
      b <- pedigree[[1L]][b]
      filt.t[,1L,1L] <- init.x[,b]
    } else {
      filt.t <- filt.t[,,-1L,drop=FALSE]
    }
  }

  if (!save.states) xparticles <- list()

  if (nfail>0)
    warning(
      ep,nfail,
      ngettext(
        nfail,
        msg1=" filtering failure occurred.",
        msg2=" filtering failures occurred."
      ),
      call.=FALSE
    )

  pompUnload(object,verbose=verbose)

  new(
    "psmooth.pomp",
    object,
    pred.mean=pred.m,
    pred.var=pred.v,
    filter.mean=filt.m,
    filter.traj=filt.t,
    paramMatrix=array(data=numeric(0),dim=c(0,0)),
    eff.sample.size=eff.sample.size,
    cond.loglik=loglik,
    saved.states=xparticles,
    saved.params=pparticles,
    Np=as.integer(Np),
    tol=tol,
    nfail=as.integer(nfail),
    loglik=sum(loglik),
    phats=phats,
    covhats=covhats,
    pcovhats=pcovhats,
    lag=lag

  )
}

setMethod(
  "psmooth",
  signature=signature(object="pomp"),
  function (object, params, Np,
            tol = 1e-17,
            max.fail = Inf,
            pred.mean = FALSE,
            pred.var = FALSE,
            filter.mean = FALSE,
            filter.traj = FALSE,
            save.states = FALSE,
            save.params = FALSE,
            lag =0,
            verbose = getOption("verbose"),
            ...) {
    if (missing(params)) params <- coef(object)
    psmooth.internal(
      object=object,
      params=params,
      Np=Np,
      tol=tol,
      max.fail=max.fail,
      pred.mean=pred.mean,
      pred.var=pred.var,
      filter.mean=filter.mean,
      filter.traj=filter.traj,
      save.states=save.states,
      save.params=save.params,
      lag=lag,
      verbose=verbose,
      ...
    )
  }
)

setMethod(
  "psmooth",
  signature=signature(object="psmooth.pomp"),
  function (object, params, Np, tol, ...) {
    if (missing(params)) params <- coef(object)
    if (missing(Np)) Np <- object@Np
    if (missing(tol)) tol <- object@tol
    f <- selectMethod("psmooth","pomp")
    f(
      object=object,
      params=params,
      Np=Np,
      tol=tol,
      ...
    )
  }
)

nonfinite_dmeasure_error <- function (time, lik, datvals, states, params) {
  showvals <- c(time=time,lik=lik,datvals,states,params)
  m1 <- formatC(names(showvals),preserve.width="common")
  m2 <- formatC(showvals,digits=6,width=12,format="g",
                preserve.width="common")
  paste0(
    sQuote("dmeasure")," returns non-finite value.\n",
    "likelihood, data, states, and parameters are:\n",
    paste0(m1,": ",m2,collapse="\n")
  )
}
