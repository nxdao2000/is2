## define the pmif class
setClass(
    'pmif',
    contains='pfilterd.pomp',
    slots=c(
        pars = 'character',
        Nmif = 'integer',
        accepts = 'integer',
        proposal = 'function',
        conv.rec = 'array',
        log.prior = 'numeric'
    ),
    prototype=prototype(
        pars = character(0),
        Nmif = 0L,
        accepts = 0L,
        proposal = function (...)
            stop("in ",sQuote("pmif"),": proposal not specified",call.=FALSE),
        conv.rec=array(dim=c(0,0)),
        log.prior=numeric(0)
    )
)

pmif.internal <- function (object, Nmif,
                            start, proposal,
                            Np, tol, max.fail,
                            verbose,
                            .ndone = 0L,
                            .accepts = 0L,
                            .prev.pfp = NULL, .prev.log.prior = NULL,
                            .getnativesymbolinfo = TRUE) {

    object <- as(object,"pomp")
    gnsi <- as.logical(.getnativesymbolinfo)
    verbose <- as.logical(verbose)
    .ndone <- as.integer(.ndone)
    .accepts <- as.integer(.accepts)
    
    ep <- paste0("in ",sQuote("pmif"),": ")

    pompLoad(object,verbose=verbose)

    if (missing(start))
        stop(ep,sQuote("start")," must be specified",call.=FALSE)
    if (length(start)==0)
        stop(ep,sQuote("start")," must be specified if ",
             sQuote("coef(object)")," is NULL",call.=FALSE)
    if (is.null(names(start)))
        stop(ep,sQuote("start")," must be a named vector",call.=FALSE)

    if (!is.function(proposal))
        stop(ep,sQuote("proposal")," must be a function",call.=FALSE)

    ## test proposal distribution
    theta <- tryCatch(
        proposal(start,.n=0),
        error = function (e) {
            stop(ep,"error in proposal function: ",conditionMessage(e),call.=FALSE)
        }
    )
    if (is.null(names(theta)) || !is.numeric(theta) || any(names(theta)==""))
        stop(ep,sQuote("proposal")," must return a named numeric vector",call.=FALSE)

    ntimes <- length(time(object))
    if (missing(Np))
        stop(ep,sQuote("Np")," must be specified",call.=FALSE)
    if (is.function(Np)) {
        Np <- tryCatch(
            vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
            error = function (e) {
                stop(ep,"if ",sQuote("Np")," is a function, it must return a single positive integer",call.=FALSE)
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

    if (missing(Nmif))
        stop(ep,sQuote("Nmif")," must be specified",call.=FALSE)
    Nmif <- as.integer(Nmif)
    if (Nmif<0)
        stop(ep,sQuote("Nmif")," must be a positive integer",call.=FALSE)

    if (verbose) {
        cat("performing",Nmif,"pmif iteration(s) using",Np[1L],"particles\n")
    }

    conv.rec <- matrix(
        data=NA,
        nrow=Nmif+1,
        ncol=length(theta)+3,
        dimnames=list(
            iteration=seq(from=0,to=Nmif,by=1),
            variable=c('loglik','log.prior','nfail',names(theta))
        )
    )

    if (.ndone==0L) { ## compute prior and likelihood on initial parameter vector
        pfp <- tryCatch(
            pfilter1.internal(
                object=object,
                params=theta,
                Np=Np,
                tol=tol,
                max.fail=max.fail,
                pred.mean=FALSE,
                pred.var=FALSE,
                filter.mean=TRUE,
                filter.traj=TRUE,
                save.states=FALSE,
                save.params=FALSE,
                verbose=FALSE,
                .getnativesymbolinfo=gnsi
            ),
            error = function (e) {
                stop(ep,conditionMessage(e),call.=FALSE) # nocov
            }
        )
        log.prior <- tryCatch(
            dprior(object,params=theta,log=TRUE,.getnativesymbolinfo=gnsi),
            error = function (e) {
                stop(ep,sQuote("dprior")," error: ",conditionMessage(e),call.=FALSE)
            }
        )
        gnsi <- FALSE
    } else { ## has been computed previously
        pfp <- .prev.pfp
        log.prior <- .prev.log.prior
        pfp@filter.traj <- pfp@filter.traj[,.ndone,,drop=FALSE]
    }
    conv.rec[1,names(theta)] <- theta
    conv.rec[1,c(1,2,3)] <- c(pfp@loglik,log.prior,pfp@nfail)
    
    filt.t <- array(
        data=0,
        dim=replace(dim(pfp@filter.traj),2L,Nmif),
        dimnames=replace(dimnames(pfp@filter.traj),2L,
                         list(as.character(seq_len(Nmif))))
    )

    for (n in seq_len(Nmif)) { # main loop

        theta.prop <- tryCatch(
            proposal(theta,.n=n+.ndone,.accepts=.accepts,verbose=verbose),
            error = function (e) {
                stop(ep,"error in proposal function: ",conditionMessage(e),call.=FALSE)
            }
        )
        pfp.prop <- tryCatch(
                pfilter1.internal(
                    object=pfp,
                    params=theta.prop,
                    Np=Np,
                    tol=tol,
                    max.fail=max.fail,
                    pred.mean=FALSE,
                    pred.var=FALSE,
                    filter.mean=TRUE,
                    filter.traj=TRUE,
                    save.states=FALSE,
                    save.params=TRUE,
                    verbose=FALSE,
                    .getnativesymbolinfo=gnsi
                ),
                error = function (e) {
                    stop(ep,conditionMessage(e),call.=FALSE) # nocov
                }
            )
            
            
        theta.prop <- tryCatch(
            proposal(pfp.prop@params,.n=n+.ndone,.accepts=.accepts,verbose=verbose),
            error = function (e) {
                stop(ep,"error in proposal function: ",conditionMessage(e),call.=FALSE)
            }
        )
        
    
            #theta.prop <- pfp.prop@params
#            tryCatch(
 #           proposal(pfp.prop@params,.n=n+.ndone,.accepts=.accepts,verbose=verbose),
 #           error = function (e) {
 #               stop(ep,"error in proposal function: ",conditionMessage(e),call.=FALSE)
 #           }
 #       )
        
            
        ## compute log prior
        log.prior.prop <- tryCatch(
            dprior(object,params=theta.prop,log=TRUE,.getnativesymbolinfo=gnsi),
            error = function (e) {
                stop(ep,sQuote("dprior")," error: ",conditionMessage(e),call.=FALSE)
            }
        )

        if (is.finite(log.prior.prop)) {
            
            ## run the particle filter on the proposed new parameter values
            pfp.prop <- tryCatch(
                pfilter1.internal(
                    object=pfp,
                    params=theta.prop,
                    Np=Np,
                    tol=tol,
                    max.fail=max.fail,
                    pred.mean=FALSE,
                    pred.var=FALSE,
                    filter.mean=TRUE,
                    filter.traj=TRUE,
                    save.states=FALSE,
                    save.params=FALSE,
                    verbose=FALSE,
                    .getnativesymbolinfo=gnsi
                ),
                error = function (e) {
                    stop(ep,conditionMessage(e),call.=FALSE) # nocov
                }
            )
            gnsi <- FALSE

            ## pmif update rule (OK because proposal is symmetric)
            alpha <- exp(pfp.prop@loglik+log.prior.prop-pfp@loglik-log.prior)
            if (runif(1) < alpha) {
                pfp <- pfp.prop
                theta <- theta.prop
                log.prior <- log.prior.prop
                .accepts <- .accepts+1L
            }
        }

        ## add filtered trajectory to the store
        filt.t[,n,] <- pfp@filter.traj[,1L,]

        ## store a record of this iteration
        conv.rec[n+1,names(theta)] <- theta
        conv.rec[n+1,c(1,2,3)] <- c(pfp@loglik,log.prior,pfp@nfail)

        if (verbose) cat("pmif iteration",n+.ndone,"of",Nmif+.ndone,
                         "completed\nacceptance ratio:",
                         round(.accepts/(n+.ndone),3),"\n")

    }

    pars <- apply(conv.rec,2,function(x)diff(range(x))>0)
    pars <- setdiff(names(pars[pars]),c("loglik","log.prior","nfail"))

    pompUnload(object,verbose=verbose)

    new(
        "pmif",
        pfp,
        params=theta,
        pars=pars,
        Nmif=Nmif,
        accepts=.accepts,
        proposal=proposal,
        Np=Np,
        tol=tol,
        conv.rec=conv.rec,
        log.prior=log.prior, 
        filter.traj=filt.t
    )
}

setMethod(
    "pmif",
    signature=signature(object="pomp"),
    function (object, Nmif = 1,
              start, proposal, Np,
              tol = 1e-17, max.fail = Inf,
              verbose = getOption("verbose"),
              ...) {
        
        ep <- paste0("in ",sQuote("pmif"),": ")
        
        if (missing(start)) start <- coef(object)
        if (missing(Np))
            stop(ep,sQuote("Np")," must be specified",call.=FALSE)
        
        if (missing(proposal)) proposal <- NULL

        if (is.null(proposal))
            stop(ep,sQuote("proposal")," must be specified",call.=FALSE)

        pmif.internal(
            object=object,
            Nmif=Nmif,
            start=start,
            proposal=proposal,
            Np=Np,
            tol=tol,
            max.fail=max.fail,
            verbose=verbose,
            ...
        )
    }
)

setMethod(
    "pmif",
    signature=signature(object="pfilterd.pomp"),
    function (object, Nmif = 1, Np, tol, ...) {

        if (missing(Np)) Np <- object@Np
        if (missing(tol)) tol <- object@tol
        
        f <- selectMethod("pmif","pomp")

        f(object,Nmif=Nmif,Np=Np,tol=tol,...)
    }
)

setMethod(
    "pmif",
    signature=signature(object="pmif"),
    function (object, Nmif,
              start, proposal,
              Np, tol, max.fail = Inf,
              verbose = getOption("verbose"),
              ...) {

        if (missing(Nmif)) Nmif <- object@Nmif
        if (missing(start)) start <- coef(object)
        if (missing(proposal)) proposal <- object@proposal
        if (missing(Np)) Np <- object@Np
        if (missing(tol)) tol <- object@tol

        f <- selectMethod("pmif","pomp")
        
        f(object,Nmif=Nmif,start=start,proposal=proposal,
          Np=Np,tol=tol,max.fail=max.fail,verbose=verbose,...)
    }
)

setMethod(
    'continue',
    signature=signature(object='pmif'),
    function (object, Nmif = 1, ...) {

        ndone <- object@Nmif
        accepts <- object@accepts

        obj <- pmif(
            object=object,
            Nmif=Nmif,
            ...,
            .ndone=ndone,
            .accepts=accepts,
            .prev.pfp=as(object,"pfilterd.pomp"),
            .prev.log.prior=object@log.prior
        )
        
        obj@conv.rec <- rbind(
            object@conv.rec[,colnames(obj@conv.rec)],
            obj@conv.rec[-1,]
        )
        names(dimnames(obj@conv.rec)) <- c("iteration","variable")
        ft <- array(dim=replace(dim(obj@filter.traj),2L,ndone+Nmif),
                    dimnames=replace(dimnames(obj@filter.traj),2L,
                                     list(seq_len(ndone+Nmif))))
        ft[,seq_len(ndone),] <- object@filter.traj
        ft[,ndone+seq_len(Nmif),] <- obj@filter.traj
        obj@filter.traj <- ft
        obj@Nmif <- as.integer(ndone+Nmif)
        obj@accepts <- as.integer(accepts+obj@accepts)

        obj
    }
)
