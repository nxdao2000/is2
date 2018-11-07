## this file contains short definitions of methods for the 'mif1' class

## extract the convergence record
conv.rec.internal <- function (object, pars, transform = FALSE, ...) {
    if (transform) {
        pars.proper <- names(coef(object))
        pars.improper <- setdiff(colnames(object@conv.rec),pars.proper)
        retval <- cbind(
            t(
                partrans(
                    object,
                    params=t(object@conv.rec[,pars.proper]),
                    dir="fromEstimationScale"
                )
            ),
            object@conv.rec[,pars.improper]
        )
    } else {
        retval <- object@conv.rec
    }
    if (missing(pars))
        retval
    else {
        bad.pars <- setdiff(pars,colnames(retval))
        if (length(bad.pars)>0)
            stop(
                "in ",sQuote("conv.rec"),": name(s) ",
                paste(sQuote(bad.pars),collapse=","),
                " correspond to no parameter(s) in ",
                if (transform) sQuote("conv.rec(object,transform=TRUE)")
                else sQuote("conv.rec(object,transform=FALSE)"),
                call.=FALSE
            )
        retval[,pars]
    }
}

setMethod('conv.rec','mif1',
          function (object, pars, transform = FALSE, ...) {
              conv.rec.internal(object=object,pars=pars,transform=transform,...)
          }
          )

## plot mif1 object
setMethod(
    "plot",
    "mif1",
    function (x, y, ...) {
        if (!missing(y))
            warning("in ",sQuote("plot-mif1"),": ",
                    sQuote("y")," is ignored",call.=FALSE)
        mif1.diagnostics(list(x))
    }
)

## mif1List class
setClass(
    'mif1List',
    contains='list',
    validity=function (object) {
        if (!all(sapply(object,is,'mif1'))) {
            retval <- paste0(
                "error in ",sQuote("c"),
                ": dissimilar objects cannot be combined"
            )
            return(retval)
        }
        d <- sapply(object,function(x)dim(x@conv.rec))
        if (!all(apply(d,1,diff)==0)) {
            retval <- paste0(
                "error in ",sQuote("c"),
                ": to be combined, ",sQuote("mif1"),
                " objects must equal numbers of iterations"
            )
            return(retval)
        }
        TRUE
    }
)

setMethod(
    'c',
    signature=signature(x='mif1'),
    definition=function (x, ...) {
        y <- list(...)
        if (length(y)==0) {
            new("mif1List",list(x))
        } else {
            p <- sapply(y,is,'mif1')
            pl <- sapply(y,is,'mif1List')
            if (!all(p||pl))
                stop("in ",sQuote("c"),": cannot mix ",sQuote("mif1"),
                     " and non-",sQuote("mif1")," objects",call.=FALSE)
            y[p] <- lapply(y[p],list)
            y[pl] <- lapply(y[pl],as,"list")
            new("mif1List",c(list(x),y,recursive=TRUE))
        }
    }
)

setMethod(
    'c',
    signature=signature(x='mif1List'),
    definition=function (x, ...) {
        y <- list(...)
        if (length(y)==0) {
            x
        } else {
            p <- sapply(y,is,'mif1')
            pl <- sapply(y,is,'mif1List')
            if (!all(p||pl))
                stop("in ",sQuote("c"),": cannot mix ",sQuote("mif1"),
                     " and non-",sQuote("mif1")," objects",call.=FALSE)
            y[p] <- lapply(y[p],list)
            y[pl] <- lapply(y[pl],as,"list")
            new("mif1List",c(as(x,"list"),y,recursive=TRUE))
        }
    }
)

setMethod(
    "[",
    signature=signature(x="mif1List"),
    definition=function(x, i, ...) {
        new('mif1List',as(x,"list")[i])
    }
)

setMethod(
    'conv.rec',
    signature=signature(object='mif1List'),
    definition=function (object, ...) {
        lapply(object,conv.rec,...)
    }
)

setMethod(
    'coef',
    signature=signature(object='mif1List'),
    definition=function (object, ...) {
        do.call(rbind,lapply(object,coef,...))
    }
)

setMethod(
    "plot",
    signature=signature(x='mif1List'),
    definition=function (x, y, ...) {
        if (!missing(y))
            warning("in ",sQuote("plot-mif1"),": ",
                    sQuote("y")," is ignored",call.=FALSE)
        mif1.diagnostics(x)
    }
)

## predvarplot.mif1 <- function (object, pars, type = 'l', mean = FALSE, ...) {
##   if (!is(object,'mif1'))
##     stop("predvarplot error: ",sQuote("object")," must be of class ",sQuote("mif1"),call.=FALSE)
##   if (missing(pars))
##     pars <- object@pars
##   npv <- pred.var(object,pars)/(object@random.walk.sd[pars]^2)
##   if (!is.null(dim(npv))) npv <- t(npv)
##   if (mean && !is.null(dim(npv)))
##     npv <- rowMeans(npv)
##   if (!is.null(dim(npv))) {
##     matplot(time(object),npv,type=type,ylab='prediction variance',xlab='time',...)
##     legend(x='topright',legend=pars,col=1:length(pars),lty=1:length(pars),bty='n')
##   } else {
##     plot(time(object),npv,type=type,ylab='prediction variance',xlab='time',...)
##   }
## }

mif1.diagnostics <- function (z) {
    ## assumes that z is a list of mif1s with identical structure
    mar.multi <- c(0,5.1,0,2.1)
    oma.multi <- c(6,0,5,0)
    xx <- z[[1]]
    ivpnames <- xx@ivps
    estnames <- c(xx@pars,ivpnames)
    parnames <- names(coef(xx,transform=xx@transform))
    unestnames <- parnames[-match(estnames,parnames)]

    ## plot filter means
    filt.diag <- rbind("eff. sample size"=xx@eff.sample.size,filter.mean(xx))
    filtnames <- rownames(filt.diag)
    plotnames <- if(length(unestnames)>0) filtnames[-match(unestnames,filtnames)] else filtnames
    lognames <- filtnames[1] # eff. sample size
    nplots <- length(plotnames)
    n.per.page <- min(nplots,10)
    if(n.per.page<=4) nc <- 1 else nc <- 2
    nr <- ceiling(n.per.page/nc)
    oldpar <- par(mar=mar.multi,oma=oma.multi,mfcol=c(nr,nc),ask=dev.interactive(orNone=TRUE))
    on.exit(par(oldpar)) 
    low <- 1
    hi <- 0
    time <- time(xx)
    while (hi<nplots) {
        hi <- min(low+n.per.page-1,nplots)
        for (i in seq(from=low,to=hi,by=1)) {
            n <- i-low+1
            logplot <- if (plotnames[i]%in%lognames) "y" else ""
            dat <- sapply(
                z,
                function(po, label) {
                    if (label=="eff. sample size") 
                        po@eff.sample.size
                    else
                        filter.mean(po,label)
                },
                label=plotnames[i]
            )
            matplot(
                y=dat, 
                x=time,
                axes = FALSE,
                xlab = "",
                log=logplot,
                ylab = "",
                type = "l"
            )
            box()
            y.side <- 2
            axis(y.side, xpd = NA)
            mtext(plotnames[i], y.side, line = 3)
            do.xax <- (n%%nr==0||n==n.per.page)
            if (do.xax) axis(1,xpd=NA)
            if (do.xax) mtext("time",side=1,line=3)
        }  
        low <- hi+1
        mtext("Filter diagnostics (last iteration)",3,line=2,outer=TRUE)
    } 

    ## plot mif1 convergence diagnostics
    other.diagnostics <- c("loglik", "nfail")
    plotnames <- c(other.diagnostics,estnames)
    nplots <- length(plotnames)
    n.per.page <- min(nplots,10)
    nc <- if (n.per.page<=4) 1 else 2
    nr <- ceiling(n.per.page/nc)
    par(mar=mar.multi,oma=oma.multi,mfcol=c(nr,nc))
    ## on.exit(par(oldpar)) 
    low <- 1
    hi <- 0
    iteration <- seq(0,xx@Nmif)
    while (hi<nplots) {
        hi <- min(low+n.per.page-1,nplots)
        for (i in seq(from=low,to=hi,by=1)) {
            n <- i-low+1
            dat <- sapply(z,function(po,label) conv.rec(po,label),label=plotnames[i])
            matplot(
                y=dat, 
                x=iteration,
                axes = FALSE,
                xlab = "",
                ylab = "",
                type = "l"
            )
            box()
            y.side <- 2
            axis(y.side,xpd=NA)
            mtext(plotnames[i],y.side,line=3)
            do.xax <- (n%%nr==0||n==n.per.page)
            if (do.xax) axis(1,xpd=NA)
            if (do.xax) mtext("MIF iteration",side=1,line=3)
        }  
        low <- hi+1
        mtext("MIF convergence diagnostics",3,line=2,outer=TRUE)
    }
    invisible(NULL)
}
