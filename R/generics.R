

## iterated smoothing
setGeneric('is2',function(object,...)standardGeneric("is2"))

## accelerate iterated smoothing
setGeneric('is3',function(object,...)standardGeneric("is3"))

## SMC (particle filter) and use for fixed lag smoothing
setGeneric("pfilter1",function(object,...)standardGeneric("pfilter1"))

setMethod(
    "pfilter1",
    signature=signature(object="missing"),
    definition=function (...) {
        stop("in ",sQuote("pfilter1"),": ",sQuote("object")," is a required argument",call.=FALSE)
    }
)

setMethod(
    "pfilter1",
    signature=signature(object="ANY"),
    definition=function (object, ...) {
        stop(sQuote("pfilter1")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
    }
)


## SMC (particle filter) and use for fixed lag smoothing
setGeneric("pfilter2",function(object,...)standardGeneric("pfilter2"))

setMethod(
    "pfilter2",
    signature=signature(object="missing"),
    definition=function (...) {
        stop("in ",sQuote("pfilter2"),": ",sQuote("object")," is a required argument",call.=FALSE)
    }
)

setMethod(
    "pfilter2",
    signature=signature(object="ANY"),
    definition=function (object, ...) {
        stop(sQuote("pfilter2")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
    }
)


## SMC (particle filter) and use for fixed lag smoothing
setGeneric("pfilter3",function(object,...)standardGeneric("pfilter3"))

setMethod(
    "pfilter3",
    signature=signature(object="missing"),
    definition=function (...) {
        stop("in ",sQuote("pfilter3"),": ",sQuote("object")," is a required argument",call.=FALSE)
    }
)

setMethod(
    "pfilter3",
    signature=signature(object="ANY"),
    definition=function (object, ...) {
        stop(sQuote("pfilter3")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
    }
)



## SMC (particle filter) and use for fixed lag smoothing
setGeneric("pfilter4",function(object,...)standardGeneric("pfilter4"))

setMethod(
    "pfilter4",
    signature=signature(object="missing"),
    definition=function (...) {
        stop("in ",sQuote("pfilter4"),": ",sQuote("object")," is a required argument",call.=FALSE)
    }
)

setMethod(
    "pfilter4",
    signature=signature(object="ANY"),
    definition=function (object, ...) {
        stop(sQuote("pfilter4")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
    }
)

setMethod(
    "is2",
    signature=signature(object="missing"),
    definition=function (...) {
        stop("in ",sQuote("is2"),": ",sQuote("object")," is a required argument",call.=FALSE)
    }
)

setMethod(
    "is2",
    signature=signature(object="ANY"),
    definition=function (object, ...) {
        stop(sQuote("is2")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
    }
)

setMethod(
    "is3",
    signature=signature(object="missing"),
    definition=function (...) {
        stop("in ",sQuote("is3"),": ",sQuote("object")," is a required argument",call.=FALSE)
    }
)

setMethod(
    "is3",
    signature=signature(object="ANY"),
    definition=function (object, ...) {
        stop(sQuote("is3")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
    }
)

## iterated smoothing
setGeneric('mif1',function(object,...)standardGeneric("mif1"))


setMethod(
    "mif1",
    signature=signature(object="missing"),
    definition=function (...) {
        stop("in ",sQuote("mif1"),": ",sQuote("object")," is a required argument",call.=FALSE)
    }
)

setMethod(
    "mif1",
    signature=signature(object="ANY"),
    definition=function (object, ...) {
        stop(sQuote("mif1")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
    }
)

## accelerate iterated smoothing
setGeneric('aif',function(object,...)standardGeneric("aif"))


setMethod(
    "aif",
    signature=signature(object="missing"),
    definition=function (...) {
        stop("in ",sQuote("aif"),": ",sQuote("object")," is a required argument",call.=FALSE)
    }
)

setMethod(
    "aif",
    signature=signature(object="ANY"),
    definition=function (object, ...) {
        stop(sQuote("aif")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
    }
)

## accelerate iterated smoothing
setGeneric('avif',function(object,...)standardGeneric("avif"))


setMethod(
    "avif",
    signature=signature(object="missing"),
    definition=function (...) {
        stop("in ",sQuote("avif"),": ",sQuote("object")," is a required argument",call.=FALSE)
    }
)

setMethod(
    "avif",
    signature=signature(object="ANY"),
    definition=function (object, ...) {
        stop(sQuote("avif")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
    }
)

## iterated filtering
setGeneric('mif3',function(object,...)standardGeneric("mif3"))


setMethod(
    "mif3",
    signature=signature(object="missing"),
    definition=function (...) {
        stop("in ",sQuote("mif3"),": ",sQuote("object")," is a required argument",call.=FALSE)
    }
)

setMethod(
    "mif3",
    signature=signature(object="ANY"),
    definition=function (object, ...) {
        stop(sQuote("mif3")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
    }
)

setGeneric('mifMomentum',function(object,...)standardGeneric("mifMomentum"))

setMethod(
  "mifMomentum",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("mifMomentum"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "mifMomentum",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("mifMomentum")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

## particle Markov chain Monte Carlo (Pmif)
setGeneric('pmif',function(object,...)standardGeneric("pmif"))

setMethod(
    "pmif",
    signature=signature(object="missing"),
    definition=function (...) {
        stop("in ",sQuote("pmif"),": ",sQuote("object")," is a required argument",call.=FALSE)
    }
)

setMethod(
    "pmif",
    signature=signature(object="ANY"),
    definition=function (object, ...) {
        stop(sQuote("pmif")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
    }
)

## basic SMC (particle filter)
setGeneric("psmooth",function(object,...)standardGeneric("psmooth"))

setMethod(
    "psmooth",
    signature=signature(object="missing"),
    definition=function (...) {
        stop("in ",sQuote("psmooth"),": ",sQuote("object")," is a required argument",call.=FALSE)
    }
)

setMethod(
    "psmooth",
    signature=signature(object="ANY"),
    definition=function (object, ...) {
        stop(sQuote("psmooth")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
    }
)


