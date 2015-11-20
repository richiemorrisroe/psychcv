##' {Return an averaged factor analysis solution over a number of cross-validated splits}
##' {Takes two fa objects, and aggregates them into one set of coefficients. Can be useful either for combining CV selected splits of the correct size, and for bootstrapped FA results}
##' @title factor_average
##' @param sols a list of factor solutions to average over
##' @param mynames who the hell knows?
##' @param FUN function to aggregate by
##' @param ... other arguments passed to fun
##' @return a matrix containing the averaged fa solutions
##' @author Richie Morrisroe
factor_average <- function (sols=list(), mynames=NULL, FUN=mean, ...) {

    sols_coeff_list <- list()

    for(i in 1:length(sols)) {
        coeff <- as.data.frame(factor_coeff(sols[[i]]))
        if(!is.null(mynames)) {
            coeff_ord <- coeff[,mynames]
        }
        else {
            coeff_ord <- coeff
        }
        sols_coeff_list[[i]] <- coeff_ord
    }
    sols_coeff_list
    sols_list <- lapply(sols_coeff_list, as.matrix)

    resmat <- apply(simplify2array(sols_list), c(1,2), FUN)
    return(resmat)

    
}
##' {Fit a series of factor solutions from 1:k}
##' {Fit a series of factor solutions}
##' @title fit_factor_series
##' @inheritParams psych::fa
##' @param data a matrix of numeric data
##' @param factors a sequence of numbers
##' @param meth method of rotation (see \link[psych]{fa}), defaults to minres 
##' @param rotation a method of rotation see \link[psych]{fa}. defaults to oblimin
##' @param scores which kind of factor scores to compute
##' @param ... further arguments passed through to psych::fa
##' @return a factor_series object
##' @author Richard Morrisroe
fit_factor_series <- function(data, factors, meth, rotation, scores, ...) {
    fno <- factors
    if(length(fno) == 1) {
        stop("factors should be a range of numbers")
    }
      fnolist <-  list()
  for (k in seq_along(along.with=fno)) {
      z <- psych::fa(na.omit(data), nfactors=fno[k],
                     rotate="oblimin", fm="minres", ...)
      fnolist[[k]] <- z
      ## browser()
  }
    class(fnolist) <- c("factor_series", "fa")
    fnolist
}
get_component <- function(fs, component) {
    fac <- unlist(sapply(fs, `[`, "factors"))
    compdata <- sapply(fs, `[`, component)
    if(length(compdata[[1]]) > 1) {
        warning("aggregating results using the mean function")
        compdata <- unlist(lapply(compdata, mean, na.rm=TRUE))
    }
    data.frame(component=compdata, factors=fac)
}
extractor <- function (fs, parameter) {
    function(fs) {
    get_component(fs, component=parameter)
}
    #return a function that extracts the given parameter
}
chi <- extractor(parameter="chi")
rms <- extractor(parameter="rms")
communalities <- extractor(parameter="communality")
communality <- function(fs) {
    #deserves its own function
     sapply(fs, `[`, "communality")

}
##' {Extract loadings from a list of multiple factor solutions}
##' {As above}
##' @title get_loadings
##' @param fs a factor series object
##' @param loadings a scalar specifying the minimum threshold for items to be returned
##' @return a list containing the loadings of all fs solutiions
##' @author Richard Morrisroe
get_loadings <- function (fs, loadings=0.3) {
  ind <- lapply(fs, extract_loadings, loadings)
  ind
}
##' {This seems very similar to FactorAverage, rationalise these functions ASAP}
##' (no empty lines)
##' {}
##' @title combine_loadings
##' @param mfa a multi-factor object
##' @return a list of mean loadings averaged over all potential factor solutions
##' @author Richard Morrisroe
combine_loadings <-  function (mfa) {
  loadlist <- list()
  for (i in seq_along(along.with=mfa)) {
    loadlist[[i]] <- mfa[[i]]$loadings

  }
  loadlist
  loadings <- list()
  for (j in seq_along(along.with=loadlist)) {


    loadings[[j]] <- as.matrix(unclass(loadlist[[j]]))
}
  loadings

  lapply(loadings, function (x) Reduce("+", x))
}
fortify.mfa <- function(model, data, parameter, ...) {
    stopifnot(class(model)=="mfa")
    if(parameter != communality) {
        stop("communality only part implemented")
    }
    comms <- communality(model)
    comms.dc <- do.call("rbind", comms)
    comms.df <- sapply(comms.dc, unlist)
    comms.m <- reshape2::melt(comms.df, id.vars=rownames(mfa.df))
    comms.m[,"Item"] <- with(mfa.m, gsub("communality", "",  x=Var1))
    return(comms.m)
}
