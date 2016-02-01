##' {Takes an fa object, extracts the coefficients, and returns an xtable object of these coefficients (including communalities)}
##'
##' {Returns a table suitable for research. (Probably) doesn't meet APA, or indeed any, standards. Thin wrapper around factor_coeff}
##' @title factor_xtab
##' @param x an object of class "psych" "fa"
##' @param names names to be given to the output. Should be a vector of strings
##' @param ... further arguments passed to the FactorCoeff method
##' @return xtable object of coefficients
##' @author Richard Morrisroe
factor_xtab <-  function (x, names=NULL, ...) {
    xmatdf <- factor_coeff(x, names=names)
    fact.xtab <- xtable::xtable(xmatdf, ...)
    return(fact.xtab)
}
##' {This function takes an fa object and extracts the coefficients from it}
##'  {Extracts the loadings and communalities separately and cbinds them}
##' @title factor_coeff
##' @param x a factor solution
##' @param names names to give the output
##' @return a dataframe suitable for passing to xtable
##' @author Richie Morrisroe
factor_coeff <- function (x, names=NULL) {
   xload <- x$loadings
   xcomm <- x$communality
   xnames <- colnames(xload)
   len <- length(colnames(xload))

    xcommload <- cbind(xload, xcomm)

   xmatdf <- as.matrix.data.frame(xcommload)
   if(!is.null(names)) {
       names2 <- c(names, "Communalites")
       colnames(xmatdf) <- names2
   }
   else {
       xnames[len + 1] <- "Communalities"
       colnames(xmatdf) <- xnames

   }
   return(xmatdf)
}
##' {Extracts the correlations of factors from an fa object}
##' {gets the estimated factor correlations. Only really useful with oblique rotations}
##' @title factor_cor
##' @param x a factor solution
##' @param ... Further arguments passed to the xtable method
##' @return an xtable object containing the between-factor correlations
##' @author Richard Morrisroe
factor_cor <- function (x, ...) {
  res <- x$score.cor
  #allnames <- attr(x$loadings, "dimnames")
  factnames <- colnames(x$loadings)
  res <- as.data.frame(res)
  #names(res) <- factnames
  resx <- xtable::xtable(res, ...)
  resx
}
##'  {This function returns a list of factor names (taken from the FA object) and the items which have an absolute correlation of greater than loadings}
##' {Its not really that complicated, and honestly should be implemented generically. This could prove useful for doing CV across multiple settings, if it returned the indices as well as the names}
##' @title extract_loadings
##' @param x an fa object
##' @param loadings the cutoff point for reporting an association
##' @return the rownames of all items which have an absolute value greater than or equal to loadings
##' @author Richard Morrisroe
extract_loadings <- function (x, loadings=0.3) {
    stopifnot(class(x)==c("psych", "fa"))
  xload <- x$loadings
  xucmat <- as.data.frame(unclass(xload))
    xitemsload <- apply(xucmat, 2,
                        function (y)
                        names(y[which(abs(y) >=
                                      loadings)])) 
  xitemsload
}
##' {Takes a psych fa object, and extracts the communalities and uniquenesses}
##' {See description}
##' @title extract_h2u2
##' @param x a psych fa object
##' @return a data frame containing uniquenesses and communalities
##' @author Richard Morrisroe
extract_h2u2 <- function (x) {
  x.comm <- x$communality
  x.uniq <- x$uniquenesses
  x.ratio <- x.comm / x.uniq
  x.h2u2 <- as.data.frame(cbind(x.comm, x.uniq, x.ratio))
  x.h2u2
 }
##' {Extract TLI (NNFI), BIC and RMSEA from a fa psych object}
##' {see description}
##' @title fit_indices
##' @param x a psych fa object
##' @param labels something unused
##' @return A dataframe containing the fit indices
##' @author Richard Morrisroe
fit_indices <- function (x, labels=NULL) {
    stopifnot(class(x) %in% c("psych", "fa"))
    tli <- x$TLI
    bic <- x$BIC
    rmsea <- x$RMSEA
    rmsnames <- attr(x$RMSEA, "names")
    res <- as.data.frame(cbind(tli, bic, rmsea[1], rmsea[2], rmsea[3]))
    print(length(res))
    colnames(res) <- c(
        paste(substitute(x),"TLI", sep=""),
        paste(substitute(x),"BIC", sep=""),
        paste(substitute(x),"RMSEA", sep=""),
        paste(substitute(x),"-90CInt", sep=""),
        paste(substitute(x),"+90CInt", sep=""))
    res
}
##' {Performs a SVD based CV metric used in chemometrics}
##' {Takes either a wold or gabriel cross validation statistic, and returns a dataframe suitable for printing}
##' @title svd_cv
##' @param x a svdcv factor solution
##' ##' @return a dataframe containing the SVD rank, prediction error and SD of prediction error
##' @author Richie Morrisroe
svd_cv <- function(x) {
    if(!requireNamespace("bcv", quietly=TRUE)) {
        stop("cvsvd is required for this function")
    }
    stopifnot(inherits(x, "cvsvd"))
    msep <- x$msep
    K <- nrow(msep)
    rank <- seq(from = 0, to = x$maxrank, by = 1)
    msep.mean <- apply(x$msep, 2, mean)
    msep.se <- apply(x$msep, 2, sd) / sqrt(K)
    res <- as.data.frame(cbind(rank, msep.mean, msep.se))
    names(res) <- c("Rank", "Prediction Error", "Prediction Error SE")
    res
}
##' {Return the mean, SD, min and max of a set of variables stored in a dataframe}
##'{return an APA standard mean, SE, min and max table for summary statistics}
##' @title apa_demo_tables
##' @param data a dataframe containing all numeric variables
##' @param FUN function to use (currently not implemented)
##' @param xtable Boolean. Return an xtable formatted object or not
##' @param ... other arguments passed to xtable
##' @return a dataframe containing the summary statistics
##' @author Richie Morrisroe
apa_demo_tables <- function(data, FUN=mean, xtable=FALSE, ...) {
    ## stopifnot(require(reshape2), require(plyr))
    fun <- match.call(FUN)
    data.m <- reshape2::melt(data)
    data.tab <- plyr::ddply(data.m,
                            plyr::.(variable),
                            plyr::summarise,
                            Mean=mean(value, na.rm=TRUE),
                            SD=sd(value, na.rm=TRUE),
                            Min=min(value, na.rm=TRUE),
                            Max=max(value, na.rm=TRUE))
    names(data.tab)[1] <- ""
    if(xtable == TRUE) {
        data.tab <- xtable::xtable(data.tab)
    }
    return(data.tab)
}


##' {Give meaningful names to an fa solution}
##' {As description}
##' @title factor_names
##' @param fac a psych fa object
##' @param names a character vector of names to give to each factor (in order)
##' @return a psych fa object with the names applied to the coefficient matrix
##' @author Richie Morrisroe
factor_names <- function(fac, names=NULL) {
    if(is.null(names)) {
        stop(
            "Calling a function based on factor names
with no names seems like a bad idea, don't you think?")
    }
    colnames(fac$loadings) <- names
    return(fac)
}
