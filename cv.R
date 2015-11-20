##' {Yet another split sample function}
##' @title split_sample
##' @param x the data
##' @param split the number of splits to make
##' @return a list containing the splits
##' @author Richie Morrisroe
split_sample <- function(x, split) {
    xlen <- nrow(x)
    indices <- sample(1:xlen, xlen, replace=FALSE)
    splitlen <- xlen/split
    splits <- cut(indices, split, labels=FALSE)
    samplist <- list()
    for(i in 1:max(split)) {
        samp <- x[splits == i,]
        samplist[[i]] <- samp
    }
    samplist
}
##' {Something}
##'
##' {more things} 
##' @title create_combinations
##' @param splits a list containing splits from split_sample
##' @return a list containing the splits into test and train
##' @author Richard Morrisroe
create_combinations <-  function(splits) {
    stopifnot(class(splits) == "list")
    splitnumbers <- length(splits)
    facsplits <- choose(splitnumbers, k = (splitnumbers - 1))
    reslist <- list()
    for (i in 1:facsplits) {
        samples <- sample(
            1:splitnumbers,
            size=1,
            replace=FALSE)
        train <- do.call("rbind", splits[-c(samples)])
        test <- as.data.frame(splits[c(samples)])
        reslist[[i]] <- list(train=train, test=test)
    }
    reslist
}
