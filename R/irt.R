##' compare Z-scores for two IRT models
##' {Examine the differences in z-scores between two (and only two) IRT models}
##' @title compare_irt_scores
##' @param x an IRT model
##' @param y an IRT model
##' @return a list containing the (Pearson) correlations between the two z-scores, and the squared differences between the two sets of scores
##' @author Richie Morrisroe
compare_irt_scores <- function (x, y) {
    scores_x <- x$z1
    scores_y <- y$z1
    cor_xy <- cor(scores_x, scores_y,
                  method="pearson", use="pairwise.complete.obs")
    diff_xy <- (scores_x - scores_y) ^ 2
    res <- list(cor=cor_xy, differences=diff_xy)
    res
}
##' {Actual implemented IRT test on new data function}
##' {Examines the difference in accuracy between a model estimated on the new data, versus the predictions from the old model on the new data} 
##' @title test_irt_models
##' @param oldmodel the original IRT model
##' @param newdata the new data
##' @param gpcmconstraint the constraint if the model is gpcm
##' @param grmconstraint the constraint if the model is grm
##' @param ... other arguments passed through
##' @return A dataframe containing two columns, ErrorApproximation and Correlation between models
##' @author Richie Morrisroe
test_irt_models <- function(oldmodel,
                          newdata,
                          gpcmconstraint=c("rasch", "1PL", "gpcm",),
                          grmconstraint= c(TRUE, FALSE), ...) {
    if(class(oldmodel)=="gpcm") {
        constraint <- gpcmconstraint
    }
    else {
        constraint <- grmconstraint
    }

    comp.para <- length(unique(as.vector(coef(oldmodel))))
    predscores <- get_irt_estimates(
        ltm::factor.scores(oldmodel, resp.patterns=newdata))
    if(class(oldmodel) == "gpcm") {
        newmodel <- ltm::gpcm(newdata, constraint=constraint)
    }
    else {
        newmodel <- ltm::grm(newdata, constrained=constraint)
    }
    newscores <- get_irt_estimates(
        ltm::factor.scores(newmodel, resp.patterns=newdata))
    diffscores <- mapply("-", predscores[,1], newscores[,1])
    rea <- sqrt(sum(diffscores ^ 2)) * log(comp.para)
    scorescor <- cor(predscores[,1], newscores[,1], ...)
    res <- data.frame(error_approximation=rea, correlation=scorescor)
    return(res)
}
##' {Average a set of IRT Models fit on different cross-validation splits} 
##' {Simple averaging of the coefficients}
##' @title irt_average
##' @param sols IRT Models
##' @return a dataframe containing the averaged coefficients
##' @author Richie Morrisroe
irt_average <- function(sols=list()) {
    coef <- lapply(sols, coef)
    res <- Reduce(`+`, x=coef)/length(coef)
    return(res)
}
##' {A function that returns the coefficients from an ltm model} 
##'
##' {Additionally handles averaged coefficients and standard errors. This needs to be refactored into calculation vs display functions and is waaaayyyy too long}
##' @title coef_irt
##' @param grm model
##' @param se return standard errors?
##' @param averaged return averaged errors
##' @return a set of IRT coefficients in
##' @author Richard Morrisroe
coef_irt <- function(grm, se=FALSE, averaged=FALSE) {
    ##TODO: Refactor this into a coef method, and a display method
    if(class(grm)=="grm" || class(grm)=="gpcm") {
        dat <- coef(grm)
        itemnames <- rownames(coef(grm))
    }
    else {
        dat <- grm
    }
    if(averaged) {
        if(se) {
            dat <- grm$coef
        } else {
            dat <- grm
        }
        dat <- round(dat, 2)
        itemnames <- rownames(dat)
        if(se) {
        standerr <- t(grm$se)
        standerr <- round(standerr, 2)
    }
    }


    if(se & !averaged ) {
        standard.errors <- sapply(coef(summary(grm)), '[[', "std.err")
        standerr <- round(standard.errors, 2)
    }
    dims <- dim(dat)
    betas <- dims[2]-1
    if(se) {
        dat <- as.matrix(dat)
        dat.se <- matrix(paste(dat,
                               " (",
                               standerr,
                               ")",
                               sep=""),
                         nrow=dims[1],
                         ncol=dims[2])

        rownames(dat.se) <- itemnames
        dat <- dat.se
    }
    beta.names <- paste("$", "\\beta",  "^", 1:betas, " (se) ", "$", sep="")
    alpha.name <- "$\\alpha$"
    allnames <- c(beta.names, alpha.name)
    colnames(dat) <- allnames
    dat
}
##' {Return an xtable object of an IRT easiness/difficulty parameters}
##' {Takes a GRM or GPCM object and returns a decidely non-standard table} 
##' @title irt_xtab
##' @param x an IRT model object of <some_bunch> of classes
##' @param ... arguments based to xtable function
##' @return an xtable representation of the difficulty and/or discrimination parameters
##' @author Richard Morrisroe
irt_xtab <- function (x, ...) {
    eta <- x$etapar #$
    se <- x$se.eta #$
    eta_mat <- as.matrix(eta)
    se_eta_mat <- as.matrix(se)
    eta_par_mat <- cbind(eta_mat, se_eta_mat)
    colnames(eta_par_mat) <- c("Ability Estimate", "Standard Error")
    coef_xtab <- xtable::xtable(eta_par_mat, ...)
    coef_xtab
}
##' {Wrapper around ggplot for a person-item difficulty plot}
##' {Not much, really}
##' @title ggplot_grm
##' @param grm an IRT GRM model object
##' @param ... other methods passed to plotting function
##' @return a ggplot object
##' @author Richard Morrisroe
ggplot_grm <- function (grm, ...) {
    stopifnot(class(grm) == "grm")
    x <- coef(grm)
    x <- as.matrix(x)
    x <- x[,-ncol(x)]
    xt <- t(x)
    response <- 1:nrow(xt)
    respind <- ncol(xt)+1
    xt <- as.data.frame(xt)
    xt$response <- response
    xtm <- reshape2::melt(xt, id="response")
    names(xtm) <- c("threshold", "item", "ability")
    plot1 <- ggplot2::ggplot(xtm,
                             aes(x=ability,
                                 y=item,
                                 shape=as.factor(threshold),
                                 colour=as.factor(threshold)), ...)
    plot2 <- plot1 + ggplot2::geom_point() + ggplot2::geom_rug()
    plot2
}
##' {Convert a GPCM object to a matrix for turning into a table} 
##' {In some cases, the output of a gpcm will have a different number of threshold parameters for different items.
##' This function extracts the coefficients from an object of class gpcm, and solves this problem so that the coefficients can be coerced to a data.frame or matrix and the tables reported easily}
##' @title coef2mat
##' @param gpcm a gpcm object (not relevant for GRM)
##' @return a matrix containing the estimated parameters from the gpcm model
##' @author Richard Morrisroe
coef2mat <- function (gpcm) {
    if(is.matrix(gpcm)) {
        return(gpcm)
    }
    else {

        len <- lapply(gpcm, length)
        probelem <- which.min(as.matrix(unlist(len)))
        dimcols <- max(as.matrix(unlist(len)))
        dimrows <- length(names(gpcm))
        mat_res <- matrix(NA, nrow=dimrows, ncol=dimcols)
        modlength <- lapply(gpcm, length)
        maxlength <- max(as.matrix(unlist(modlength)))
        for (i in 1:maxlength) {
            column <- lapply(gpcm, "[", i)
            column <- as.matrix(unlist(column))
            mat_res[1:length(column),i] <- column
            mat_res
        }
        rownames(mat_res) <- names(gpcm)
        probelemlength <- length(gpcm[[probelem]])
        ##this gives a scalar, as internally matrices are stored as vectors
        missingvalue <- which(is.na(mat.res))
        #get the element where is the discrimination parameter has ended up
        wrongvalue <- missingvalue - nrow(mat.res) 
        mat_res[missingvalue] <- mat_res[wrongvalue]
        mat_res[wrongvalue] <- NA

        categories <- lapply(gpcm, names)
        categorynames <- categories[[which.max(sapply(categories, length))]]
        colnames(mat_res) <- categorynames
        return(mat_res)
    }
    mat_res
}
##' Extract the predictions from an IRT fascore object
##'
##' {Selects useful information from an object inheriting from grm or gpcm}
##' @title get_irt_preds
##' @param x  an IRT model object
##' @return a dataframe containing observed scores, expected scores, the results of a z-test, and the se of the z-test
##' @author Richie Morrisroe
get_irt_preds <- function (x) {
    res <- x$score.dat[,c("Obs", "Exp", "z1","se.z1")]
    res
}

##' {Unfinished function used to perform cross-validation over IRT models}
##' {See Description}
##' @title irt_cv
##' @param data a dataframe containing the data to be used
##' @param model the kind of model (either grm or gpcm)
##' @param constraint the constraint to use - see documentation for grm and gpcm objects
##' @param splits the number of splits to use
##' @param ... extra arguments passed through to other methods 
##' @return a test and train set
##' @author Richie Morrisroe
irt_cv <- function (data,
                    model=c("grm", "gpcm"),
                    constraint = c(TRUE, FALSE, "rasch", "1PL", "gpcm"),
                    splits = 10, ...) {
    if(is.dataframe(data) ||is.matrix(data))
        stop("this function needs matrix or dataframe input")
    splittedsamples <- split_sample(data, splits)
    for (i in 1:length(splittedsamples)) {
        testset <- splittedsamples[i]
        trainset <- splittedsamples[!i]
        return(list(train=trainset, test=testset))
    }
}
##' Yet another unsuccessful IRT CV function (look at this one, there were good ideas in there)
##'
##' {See description}
##' @title irt_cross_validate
##' @param x a dataframe containing all numeric variables
##' @return a dataframe containing observed and expected scores
##' @author Richie Morrisroe
irt_cross_validate <- function(x) {
#get observed frequencies from display command in package ltm
    obs <- ltm::descript(x)$perc
    totscores <- grm::descript(x)$items
    totscores[totscores == 0] <- NA
    model <- ltm::grm(x)
    model.scores <- ltm::factor.scores(model, resp.patterns=x)
    abilities <- model.scores$score.dat["z1"]
    pointsweights <- model$GH
    cutpoints <- pointsweights[[1]]
    weights <- pointsweights[[2]]
    q <- seq(from=0, to=1, by=0.05) #create 21 points
    quadnorm <- qnorm(q) # map 21 points to the normal quantiles
    totscores2 <- rowSums(x, na.rm=TRUE)
    totscores2[totscores2 == 0] <- NA
    ab.scores <- as.matrix(cbind(totscores2, abilities))
    res <- list(obsscores=obs,
                totscores=totscores2,
                abscores=ab.scores,
                model=model,
                scores=model.scores,
                abilities=abilities,
                weights=weights)
}
##' {Average a set of IRT Factor Scores across Cross-validation Splits}
##' {Thought I did this above?}
##' @title irt_average_factor_scores
##' @param scores IRT scores 
##' @return the average across all splits
##' @author Richie Morrisroe
irt_average_factor_scores <- function (scores=list) {
    abilities <- sapply(scores, `[`, 1)
    ab_average <- Reduce(`+`, abilities) / length(abilities)
    names(ab_average) <- "ability_estimation"
    return(ab_average)
}
