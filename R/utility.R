##' Find index of maximum, breaking ties at random and ignoring NAs.
##' @title Index of maximum with random ties
##' @param x Input vector.
##' @return Index of the maximum value.
##' @author Marvin N. Wright
which.max.random <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  which(rank(x, ties.method = "random", na.last = FALSE) == length(x))
}

##' Find index of minimum, breaking ties at random and ignoring NAs.
##' @title Index of minimum with random ties
##' @param x Input vector.
##' @return Index of the minimum value.
##' @author Marvin N. Wright
which.min.random <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  which(rank(x, ties.method = "random", na.last = TRUE) == 1)
}

##' Order factor levels with correlation approach
##' @title Order factor levels with correlation approach
##' @param y Response factor.
##' @param x Covariate factor.
##' @return Ordered factor levels
##' @author Marvin N. Wright
cor.order <- function(y, x) {
  ## Create contingency table of the nominal outcome with the nominal covariate
  tab <- table(droplevels(y), droplevels(x))
  
  ## Compute correlation matrix of the contingency table (correlation of the covariate levels w.r.t outcome)
  cr <- suppressWarnings(cor(tab))
  cr[is.na(cr)] <- 0                      
  diag(cr) <- NA
  
  ## Start with a random level and select as next level the level with highest correlation to the current level (excluding already selected levels)
  num_levels <- nlevels(droplevels(x))
  next_level <- sample(num_levels, 1)
  res <- c(next_level, rep(NA, num_levels - 1))
  for (i in 2:num_levels) {
    cr[, next_level] <- NA
    next_level <- which.max.random(cr[next_level, ])
    res[i] <- next_level
  }
  
  ## Return ordered factor levels
  as.character(levels(droplevels(x))[res])
}

##' Order factor levels with PCA approach
##' @title Order factor levels with PCA approach
##' @param y Response factor.
##' @param x Covariate factor.
##' @return Ordered factor levels
##' @author Marvin N. Wright
##' @references Coppersmith, D., Hong, S.J. & Hosking, J.R. (1999) Partitioning Nominal Attributes in Decision Trees. Data Min Knowl Discov 3:197. \url{https://doi.org/10.1023/A:1009869804967}.
pc.order <- function(y, x) {
  if (nlevels(droplevels(x)) < 2) {
    return(as.character(levels(droplevels(x))))
  }
  
  ## Create contingency table of the nominal outcome with the nominal covariate
  N <- table(droplevels(x), droplevels(y))

  ## PCA of weighted covariance matrix of class probabilites
  P <- N/rowSums(N)
  S <- cov.wt(P, wt = rowSums(N))$cov
  pc1 <- prcomp(S, rank. = 1)$rotation
  score <- P %*% pc1
  
  ## Return ordered factor levels
  as.character(levels(droplevels(x))[order(score)])
}

##' Reorder factor columns. Use mean for continuous response, class counts for factors and mean survival for survival response.
##' @title Reorder factor columns
##' @param data Data with factor columns.
##' @return Data with reordered factor columns.
##' @importFrom coin logrank_trafo
##' @author Marvin N. Wright
reorder.factor.columns <- function(data) {
  ## Recode characters and unordered factors
  character.idx <- sapply(data[, -1], is.character)
  ordered.idx <- sapply(data[, -1], is.ordered)
  factor.idx <- sapply(data[, -1], is.factor)
  recode.idx <- character.idx | (factor.idx & !ordered.idx)
  
  ## Numeric response
  response <- data[, 1]
  if (is.factor(response)) {
    num.response <- as.numeric(response)
  } else if ("Surv" %in% class(response)) {
    num.response <- coin::logrank_trafo(response, ties.method = "Hothorn-Lausen")
  } else {
    num.response <- response
  }
  
  ## Recode each column
  data[, -1][, recode.idx] <- lapply(data[, -1][, recode.idx, drop = FALSE], function(x) {
    if (is.factor(response) & nlevels(response) > 2) {
      levels.ordered <- pc.order(y = response, x = x)
    } else {
      ## Order factor levels by num.response
      means <- sapply(levels(x), function(y) {
        mean(num.response[x == y])
      })
      levels.ordered <- as.character(levels(x)[order(means)])
    }
    
    ## Return reordered factor
    factor(x, levels = levels.ordered, ordered = TRUE)
  })
  
  ## Return data
  data
}

##' Convert number to bit vector.
##' @title Number to bit vector
##' @param x Input number.
##' @param length Length of output bit vector.
##' @return Bit vector.
##' @author Marvin N. Wright
as.bitvect <- function(x, length = 32) {
  i <- 1
  string <- numeric(length)
  while(x > 0) {
    string[i] <- x %% 2
    x <- x %/% 2
    i <- i + 1 
  }
  as.logical(string)
}


