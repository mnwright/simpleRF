##' Implements Random Forests (Breiman 2001) with emphasis on simplicity. 
##' Uses reference classes and only plain \code{R}. 
##' Not optimized for computation speed. 
##' Allows rapid prototyping of RF-type algorithms.
##' 
##' Unordered factor variables can be handled in different ways. 
##' Use "ignore" to treat them as ordered in the order of the factor levels. 
##' With "order_once" and "order_split" they are ordered by their response values. For "order_once" this is done once before the analysis, for "order_split" this is done in each split.
##' With "partition" all 2-partitions of the factor levels are considered for splitting.
##' 
##' @title simpleRF
##' @param formula Object of class \code{formula} or \code{character} describing the model to fit.
##' @param data Training data of class \code{data.frame}.
##' @param num_trees Number of trees.
##' @param mtry Number of variables to possibly split at in each node.
##' @param min_node_size Minimal node size. Default 1 for classification, 5 for regression and 3 for survival.
##' @param replace Sample with replacement. Default TRUE.
##' @param probability Grow a probability forest. Default FALSE.
##' @param splitrule Splitrule to use in trees. Default "Gini" for classification forests, "Variance" for regression and probability forests and "Logrank" for survival forests.
##' @param unordered_factors How to handle unordered factor variables. One of "ignore", "order_once", "order_split" and "partition" with default "ignore".
##' @examples 
##' \donttest{
##' library(simpleRF)
##' 
##' # Classification
##' simpleRF(Species ~ ., iris)
##' 
##' # Prediction
##' train_idx <- sample(nrow(iris), 2/3 * nrow(iris))
##' iris_train <- iris[train_idx, ]
##' iris_test <- iris[-train_idx, ]
##' rf_iris <- simpleRF(Species ~ ., data = iris_train)
##' pred_iris <- rf_iris$predict(iris_test)
##' table(iris_test$Species, pred_iris)
##' }
##' 
##' @author Marvin N. Wright
##' @references
##' Breiman, L. (2001). Random forests. Mach Learn, 45(1), 5-32. \cr
##' @export
simpleRF <- function(formula, data, num_trees = 50, mtry = NULL, 
                     min_node_size = NULL, replace = TRUE, probability = FALSE, 
                     splitrule = NULL, unordered_factors = "ignore") {
  
  model.data <- model.frame(formula, data)
  
  if (class(model.data[, 1]) == "factor") {
    if (probability) {
      treetype <- "Probability" 
    } else {
      treetype <- "Classification"
    }
  } else if (class(model.data[, 1]) == "numeric") {
    treetype <- "Regression"
  } else if (class(model.data[, 1]) == "Surv") {
    treetype <- "Survival"
  } else {
    stop("Unkown response type.")
  }
  
  ## Check parameters
  if (is.null(mtry)) {
    mtry <- sqrt(ncol(model.data)-1)
  } else if (mtry > ncol(model.data)-1) {
    stop("Mtry cannot be larger than number of independent variables.")
  }
  if (is.null(min_node_size)) {
    if (treetype == "Classification") {
      min_node_size <- 1
    } else if (treetype == "Probability") {
      min_node_size <- 1
    } else if (treetype == "Regression") {
      min_node_size <- 5
    } else if (treetype == "Survival") {
      min_node_size <- 3
    }
  }
  
  ## Splitrule
  if (is.null(splitrule)) {
    if (treetype == "Classification") {
      splitrule <- "Gini"
    } else if (treetype == "Probability") {
      splitrule <- "Variance"
    } else if (treetype == "Regression") {
      splitrule <- "Variance"
    } else if (treetype == "Survival") {
      splitrule <- "Logrank"
    }
  }
  
  ## Unordered factors
  if (!(unordered_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for unordered_factors.")
  }
  if (unordered_factors == "order_once") {
    ## TODO: Move to function!
    ## Recode characters and unordered factors
    character.idx <- sapply(model.data[, -1], is.character)
    ordered.idx <- sapply(model.data[, -1], is.ordered)
    factor.idx <- sapply(model.data[, -1], is.factor)
    recode.idx <- character.idx | (factor.idx & !ordered.idx)
    
    ## Numeric response
    response <- model.data[, 1]
    if (is.factor(response)) {
      num.response <- as.numeric(response)
    } else if ("Surv" %in% class(response)) {
      num.response <- response[, 1]
    } else {
      num.response <- response
    }
    
    ## Recode each column
    model.data[, -1][, recode.idx] <- lapply(model.data[, -1][, recode.idx], function(x) {
      ## Order factor levels
      means <- aggregate(num.response~x, FUN=mean)
      levels.ordered <- means$x[order(means$num.response)]
      
      ## Return reordered factor
      factor(x, levels = levels.ordered, ordered = TRUE)
    })
    
    ## TODO: Where to save levels? Or just save left child values in all cases?
    ## Save levels
    covariate.levels <- lapply(model.data[, -1], levels)
  }
  if (unordered_factors == "ignore") {
    ## Just set to ordered if "ignore"
    model.data[, -1] <- lapply(model.data[, -1], as.ordered)
  }
    
  ## Create forest object
  if (treetype == "Classification") {
    forest <- ForestClassification$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                       min_node_size = as.integer(min_node_size), 
                                       replace = replace, splitrule = splitrule,
                                       data = Data$new(data = model.data), 
                                       formula = formula,  unordered_factors = unordered_factors, 
                                       response_levels = levels(model.data[, 1]))
  } else if (treetype == "Probability") {
    forest <- ForestProbability$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                   min_node_size = as.integer(min_node_size), 
                                   replace = replace, splitrule = splitrule,
                                   data = Data$new(data = model.data), 
                                   formula = formula, unordered_factors = unordered_factors,
                                   response_levels = levels(model.data[, 1]))
  } else if (treetype == "Regression") {
    forest <- ForestRegression$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                   min_node_size = as.integer(min_node_size), 
                                   replace = replace, splitrule = splitrule,
                                   data = Data$new(data = model.data), 
                                   formula = formula, unordered_factors = unordered_factors)
  } else if (treetype == "Survival") {
    idx.death <- model.data[, 1][, 2] == 1
    timepoints <- sort(unique(model.data[idx.death, 1][, 1]))
    forest <- ForestSurvival$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                 min_node_size = as.integer(min_node_size), 
                                 replace = replace, splitrule = splitrule,
                                 data = Data$new(data = model.data), 
                                 formula = formula, unordered_factors = unordered_factors, 
                                 timepoints = timepoints)
  } else {
    stop("Unkown tree type.")
  }

  ## Grow forest
  forest$grow()

  ## Return forest
  return(forest) 
}