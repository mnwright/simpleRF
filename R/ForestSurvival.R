
##' @title Survival forest class
##' @description Subclass for survival forest. 
##' Contains all fields and methods used special for survival forests.
##' @importFrom ipred sbrier
ForestSurvival <- setRefClass("ForestSurvival", 
  contains = "Forest", 
  fields = list(
     timepoints = "numeric"), 
  methods = list(
    
    grow = function(num_threads) {      
      treetype <<- "Survival"
      
      ## Create trees
      trees <<- replicate(num_trees, TreeSurvival$new(timepoints = timepoints))
      
      ## Call parent method
      callSuper(num_threads)
    }, 
    
    predict = function(newdata) {
      callSuper(newdata)
    },
    
    aggregatePredictions = function(predictions) {
      ## For each person and timepoint mean over trees      
      result <- apply(predictions, c(2,1), mean)      
      return(result)
    }, 
    
    predictionError = function() {      
      ## For each sample mean over trees where sample is OOB
      tree_predictions <- simplify2array(lapply(trees, function(x) {
        oob_samples <- x$oob_sampleIDs
        result <- matrix(NA, length(timepoints), data$nrow)
        result[, oob_samples] <- x$predictOOB()        
        return(result)
      }))
      oob_predictions <- apply(tree_predictions, c(2,1), mean, na.rm = TRUE)
      
      ## Return brier score for OOB predictions
      idx <- rowSums(is.na(oob_predictions)) == 0
      bs <- ipred::sbrier(data$subset(idx, 1), t(oob_predictions[idx, ]), timepoints)
      return(bs[1])
    })
)


