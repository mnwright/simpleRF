
##' @title Probability forest class
##' @description Subclass for probability forest. 
##' Contains all fields and methods used special for probability forests.
ForestProbability <- setRefClass("ForestProbability", 
  contains = "Forest", 
  fields = list(
    response_levels = "character"), 
  methods = list(
    
    grow = function(num_threads) {      
      treetype <<- "Probability"
      
      ## Create trees
      trees <<- replicate(num_trees, TreeProbability$new())
      
      ## Call parent method
      callSuper(num_threads)
    }, 
    
    predict = function(newdata) {
      callSuper(newdata)
    },
    
    aggregatePredictions = function(predictions) {
      ## Return class probabilities per sample
      apply(predictions, 1, rowMeans, na.rm = TRUE)
    }, 
    
    predictionError = function() {
        
      ## For each tree loop over OOB samples and get predictions
      tree_predictions <- simplify2array(lapply(trees, function(x) {
        oob_samples <- x$oob_sampleIDs
        result <- matrix(NA, nrow = length(response_levels), ncol = data$nrow)
        result[, oob_samples] <- x$predictOOB()
        return(result)
      }))
      
      ## Compute prediction for each sample
      probabilities <- apply(tree_predictions, 1, rowMeans, na.rm = TRUE)
      
      ## Get probabilities for true classes
      response <- as.numeric(data$column(1))
      num_samples <- nrow(probabilities)
      num_classes <- ncol(probabilities)
      idx <- response == matrix(1:num_classes, byrow = TRUE, 
                                nrow = num_samples, ncol = num_classes)
      true_class_probabilities <- probabilities[idx]
      
      ## Compute mse
      sum((true_class_probabilities - 1)^2, na.rm = TRUE) / data$nrow
    })
)

