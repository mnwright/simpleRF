
##' @title Survival tree class
##' @description Subclass for survival tree.
##' Contains all fields and methods used special for survival trees.
##' @importFrom coin logrank_trafo
TreeSurvival <- setRefClass("TreeSurvival", 
  contains = "Tree",
  fields = list(
     timepoints = "numeric",
     survival = "list"),
  methods = list(
    
    splitNodeInternal = function(nodeID, possible_split_varIDs) {
      response <- data$subset(sampleIDs[[nodeID]], 1)
      
      ## Check node size, stop if maximum reached
      if (length(sampleIDs[[nodeID]]) <= min_node_size) {
        computeSurvival(nodeID, response)
        return(NULL)
      }
      
      ## Stop if node is pure      
      if (length(unique(response[, 1])) < 2 & length(unique(response[, 2])) < 2) {
        computeSurvival(nodeID, response)
        return(NULL)
      }
      
      ## Stop if no events      
      if (all(response[, 2] == 0)) {
        computeSurvival(nodeID, response)
        return(NULL)
      }
            
      ## Find best split, stop if no decrease of impurity
      return(findBestSplit(nodeID, possible_split_varIDs))
    }, 
    
    findBestSplit = function(nodeID, possible_split_varIDs) {   
      ## Initialize
      best_split <- NULL
      best_split$teststat <- -1
      best_split$varID <- -1
      best_split$value <- -1
      
      ## Get response
      response <- data$subset(sampleIDs[[nodeID]], 1)
      
      ## For all possible variables
      for (i in 1:length(possible_split_varIDs)) {
        split_varID <- possible_split_varIDs[i]
        data_values <- data$subset(sampleIDs[[nodeID]], split_varID)
        
        ## Try next variable if all equal for this
        if (length(unique(data_values)) < 2) {
          next()
        }
        
        ## Handle ordered factors
        if (!is.numeric(data_values) & !is.ordered(data_values) & unordered_factors == "order_split") {
          ## Order factor levels
          num.response <- coin::logrank_trafo(response, ties.method = "Hothorn-Lausen")
          means <- sapply(levels(data_values), function(x) {
            mean(num.response[data_values == x])
          })
          levels.ordered <- as.character(levels(data_values)[order(means)])
          
          ## Get all levels not in node
          levels.missing <- setdiff(levels(data_values), levels.ordered)
          levels.ordered <- c(levels.missing, levels.ordered)

          ## Return reordered factor
          data_values <- factor(data_values, levels = levels.ordered, ordered = TRUE)
        }
        
        ## If still not ordered, use partition splitting
        if (!is.numeric(data_values) & !is.ordered(data_values)) {
          best_split = findBestSplitValuePartition(split_varID, data_values, best_split, response)
          
          ## Set split levels left
          if (best_split$varID == split_varID) {
            split_levels_left[[nodeID]] <<- best_split$values_left
          }
        } else {
          best_split = findBestSplitValueOrdered(split_varID, data_values, best_split, response)
          
          ## Set split levels left (empty if ordered splitting)
          if (unordered_factors == "order_split") {
            if (best_split$varID == split_varID) {
              split_levels_left[[nodeID]] <<- unique(data_values[data_values <= best_split$value])
              
              if (is.factor(data_values)) {
                ## Use same splits as in partition
                ints <- as.integer(factor(split_levels_left[[nodeID]], levels = levels(data$subset(sampleIDs[[nodeID]], split_varID))))
                if (sum(2^(ints-1)) >= 2^(max(as.numeric(data$subset(sampleIDs[[nodeID]], split_varID))) - 1)) {
                  split_levels_left[[nodeID]] <<- unique(data_values[data_values > best_split$value])
                }
              }
            }
          } else {
            if (best_split$varID == split_varID) {
              split_levels_left[[nodeID]] <<- list()
            }
          }
        }
      }
      
      if (best_split$varID < 0) {
        ## Stop if no good split found
        computeSurvival(nodeID, response)
        return(NULL)
      } else {
        ## Return best split
        result <- NULL
        result$varID <- as.integer(best_split$varID)
        result$value <- best_split$value
        return(result)
      }      
    }, 
    
    findBestSplitValueOrdered = function(split_varID, data_values, best_split, response) {
      ## For all possible splits
      possible_split_values <- unique(data_values)
      for (j in 1:length(possible_split_values)) {
        split_value <- possible_split_values[j]
        idx <- data_values <= split_value
        
        ## Skip if one child empty
        if (sum(idx) == 0 | sum(!idx) == 0) {
          next
        }
        
        ## Compute test statistic depending on splitrule
        if (splitrule == "Logrank") {
          ## Compute logrank test (skip split point on error)
          teststat <- tryCatch(survdiff(Surv(response[, 1], response[, 2]) ~ idx)$chisq, 
                   error = function(e) e)
          if (inherits(teststat, "error")) {
            next
          }
        } else if (splitrule == "mean") {
          time <- response[, 1]
          time_child1 <- time[idx]
          time_child2 <- time[!idx]
          
          mean_child1 <- mean(time_child1)
          mean_child2 <- mean(time_child2)
          
          mean_diff_squared <- (mean_child1 - mean_child2)^2
          teststat <- mean_diff_squared
        } else {
          stop("Unknown splitrule.")
        }
        
        ## Use this split if better than before
        if (teststat > best_split$teststat) {
          best_split$value <- split_value
          best_split$varID <- split_varID
          best_split$teststat <- teststat
        }
      }
      return(best_split)
    },
    
    findBestSplitValuePartition = function(split_varID, data_values, best_split, response) {
      ## For all possible splits
      possible_split_values <- sort(unique(data_values))
      
      ## For all 2^(n-1)-1 2-partitions
      num_partitions <- 2^(length(possible_split_values) - 1) - 1
      for (j in 1:num_partitions) {
        ## Convert number to logic vector
        left_idx <- as.bitvect(j, length = length(possible_split_values))
        values_left <- possible_split_values[left_idx]
        idx <- data_values %in% values_left
        
        ## Skip if one child empty
        if (sum(idx) == 0 | sum(!idx) == 0) {
          next
        }
        
        ## Compute test statistic depending on splitrule
        if (splitrule == "Logrank") {
          ## Compute logrank test (skip split point on error)
          teststat <- tryCatch(survdiff(Surv(response[, 1], response[, 2]) ~ idx)$chisq, 
                               error = function(e) e)
          if (inherits(teststat, "error")) {
            next
          }
        } else if (splitrule == "mean") {
          time <- response[, 1]
          time_child1 <- time[idx]
          time_child2 <- time[!idx]
          
          mean_child1 <- mean(time_child1)
          mean_child2 <- mean(time_child2)
          
          mean_diff_squared <- (mean_child1 - mean_child2)^2
          teststat <- mean_diff_squared
        } else {
          stop("Unknown splitrule.")
        }
        
        ## Use this split if better than before
        if (teststat > best_split$teststat) {
          best_split$values_left <- values_left
          best_split$varID <- split_varID
          best_split$teststat <- teststat
        }
      }
      return(best_split)
    },
    
    estimate = function(nodeID) {      
      ## Return only NA, value is not used later
      return(NA)
    }, 

    getNodePrediction = function(nodeID) {
      return(survival[[nodeID]])
    }, 
    
    computeSurvival = function(nodeID, response) {      
      survival[[nodeID]] <<- rep(NA, length(timepoints))
      surv_value <- 1    
      
      ## Compute Kaplan-Meier estimator for each timepoint
      for (i in 1:length(timepoints)) {
        num_death <- sum(response[, 2] == 1 & response[, 1] == timepoints[i])
        num_samples_at_risk <- sum(response[, 1] > timepoints[i]) + num_death
        
        if (num_samples_at_risk > 0) {
          surv_value <- surv_value * (num_samples_at_risk - num_death) / num_samples_at_risk
        } 
        survival[[nodeID]][i] <<- surv_value
      }  
    }, 
    
    predictionError = function(pred = NULL) {
      if (is.null(pred)) {
        pred <- predictOOB()
      }
      sum(pred != as.numeric(data$subset(oob_sampleIDs, 1)), na.rm = TRUE) / data$nrow
      
      ## Return brier score for OOB predictions
      bs <- ipred::sbrier(data$subset(oob_sampleIDs, 1), t(pred), timepoints)
      bs[1]
    })
)