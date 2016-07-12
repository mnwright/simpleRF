
##' @title Survival tree class
##' @description Subclass for survival tree.
##' Contains all fields and methods used special for survival trees.
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
            
      ## Find best split, stop if no decrease of impurity
      return(findBestSplit(nodeID, possible_split_varIDs))
    }, 
    
    findBestSplit = function(nodeID, possible_split_varIDs) {   
      ## Initialize
      best_teststat <- -1;
      best_varID <- -1;
      best_value <- -1;
      response <- data$subset(sampleIDs[[nodeID]], 1)
      
      ## For all possible variables
      for (i in 1:length(possible_split_varIDs)) {
        split_varID <- possible_split_varIDs[i]
        data_values <- data$subset(sampleIDs[[nodeID]], split_varID)
        
        ## For all possible splits
        possible_split_values <- unique(data_values)
        
        ## Try next variable if all equal for this
        if (length(possible_split_values) < 2) {
          next
        }
        
        for (j in 1:length(possible_split_values)) {
          split_value <- possible_split_values[j]
          
          ## Count classes in childs
          idx <- data_values <= split_value
          
          ## Skip if one child empty
          if (sum(idx) == 0 | sum(!idx) == 0) {
            next
          }
          
          ## Compute test statistic depending on splitrule
          if (splitrule == "Logrank") {
            ## Compute logrank test
            teststat <- survdiff(Surv(response[, 1], response[, 2]) ~ idx)$chisq
          } else if (splitrule == "mean") {
            time <- response[, 1]
            time_child1 <- time[idx]
            time_child2 <- time[!idx]
            
            mean_child1 <- mean(time_child1)
            mean_child2 <- mean(time_child2)
            
            mean_diff_squared <- (mean_child1 - mean_child2)^2
            teststat <- mean_diff_squared
          } else if (splitrule == "Likelihood") {
            stop("Not implemented yet.")
          } else {
            stop("Unknown splitrule.")
          }
          
          ## Use this split if better than before
          if (teststat > best_teststat) {
            best_value <- split_value
            best_varID <- split_varID
            best_teststat <- teststat
          }
        }
      }
      
      if (best_varID < 0) {
        ## Stop if no good split found
        computeSurvival(nodeID, response)
        return(NULL)
      } else {
        ## Return best split
        result <- NULL
        result$varID <- as.integer(best_varID)
        result$value <- best_value
        return(result)
      }      
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
    })
)