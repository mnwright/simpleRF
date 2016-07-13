
##' @title Regression tree class
##' @description Subclass for regression tree.
##' Contains all fields and methods used special for regression trees.
TreeRegression <- setRefClass("TreeRegression", 
  contains = "Tree",
  fields = list(),
  methods = list(
 
    splitNodeInternal = function(nodeID, possible_split_varIDs) {
      ## Check node size, stop if maximum reached
      if (length(sampleIDs[[nodeID]]) <= min_node_size) {
        return(NULL)
      }
      
      ## Find best split, stop if no decrease of impurity
      return(findBestSplit(nodeID, possible_split_varIDs))
    }, 
    
    findBestSplit = function(nodeID, possible_split_varIDs) {
      
      ## Initialize
      best_decrease <- -1;
      best_varID <- -1;
      best_value <- -1;
      response <- data$subset(sampleIDs[[nodeID]], 1)
      
      ## For all possible variables
      for (i in 1:length(possible_split_varIDs)) {
        split_varID <- possible_split_varIDs[i]
        data_values <- data$subset(sampleIDs[[nodeID]], split_varID)
        
        ## For all possible splits
        possible_split_values <- unique(data_values)
        for (j in 1:length(possible_split_values)) {
          split_value <- possible_split_values[j]
          
          ## TODO: Handle unordered factors
          ## Sum responses in childs
          idx <- data_values <= split_value
          response_left <- response[idx]
          response_right <- response[!idx]
          
          ## Skip if one child empty
          if (length(response_left) == 0 | length(response_right) == 0) {
            next
          }
          
          ## TODO: Use splitrule parameter
          ## Decrease of impurity
          decrease <- sum(response_left)^2/length(response_left) + 
            sum(response_right)^2/length(response_right)
          
          ## Use this split if better than before
          if (decrease > best_decrease) {
            best_value <- split_value
            best_varID <- split_varID
            best_decrease <- decrease
          }
        }
      }
      
      if (best_varID < 0) {
        ## Stop if no good split found
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
      ## Return mean response
      return(mean(data$subset(sampleIDs[[nodeID]], 1)))
    }, 
    
    getNodePrediction = function(nodeID) {
      return(split_values[nodeID])
    })
    
)