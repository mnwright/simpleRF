
##' @title Probability tree class
##' @description Subclass for probability tree.
##' Contains all fields and methods used special for probability trees.
TreeProbability <- setRefClass("TreeProbability", 
  contains = "Tree",
  fields = list(),
  methods = list(
    
    splitNodeInternal = function(nodeID, possible_split_varIDs) {
      ## Check node size, stop if maximum reached
      if (length(sampleIDs[[nodeID]]) <= min_node_size) {
        return(NULL)
      }
      
      ## Stop if node is pure      
      unique_response <- unique(data$subset(sampleIDs[[nodeID]], 1))
      if (length(unique_response) == 1) {
        return(NULL)
      }
      
      ## Find best split, stop if no decrease of impurity
      return(findBestSplit(nodeID, possible_split_varIDs))
    }, 
    
    findBestSplit = function(nodeID, possible_split_varIDs) {
      
      ## Initialize
      best_decrease <- -1
      best_varID <- -1
      best_value <- -1
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
          ## Count classes in childs
          idx <- data_values <= split_value
          class_counts_left <- tabulate(response[idx])
          class_counts_right <- tabulate(response[!idx])
          
          ## Skip if one child empty
          if (sum(class_counts_left) == 0 | sum(class_counts_right) == 0) {
            next
          }
          
          ## TODO: Use splitrule parameter
          ## Decrease of impurity
          decrease <- sum(class_counts_left^2)/sum(class_counts_left) + 
            sum(class_counts_right^2)/sum(class_counts_right)
          
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
      ## Return only NA, value is not used later
      return(NA)
    }, 
    
    getNodePrediction = function(nodeID) {
      ## Return class fractions
      node_samples <- sampleIDs[[nodeID]]
      table(data$subset(node_samples, 1))/length(node_samples)
    })
)