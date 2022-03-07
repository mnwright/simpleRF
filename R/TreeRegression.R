
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
      best_split <- NULL
      best_split$decrease <- -1
      best_split$varID <- -1
      best_split$value <- -1
      
      ## Get response
      response <- data$subset(sampleIDs[[nodeID]], 1)
      
      ## For all possible variables
      for (i in 1:length(possible_split_varIDs)) {
        split_varID <- possible_split_varIDs[i]
        data_values <- data$subset(sampleIDs[[nodeID]], split_varID)
        
        ## Handle ordered factors
        if (!is.numeric(data_values) & !is.ordered(data_values) & unordered_factors == "order_split") {
          ## Order factor levels
          num.response <- as.numeric(response)
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

        ## Sum responses in childs
        idx <- data_values <= split_value
        response_left <- response[idx]
        response_right <- response[!idx]
        
        ## Skip if one child empty
        if (length(response_left) == 0 | length(response_right) == 0) {
          next
        }
        
        if (splitrule == "Variance") {
          ## Decrease of impurity
          decrease <- sum(response_left)^2/length(response_left) + 
            sum(response_right)^2/length(response_right)
        } else {
          stop("Unknown splitrule.")
        }

        ## Use this split if better than before
        if (decrease > best_split$decrease) {
          best_split$value <- split_value
          best_split$varID <- split_varID
          best_split$decrease <- decrease
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
        
        ## Sum responses in childs
        idx <- data_values %in% values_left
        response_left <- response[idx]
        response_right <- response[!idx]
        
        ## Skip if one child empty
        if (length(response_left) == 0 | length(response_right) == 0) {
          next
        }
        
        if (splitrule == "Variance") {
          ## Decrease of impurity
          decrease <- sum(response_left)^2/length(response_left) + 
            sum(response_right)^2/length(response_right)
        } else {
          stop("Unknown splitrule.")
        }
        
        ## Use this split if better than before
        if (decrease > best_split$decrease) {
          best_split$values_left <- values_left
          best_split$varID <- split_varID
          best_split$decrease <- decrease
        }
      }
      return(best_split)
    },
    
    estimate = function(nodeID) {      
      ## Return mean response
      return(mean(data$subset(sampleIDs[[nodeID]], 1)))
    }, 
    
    getNodePrediction = function(nodeID) {
      return(split_values[nodeID])
    }, 
    
    predictionError = function(pred = NULL) {
      if (is.null(pred)) {
        pred <- predictOOB()
      }
      sum((pred - data$subset(oob_sampleIDs, 1))^2, na.rm = TRUE) / length(oob_sampleIDs)
    })
    
)