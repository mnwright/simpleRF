
##' @title Classification tree class
##' @description Subclass for classification tree.
##' Contains all fields and methods used special for classification trees.
TreeClassification <- setRefClass("TreeClassification", 
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
      best_split$decrease <- -1;
      best_split$varID <- -1;
      best_split$value <- -1;
      
      ## Get response
      response <- data$subset(sampleIDs[[nodeID]], 1)
      
      ## For all possible variables
      for (i in 1:length(possible_split_varIDs)) {
        split_varID <- possible_split_varIDs[i]
        data_values <- data$subset(sampleIDs[[nodeID]], split_varID)
        
        ## Handle ordered factors
        if (!is.ordered(data_values) & unordered_factors == "order_split") {
          ## TODO: Move to function?
          ## Numeric response
          if (is.factor(response)) {
            num.response <- as.numeric(response)
          } else if ("Surv" %in% class(response)) {
            num.response <- response[, 1]
          } else {
            num.response <- response
          }
          
          ## Order factor levels
          means <- aggregate(num.response ~ data_values, FUN=mean)
          levels.ordered <- means$data_values[order(means$num.response)]
          
          ## Return reordered factor
          data_values <- factor(data_values, levels = levels.ordered, ordered = TRUE)
        }
        
        ## If still not ordered, use partition splitting
        if (!is.ordered(data_values)) {
          ## TODO: set split_levels_left somewhere
          best_split = findBestSplitValuePartition(data_values, best_split)
        } else {
          best_split = findBestSplitValueOrdered(split_varID, data_values, best_split, response)
          
          ## Set split levels left (empty if ordered splitting)
      #    if (unordered_factors == "order_split") {
          #browser()
          if (best_split$varID == split_varID) {
            split_levels_left[[nodeID]] <<- unique(data_values[data_values <= best_split$value])
          }
            #browser()
      #    } else {
      #      split_levels_left[[nodeID]] <<- list()
      #    }
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
      ## Convert to numeric if necessary
      # if(!is.numeric(data_values)) {
      #   data_values <- as.numeric(data_values)
      # }

      ## For all possible splits
      possible_split_values <- unique(data_values)
      for (j in 1:length(possible_split_values)) {
        split_value <- possible_split_values[j]
        
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
        if (decrease > best_split$decrease) {
          best_split$value <- split_value
          best_split$varID <- split_varID
          best_split$decrease <- decrease
        }
      }
      return(best_split)
    },
    
    findBestSplitValuePartition = function(data_values, best_split) {
      
      ## TODO: Implement
      stop("Partition splitting not implemented yet")
      return(best_split)
    },
    
    estimate = function(nodeID) {      
      ## Return class with maximal count, random at ties
      class_counts <- table(data$subset(sampleIDs[[nodeID]], 1))
      which.max.random(class_counts)
    }, 
    
    getNodePrediction = function(nodeID) {
      return(split_values[nodeID])
    })
)