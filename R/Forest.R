
##' @title Forest class
##' @description Virtual class for Random forest. 
##' Contains all fields and methods used in all Forest subclasses.
##' @importFrom parallel mclapply
##' @import methods
Forest <- setRefClass("Forest", 
  fields = list(
    num_trees = "integer", 
    mtry = "integer", 
    min_node_size = "integer", 
    splitrule = "character",
    unordered_factors = "character",
    data = "Data",
    predict_data = "Data",
    formula = "formula",
    trees = "list",
    treetype = "character",
    replace = "logical", 
    covariate_levels = "list"),
  methods = list(
    
    grow = function(num_threads) { 
      
      ## Init trees
      temp <- lapply(trees, function(x) {
        x$mtry <- mtry
        x$min_node_size <- min_node_size
        x$splitrule <- splitrule
        x$unordered_factors <- unordered_factors
        x$data <- data
      })
      
      ## Grow trees
      trees <<- mclapply(trees, function(x) {
        x$grow(replace)
        x
      }, mc.cores = num_threads)
    }, 
    
    predict = function(newdata) {
      model.data <- model.frame(formula, newdata)

      ## Recode factors if forest grown 'order_once' mode
      if (unordered_factors == "order_once" & length(covariate_levels) > 0) {
        model.data[, -1] <- mapply(function(x, y) {
          if(is.null(y)) {
            x
          } else {
            new.levels <- setdiff(levels(x), y)
            factor(x, levels = c(y, new.levels), ordered = TRUE)
          }
        }, model.data[, -1], covariate_levels, SIMPLIFY = FALSE)
      }

      ## Save prediction data in model
      predict_data <<- Data$new(data = model.data)
      
      ## Predict in trees
      predictions <- simplify2array(lapply(trees, function(x) {
        x$predict(predict_data)
      }))
      
      ## Aggregate predictions
      return(aggregatePredictions(predictions))
    }, 
    
    aggregatePredictions = function(predictions) {
      ## Empty virtual function
    }, 
    
    predictionError = function() {
      ## Empty virtual function
    },
    
    show = function() {
      cat("simpleRF Forest\n")
      cat("Type:                            ", treetype, "\n")
      cat("Splitrule:                       ", splitrule, "\n")
      cat("Number of trees:                 ", num_trees, "\n")
      cat("Sample size:                     ", data$nrow, "\n")
      cat("Number of independent variables: ", data$ncol-1, "\n")
      cat("Mtry:                            ", mtry, "\n")
      cat("Target node size:                ", min_node_size, "\n")
      cat("Replace                          ", replace, "\n")
      cat("Unordered factor handling        ", unordered_factors, "\n")
      cat("OOB prediction error:            ", predictionError(), "\n")
    }, 
    
    print = function() {
      show()
    })
)


