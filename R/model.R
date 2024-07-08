setOldClass("cv.glmnet")
#' ROCKET model object
#'
#' @slot coefficients numeric.
#' @slot model cv.glmnet.
#' @slot kernels list.
#' @slot features matrix.
#' @slot rel_features matrix.
#' @slot call call.
#'
#' @return none
#' @export
#'
setClass(
  Class = "rROCKET",
  slots = c(
    coefficients = "numeric",
    model = "cv.glmnet",
    kernels = "list",
    features = "matrix",
    rel_features = "matrix",
    call = "call"
  )
)

#' Fit a ROCKET model
#'
#' @param data a matrix (samples, timepoints) or array (samples, channels, timepoints)
#' @param targets vector of regression dependent values
#' @param n number of kernels to generate
#' @param type "classification" or "regression"
#' @param separate_mv method to use for separating channels, keep at FALSE
#' @param features_to_use which features to use for the fit. If all, all features are retained. Otherwise, a character or character vector with feature types
#'
#' @return A rROCKET model object
#' @export
#'
fit <- function(data, targets, n, type = "classification", separate_mv = FALSE, features_to_use = "all") {
  if (length(class(data)) == 2){
    data <- array(data, dim = c(nrow(data), 1, ncol(data)))
  }
  y_train <- targets

  kernel_list <- NULL
  kernel_list <- generate_kernels_for_data(data[1, , ], n)
  rocket_features <- rocket_transform(data, kernel_list)
  features <- NULL
  if (features_to_use == "all"){
    features <- rocket_features[[1]]
    for (i in 2:length(rocket_features)){
      features <- cbind(features, rocket_features[[i]])
    }
  }
  else{
    features_temp <- rocket_features[features_to_use]
    features <- features_temp[[1]]
    for (i in 2:length(features_temp)){
      features <- cbind(features, features_temp[[i]])
    }
  }

  #implement logistic regression
  if (type == "regression"){
    cv_model <- cv.glmnet(features, y_train, nfolds = 5, alpha = 0.5)
  }
  else if (type == "classification"){
    #logistic regression
    cv_model <- NULL
    if (length(levels(factor(targets))) == 2){
      #binary targets, classify
      cv_model <- cv.glmnet(features, y_train, nfolds = 5, alpha = 0.5, family = "binomial")
    }
    else{
      cv_model <- cv.glmnet(features, y_train, nfolds = 5, alpha = 0.5, family = "multinomial")
    }
  }
  best_lambda <- cv_model$lambda.min
  coefs <- coef(cv_model, s = "lambda.min")
  print(coefs)
  print(colnames(features))
  #select kernels, but return all features.
  rel_features <- features[,as.matrix(coefs[-1] > 0)]
  print(class(cv_model))

  #fitted kernels
  model <- new("rROCKET",
               coefficients = as.vector(coefs),
               model = cv_model,
               kernels = kernel_list,
               features = features,
               rel_features = rel_features,
               call = match.call()
  )
  return(model)
}

#' Summary
#'
#' @param rROCKET the rROCKET model to summarize
#'
#' @return none
#' @export
#'
setMethod(
  f = "summary",
  signature = "rROCKET",
  definition = function(object) {
    cat("Call:\n")
    print(object@call)

    cat("\nCoefficients:\n")
    print(object@coefficients)

    cat("\nResiduals:\n")
    print(object@lambda)

  }
)

#' Predict rROCKET
#'
#' @param rROCKET the rROCKET model to run a prediction on.
#' @param data the data to predict from
#'
#' @return a list containing $predict, predicted values and $features, rocket features generated from data
#' @export
#'
setMethod(
  f = "predict",
  signature = "rROCKET",
  definition = function(object, data) {
    if (length(class(data)) == 2){
      data <- array(data, dim = c(nrow(data), 1, ncol(data)))
    }
    kernel_list <- object@kernels
    rocket_features <- rocket_transform(data, kernel_list)
    features <- rocket_features[[1]]
    for (i in 2:length(rocket_features)){
      features <- cbind(features, rocket_features[[i]])
    }
    predicted.values <- predict(object@model, newx = features,  s = "lambda.min")
    return(list(predict = predicted.values, feat = features))
  }
)

