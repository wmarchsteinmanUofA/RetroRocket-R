#need to pass in a list of kernels, list of ts
#data should be a 3d array
#' Rocket_transform
#' A custom implementation of a ROCKET transform
#'
#' @param data  a 3D array of data, of the form (samples, channels, timepoints)
#' @param kernels a list of kernels generated from generate_kernels_from_data
#' @param ... other arguments to pass through
#'
#' @return a list of transformed data, one matrix per feature type.
#' @export
#'

rocket_transform <- function(data, kernels, ...){
  transform_list <- list()
  #add more features here.
  ns <- unlist(kernels)
  ker_names <- ns[names(ns) == "name"]

  transform_list$max <- matrix(rep(0, dim(data)[1]*length(kernels)), nrow = dim(data)[1],
                               dimnames = list(paste0("ts",1:dim(data)[1]),paste0(ker_names, "_max")))
  transform_list$av_max <- matrix(rep(0, dim(data)[1]*length(kernels)), nrow = dim(data)[1],
                               dimnames = list(paste0("ts",1:dim(data)[1]),paste0(ker_names, "_avmax")))
  transform_list$ppv <- matrix(rep(0, dim(data)[1]*length(kernels)), nrow = dim(data)[1],
                               dimnames = list(paste0("ts",1:dim(data)[1]),paste0(ker_names, "_ppv")))
  transform_list$max_ind <- matrix(rep(0, dim(data)[1]*length(kernels)), nrow = dim(data)[1],
                                   dimnames = list(paste0("ts",1:dim(data)[1]),paste0(ker_names, "_max_ind")))
  transform_list$max_row <- matrix(rep(0, dim(data)[1]*length(kernels)), nrow = dim(data)[1],
                                   dimnames = list(paste0("ts",1:dim(data)[1]),paste0(ker_names, "_max_row")))
  transform_list$lspv <- matrix(rep(0, dim(data)[1]*length(kernels)), nrow = dim(data)[1],
                                dimnames = list(paste0("ts",1:dim(data)[1]),paste0(ker_names, "_lspv")))
  transform_list$lspv_ind <- matrix(rep(0, dim(data)[1]*length(kernels)), nrow = dim(data)[1],
                                    dimnames = list(paste0("ts",1:dim(data)[1]),paste0(ker_names, "_lspv_ind")))
  transform_list$lspv_row <- matrix(rep(0, dim(data)[1]*length(kernels)), nrow = dim(data)[1],
                                    dimnames = list(paste0("ts",1:dim(data)[1]),paste0(ker_names, "_lspv_row")))
  transform_list$mean_ipv <- matrix(rep(0, dim(data)[1]*length(kernels)), nrow = dim(data)[1],
                                    dimnames = list(paste0("ts",1:dim(data)[1]),paste0(ker_names, "_mean_ipv")))
  transform_list$mean_pv <- matrix(rep(0, dim(data)[1]*length(kernels)), nrow = dim(data)[1],
                                   dimnames = list(paste0("ts",1:dim(data)[1]),paste0(ker_names, "_mean_pv")))

  for (i in 1:dim(data)[1]){
    for (j in 1:length(kernels)){
      new_ts <- pad_ts(data[i, , ], kernels[[j]]$padding)
      ker <- create_dilation(kernels[[j]])
      num_dots <- (dim(new_ts)[2] - dim(ker)[2])
      transform_vals <- matrix(rep(0, num_dots * dim(ker)[1]), ncol = num_dots)
      for (k in 1:num_dots){
        transform_vals[,k] <- rowSums(ker * new_ts[,k:k + dim(ker)[2]]) + kernels[[j]]$bias
      }
      v <- ramify::flatten(transform_vals)
      transform_list$max[i, j] <- max(transform_vals)
      wher <- which(v == max(transform_vals))
      if (length(wher) > 1){
        wher <- wher[1]
      }
      if (length(wher) == 1){
        tem <- wher%%ncol(transform_vals)
        if (tem == 0) tem <- 12
        transform_list$max_ind[i, j] <- tem
        transform_list$max_row[i, j] <- floor(wher/ncol(transform_vals) + 1)
      }
      transform_list$av_max[i,j] <- mean(apply(transform_vals, 1, max, na.rm=TRUE))
      m <- mean(v[v>0])
      transform_list$mean_pv[i, j] <- ifelse(is.nan(m),0, m)
      transform_list$ppv[i, j] <- sum(v > 0)/length(v)
      lspv_t <- calculate_lspv_mat(transform_vals)
      lspv<-max(lspv_t[,1])
      lspv_row <- which(lspv_t[,1] == lspv)
      if (length(lspv_row) > 1){
        lspv_row <- lspv_row[1]
      }
      if (length(lspv_row) == 1){
        tem <- wher%%ncol(transform_vals)
        transform_list$lspv_row[i, j] <- lspv_row
        transform_list$lspv_ind[i, j] <- lspv_t[lspv_row, 2]
      }
      transform_list$lspv[i, j] <- lspv

    }
  }
  return(transform_list)
}
calculate_lspv_mat <- function(transform){
  ans <- matrix(0, nrow = nrow(transform), ncol = 2)
  for (i in 1:nrow(transform)){
    ans[i,] <- calculate_lspv(transform[i,])
  }
  return(ans)
}
calculate_lspv <- function(vec){
  tf <- vec > 0
  flg <- FALSE
  ls <- 0
  ls_temp <- 0
  loc <- -1
  for (i in 1:length(tf)){
    if (flg == TRUE){
      if (tf[i]){
        ls_temp <- ls_temp + 1
        if (ls_temp > ls){
          ls <- ls_temp
          loc <- i - ls_temp + 1
        }
      }
      else{
        flg <- FALSE
        if (ls_temp > ls){
          ls <- ls_temp
          loc <- i - ls_temp
        }
        ls_temp <- 0
      }
    }
    else{
      if (tf[i]){
        ls_temp <- 1
        flg <- TRUE
      }
    }
  }
  return(c(ls, loc))
}

pad_ts <- function(ts, padding){
  tr <- as.matrix(ts)
  if (dim(tr)[2] == 1) tr <- t(tr)
  new_t <- matrix(rep(0,dim(tr)[1]*(2*padding + dim(tr)[2])), nrow = dim(tr)[1])
  for (row in 1:dim(tr)[1]){
    new_t[row,] <- c(rep(0, padding), tr[row, ], rep(0, padding))
  }
  return(new_t)
}

create_dilation <- function(kernel){
  weights <- kernel$weights
  dilation <- kernel$dilation
  new_mat <- list()
  for (row in 1:dim(weights)[1]){
    new_mat[[row]] <- insert_zeros(weights[row,], dilation-1)
  }
  return(t(matrix(unlist(new_mat), ncol = dim(kernel$weights)[1])))
}

insert_zeros <- function(vec, num_zeros) {
  # Create a vector of zeros with the specified length
  zeros <- rep(0, num_zeros)
  # Create an empty vector to store the result
  result <- numeric(0)

  # Loop through the input vector and insert zeros between elements
  for (i in seq_along(vec)) {
    result <- c(result, vec[i])
    if (i < length(vec)) {
      result <- c(result, zeros)
    }
  }

  return(result)
}


generate_kernel <- function(dims, dilation, padding, len, bias, name, weight_values = c(1, 0, -1), seed = NA){
  if (!is.na(seed)){
    set.seed(seed)
  }
  ker <- list()
  ker$dilation <- dilation
  ker$padding <- padding
  ker$length <- len
  ker$weights <- NULL
  ker$name <- name
  ker$bias <- bias
  if (typeof(weight_values) == "character" && weight_values == "unif"){
    ker$weights <- matrix(runif(dims[1]*len, min = -1, max = 1), nrow = dims[1])
  }
  else if (typeof(weight_values) == "character" && weight_values == "norm"){
    ker$weights <- matrix(rnorm(dims[1]*len, min = -1, max = 1), nrow = dims[1])
  }
  else if (class(weight_values) == "numeric"){
    ker$weights <- matrix(sample(weight_values, dims[1]*len, replace = TRUE), nrow = dims[1])
  }
  return(ker)
}


generate_kernels_for_data <- function(data_entry, n, lengths = c(7, 9, 11, 13, 15), seed = "none", tag = ""){
  kernel_list <- list(NULL)
  dims <- dim(data_entry)
  if (is.null(dims)){
    dims <- c(1, length(data_entry))
  }
  if (seed != "none") {
    set.seed(seed)
  }
  for (i in 1:n){
    len <- sample(lengths, 1)
    if (len > dims[2]) len <- floor(dims[2]/2)
    dilation <- floor(2^(runif(1, min = 0, max = 1) * log2( (dims[2]-1) / (len-1))))
    padding <- 0
    if (runif(1) < 0.5){
      padding= floor( (len-1) * dilation/2)
    }

    ker <- generate_kernel(dims, dilation, padding, len, runif(1, min = -1, max = 1), paste0("ker", i, tag),weight_values = "norm")
    kernel_list <- append(kernel_list, list(ker))
  }
  return(kernel_list[-1])
}

