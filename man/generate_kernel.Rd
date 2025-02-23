% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rocket.R
\name{generate_kernel}
\alias{generate_kernel}
\title{Generate Kernel
Generates a single kernel}
\usage{
generate_kernel(
  dims,
  dilation,
  padding,
  len,
  bias,
  name,
  weight_values = c(1, 0, -1),
  seed = NA
)
}
\arguments{
\item{dims}{dimensions for weights (should be a two-element vector)}

\item{dilation}{dilation value, should be generated from data}

\item{padding}{padding value to pre- and post-pend to ts samples}

\item{len}{kernel length, integer}

\item{bias}{bias value, numeric}

\item{name}{kernel name}

\item{weight_values}{weight values from which to sample. "unif" and "norm" are options, as are a vector of numeric weights.}

\item{seed}{random seed for reproducibility}
}
\value{
A ROCKET kernel
}
\description{
Generate Kernel
Generates a single kernel
}
