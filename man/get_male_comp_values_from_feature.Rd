% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_components.R
\name{get_male_comp_values_from_feature}
\alias{get_male_comp_values_from_feature}
\title{Draw male competitive value relative to another feature (e.g. pollen set size)}
\usage{
get_male_comp_values_from_feature(
  mean_comp_value,
  sd_comp_value,
  rho,
  feature,
  translation = TRUE,
  plots = FALSE
)
}
\arguments{
\item{mean_comp_value}{Mean competitive value for males}

\item{sd_comp_value}{Standard deviation of competitive values for males}

\item{rho}{Correlation between competitive values and provided feature}

\item{feature}{Vector of size n_males with}

\item{translation}{Should the resulting comp. value be translated to be > 0 ? (default TRUE)}

\item{plots}{Should the comp. value histogram and link between feature and comp. values be plotted ? (default FALSE)}
}
\value{
Vector of competitive values
}
\description{
Produce a vector of competitive values from conditional normal distribution
}
\details{
Draw males competitive values from conditional distribution of a bivariate normal distribution. (to detailed; could probably be easily extended to multivariate cases)

After drawing, values are transformed such that all values are strictly above zero (i.e. x' = x + min(x) + 0.01).
This is needed because these value are used later as probabilities. (that default behavior can be modified with the translation arg.)
Note however this may change mean and variance of specified distribution !
}
