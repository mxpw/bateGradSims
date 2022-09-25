#' Sampling groundtruth
#'
#' Not sampling per se, compute exact RS/MS (genetic - i.e. in the sens of 'from fertilized eggs')
#'
#' @param fertilized_eggs List of size n_females containing father for each fertilized eggs (from e.g., pollen_competition() function) (no default)
#' @param n_males Number of males (no default)
#'
#' @details Compute exact mating success (i.e. number of 'genetic' mates) and reproductive success (i.e. number of offspring) for males and females
#'
#' @return List containing MS & RS for females and for males
#'
#' @export
#'

sampling_groundtruth = function( fertilized_eggs, n_males){
  msg_female = unlist(lapply(fertilized_eggs, FUN = function(x) length(unique(x)) ))
  rsg_female = unlist(lapply(fertilized_eggs, length))

  msg_male = rep(0, n_males)
  rsg_male = rep(0, n_males)

  for(m in 1:n_males){
    msg_male[m] = sum(unlist(lapply(fertilized_eggs, FUN = function(x) m %in% x )))
    rsg_male[m] = sum(unlist(fertilized_eggs) == m)
  }

  return(list(msg_female = msg_female,
              rsg_female = rsg_female,
              msg_male = msg_male,
              rsg_male = rsg_male))
}
