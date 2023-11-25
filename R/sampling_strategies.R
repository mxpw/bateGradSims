#' Sampling groundtruth
#'
#' Not sampling per se, compute exact RS/MS (genetic - i.e. in the sens of 'from fertilized eggs')
#'
#' (Should be kept as an internal function ?)
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

#' Fixed number sampling
#'
#' Samples females eggs using a fixed number (So call 'Fixed Strategy', see details)
#'
#' (Should be kept as an internal function ?)
#'
#' @param fertilized_eggs List of size n_females containing father for each fertilized eggs (from e.g., pollen_competition() function) (no default)
#' @param n_males Number of males (no default)
#' @param total_samples Total number of eggs to be sample (default 1000)
#' @param by_female_samples Number of eggs to be sample by female (no default)
#' @param undercount_female Which strategy to use when some females doesn't have enough eggs to be sample ? (see details, default 'remove')
#'
#' @details Compute mating success (i.e. number of 'genetic' mates) and reproductive success (i.e. number of offspring) for males and females
#' using a 'fixed number' strategy (i.e. a fixed number of eggs by female is sampled). This number is determined using either total_samples or by_female_samples argument
#' (thus at least one should be provided if not using the default value. If both are provided, the by_female_samples will be used).
#'
#' Note that in all cases, the RS for females is computed using all eggs ! (should add another option here to allow this behavior to be changed ?)
#'
#' In some cases, some females might not have enough eggs to sample. Two strategies can be used to consider these cases (i) 'remove' strategy,
#' (ii) 'keep' strategy. The default 'remove' strategy doesn't consider females w/o enough eggs (will be NA in the output), the 'keep' strategy
#' consider these females by sampling all their eggs.
#'
#' @return List containing MS & RS for females and for males
#'
#' @export
#'

sampling_fixed = function( fertilized_eggs, n_males, total_samples = 1000, by_female_samples = NULL, undercount_female = 'remove' ){

  n_females = length(fertilized_eggs)

  if(is.null(by_female_samples))
    sample_by_female = floor( total_samples / n_females )
  else
    sample_by_female = by_female_samples

  if(sample_by_female < 1)
    stop("Samples by female should be >= 1")

  print(paste0(sample_by_female, " eggs will be sampled by female"))


  total_eggs = sum(unlist(lapply(fertilized_eggs, length)))

  # Keep track of females w/o enough eggs
  female_to_remove = rep(F, n_females)

  if( any( lapply(fertilized_eggs, length) < sample_by_female ) ){
    print(paste0("Some female(s) doesn't have enough fruits! Use strategy '",undercount_female,"'"))

    if(undercount_female == 'keep'){
      fertilized_eggs_reduced = lapply(fertilized_eggs,
                                   FUN = function(x) sample_handmade(x, size = sample_by_female) )

    }else if(undercount_female == 'remove'){
      fertilized_eggs_reduced = lapply(fertilized_eggs,
                                   FUN = function(x) {
                                     if(length(x) >= sample_by_female)
                                       sample_handmade(x, size = sample_by_female)
                                     else
                                       NULL
                                   })

      female_to_remove = unlist(lapply(fertilized_eggs, length)) < sample_by_female
      print(paste0("Female(s) removed from dataset : ",sum(female_to_remove)," idv(s) (",sum(female_to_remove)/length(female_to_remove),"%)"))

    }else{
      print(paste0("Strategy ", undercount_female, "doesn't exists !"))
    }
  }else{
    fertilized_eggs_reduced = lapply(fertilized_eggs,
                                 FUN = function(x) sample(x, size = sample_by_female, replace = FALSE) )
  }

  msg_female = unlist(lapply(fertilized_eggs_reduced, FUN = function(x) length(unique(x)) ))
  rsg_female = unlist(lapply(fertilized_eggs, length)) # Compute on non-reduced dataset

  msg_female[female_to_remove] = NA
  rsg_female[female_to_remove] = NA

  msg_male = rep(0, n_males)
  rsg_male = rep(0, n_males)
  for(m in 1:n_males){
    msg_male[m] = sum(unlist(lapply(fertilized_eggs_reduced, FUN = function(x) m %in% x )))
    rsg_male[m] = sum(unlist(fertilized_eggs_reduced) == m)
  }
  return(list(msg_female = msg_female,
              rsg_female = rsg_female,
              msg_male = msg_male,
              rsg_male = rsg_male))
}

#' Prorata sampling
#'
#' Samples females eggs proportionally to their total eggs count (So call 'Prorata Strategy', see details)
#'
#' (Should be kept as an internal function ?)
#'
#' @param fertilized_eggs List of size n_females containing father for each fertilized eggs (from e.g., pollen_competition() function) (no default)
#' @param n_males Number of males (no default)
#' @param total_samples Total number of eggs to be sample (default 1000)
#' @param by_female_prop Proportion of eggs to be sample by female (no default)
#' @param min_threshold Minimal threshold for eggs by females (see details, default : 0)
#' @param undercount_female Which strategy to use when some females doesn't have enough eggs to be sample ? (see details, default 'remove')
#' @param upsample_strategy Which strategy to use ........
#' @param upsampling_plot Should upsampling effects (i.e. deviations from exact proportions) be plotted ? (default FALSE)
#'
#' @details Compute mating success (i.e. number of 'genetic' mates) and reproductive success (i.e. number of offspring) for males and females
#' using a 'Prorata' strategy (i.e. sample eggs by female proportionally to their total eggs count).
#' This proportion is determined using either total_samples or directly by_female_prop argument
#' (thus at least one should be provided if not using the default value. If both are provided, the by_female_prop will be used).
#' When computing the number of eggs to be sampled, the value is rounded.
#'
#' Note that in all cases, the RS for females is computed using all eggs ! (should add another option here to allow this behavior to be changed ?)
#'
#' In some cases, some females might not have enough eggs to sample.
#' Two cases can be distinguished: (a) total eggs * proportion to be sampled < min_threshold and total eggs < min_threshold,
#' (b) total eggs * proportion to be sampled < min_threshold and total eggs >= min_threshold.
#'
#' Three strategies can be used in these cases : (i) 'keep' strategy, (ii) 'remove' strategy, and (iii) 'remove_and_upsample' strategy.
#' The 'keep' strategy actually does nothing (some females might then have a sampled number < min_threshold),
#' the 'remove' strategy doesn't consider females w/o enough eggs (i.e., in cases (a) and (b)), and will be NA in the output,
#' the 'remove_and_upsample' doesn't consider females in cases (a), but upsample females (to min_threshold) in (b) case.
#'
#' (Should add possibility to just upsample females to their total eggs without removing from others - meaning the total samples increases)
#' For now, two upsampling strategies (which keep total eggs to be sampled constant) are implemented.
#' (i) 's1' and (ii) 's2' strategies. 's1' strategy removes recursively one eggs from the female with the larger eggs count to allow
#' one more eggs to be sampled from females that are below min_threshold. The 's2' strategy operate differently; when n eggs have to be redistributed
#' then one eggs is removed from the n first females (if possible) with the larger eggs counts, if the number of females above the threshold
#' is less than n, then, the operation is done several times.
#'
#' Each strategy leads to a different re-repartion (either 'truncating' each side of the distribution - s1,
#'  or 'truncating' at left and homogeneously reduced sampled eggs at right - s2) - the redistribution can be plotted by setting upsampling_plot arg to TRUE.
#'
#' Note that upsampling might not always be possible (while keeping constant total), in such cases a error message will be prompted.
#'
#' Carefull when number of eggs by idv is low ... not sure about the behavior (rounding close from zero) - might need an update
#'
#' What to do when theo_perc * total eggs ~ e.g., 0.3 (e.g. an idv with 4 eggs, and theo_prec = 1\%) - take at least one, or not ?
#'
#'
#' @return List containing MS & RS for females and for males
#'
#' @importFrom graphics points
#' @export
#'

sampling_prorata = function( fertilized_eggs, n_males, total_samples = 1000, by_female_prop = NULL,
                           min_threshold = 0, undercount_female = 'remove', upsample_strategy = 's1', upsampling_plot = FALSE ){

  n_females = length(fertilized_eggs)

  if(is.null(by_female_prop)){
    theoretical_perc = total_samples / length(unlist(fertilized_eggs))
  }else{
    theoretical_perc = by_female_prop
  }

  print(paste0("Initial theoretical % to be sampled : ", round(100*theoretical_perc, 2), " %"))


  # Other check might be needed
  if(is.null(by_female_prop) & (min_threshold*n_females>total_samples))
    stop("Required minimal threshold is not compatible with total sample expected !!!")

  # Check eggs by female
  condition_one = unlist(lapply(fertilized_eggs, FUN = function(x) round(theoretical_perc * length(x)) < min_threshold ))
  condition_two = unlist(lapply(fertilized_eggs, FUN = function(x) length(x) < min_threshold ))

  if( any(condition_two) )
    print(paste0("Some female(s) (",sum(condition_two),", ",mean(100*condition_two)," %) doesn't have enough fruits at all."))

  if( any(condition_one) )
    print(paste0("Some female(s) (",sum(condition_one),", ",round(100*mean(condition_one),2)," %) will have less than expected fruits after sampling.",
                 " (",sum(condition_one & !condition_two),", ", round(100*mean(condition_one & !condition_two), 2)," %) of them can be upsampled."))

  print(paste0("Sampling with strategy '",undercount_female,"'"))

  female_to_remove = rep(F, n_females)

  # In all cases, will be update later
  fertilized_eggs_reduced = lapply(fertilized_eggs, FUN = function(x) sample_handmade(x, size = round(theoretical_perc * length(x))) )

  # Size before upscaling
  # print(paste("(debug) Size before upscaling ::: > ", length(unlist(fertilized_eggs_reduced))))

  if(undercount_female == 'keep'){
  }
  else if(undercount_female == 'remove'){
    female_to_remove = condition_one
    print(paste0("--> ",sum(condition_one)," (",round(100*mean(condition_one),2)," %) female(s) will be removed"))

    # if the total expected was provided, then new % have to be computed
    if(is.null(by_female_prop)){
      theoretical_perc = total_samples / length(unlist(fertilized_eggs[!female_to_remove]))
      print(paste0("Updated theoretical % to be sample : ", round(100*theoretical_perc, 2), " %"))
    }

    fertilized_eggs_reduced = lapply(fertilized_eggs, FUN = function(x) sample_handmade(x, size = round(theoretical_perc * length(x))) )

  }else if(undercount_female == 'remove_and_upsample'){

    female_to_remove = condition_two
    print(paste0("--> ",sum(condition_two)," (",round(100*mean(condition_two),2)," %) female(s) will be removed"))

    # if the total expected was provided, then new % have to be computed
    if(is.null(by_female_prop)){
      theoretical_perc = total_samples / length(unlist(fertilized_eggs[!female_to_remove]))
      print(paste0("Updated theoretical % to be sample : ", round(100*theoretical_perc, 2), " %"))
    }

    # Eggs to consider
    eggs = unlist(lapply(fertilized_eggs, FUN = function(x) round(theoretical_perc * length(x))))
    eggs[female_to_remove] = NA
    eggs = eggs - min_threshold

    if(upsampling_plot)
      plot(sort(eggs + min_threshold))

    if( sum(eggs, na.rm=T) < 0 ){
      stop("Upsampling is not possible !")
    }

    sorted_idx = sort(eggs, decreasing = T, na.last = T, index.return = T)$ix

    print(paste0("Upsampling with strategy ", upsample_strategy))

    if(upsample_strategy == 's1'){
      n_transfert = -sum(eggs[eggs < 0], na.rm = T)
      eggs[eggs < 0] = 0
      while(n_transfert > 0){
        eggs[sorted_idx][1] = eggs[sorted_idx][1] - 1
        sorted_idx = sort(eggs, decreasing = T, na.last = T, index.return = T)$ix
        n_transfert = n_transfert - 1
      }
    }else if(upsample_strategy == "s2"){
      # Method 2
      n_transfert = -sum(eggs[eggs < 0], na.rm = T)
      eggs[eggs < 0] = 0
      while(n_transfert > 0){
        aside = 0
        local_transfert = n_transfert
        if(n_transfert > length(eggs[!is.na(eggs)])){
          aside = n_transfert - length(eggs[!is.na(eggs)])
          local_transfert = length(eggs[!is.na(eggs)])
        }
        eggs[sorted_idx][1:local_transfert] = eggs[sorted_idx][1:local_transfert] - 1
        n_transfert = -sum(eggs[eggs < 0], na.rm = T) + aside
        eggs[eggs < 0] = 0
      }
    }else{
      print(paste0("Strategy ", upsample_strategy, " doesn't exists !"))
    }


    # After upsampling
    if(upsampling_plot)
      points(sort(eggs + min_threshold), col = 'red')

    # Recompute Female reduced
    eggs = eggs + min_threshold

    for(i in 1:length(fertilized_eggs)){
      fertilized_eggs_reduced[[i]] = sample_handmade(fertilized_eggs[[i]], size = eggs[i])
    }

  }else{
    print(paste0("Strategy ", undercount_female, " doesn't exists !"))
  }

  fertilized_eggs_reduced[ female_to_remove ] = list(NULL)

  print(paste0("Number of fruits to be genotyped = ", length(unlist(fertilized_eggs_reduced))))

  msg_female = unlist(lapply(fertilized_eggs_reduced, FUN = function(x) length(unique(x)) ))
  rsg_female = unlist(lapply(fertilized_eggs, length)) # Compute on non-reduced dataset !

  msg_female[female_to_remove] = NA
  rsg_female[female_to_remove] = NA

  msg_male = rep(0, n_males)
  rsg_male = rep(0, n_males)
  for(m in 1:n_males){
    msg_male[m] = sum(unlist(lapply(fertilized_eggs_reduced, FUN = function(x) m %in% x )))
    rsg_male[m] = sum(unlist(fertilized_eggs_reduced) == m)
  }
  return(list(msg_female = msg_female,
              rsg_female = rsg_female,
              msg_male = msg_male,
              rsg_male = rsg_male))
}

#' Random sampling
#'
#' Samples eggs randomly (i.e. independently of female identity, see details)
#'
#' (Should be kept as an internal function ?)
#'
#' @param fertilized_eggs List of size n_females containing father for each fertilized eggs (from e.g., pollen_competition() function) (no default)
#' @param n_males Number of males (no default)
#' @param total_samples Total number of eggs to be sample (default 1000)
#'
#' @details Compute mating success (i.e. number of 'genetic' mates) and reproductive success (i.e. number of offspring) for males and females
#' using a fully random strategy. A given number of eggs (total_samples) is randomly chosen, then parents identity is checked to compute MS/RS for males
#' et MS for females.
#'
#' Note that in all cases, the RS for females is computed using all eggs ! (should add another option here to allow this behavior to be changed ?)
#'
#' @return List containing MS & RS for females and for males
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble as_tibble rownames_to_column
#' @importFrom dplyr n_distinct n summarise group_by
#'
#' @export
#'

sampling_random = function( fertilized_eggs, n_males, total_samples = 1000 ){

  n_females = length(fertilized_eggs)

  # Make descendant centered results
  males_ = unlist(fertilized_eggs)
  females_ = c()
  for(f in 1:n_females)
    females_ = c(females_, rep(f, length(fertilized_eggs[[f]])) )

  spl = cbind(females_, males_)[sample(1:length(males_), size = total_samples),] %>% as_tibble()

  # Female RSg (from # of descendant)
  rsg_female_rd = unlist(lapply(fertilized_eggs, length)) # Compute on non-reduced dataset !

  # Female MSg
  female_gb = spl %>% group_by(.data$females_) %>% summarise(ms = n_distinct(.data$males_))

  msg_female_rd = rep(NA, n_females)
  for(i in 1:nrow(female_gb))
    msg_female_rd[ female_gb[[i, 1]] ] = female_gb[[i, 'ms']]

  # Male MSg/RSg
  male_gb = spl %>% group_by(.data$males_) %>% summarise(ms = n_distinct(.data$females_), rs = n())
  msg_male_rd = rep(NA, n_males)
  rsg_male_rd = rep(NA, n_males)
  for(i in 1:nrow(male_gb)){
    msg_male_rd[ male_gb[[i, 1]] ] = male_gb[[i, 'ms']]
    rsg_male_rd[ male_gb[[i, 1]] ] = male_gb[[i, 'rs']]
  }

  return(list(msg_female = msg_female_rd,
              rsg_female = rsg_female_rd,
              msg_male = msg_male_rd,
              rsg_male = rsg_male_rd))

}

#' Sampling
#'
#' Convenient wrapper around sampling strategies with convenient output format
#'
#'
#' @param fertilized_eggs List of size n_females containing father for each fertilized eggs (from e.g., pollen_competition() function) (no default)
#' @param n_males Number of males (no default)
#' @param methods List of lists containing methods to use and associated parameters (see examples)
#' @param scaled Should the MS/RS should be scaled ? (default TRUE)
#' @param mso MS obtained directly from observations (i.e. pollen repartition over females before pollen competition) (from ms_obs() function, not required)
#' @param gametes Number of gametes for each females/males (from gametes_drawing() function, not required)
#' @param n_rep Number of time each sampling method shoud be replicated (default 1)
#'
#' @details Run all required sampling methods with associated parameters n_rep time each. Results can be scaled by sex, methods & methods parameters.
#' User defined function can also be used here ; two constraints : (i) first two arguments must be 'fertilized_eggs' and 'n_males', others can be
#' anything defined by the user, (ii) the function must return list with four elements : msg_female, rsg_female, msg_male, rsg_male.
#'
#' @return MS (obs), MS (gen), RS (gen), gamete counts, sex, sampling_method, parameters used for the sampling method, and replicate IDs
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom dplyr group_by ungroup mutate across
#'
#' @examples
#' \dontrun{
#' methods = list(base = list(method = "sampling_groundtruth",
#'                            params = list()),
#'                fixed = list(method = "sampling_fixed",
#'                             params = list(by_female_samples = 50, undercount_female = 'keep')),
#'                prorata = list(method = "sampling_prorata",
#'                               params = list(by_female_prop = 0.10, min_threshold = 0,
#'                                             undercount_female = 'keep', upsample_strategy = 's2')),
#'                random = list(method = "sampling_random",
#'                              params = list()))
#'
#' sampling(fertilized_eggs, n_males, methods = methods, mso = mso, gametes = gametes, scaled = T)
#' }
#' @export
#'
#'
sampling = function(fertilized_eggs, n_males, methods = NULL, scaled = TRUE, mso = NULL, gametes = NULL, n_rep = 1){

  n_females = length(fertilized_eggs)

  output = tibble(mso = numeric(), msg = numeric(), rsg = numeric(), n_gam = numeric(),
                  sex = character(), sampling_method = character(), parameters = list(), parameters_string = character())

  if(!is.null(mso)){
    mso_vec = c(mso$mso_female, mso$mso_male)
  }else{
    mso_vec = NA
  }

  if(!is.null(gametes)){
    gam_vec = c(gametes$gam_female, gametes$gam_male)
  }else{
    gam_vec = NA
  }

  for(m in 1:length(methods)){
    for(r in 1:n_rep){
      spl = do.call( methods[[m]]$method, c(list(fertilized_eggs, n_males), methods[[m]]$params) )
      tmp = tibble(mso = mso_vec,
                   msg = c(spl$msg_female, spl$msg_male),
                   rsg = c(spl$rsg_female, spl$rsg_male),
                   n_gam = gam_vec,
                   sex = c(rep("F", n_females), rep("M", n_males)),
                   sampling_method = names(methods)[m],
                   parameters = list(methods[[m]]$params),
                   parameters_string = paste(names(methods[[m]]$params),methods[[m]]$params,sep="=",collapse=";" ),
                   replicate = r)
      output = rbind(output, tmp)
    }
  }


  print("Scaling is weird !! Should it not be true scaling instead ?! here, when zero => stay zero")
  if(scaled)
    return(output %>% group_by(.data$sex, .data$sampling_method, .data$parameters_string, .data$replicate) %>%
             mutate(across(.data$mso:.data$rsg, ~ .x / mean(.x, na.rm = T))) %>%
             ungroup())

  return(output)
}
