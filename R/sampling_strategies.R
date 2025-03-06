#' Mating and reproductive success from sampled seeds
#'
#' Compute RS/MS (genetic - i.e. in the sens of 'from fertilized eggs')
#'
#'
#' @param fertilized_eggs List of size n_females containing father for each fertilized eggs (from e.g., pollen_competition() function)
#' @param sampled_fertilized_eggs List of size n_females containing father for each fertilized eggs after samping (from any sampling functions)
#' @param n_males Number of males (no default)
#' @param paternity_share Should males RS be rescaled by total number of females eggs (default TRUE)
#'
#' @details Compute mating success and reproductive success (i.e. number of offspring) for males and females.
#' Paternity share is considered by default (see paternity_share argument), and females RS is always compute from all eggs (not sampled ones).
#'
#' @return List containing MSg & RS for females and males
#'
#' @export
#'

get_sexual_selection_components = function( fertilized_eggs, sampled_fertilized_eggs, n_males, paternity_share = TRUE){
  msg_female = unlist(lapply(sampled_fertilized_eggs, FUN = function(x) length(unique(x)) ))
  rsg_female = unlist(lapply(fertilized_eggs, length))

  msg_male = rep(0, n_males)
  rsg_male = rep(0, n_males)

  if(paternity_share){
    print("Males RS is rescaled by total eggs count by female")
    prop_sampled_eggs = unlist(lapply(sampled_fertilized_eggs, length)) / unlist(lapply(fertilized_eggs, length))
    weights = 1 / prop_sampled_eggs
    weights[is.infinite(weights) | is.na(weights)] = 0
  }

  for(m in 1:n_males){
    msg_male[m] = sum(unlist(lapply(sampled_fertilized_eggs, FUN = function(x) m %in% x )))
    rsg_male_by_female = unlist(lapply(sampled_fertilized_eggs, function(x) sum(x==m) ))

    if(paternity_share){
      rsg_male_by_female = rsg_male_by_female * weights
    }
    rsg_male[m] = sum(rsg_male_by_female)
  }

  return(list(msg_female = msg_female,
              rsg_female = rsg_female,
              msg_male = msg_male,
              rsg_male = rsg_male))
}

#' Compute Qfoc
#'
#' Compute Qfoc needed for MSgc computation
#
#' @details Product of couple production
#'
#' @return Qfoc
#'
# Il ne faut pas mettre de @export, on conserve comme fonction interne (les gens qui utiilisent le paquet ne peuvent pas l'utiliser, et ne savent pas qu'elle existe.)
#'
#'
#'

# Function to calculate Qfoc for a given group (mothers or fathers)
compute_Qfoc <- function(focal_list) {
  Qfoc <- numeric(length(focal_list))

  for (focal_it in seq_along(focal_list)) {
    pairs <- table(focal_list[[focal_it]])  # Count occurrences of each father for this mother
    prod_couple <- sum(pairs * (pairs - 1L))  # Product of pairs (count * (count - 1))
    prod_foc <- length(focal_list[[focal_it]])
    Qfoc[focal_it] <- prod_couple/(prod_foc * (prod_foc - 1L))
  }

  return(Qfoc)
}

#' Compute Qmin or Qmax
#'
#' Compute Qmin and Qmax needed for MSgc computation
#
#' @details Qmin retransmits a distribution of seeds that is as uniform as
#' possible among potential mates, while Qmax retransmits a distribution
#' of seeds that is as heterogeneous as possible, always starting with the
#' potential mate who has produced the most seeds in total.
#'
#' @return Qmin or Qmax
#'
compute_Q <- function(focal_list, method = "min") {
  Q <- numeric(length(focal_list))  # Initialize result vector

  # Step 1: Calculate the total seeds produced by each parent
  mate_total_counts <- table(unlist(focal_list))
  mate_total_counts_df <- setNames(as.data.frame(mate_total_counts), c("mate", "total"))
  rownames(mate_total_counts_df) <- mate_total_counts_df$mate

  for (focal_it in seq_along(focal_list)) {
    # Total production for the focal individual
    prod_foc <- length(focal_list[[focal_it]])

    # Step 2: Count the seeds produced specifically with the focal individual
    mate_counts <- table(focal_list[[focal_it]])
    mate_df <- setNames(as.data.frame(mate_counts), c("mate", "count"))

    # Step 3: Add parents absent from the relationship but present in total_counts_df
    mate_df <- merge(mate_total_counts_df, mate_df, by = "mate", all.x = TRUE)
    mate_df$count[is.na(mate_df$count)] <- 0L  # Replace NA with 0
    rownames(mate_df) <- mate_df$mate

    mate_df$allocated_seeds <- integer(nrow(mate_df))
    remaining_seeds_on_focal <- sum(mate_df$count)

    # Step 4: Distribution based on the chosen method
    if (method == "min") {
      # Most equitable distribution possible
      mate_df <- mate_df[order(mate_df$total), ]
      nmates <- nrow(mate_df)
      allocated_seeds <- integer(nmates)
      remaining_seeds_per_mate <- mate_df$total

      while (remaining_seeds_on_focal > 0L) {
        pos_in_mate_df <- which.max(remaining_seeds_per_mate > 0L)
        max_seeds <- remaining_seeds_per_mate[pos_in_mate_df]

        if (max_seeds > 0L) {
          mate_range <- pos_in_mate_df - 1L + seq(min(nmates + 1L - pos_in_mate_df, remaining_seeds_on_focal))
          allocated_seeds[mate_range] <- allocated_seeds[mate_range] + 1L
          remaining_seeds_per_mate[mate_range] <- remaining_seeds_per_mate[mate_range] - 1L
          remaining_seeds_on_focal <- remaining_seeds_on_focal - length(mate_range)
        } else {
          stop("Cannot allocate all seeds")
        }
      }
    } else if (method == "max") {
      # Most imbalanced distribution possible
      mate_df <- mate_df[order(-mate_df$total), ]
      allocated_seeds <- integer(nrow(mate_df))
      pos_in_mate_df <- 1L

      while (remaining_seeds_on_focal > 0L) {
        max_seeds <- min(remaining_seeds_on_focal, mate_total_counts_df[mate_df$mate[pos_in_mate_df], "total"])

        if (max_seeds > 0L) {
          allocated_seeds[pos_in_mate_df] <- allocated_seeds[pos_in_mate_df] + max_seeds
          remaining_seeds_on_focal <- remaining_seeds_on_focal - max_seeds
        }
        pos_in_mate_df <- pos_in_mate_df + 1L
      }
    } else {
      stop("Invalid method. Use 'min' or 'max'.")
    }

    # Step 6: Calculate Q
    Q[focal_it] <- sum(allocated_seeds * (allocated_seeds - 1)) / (prod_foc * (prod_foc - 1L))
  }

  return(Q)
}


#' Corrected genetic mating success
#'
#' Compute corrected genetic mating success MSgc using F. Rousset's metric
#' For the moment, this function is only applicable to all the seeds produced (i.e. sampling_groundtruth function, no sampling).
#'
#'
#' @param fertilized_eggs List of size n_females containing father for each fertilized eggs (from e.g., pollen_competition() function)
#'
#' @details Compute genetic mating success for males and females with a correction to avoid a possible artificial correlation between RS/MSg.
#'
#' @return List containing corrected MSg (MSgc) for females and males
#'
#' @export
#'

compute_msgc <- function(fertilized_eggs,n_males) {

  # List of mothers and fathers
  mothers <- seq_along(fertilized_eggs)  # IDs of mothers
  fathers <- sort(unique(unlist(fertilized_eggs)))  # IDs of unique fathers

  # Calculate Qfoc for females and males
  Qfoc_female <- compute_Qfoc(fertilized_eggs)

  # Reverse the list to obtain the fathers
  list_fathers <- split(rep(mothers, sapply(fertilized_eggs, length)), unlist(fertilized_eggs))
  Qfoc_male <- compute_Qfoc(list_fathers)

  # Calculate Qmin and Qmax for females
  Qmin_female <- compute_Q(fertilized_eggs, method = "min")
  Qmax_female <- compute_Q(fertilized_eggs, method = "max")

  # Calculate Qmin and Qmax for males
  Qmin_male <- compute_Q(list_fathers, method = "min")
  Qmax_male <- compute_Q(list_fathers, method = "max")

  # Calculate msgc_female and msgc_male
  msgc_female <- 1 - ((Qfoc_female - Qmin_female) / (Qmax_female - Qmin_female))
  msgc_male <- 1 - ((Qfoc_male - Qmin_male) / (Qmax_male - Qmin_male))

  # Add males that did not reproduce
  complete_male_df <- data.frame(fathers = 1:n_males) %>%
    left_join(data.frame(fathers,msgc_male), by= "fathers")

  # Return the results
  return(list(msgc_female = msgc_female,
              msgc_male = complete_male_df$msgc_male))
}

#' Sampling ground-truth
#'
#' Dummy/Convenience function - just return the input
#'
#' @param fertilized_eggs List of size n_females containing father for each fertilized eggs (from e.g., pollen_competition() function) (no default)
#' @param ... added for convenience
#'
#' @details Return the input list - use for generality purpose, same pipeline can be use than with other sampling methods.
#'
#' @return Sampled fertilized eggs
#'
#' @export
#'
sampling_groundtruth=function(fertilized_eggs, ...){
  print("=== Sampling : Groundtruth ===")
  fertilized_eggs
}


#' Fixed number sampling
#'
#' Samples females eggs using a fixed number (So call 'Fixed Strategy', see details)
#'
#' @param fertilized_eggs List of size n_females containing father for each fertilized eggs (from e.g., pollen_competition() function) (no default)
#' @param n_males Number of males (no default)
#' @param total_samples Total number of eggs to be sample (default 1000)
#' @param by_female_samples Number of eggs to be sample by female (no default)
#' @param undercount_female Which strategy to use when some females doesn't have enough eggs to be sample ? (see details, default 'remove')
#'
#' @details Sample eggs using a 'fixed number' strategy (i.e. a fixed number of eggs by female is sampled). This number is determined using either
#' total_samples or by_female_samples argument (thus at least one should be provided if not using the default value. If both are provided,
#' the by_female_samples will be used).
#'
#' In some cases, some females might not have enough eggs to sample. Two strategies can be used to consider these cases (i) 'remove' strategy,
#' (ii) 'keep' strategy. The default 'remove' strategy doesn't consider females w/o enough eggs (will be NA in the output), the 'keep' strategy
#' consider these females by sampling all their eggs.
#'
#' In all cases, if target number of samples is not reached because some females haven't enough eggs, additional eggs are selected from others females (in a
#' uniform and random way).
#'
#' @return Sampled fertilized eggs
#'
#' @export
#'

sampling_fixed = function( fertilized_eggs, n_males, total_samples = 1000,
                           by_female_samples = NULL, undercount_female = 'remove'){

  print("=== Sampling : Fixed ===")
  n_females = length(fertilized_eggs)

  if(is.null(by_female_samples)){
    sample_by_female = floor( total_samples / n_females )
  }else{
    sample_by_female = by_female_samples
  }

  if(sample_by_female < 1)
    stop("Samples by female should be >= 1")

  total_to_sample = sample_by_female * n_females

  print(paste0(sample_by_female, " eggs will (theoretically) be sampled by female (for a total of ",total_to_sample," eggs)"))

  total_eggs = sum(unlist(lapply(fertilized_eggs, length)))

  if(total_eggs < total_to_sample)
    stop(paste0("Target sample (",total_to_sample,") is not feasible (total eggs is ", total_eggs ,")"))

  # Keep track of females w/o enough eggs
  female_to_remove = rep(F, n_females)

  if( any( lapply(fertilized_eggs, length) < sample_by_female ) ){
    print(paste0("Some female(s) don't have enough eggs ! Use strategy '",undercount_female,"'"))

    if(undercount_female == 'keep'){

      undercounts = sapply(fertilized_eggs, length) - sample_by_female
      n_undercounts = sum(undercounts < 0)
      additional_eggs = -sum(undercounts[undercounts < 0])
      excess_sample_by_female = ceiling( additional_eggs / (n_females - n_undercounts) )

      updated_sample_by_female = sample_by_female + excess_sample_by_female

      fertilized_eggs_reduced = lapply(fertilized_eggs,
                                   FUN = function(x) sample_handmade(x, size = updated_sample_by_female) )

    }else if(undercount_female == 'remove'){

      female_to_remove = unlist(lapply(fertilized_eggs, length)) < sample_by_female
      print(paste0("Female(s) removed from dataset : ",sum(female_to_remove)," idv(s) (",100*sum(female_to_remove)/length(female_to_remove),"%)"))

      n_undercounts = sum(female_to_remove)
      additional_eggs = n_undercounts * sample_by_female
      excess_sample_by_female = ceiling( additional_eggs / (n_females - n_undercounts) )
      updated_sample_by_female = sample_by_female + excess_sample_by_female

      total_eggs = sum(unlist(lapply(fertilized_eggs[!female_to_remove], length)))

      if(total_eggs < total_to_sample)
        stop(paste0("Target sample (",total_to_sample,") is not feasible (total eggs after removing some females is ", total_eggs ,")"))

      fertilized_eggs_reduced = lapply(fertilized_eggs,
                                   FUN = function(x) {
                                     if(length(x) >= sample_by_female)
                                       sample_handmade(x, size = updated_sample_by_female)
                                     else
                                       NULL
                                   })

    }else{
      stop(paste0("Strategy ", undercount_female, "doesn't exists !"))
    }

    # Clear excess samples
    by_females_sampled_eggs = sapply(fertilized_eggs_reduced, length)
    while( sum(by_females_sampled_eggs) > total_to_sample ){
      max_ = max(by_females_sampled_eggs)
      random_idv = sample(which(by_females_sampled_eggs == max_), 1)
      fertilized_eggs_reduced[[random_idv]] = sample(fertilized_eggs_reduced[[random_idv]], size = max_ - 1)
      by_females_sampled_eggs = sapply(fertilized_eggs_reduced, length)
    }

  }else{
    fertilized_eggs_reduced = lapply(fertilized_eggs,
                                 FUN = function(x) sample(x, size = sample_by_female, replace = FALSE) )
  }

 fertilized_eggs_reduced
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
#' @details Sample eggs using a 'Prorata' strategy (i.e. sample eggs by female proportionally to their total eggs count).
#' This proportion is determined using either total_samples or directly by_female_prop argument
#' (thus at least one should be provided if not using the default value. If both are provided, the by_female_prop will be used).
#' When computing the number of eggs to be sampled, the value is rounded.
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
#' @return Sampled fertilized eggs
#'
#' @importFrom graphics points
#' @export
#'

sampling_prorata = function(fertilized_eggs, n_males, total_samples = 1000, by_female_prop = NULL,
                           min_threshold = 0, undercount_female = 'remove', upsample_strategy = 's1',
                           upsampling_plot = FALSE){

  print("=== Sampling : Prorata ===")

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

  if(undercount_female == 'keep'){}
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
      ordering_ = sort(eggs + min_threshold, index.return = T, na.last = T)$ix
      y = (eggs + min_threshold)[ordering_]
      removed_females = is.na(y)
      y[ removed_females ] = 0
      g1 = ggplot()+
        geom_point(aes(x = 1:length(ordering_),
                       y = y, shape = removed_females), show.legend = F)+
        scale_shape_manual(values=c(19, 4))+
        geom_hline(yintercept = min_threshold, linetype = "dashed")+
        theme_bw()+
        xlab("Females ID")+
        ylab("Eggs to sample")+
        expand_limits(y=0)

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
      g1 = g1 + geom_point(aes(x = 1:length(eggs), y = (eggs + min_threshold)[ordering_]), color = "red", inherit.aes = F)
      print(g1)

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

  fertilized_eggs_reduced
}

#' Random sampling
#'
#' Samples eggs randomly (i.e. independently of female identity, see details)
#'
#'
#' @param fertilized_eggs List of size n_females containing father for each fertilized eggs (from e.g., pollen_competition() function) (no default)
#' @param n_males Number of males (no default)
#' @param total_samples Total number of eggs to be sample (default 1000)
#'
#' @details Sample eggs using a fully random strategy. A given number of eggs (total_samples) is randomly chosen
#'
#' @return Sampled fertilized eggs
#'
#' @export
#'

sampling_random = function( fertilized_eggs, n_males, total_samples = 1000){

  print("=== Sampling : Random ===")

  n_females = length(fertilized_eggs)
  weights_ = sapply(fertilized_eggs, length)

  if(total_samples > sum(weights_))
    stop(paste0("Target sample (",total_samples,") is not feasible (total eggs is ", sum(weights_) ,")"))

  fertilized_eggs_reduced = vector("list", n_females)

  for(drawn in 1:total_samples){
    random_female = sample(1:n_females, size = 1, prob = weights_)
    random_idx_from_female = sample(1:length(fertilized_eggs[[random_female]]), 1)
    fertilized_eggs_reduced[[random_female]] = c(fertilized_eggs_reduced[[random_female]],
                                                 fertilized_eggs[[random_female]][random_idx_from_female])
    fertilized_eggs[[random_female]] = fertilized_eggs[[random_female]][-random_idx_from_female]
    weights_[random_female] = weights_[random_female] - 1
  }

  fertilized_eggs_reduced

}

#' Sampling
#'
#' Convenient wrapper around sampling strategies with convenient output format
#'
#'
#' @param fertilized_eggs List of size n_females containing father for each fertilized eggs (from e.g., pollen_competition() function) (no default)
#' @param n_males Number of males (no default)
#' @param methods List of lists containing methods to use and associated parameters (see examples)
#' @param mso MS obtained directly from observations (i.e. pollen repartition over females before pollen competition) (from ms_obs() function, not required)
#' @param gametes Number of gametes for each females/males (from gametes_drawing() function, not required)
#' @param n_rep Number of time each sampling method shoud be replicated (default 1)
#'
#' @details Run all required sampling methods with associated parameters n_rep time each.
#'
#' User defined function can also be used here ; two constraints : (i) first two arguments must be 'fertilized_eggs' and 'n_males' (can be made optional, see the internal sampling_groundtruth() for an example), others can be
#' anything defined by the user, (ii) the function must return list with four elements : msg_female, rsg_female, msg_male, rsg_male.
#'
#' @return Dataframe with MS (obs), MS (gen), RS (gen), gamete counts, sex, sampling_method, parameters used for the sampling method, and replicate IDs
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom dplyr group_by ungroup mutate across
#' @importFrom rlang .data
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
#' sampling(fertilized_eggs, n_males, methods = methods, mso = mso, gametes = gametes)
#' }
#' @export
#'
#'
sampling = function(fertilized_eggs, n_males, methods = NULL, mso = NULL, gametes = NULL, n_rep = 1){

  n_females = length(fertilized_eggs)

  output = tibble(mso = numeric(), msg = numeric(), msgc = numeric(), rsg = numeric(), n_gam = numeric(),
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

  # Get true MS/RS (no sampling)
  sampled_fertilized_eggs = sampling_groundtruth(fertilized_eggs)
  spl = get_sexual_selection_components(fertilized_eggs, sampled_fertilized_eggs, n_males)
  msgc_result = compute_msgc(fertilized_eggs)
  tmp = tibble(mso = mso_vec,
               msg = c(spl$msg_female, spl$msg_male),
               msgc = c(msgc_result$msgc_female, msgc_result$msgc_male),
               rsg = c(spl$rsg_female, spl$rsg_male),
               n_gam = gam_vec,
               sex = c(rep("F", n_females), rep("M", n_males)),
               sampling_method = "base",
               parameters = list(list()),
               parameters_string = "",
               replicate = 1)
  output = rbind(output, tmp)

  for(m in 1:length(methods)){
    for(r in 1:n_rep){
      sampled_fertilized_eggs = do.call( methods[[m]]$method, c(list(fertilized_eggs, n_males), methods[[m]]$params) )
      spl = get_sexual_selection_components(fertilized_eggs, sampled_fertilized_eggs, n_males)
      tmp = tibble(mso = mso_vec,
                   msg = c(spl$msg_female, spl$msg_male),
                   msgc = NA,
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

  output
}


#' Descriptive stats on MS/RS
#'
#' Compute descriptive statistics on MS/RS for all sampling method (including groundtruth) for each replicate
#'
#' @param ms_rs_from_sampling Dataframe with detailled MS/RS (from sampling() function)
#'
#' @details Return per sex, per sampling method and per replicate mean and std for MSo, MSg and RSg.
#'
#' @return Vector of size n_female with $r_p$ for each female.
#'
#' @importFrom dplyr group_by ungroup summarise across
#' @importFrom rlang .data
#'
#' @export
#'
descriptive_stats = function(ms_rs_from_sampling){
  ms_rs_from_sampling %>% group_by(.data$sex, .data$sampling_method, .data$parameters_string, .data$replicate) %>%
    summarise(across(.data$mso:.data$rsg,  list(mean = ~ mean(.x, na.rm = T),
                                           sd = ~ sd(.x, na.rm = T))))
}

#' Standardized MS/RS
#'
#' Standardized MS/RS by sex, sampling method, sampling parameters and replicate
#'
#' @param samples Dataframe similar to the one produce by sampling()
#'
#' @details Standardize MS/RS - using x / mean(x) (see Lande ??)
#'
#' @return Standardize MS/RS
#'
#' @importFrom dplyr group_by ungroup summarise across
#' @importFrom rlang .data
#'
#'

ms_rs_scaling = function(samples){
  return(samples %>% group_by(.data$sex, .data$sampling_method, .data$parameters_string, .data$replicate) %>%
           mutate(across(.data$mso:.data$rsg, ~ .x / mean(.x, na.rm = T))) %>%
           ungroup())
}


#' Correlated paternity
#'
#' Compute detailed paternity correlations (i.e. by female $r_p$)
#'
#' @param fertilized_eggs List of size n_females containing father for each fertilized eggs (from e.g., pollen_competition() function or any sampling_XXX() method) (no default)
#'
#' @details Return per female correlated paternity. See companion function pop_average_rp() for further processing.
#'
#' @return Vector of size n_female with $r_p$ for each female.
#'
#' @export
#'
correlated_paternity = function(fertilized_eggs){
  unlist(lapply(fertilized_eggs, rp_calculation))
}

#' Correlated paternity (internal function)
#'
#' @param female_seed_set female seed set

rp_calculation = function(female_seed_set){
  count_by_male = table(female_seed_set)
  numerator = sum(count_by_male^2-count_by_male)
  denominator = sum(count_by_male)^2 - sum(count_by_male)
  return( numerator / denominator )
}

#' Pop. average correlated paternity
#'
#' Return population average correlated paternity - using Dorken & Perry method and potential correction.
#'
#' @param correlated_paternity Vector of size n_females containing by female correlated paternity (no default)
#' @param n_males Number of male in the population
#'
#' @details Compute pop. average correlated paternity (following Dorken & Perry, 2017) and corrected version (pop average * (n_male - 1)).
#'
#' @return List of two elements, the uncorrected pop. average correlated paternity and the corrected version.
#'
#' @export
#'
pop_average_rp = function(correlated_paternity, n_males){
  r = list()
  if(any(is.nan(correlated_paternity)|is.na(correlated_paternity))){
    print("Warning: some female correlated paternity is NaN or NA, will be excluded from the average")
  }
  r[['pop_average_correlated_paternity']] = mean(correlated_paternity, na.rm = TRUE)
  r[['pop_average_correlated_paternity_corrected']] = r[['pop_average_correlated_paternity']] * (n_males - 1)
  return(r)
}
