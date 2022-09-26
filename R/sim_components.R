# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Get gametes per individual
#'
#' Draw the number of gametes per individual.
#'
#' @param n_females Number of females (default value : 100)
#' @param n_males   Number of males (default value : 100)
#' @param mean_gamete_female Mean number of gametes by female (default value : 50)
#' @param ratio_gamete Ratio between mean number of gametes by male over mean number of gametes by male (e.g., if ratio_gamete = 10, males have 10x more gametes than females) (default value : 10)
#' @param female_distrib Distribution of gametes by female (default value : poisson)
#' @param male_distrib Distribution of gametes by male (default value : normal)
#' @param female_distrib_params If distribution for female needs more than the mean, provide additional arguments
#' @param male_distrib_params If distribution for female needs more than the mean, provide additional arguments (here, Coefficient of variation)
#'
#' @details Implemented distributions for males/females are 'poisson' and 'normal' - others might be possible in the future.
#' Normal distribution allows user to specify larger variance (by specifying the coefficient of variation) than expected with the more constrained Poisson distribution. The standard deviation
#' of the distribution is computed as 'cv x mean', cv being the X_distrib_params parameters, and the obtained values are then rounded.
#'
#' Note: if used distribution can lead to negatives values, these will be set to zero before being returned.
#'
#' Future: other distributions might be implemented
#'
#' @return List with two vectors (gam_female, gam_male)
#'
#' @importFrom stats rnorm rpois
#'
#' @export
#'
gametes_drawing = function(n_females = 100, n_males = 100, mean_gamete_female = 50, ratio_gamete = 10,
                           female_distrib = 'poisson', male_distrib = 'normal',
                           female_distrib_params = NULL, male_distrib_params = 0.1){

  # Females
  if(female_distrib == 'poisson'){
    gam_female = rpois(n_females, lambda = mean_gamete_female)
  }else if(female_distrib == "normal"){
    gam_female = round(rnorm(n = n_females,
                             mean = mean_gamete_female,
                             sd = (female_distrib_params * mean_gamete_female) ))
  }else{
    stop("Provided distribution is not implemented !")
  }

  # Males
  if(male_distrib == 'poisson'){
    gam_male = rpois(n_males, lambda = mean_gamete_female * ratio_gamete)
  }else if(male_distrib == "normal"){
    gam_male = round(rnorm(n = n_males,
                           mean = mean_gamete_female * ratio_gamete,
                           sd = (male_distrib_params * mean_gamete_female * ratio_gamete) ))
  }else{
    stop("Provided distribution is not implemented !")
  }

  return(list(gam_female = pmax(0, gam_female), gam_male = pmax(0, gam_male)))
}

#' Draw male competitive value
#'
#' Produce a vector of competitive values
#'
#' @param n_males   Number of males (default value : 100)
#' @param distrib   Prob. distribution to use (default value : sample, i.e. uniform distrib.)
#' @param dist_params   Distribution parameters
#' @param translation Should the resulting comp. value be translated to be > 0 ? (default TRUE)
#' @param plot      Should the comp. value histogram be plotted ? (default FALSE)
#'
#' @details Draw males competitive values from provided distribution (uniform by default).
#' After drawing, values are transformed such that all values are strictly above zero (i.e. x' = x + min(x) + 0.01).
#' This is needed because these value are used later as probabilities. (that default behavior can be modified with the translation arg.)
#'
#' @return Vector of competitive values
#'
#' @importFrom graphics hist
#'
#' @export
#'
get_male_comp_values = function(n_males = 100, distrib = sample, dist_params = list(x = 100, replace = TRUE),
                                translation = TRUE, plot = FALSE){

  male_comp = do.call(distrib, c(list(n_males), dist_params))

  if(any(male_comp < 0) & translation) male_comp = male_comp - min(male_comp, 0) + 0.01

  if(plot) hist(male_comp, main = "Comp. values")

  if(any(male_comp < 0)) warning("Some comp. values are negatives - error will arise when using pollen_competition() function is not corrected")

  return(male_comp)
}

#' Modified version of base::sample
#'
#' (internal) Draw a sample from vector, some behaviors of base::sample() are modified
#'
#' @param x      See base::sample()
#' @param size   See base::sample()
#' @param prob   See base::sample()
#'
#' @details Can handle cases for which base::sample() is not made for. E.g., if x vector is empty or size is 0 or NA, return NULL. If length(x)==1, return x (base::sample would have return from 1:x)
#' In addition, the function check whether size is greater than length(x), in that case, size is reduced to length(x).
#'
#' @return A vector of length size with elements drawn from x (if x is a vector), otherwise see details
#'
#'
sample_handmade = function(x, size, prob = NULL){
  if(length(x) == 0 || size == 0 || is.na(size)) return(NULL)
  if(length(x) == 1) return(x)
  sample(x, size = min(length(x), size), prob = prob)
}

#' Pollen (male) export distribution
#'
#' Generate a matrix (n_females x n_males) with pollen load using Dirichlet & Multinomial distributions.
#'
#' @param n_females           Number of females (default value : NULL)
#' @param baseline_alpha      Baseline alpha values for Dirichlet distr. - vector of size n_females; see details (default value : NULL)
#' @param gametes_by_male     Vector of number of gametes by males (e.g., from gametes_drawing() function)
#' @param pollen_repartition  How to handle pollen heterogeneity distribution See details. (default value : 0.01)
#' @param plot                Should pollen repartition be plotted ? (default False)
#'
#' @details Return the number of pollen grains by females (rows) exported by each male (columns).
#' The general idea being to draw one vector per-male (\eqn{p_i = Dir(\alpha_i)}, for i in \eqn{1:n_males}) of relative pollen export to females (using Dirichlet distribution),
#' then, to draw a realization of this vector using Multinomial distribution (see also Polya distribution).
#'
#' The baseline_alpha (i.e. \eqn{\alpha_i = baseline_alpha} for all i) represents the alpha vector of the Dirichlet distribution (i.e. and represent how pollen from one male
#' is, on average, distributed among females), and can be tuned to e.g. represent attractiveness of females or any component
#' making pollen load more or less likely on a given female. If not provided, the n_females variable will be used
#' to create an alpha vector containing only ones. Note that it means that at least baseline_alpha or n_females should be
#' provided, and n_females will not be used if baseline_alpha is provided.
#'
#' pollen_repartition argument allows specifying various heterogeneity schems for pollen repartition beyond among-females 'preferences' (baseline_alpha).
#' Two levels of heterogenity can be considered.
#' (i) heterogeneity level among females, (ii) heterogeneity level among males for among-females heterogeneity. (Oo glurp).
#' When \eqn{\alpha} is the same for all males (no heterogeneity among males), the level of among-female heterogeneity is driven by a fixed float provided to pollen_reparition
#' (e.g., 0.01 or 10). Generally speaking, the largest is that value, the most homogeneous is the among-female pollen export. Actually, that value is used to scale the baseline_alpha
#' vector (i.e. newalpha = baseline_alpha * pollen_repartition), thus be careful as effect of this value on heterogeneity will depends on the number
#' of females as well as on the baseline_alpha values ! (see Dirichlet distribution for expressions on variance). A more convenient implementation might be attemp in the future
#' to make this dependency less embarassing. (could use that value as sum(alpha) but will still have that dependency on n_female ? or
#' trying to have a given quantity constant somewhere, but which one ? nor clear whether it is possible nor how to do it).
#'
#'
#' Alternatively, heterogeneity levels might be different among males. In these cases, pollen_repartition could be a vector of size n_males, and
#' baseline_alpha will be scaled for each male depending on its pollen_repartition value. A last option, is to specify a list(mean=XX,sd=XX) for pollen_repartition,
#' values provided here will be used to draw a value for each male from a 10^N(mean, sd) distribution.
#'
#' ## Might have some improvements to do here...
#'
#' Because details might appear unclear, it's a good practice to look at plot displayed by the function. One of them represents the among-females
#' pollen repartition for a subsample of males (here, 15 - might be added as a argument of the function...), and the other, the average number of females
#' on which pollens from one male land.
#'
#' @return A matrix (n_females x n_males) with pollen load
#'
#' @importFrom gtools rdirichlet
#' @importFrom stats rmultinom
#' @importFrom cowplot plot_grid
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom tibble tibble as_tibble rownames_to_column
#' @importFrom tidyr pivot_longer
# @importFrom rlang .data
#'
#' @export

pollen_export = function(n_females = NULL, baseline_alpha = NULL,
                         pollen_repartition = c(0.01), plot = FALSE, gametes_by_male = NULL){

  if(is.null(gametes_by_male))
    stop("agument 'gametes_by_male' is missing, with no default")

  n_males = length(gametes_by_male)

  if(is.null(baseline_alpha) & is.null(n_females))
    stop("One of 'baseline_alpha' or 'n_females' must be provided")

  if(is.null(baseline_alpha))
    baseline_alpha = rep(1, n_females)

  if(is.null(n_females))
    n_females = length(baseline_alpha)

  # Check alpha_0 for Dirichlet
  if(length(pollen_repartition) == 1){
    male_gamete_export_heterogeneity = rep(pollen_repartition, n_males)
  }else if(length(pollen_repartition) == n_males){
    male_gamete_export_heterogeneity = 10^pollen_repartition
  }else if( is.list(pollen_repartition) & all(c('mean', 'sd') %in% names(pollen_repartition)) ){
    male_gamete_export_heterogeneity = 10^do.call(rnorm, c(list(n_males), pollen_repartition))
  }else{
    stop('Argument pollen_repartition is badly formatted, check help("pollen_export")')
  }

  dirichlet_draw = sapply(male_gamete_export_heterogeneity, FUN = function(x) rdirichlet(1, alpha = baseline_alpha * x) )
  male_gamete_repartition = matrix(nrow = n_females, ncol = n_males)
  for(i in 1:n_males)
    male_gamete_repartition[,i] = rmultinom(1, gametes_by_male[i], dirichlet_draw[,i])

  if(plot){
    tp = dirichlet_draw[,sample(dim(dirichlet_draw)[2], 15)] %>%
      as_tibble() %>%
      rownames_to_column() %>%
      pivot_longer(-.data$rowname)
    distrib_ = ggplot(tp, aes(x = .data$name, y = .data$value, fill = .data$rowname))+
      geom_bar(stat = 'identity', color="black")+
      coord_flip()+
      theme(legend.position = 'none')+
      labs(x = "Males", y = "Pollen repartition among females")

    hist_ = ggplot(tibble(cnt = apply(male_gamete_repartition, 2, FUN = function(x) sum(x>0))),
                   aes(x=.data$cnt))+
      geom_histogram()+
      labs(x = "Number of females", y = "")

    print(plot_grid(distrib_, hist_, ncol = 2))
  }

  return(male_gamete_repartition)
}


#' Pollen competition
#'
#' Pollen competition to access 'eggs'
#'
#' @param pollen_repartition      Matrix females X males with number of male gametes by females (e.g., from pollen_export() function) (no default)
#' @param males_comp_values       Males competitive values (e.g., from get_male_comp_values() function ) (no default)
#' @param gametes_by_female       Vector of number of gametes by females (e.g., from gametes_drawing() function) (no default)
#' @param pollen_limitation_stats (logical) Should statistics on pollen limitation be printed ? (default False)
#'
#' @details Draw males identity for each fertilized eggs according to male comp. values (draw w/o replacement).
#' Pollen limitation (i.e. all eggs are not necessary fertilized - if pollen load is to low) is handled by the function - statistics can be displayed
#' by setting pollen_limitation_stats arguments to TRUE.
#'
#' Note: more complex competition schemes could be considered in the future (e.g. (in)compatibilities, phenomenological equivalent of (dis)assortative mating,
#' equivalent of phenological (dis)matches or kind of preemption effects(?)).
#'
#' @return List of size n_females, each element being a vector of males IDs corresponding to fertilized eggs
#'
#' @export
#'

pollen_competition = function(pollen_repartition, males_comp_values, gametes_by_female, pollen_limitation_stats = FALSE){
  # Wonder whether it might be possible to vectorized that one
  female_desc = list()
  n_female = dim(pollen_repartition)[1]
  n_male = dim(pollen_repartition)[2]

  for(f in 1:n_female){
    vector_to_sample = c()
    vector_competition = c()
    for(m in 1:n_male){
      vector_to_sample = c(vector_to_sample, rep(m, pollen_repartition[f, m]))
      vector_competition = c(vector_competition, rep(males_comp_values[m], pollen_repartition[f, m]) )
    }

    if(length(vector_to_sample) == 0){
      female_desc[f] = list(NULL)
    }else{
      female_desc[[f]] = sample_handmade(vector_to_sample,
                                         size = gametes_by_female[f], # Carefull, some pollen limitation might exists ! Meaning true draw size can be lower than gamete_by_female
                                         prob = vector_competition / sum(vector_competition))
    }
  }

  if(pollen_limitation_stats){
    print("Pollen limitation : XXXXXXXX TO DO XXXXXXXXXX")
  }

  female_desc
}

