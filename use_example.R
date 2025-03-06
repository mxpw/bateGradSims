library(bateGradSims)

n_females = 50
n_males = 50
n_gamete_fem = 250
ratio_gamete = 1
cv_normal_male = 0.1

# set.seed(42)

# Get gametes
gametes = gametes_drawing(n_females = n_females, n_males = n_males,
                          mean_gamete_female = n_gamete_fem, ratio_gamete = ratio_gamete,
                          male_distrib_params = cv_normal_male)

# Get males comp. values (either using uniform, normal, beta or any distributions)
males_comp_values = get_male_comp_values(n_males = n_males, distrib = sample, dist_params = list(x=1:100, replace = TRUE), plot = T)
males_comp_values = get_male_comp_values(n_males = n_males, distrib = sample, dist_params = list(x=c(49:50), replace = TRUE), plot = T)
# males_comp_values = get_male_comp_values(n_males = n_males, distrib = rnorm, dist_params = list(mean = 100, sd = 0), plot = T)
# males_comp_values = get_male_comp_values(n_males = n_males, distrib = rbeta, dist_params = list(shape1 = 5, shape2 = 5), plot = T)


# from a conditional distrib. e.g., depending on gametes counts
# males_comp_values = get_male_comp_values_from_feature(mean_comp_value = 10,
#                                                       sd_comp_value = 20,
#                                                       rho = -0.5,
#                                                       feature = gametes$gam_male,
#                                                       plots=T)

# Get pollen repartition
pollen_repartition = pollen_export(n_females = n_females,
                                    gametes_by_male = gametes$gam_male,
                                    pollen_repartition = c(1),
                                    plot = T)

# Mating Success (observed - exact) can be computed from the pollen_repartition
mso = ms_obs(pollen_repartition)

# Pollen competition, i.e., who fertilizes who
fertilized_eggs = pollen_competition(pollen_repartition, males_comp_values, gametes$gam_female)

# Pollen limitations stats
pollen_limitation_stats = pollen_limitation(gametes$gam_female,fertilized_eggs)

# Aborted eggs
fertilized_eggs = eggs_abortion(fertilized_eggs, aborded_fraction = 0)

# Correlated paternity (before sampling)
pop_average_rp(correlated_paternity = correlated_paternity(fertilized_eggs),
                                 n_males = n_males)


# Using any sampling method ...
# Revoir à partir d'ici
sampled_fertilized_eggs = sampling_groundtruth(fertilized_eggs)
get_sexual_selection_components(fertilized_eggs = fertilized_eggs,
                                sampled_fertilized_eggs = sampled_fertilized_eggs,
                                n_males = n_males)
# The next function is only applicable for sampling groundtruth at the moment
compute_msgc(fertilized_eggs, n_males)

sampled_fertilized_eggs = sampling_fixed(fertilized_eggs, n_males, by_female_samples = 50, undercount_female = 'keep')
get_sexual_selection_components(fertilized_eggs = fertilized_eggs,
                                sampled_fertilized_eggs = sampled_fertilized_eggs,
                                n_males = n_males)

sampled_fertilized_eggs = sampling_prorata(fertilized_eggs, n_males, by_female_prop = 0.5, min_threshold = 30,
                                           undercount_female = 'remove_and_upsample', upsample_strategy = 's1', upsampling_plot = T)
get_sexual_selection_components(fertilized_eggs = fertilized_eggs,
                                sampled_fertilized_eggs = sampled_fertilized_eggs,
                                n_males = n_males)

sampled_fertilized_eggs = sampling_random(fertilized_eggs, n_males)
get_sexual_selection_components(fertilized_eggs = fertilized_eggs,
                                sampled_fertilized_eggs = sampled_fertilized_eggs,
                                n_males = n_males)

sampled_fertilized_eggs = sampling_groundtruth(fertilized_eggs)
selection_components = get_sexual_selection_components(fertilized_eggs = fertilized_eggs,
                                sampled_fertilized_eggs = sampled_fertilized_eggs,
                                n_males = n_males)

library(tidyverse)
tab = tibble::tibble(msg = c(selection_components$msg_female, selection_components$msg_male),
       rsg = c(selection_components$rsg_female, selection_components$rsg_male),
       n_gam = c(gametes$gam_female, gametes$gam_male),
       sex = c(rep("F", n_females), rep("M", n_males)))

# If scaling is needed
tab = tab %>%
  group_by(sex) %>%
  mutate(across(msg:rsg, ~ .x / mean(.x, na.rm = T))) %>%
  ungroup()

# Bateman
summary(glm(data = tab, "rsg ~ msg * sex", family = gaussian))
# Partial Bateman
summary(glm(data = tab, "rsg ~ (msg + n_gam) * sex", family = gaussian))


# To avoid repetition and extra-gathering work, one can use the sampling() function to use multiple sampling methods (given as list as follow)
# Better to specify all parameters (even default ones) if one wants to keep tracks of all of them in outputs
n_eggs_sampled = 2000
methods = list(fixed = list(method = "sampling_fixed", params = list(total_samples = n_eggs_sampled, undercount_female = 'remove')),
               prorata = list(method = "sampling_prorata", params = list(total_samples = n_eggs_sampled, min_threshold = 10,
                                                                         undercount_female = 'remove', upsample_strategy = 's1')),
               prorata = list(method = "sampling_prorata", params = list(total_samples = n_eggs_sampled, min_threshold = 10,
                                                                         undercount_female = 'remove', upsample_strategy = 's2')),
               random = list(method = "sampling_random", params = list(total_samples = n_eggs_sampled)))

# The n_rep apply to all sampling methods
samples = sampling(fertilized_eggs, n_males, methods = methods, mso = mso, gametes = gametes, n_rep = 2)


# Descriptive stats can be obtained from samples
descriptive_stats(samples)

# Then gradients can be obtained (MS/RS can be scale before fitting, see scaled argument)
gradients = fit_gradients(samples, scaled = TRUE)
gradients$gradients

gradients$gradients %>% group_by(sampling_method, parameters_string) %>% summarise(across(abs_female_noGamControl:delta_male_gamControl, ~ mean(.x)))

library(tidyverse)
gradients$gradients %>%
  pivot_longer( cols = starts_with(c('abs', 'delta')),
                names_to = c("type", "sex", "gamControl"),
                names_pattern = "(.*)_(.*)_(.*)" )

plot_list = list()
for(i in 1:nrow(gradients$glms)){
  print(paste0( gradients$glms$sampling_method[[i]], '::' , gradients$glms$parameters_string[[i]] ))
  plot_list[[i]] = ggplot( (gradients$glms %>% pull(data))[[i]], aes(x = msg, y = rsg, color = sex))+
    geom_point()+
    geom_smooth(method='lm')+
    theme_bw()
}
cowplot::plot_grid(plotlist = plot_list)
