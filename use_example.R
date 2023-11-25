library(bateGradSims)

n_females = 100
n_males = 90
n_gamete_fem = 200
ratio_gamete = 0.8
cv_normal_male = 0.1

# Get gametes
gametes = gametes_drawing(n_females = n_females, n_males = n_males,
                          mean_gamete_female = n_gamete_fem, ratio_gamete = ratio_gamete,
                          male_distrib_params = cv_normal_male)

# Get males comp. values
males_comp_values = get_male_comp_values(n_males = n_males, dist_params = list(x=1:100, replace = TRUE), plot = T)

# Get pollen repartition
pollen_repartition = pollen_export(n_females = n_females,
                                    gametes_by_male = gametes$gam_male,
                                    pollen_repartition = c(0.01),
                                    plot = T)

# Mating Success (observed - exact) can then be computed from the pollen_repartition
mso = ms_obs(pollen_repartition)

# Pollen competition, i.e., who fertilizes who
fertilized_eggs = pollen_competition(pollen_repartition, males_comp_values, gametes$gam_female)

# Aborted eggs
fertilized_eggs = eggs_abortion(fertilized_eggs, aborded_fraction = 0)


# Get samples 'by hand' - not the better strategy in order to compare method afterward
sampling_groundtruth(fertilized_eggs, n_males)
sampling_fixed(fertilized_eggs, n_males, by_female_samples = 50, undercount_female = 'keep')
sampling_prorata(fertilized_eggs, n_males, by_female_prop = 0.10, min_threshold = 0,
                 undercount_female = 'keep', upsample_strategy = 's2', upsampling_plot = T)
sampling_random(fertilized_eggs, n_males)


# Get samples with wrapper (groundtruth should always be estimated - that's the baseline, might be added directly within the sampling() function)
# Better to specify all parameters (even default ones) if one wants to keep tracks of all of them
methods = list(base = list(method = "sampling_groundtruth", params = list()),
               fixed = list(method = "sampling_fixed", params = list(total_samples = 4000, undercount_female = 'remove')),
               prorata = list(method = "sampling_prorata", params = list(total_samples = 4000, min_threshold = 5,
                                                                         undercount_female = 'remove', upsample_strategy = 's1')),
               prorata = list(method = "sampling_prorata", params = list(total_samples = 4000, min_threshold = 5,
                                                                         undercount_female = 'remove', upsample_strategy = 's2')),
               random = list(method = "sampling_random", params = list(total_samples = 4000)))

samples = sampling(fertilized_eggs, n_males, methods = methods, mso = mso, gametes = gametes, scaled = T)
gradients = fit_gradients(samples)
gradients$gradients

library(tidyverse)
gradients$gradients %>%
  pivot_longer( cols = starts_with(c('abs', 'delta')),
                names_to = c("type", "sex", "gamControl"),
                names_pattern = "(.*)_(.*)_(.*)" )

plot_list = list()
for(i in 1:nrow(gradients$glms)){
  print(paste0( gradients$glms$sampling_method[[i]], '::' , gradients$glms$parameters_string[[i]] ))
  plot_list[[i]] = ggplot( (gradients$glms %>% pull(data))[[i]], aes(x = msg, y = rsg, color = sex))+geom_point()+geom_smooth(method='lm')
}
cowplot::plot_grid(plotlist = plot_list)
