load('environments/data.RData')

library(ggplot2)


# Heatmap


data %>%
    pivot_wider(names_from = treatment, values_from = value) %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac", "fox_abs")),
                 names_to = "treatment") %>%
    filter(!grepl("h3k27ac|dox_foxa1_low", treatment)) %>%
    arrange(factor(treatment, levels = c("ed_foxa1", "", "")), value) %>%
    ggplot(aes(factor(
        treatment, c("ed_foxa1", "dox_foxa1_high", "fox_abs")
    ),
    reorder(name, value),
    fill =  value)) +
    geom_tile() +
    scale_fill_gradient2(low = "purple",
                         mid = "white",
                         high = "green")
