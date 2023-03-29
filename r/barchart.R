load('environments/data.RData')

total %>%
    mutate(activity = if_else(value > total_threshold, "active", "inactive")) %>%
    filter(!grepl("dox_foxa1_low", treatment)) %>%
    filter(!grepl("inactive", activity)) %>%
    mutate(target = if_else(grepl("foxa1", treatment), "fox", "h3k27ac")) %>%
    mutate(treatment = if_else(grepl("ed", treatment), "ed", "dox")) %>%
    ggplot(aes(
        fill = factor(position, level = c("upstream", "downstream")),
        x = factor(treatment, level = c("ed", "dox")),
        y = after_stat(count / sum(count / 2))
    )) +
    geom_bar(position = "dodge") +
    xlab("") +
    scale_fill_discrete(name = "") +
    theme_ipsum() +
    facet_grid( ~ treatment) +
    scale_y_continuous(labels = percent)

save.image(file = 'environments/barchart.RData')
