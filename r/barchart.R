load('environments/data.RData')

x <- data %>%
    filter(grepl("ed_h3k27ac|dox_h3k27ac", treatment)) %>%
    mutate(h3_status = case_when(value > threshold ~ "active", TRUE ~ "inactive")) %>%
    group_by(treatment, h3_status) %>%
    mutate(active_count = n() / 2) %>%
    ungroup() %>%
    group_by(treatment) %>%
    mutate(gene_count = n() / 2) %>%
    ungroup()



x %>% ggplot(aes(
    fill = h3_status,
    x = treatment,
    y = active_count / gene_count
)) +
    geom_bar(position = "dodge", stat = "identity") +
    xlab("") +
    scale_fill_discrete(name = "") +
    theme_ipsum() +
    scale_y_continuous(limits = c(0, 0.6),
                       labels = scales::percent)



## IGNORE !!! ACTIVITY APPLIED TO FOX FIX
##
## total %>%
mutate(activity = if_else(value > total_threshold, "active", "inactive")) %>%
    filter(!grepl("dox_foxa1_low", treatment)) %>%
    filter(!grepl("inactive", activity)) %>%
    mutate(treatment = if_else(grepl("ed", treatment), "ed", "dox")) %>%
    ggplot(aes(
        fill = factor(position, level = c("upstream", "downstream")),
        x = factor(treatment, level = c("ed", "dox")),
        y = after_stat(count / sum(count / 4))
    )) +
    geom_bar(position = "dodge") +
    xlab("") +
    scale_fill_discrete(name = "") +
    theme_ipsum() +
    facet_grid(~ target) +
    scale_y_continuous(labels = scales::percent)

save.image(file = 'environments/barchart.RData')
