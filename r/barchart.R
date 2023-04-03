load('environments/data.RData')


binding <- raw %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    mutate(dox_fox_bound = ifelse(dox_foxa1_high > 1, '1', '0')) %>%
    mutate(ed_fox_bound = ifelse(ed_foxa1 > 1, '1', '0')) %>%
    mutate(dox_h3k27ac_bound = ifelse(dox_h3k27ac > 1, '1', '0')) %>%
    mutate(ed_h3k27ac_bound = ifelse(ed_h3k27ac > 1, '1', '0')) %>%
    pivot_longer(cols = contains("bound"),
                 names_to = "binding",
                 values_to = "bound") %>%
    mutate(condition = case_when(grepl("ed", binding) ~ "ed",
                                 grepl("dox", binding) ~ "dox")) %>%
    mutate(binding = case_when(
        grepl("fox", binding) ~ "fox",
        grepl("h3k27ac", binding) ~ "h3k27ac"
    )) %>%
    group_by(condition, binding, bound) %>%
    mutate(bound_count = n()) %>%
    ungroup() %>%
    group_by(binding) %>%
    mutate(gene_count = n()) %>%
    ungroup()




binding %>% ggplot(aes(
    x = factor(condition, levels = c("ed", "dox")),
    y = bound_count,
    fill = bound
)) +
    geom_bar(position = "dodge", stat = "identity") +
    xlab("") +
    scale_fill_discrete(name = "") +
    theme_ipsum() +
    scale_y_continuous() +
    geom_text(
        position = position_dodge(width = 1),
        aes(label = paste(bound_count),
            y = bound_count),
        size = 3
    ) +
    facet_wrap( ~ binding)

save.image(file = 'environments/barchart.RData')



x <- binding %>%
    filter(binding == "fox") %>%
    filter(bound == "1") %>%
    filter(condition == "ed") %>%
    select(name)

y <- binding %>%
    filter(binding == "fox") %>%
    filter(bound == "1") %>%
    filter(condition == "dox") %>%
    select(name)

z <- anti_join(x, y)

binding %>%
    filter(binding == "fox")  %>%
    filter(activity == "gain") %>%
    count()

binding %>%
    filter(binding == "fox")  %>%
    filter(activity == "loss") %>%
    count()



x <- binding %>%
    filter(binding == "h3k27ac") %>%
    filter(bound == "1") %>%
    filter(condition == "ed") %>%
    select(name)

y <- binding %>%
    filter(binding == "h3k27ac") %>%
    filter(bound == "1") %>%
    filter(condition == "dox") %>%
    select(name)

z <- anti_join(x, y)
