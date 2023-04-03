load('environments/data.RData')


# active genes
data %>%
    filter(grepl("h3k27ac", treatment)) %>%
    mutate(h3_status = ifelse(value > activity_threshold,
                              'active', 'inactive')) %>%
    group_by(h3_status, condition) %>%
    mutate(active_count = n()) %>%
    ungroup() %>%
    group_by(condition) %>%
    mutate(gene_count = n()) %>%
    ungroup() %>%
    ggplot(aes(
        x = factor(condition, levels = c("ed", "dox")),
        y = active_count / gene_count,
        fill = h3_status
    )) +
    geom_bar(position = "dodge", stat = "identity") +
    xlab("") +
    scale_fill_discrete(name = "") +
    theme_ipsum() +
    scale_y_continuous(labels = scales::percent) +
    geom_text(
        position = position_dodge(width = 1),
        aes(label = active_count / 2,
            y = active_count / gene_count),
        size = 3
    )





binding <- raw %>%
    mutate(dox_fox_bound = ifelse(dox_foxa1_high > dox_biding_threshold, '1', '0')) %>%
    mutate(ed_fox_bound = ifelse(ed_foxa1 > ed_biding_threshold, '1', '0')) %>%
    mutate(dox_h3k27ac_bound = ifelse(dox_h3k27ac > dox_biding_threshold, '1', '0')) %>%
    mutate(ed_h3k27ac_bound = ifelse(ed_h3k27ac > ed_biding_threshold, '1', '0')) %>%
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
    y = bound_count / gene_count,
    fill = bound
)) +
    geom_bar(position = "dodge", stat = "identity") +
    xlab("") +
    scale_fill_discrete(name = "") +
    theme_ipsum() +
    scale_y_continuous(labels = scales::percent) +
    geom_text(position = position_dodge(width = 1),
              aes(
                  label = paste(bound_count),
                  y = bound_count / gene_count,
                  size = 3
              )) +
    facet_wrap( ~ binding)

save.image(file = 'environments/barchart.RData')



x <- binding %>%
    filter(binding == "fox") %>%
    filter(bound == "1") %>% 
    filter(condition == "ed")

y <- binding %>%
    filter(binding == "fox") %>%
    filter(bound == "1")%>%
    filter(condition == "dox") %>%
    select(name)

z <- inner_join(x, y)
