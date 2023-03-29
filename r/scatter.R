load('environments/data.RData')


# q
x <- data %>% pivot_wider(names_from = treatment)

up %>%
    ggplot(aes(x = log2(h3_fold),
               y = log2(fox_fold))) +
    geom_point(color = "black") +
    geom_smooth(
        method = lm ,
        color = "red",
        fill = "#69b3a2",
        se = TRUE
    ) +
    theme_classic()


save.image(file = 'environments/scatter.RData')
