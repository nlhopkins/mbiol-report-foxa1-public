load('environments/data.RData')

eulerr_options(labels = list(cex = 3))

euler(
    c(
        Fox_Gain = fox_gain,
        H3K27ac_Gain = h3_gain,
        "Fox_Gain&H3K27ac_Gain" = gain_gain
    )
) %>% plot(
    quantities = list(type = c("counts", "percent"), cex = 3),
    lty = 0,
    alpha = 0.75,
    fills = c(
        Fox_Gain = "#ACA4E1",
        H3K27ac_Gain = "#DB9D85"
    ),
    labels = c(
        Fox_Gain = "FOXA1 UP",
        H3K27ac_Gain = "H3K27ac UP"
    )
)


x <- diffexpressed %>% pivot_wider(names_from = comparison,
                                   values_from = diffexpressed) %>%
    filter(fox_diff == "GAIN")

y <- data %>% filter(grepl("dox_h3k27ac", treatment)) %>%
    filter(h3_status == "gain")

z <- inner_join(x, y, by = "name")


### Active UP
p <- data  %>%
    mutate(activity = ifelse(value > activity_threshold, '1', '0')) %>%
    merge(diffexpressed, by = "name") %>%
    pivot_wider(names_from = "comparison",
                values_from = "diffexpressed")



euler(
    c(
        GAIN = 23,
        FOXUP = 25,
        H3UP = 54,
        "GAIN&FOXUP" = 12,
        "GAIN&H3UP" = 17,
        "GAIN&FOXUP&H3UP" = 9
    ),
    shape = "ellipse"
)  %>%
    plot(
        quantities = list(type = c("counts"), cex = 3),
        lty = 0,
        alpha = 0.75,
        fills = c(
            GAIN = "#86B875",
            FOXUP = "#ACA4E1",
            H3UP = "#DB9D85"
        ),
        labels = c(GAIN = "Activity\nGAIN",
                   FOXUP = "FOXA1\nUP",
                   H3UP = "H3K27ac\nUP")
    )


data.frame(
    group = c("FOX UP",
              "FOX DN",
              "FOX NS",
              "H3 UP",
              "H3 DN",
              "H3 NS",
              "Co UP",
              "Co DN"),
    value = c(12,
              0,
              11,
              17,
              0,
              5,
              9,
              0)
) %>% ggplot(aes(x = "", y = value, fill =
                     group)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0)


save.image(file = 'environments/venn.RData')
