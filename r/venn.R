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
    merge(diffexpressed, by = "name") %>%
    pivot_wider(names_from = "comparison",
                values_from = "diffexpressed")



euler(
    c(
        GAIN = 26,
        FOXUP = 96,
        H3UP = 54,
        "GAIN&FOXUP" = 13,
        "GAIN&H3UP" = 19,
        "GAIN&FOXUP&H3UP" = 9
    ),
    shape = "circle"
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
