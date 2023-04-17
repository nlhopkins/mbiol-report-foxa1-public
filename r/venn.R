load('environments/data.RData')

eulerr_options(labels = list(cex = 3))

euler(c(
    Fox_Gain = fox_gain,
    H3K27ac_Gain = h3_gain,
    "Fox_Gain&H3K27ac_Gain" = gain_gain
)) %>% plot(
    quantities = list(type = c("counts", "percent"), cex = 3),
    lty = 0,
    alpha = 0.75,
    fills = c(
        Fox_Gain = "#ACA4E1",
        H3K27ac_Gain = "#DB9D85"
    ),
    labels = c(Fox_Gain = "FOXA1 UP",
               H3K27ac_Gain = "H3K27ac UP")
)


euler(
    c(
        GAIN = 26,
        FOXUP = 96,
        H3UP = 54,
        "GAIN&FOXUP" = 13,
        "GAIN&H3UP" = 19,
        "GAIN&FOXUP&H3UP" = 9,
        "FOXUP&H3UP" = 16
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
        labels = c(GAIN = "Activity GAIN",
                   FOXUP = "FOXA1 UP",
                   H3UP = "H3K27ac UP")
    )




euler(
    c(
        LOSS = 27,
        FOXDN = 62,
        H3DN = 124,
        "LOSS&FOXDN" = 5,
        "LOSS&H3DN" = 24,
        "LOSS&FOXDN&H3DN" = 5,
        "FOXDN&H3DN" = 33
    )
    ,
    shape = "circle"
)  %>%
    plot(
        quantities = list(type = c("counts"), cex = 3),
        lty = 0,
        alpha = 0.75,
        fills = c(
            LOSS = "grey",
            FOXDN = "#39BDB1",
            H3DN = "#ABB064"
        ),
        labels = c(LOSS = "Activity LOSS",
                   FOXDN = "FOXA1 DN",
                   H3DN = "H3K27ac DN")
    )
