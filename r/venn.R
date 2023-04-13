load('environments/data.RData')

eulerr_options(labels = list(cex = 3))

plot <- euler(
    c(
        Fox_Gain = fox_gain,
        H3K27ac_Gain = h3_gain,
        Fox_Loss = fox_loss,
        H3K27ac_Loss = h3_loss,
        "Fox_Gain&H3K27ac_Gain" = gain_gain,
        "Fox_Loss&H3K27ac_Loss" = loss_loss,
        "Fox_Gain&H3K27ac_Loss" = gain_loss,
        "Fox_Loss&H3K27ac_Gain" = loss_gain
    )
) %>% plot(
    quantities = list(type = c("counts"), cex = 3),
    lty = 0,
    alpha = 0.75,
    fills = c(
        Fox_Gain = "#ACA4E1",
        Fox_Loss = "#39BDB1",
        H3K27ac_Gain = "#DB9D85",
        H3K27ac_Loss = "#ABB064"
    ),
    labels = c(
        Fox_Gain = "FOXA1 UP",
        H3K27ac_Gain = "H3K27ac UP",
        Fox_Loss = "FOXA1 DN",
        H3K27ac_Loss = "H3K27ac DN"
    )
)

tags <- plot$children[[1]]$children[[1]]$children$tags$children

tags <- do.call(grid::gList, lapply(tags, function(x) {
    x$children[[2]]$just <- NULL
    x$children[[2]]$hjust <- 0.5
    x$children[[2]]$vjust <- 1
    x$children[[1]]$vjust <- .75
    x
}))

plot$children[[1]]$children[[1]]$children$tags$children <- tags

plot


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
