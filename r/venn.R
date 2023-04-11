load('environments/data.RData')

eulerr_options(labels = list(cex = 1.5))

euler(
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
    quantities = list(type = c("counts"), cex = 1.5),
    lty = 0,
    fills = colorspace::rainbow_hcl(n = 4, alpha = 0.75),
    labels = c(
        "FOXA1 UP",
        "H3K27ac UP",
        "FOXA1 DN",
        "H3K27ac DN"
    )
)

tags <- plot$children[[1]]$children[[1]]$children$tags$children

tags <- do.call(grid::gList, lapply(tags, function(x) {
    x$children[[2]]$just <- NULL
    x$children[[2]]$hjust <- 0.5
    x$children[[2]]$vjust <- 1
    x$children[[1]]$vjust <- .75
    x}))

plot$children[[1]]$children[[1]]$children$tags$children <- tags

plot


x <- diffexpressed %>% pivot_wider(names_from = comparison,
                                   values_from = diffexpressed) %>%
    filter(fox_diff == "GAIN")

y <- data %>% filter(grepl("dox_h3k27ac", treatment)) %>%
    filter(h3_status == "gain")

z <- inner_join(x, y, by = "name")



save.image(file = 'environments/venn.RData')
