load('environments/data.RData')


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
) %>% plot



x <- diffexpressed %>% pivot_wider(names_from = comparison,
                                   values_from = diffexpressed) %>% 
        filter(fox_diff == "GAIN" & h3_diff == "GAIN") %>%
        select(name)

y <- data %>% filter(grepl("ed_h3k27ac", treatment)) %>%
        filter(value > activity_threshold) %>%
        select(c(name, h3_status, value, h3_fold))

z <- inner_join(x, y, by = "name")






save.image(file = 'environments/venn.RData')
