load('environments/data.RData')


euler(c(Fox_Gain = fox_gain, H3K27ac_Gain = h3_gain,
        Fox_Loss = fox_loss, H3K27ac_Loss = h3_loss,
        "Fox_Gain&H3K27ac_Gain" = gain_gain,
        "Fox_Loss&H3K27ac_Loss" = loss_loss,
        "Fox_Gain&H3K27ac_Loss" = gain_loss,
        "Fox_Loss&H3K27ac_Gain" = loss_gain)) %>% plot


save.image(file = 'environments/venn.RData')

gain_gain
