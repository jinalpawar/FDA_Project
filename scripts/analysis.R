library(tidyverse)
require(withr)
source('scripts/varselect.R')

load("data/ca_weather.RData")


#Number of unique observations
N <- length(unique(ca_weather_df$station))

#### Visualization ####
#Choosing some individuals to visualize at random
sample_ids <- withr::with_seed(1234, sample(ca_weather_df$station, 5))

#Choosing some colors
sample_colors <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#AA3377")


#Adding a color and an alpha column
ca_weather_plot <- ca_weather_df |>
  mutate(
    color = as_factor(if_else(station %in% sample_ids, station, "Other"))
  )


ggplot(
  data = ca_weather_plot,
  mapping = aes(x = day, y = Temperature.C, col = color)
) +
  geom_line(aes(alpha = color), linewidth = 1) +
  scale_color_manual(
    values = setNames(
      c("gray", sample_colors),
      c("Other", as.character(sample_ids))
    )
  ) +
  scale_alpha_manual(
    values = setNames(
      c(0.7, rep(1, 5)),
      c("Other", as.character(sample_ids))
    )
  ) +
  xlab("Day") +
  ylab("Average Temperature") +
  ggtitle("(A)")

ggplot(
  data = ca_weather_plot,
  mapping = aes(x = day, y = log10precip, col = color)
) +
  geom_line(aes(alpha = color), linewidth = 1) +
  scale_color_manual(
    values = setNames(
      c("gray", sample_colors),
      c("Other", as.character(sample_ids))
    )
  ) +
  scale_alpha_manual(
    values = setNames(
      c(0.7, rep(1, 5)),
      c("Other", as.character(sample_ids))
    )
  ) +
  xlab("Day") +
  ylab("log10precip") +
  ggtitle("(B)")
