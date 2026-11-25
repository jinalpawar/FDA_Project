library(tidyverse)
source('project/varselect.R')

load('project/data/ca_weather.RData')

ca_weather <- ca_weather |>
  select(station_id, date, mean_temp, total_rain, total_snow)


ggplot(
  data = ca_weather,
  mapping = aes(x = date, y = mean_temp, col = station_id)
) +
  geom_line(alpha = 0.5)
