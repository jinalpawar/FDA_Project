#Getting Data
library(tidyverse)
ca_weather <- fda::CanadianWeather

ca_weather_df <- ca_weather$dailyAv |>
  as.data.frame.table(stringsAsFactors = FALSE) |>
  as.tibble() |>
  rename(
    day = Var1,
    station = Var2,
    variable = Var3,
    value = Freq
  ) |>
  mutate(day = match(day, rownames(CanadianWeather$dailyAv))) |>
  pivot_wider(
    names_from = variable,
    values_from = value
  ) |>
  mutate(station = as_factor(station))

save(ca_weather_df, file = "data/ca_weather.RData")
