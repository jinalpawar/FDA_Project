library(tidyverse)
library(Matrix)
library(grpreg)
library(fda)
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

#Preprocessing data
ca_weather_processed <- ca_weather_df |>
  relocate(station, day) |>
  rename(id = station, time = day) |>
  mutate(
    id = as.integer(factor(id)),
    time = as.numeric(time)
  ) |>
  arrange(id, time) |>
  as.data.frame() |>
  preprocess()

y <- ca_weather_processed$log10precip
covs <- ca_weather_processed |>
  dplyr::select(!c(log10precip, Precipitation.mm))


selection <- FLCM.select(
  y,
  covs
) #I guess this is it??

print(selection$selected)

library(fda)
library(ggplot2)

extract_flcm_coefs <- function(
  out,
  mydata,
  method = c("scad", "mcp"),
  grid_length = 200,
  include_intercept = FALSE
) {
  method <- match.arg(method)
  est <- out$estimated

  gamma <- if (method == "scad") est$gamma2 else est$gamma3
  M <- if (method == "scad") est$Mscad else est$Mmcp

  nbasis <- est$bbasis$nbasis
  nbasisint <- est$bbasisint$nbasis
  rng <- est$bbasis$rangeval
  t_grid <- seq(rng[1], rng[2], length.out = grid_length)

  cov_names <- names(mydata)[-c(1, 2)]
  p <- length(cov_names)

  # Drop grpreg's leading intercept scalar so indices align with `group`
  coefs_no_scalar <- gamma[-1]

  # Basis matrices evaluated on t_grid
  phi_cov <- eval.basis(t_grid, est$bbasis) %*% M
  phi_int <- eval.basis(t_grid, est$bbasisint)

  out_list <- list()

  if (include_intercept) {
    b0_coef <- coefs_no_scalar[seq_len(nbasisint)]
    out_list[[1]] <- data.frame(
      time = t_grid,
      estimate = as.vector(phi_int %*% b0_coef),
      term = "Intercept (beta_0)",
      selected = any(b0_coef != 0),
      stringsAsFactors = FALSE
    )
  }

  offset <- nbasisint
  for (j in seq_len(p)) {
    idx <- offset + seq_len(nbasis)
    bj_coef <- coefs_no_scalar[idx]
    offset <- offset + nbasis

    out_list[[length(out_list) + 1]] <- data.frame(
      time = t_grid,
      estimate = as.vector(phi_cov %*% bj_coef),
      term = cov_names[j],
      selected = any(bj_coef != 0),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out_list)
}

# Extract curves (switch method = "mcp" for grMCP)
df_coefs <- extract_flcm_coefs(
  selection,
  covs,
  method = "mcp",
  include_intercept = FALSE
)

# Option A: Plot all covariate curves in facets
ggplot(df_coefs, aes(x = time, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.9, color = "#1f77b4") +
  facet_wrap(~term, scales = "free_y") +
  labs(
    title = "Estimated Functional Coefficients (grSCAD)",
    x = "Time",
    y = expression(hat(beta)[j](t))
  ) +
  theme_minimal()

# Option B: Plot ONLY non-zero (selected) covariate curves
df_selected <- subset(df_coefs, selected)

ggplot(df_selected, aes(x = time, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1, color = "#d62728") +
  facet_wrap(~term, scales = "free_y") +
  labs(
    title = "Selected Functional Coefficients (grSCAD)",
    x = "Time",
    y = expression(hat(beta)[j](t))
  ) +
  theme_minimal()
