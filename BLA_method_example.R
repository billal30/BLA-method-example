
# ============================================================
# Boundary Line Analysis of Soil Organic Carbon
#
# Rainfall region: Summer 350-650 mm
# Predictors: Clay and electrical conductivity (EC)
#
# Input file: example_data.csv
# Required columns: bom_class, OC, Clay, EC
#
# Run this script from the repository's main folder.
# ============================================================

# Install packages once if required:
# install.packages(c("MASS", "ggplot2", "dplyr"))

library(MASS)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------
# 1. Read and check input data
# ------------------------------------------------------------

input_file <- "github/example_data.csv"

if (!file.exists(input_file)) {
  stop("Place example_data.csv in the working directory.")
}

df <- read.csv(
  input_file,
  stringsAsFactors = FALSE
)

required_columns <- c("bom_class", "OC", "Clay", "EC")

missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

numeric_columns <- c("OC", "Clay", "EC")

if (!all(vapply(df[numeric_columns], is.numeric, logical(1)))) {
  stop("OC, Clay and EC must be numeric.")
}

# ------------------------------------------------------------
# 2. Select one rainfall region
# ------------------------------------------------------------

selected_region <- "Summer 350-650 mm"

df <- df %>%
  select(all_of(required_columns)) %>%
  filter(bom_class == selected_region)

if (nrow(df) == 0) {
  stop("No observations found for: ", selected_region)
}

message("Rainfall region: ", selected_region)
message("Input observations: ", nrow(df))

# ------------------------------------------------------------
# 3. Helper function: estimate density at each observation
# ------------------------------------------------------------

get_density <- function(x, y, n = 100) {
  
  density_grid <- MASS::kde2d(x, y, n = n)
  
  ix <- findInterval(x, density_grid$x)
  iy <- findInterval(y, density_grid$y)
  
  # Keep indices within grid limits.
  ix <- pmax(1L, pmin(ix, length(density_grid$x)))
  iy <- pmax(1L, pmin(iy, length(density_grid$y)))
  
  density_grid$z[cbind(ix, iy)]
}

# ------------------------------------------------------------
# 4. Function: fit and plot one boundary model
# ------------------------------------------------------------

fit_boundary <- function(data, predictor, x_label) {
  
  # Retain finite predictor-SOC pairs.
  df_sub <- data.frame(
    x = data[[predictor]],
    oc = data$OC
  )
  
  df_sub <- df_sub %>%
    filter(is.finite(x), is.finite(oc))
  
  if (nrow(df_sub) < 10) {
    stop("Fewer than 10 valid observations for ", predictor, ".")
  }
  
  # ----------------------------------------------------------
  # 4a. Classical Mahalanobis screening
  # ----------------------------------------------------------
  
  subdat <- as.matrix(df_sub[, c("x", "oc")])
  
  centre <- colMeans(subdat)
  covmat <- cov(subdat)
  
  if (qr(covmat)$rank < 2) {
    stop("Singular covariance matrix for ", predictor, ".")
  }
  
  md <- mahalanobis(
    subdat,
    center = centre,
    cov = covmat
  )
  
  # Empirical 98th percentile of squared Mahalanobis distances.
  md_threshold <- quantile(
    md,
    probs = 0.98,
    na.rm = TRUE,
    names = FALSE
  )
  
  df_clean <- df_sub[md <= md_threshold, , drop = FALSE]
  
  if (nrow(df_clean) < 10) {
    stop("Fewer than 10 retained observations for ", predictor, ".")
  }
  
  if (
    length(unique(df_clean$x)) < 2 ||
    length(unique(df_clean$oc)) < 2
  ) {
    stop("Insufficient variation for density estimation: ", predictor)
  }
  
  # ----------------------------------------------------------
  # 4b. Two-dimensional kernel density estimation
  # ----------------------------------------------------------
  
  df_clean$density <- get_density(
    x = df_clean$x,
    y = df_clean$oc,
    n = 100
  )
  
  # ----------------------------------------------------------
  # 4c. Define 15 overlapping predictor windows
  # ----------------------------------------------------------
  
  predictor_range <- range(df_clean$x)
  
  half_width <- (
    predictor_range[2] - predictor_range[1]
  ) / 15
  
  window_centres <- seq(
    from = predictor_range[1],
    to = predictor_range[2],
    length.out = 15
  )
  
  # ----------------------------------------------------------
  # 4d. Select the maximum eligible SOC in each window
  # ----------------------------------------------------------
  
  selected <- lapply(window_centres, function(cx) {
    
    in_window <-
      df_clean$x >= (cx - half_width) &
      df_clean$x <= (cx + half_width)
    
    window_data <- df_clean[in_window, , drop = FALSE]
    
    if (nrow(window_data) == 0) {
      return(NULL)
    }
    
    # Exclude observations below the fifth density percentile
    # within this window from boundary-point selection.
    density_threshold <- quantile(
      window_data$density,
      probs = 0.05,
      na.rm = TRUE,
      names = FALSE
    )
    
    eligible <- window_data[
      is.finite(window_data$density) &
        window_data$density >= density_threshold,
      ,
      drop = FALSE
    ]
    
    if (nrow(eligible) == 0) {
      return(NULL)
    }
    
    # Assign maximum SOC to the window centre.
    data.frame(
      x = cx,
      OC = max(eligible$oc, na.rm = TRUE)
    )
  })
  
  boundary_points <- bind_rows(selected)
  
  if (nrow(boundary_points) < 5) {
    stop("Fewer than five boundary points for ", predictor, ".")
  }
  
  # ----------------------------------------------------------
  # 4e. Fit the LOESS boundary
  # ----------------------------------------------------------
  
  boundary_model <- loess(
    OC ~ x,
    data = boundary_points,
    span = 2
  )
  
  training_range <- range(boundary_points$x)
  
  # Predict at retained observations.
  df_clean$OC_pred <- as.numeric(
    predict(
      boundary_model,
      newdata = data.frame(x = df_clean$x)
    )
  )
  
  # Predict on an ordered grid for plotting.
  boundary_line <- data.frame(
    x = seq(
      from = training_range[1],
      to = training_range[2],
      length.out = 300
    )
  )
  
  boundary_line$OC_pred <- as.numeric(
    predict(
      boundary_model,
      newdata = boundary_line
    )
  )
  
  # ----------------------------------------------------------
  # 4f. Create the individual boundary plot
  # ----------------------------------------------------------
  
  boundary_plot <- ggplot() +
    geom_point(
      data = df_clean,
      aes(x = x, y = oc),
      colour = "green",
      alpha = 0.6,
      size = 2
    ) +
    geom_line(
      data = boundary_line,
      aes(x = x, y = OC_pred),
      colour = "black",
      linewidth = 1.2
    ) +
    labs(
      x = x_label,
      y = "SOC (%)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 14)
    )
  
  # Summary of this model.
  model_summary <- data.frame(
    rainfall_region = selected_region,
    predictor = predictor,
    n_valid = nrow(df_sub),
    n_removed_mahalanobis = nrow(df_sub) - nrow(df_clean),
    n_retained = nrow(df_clean),
    n_boundary_points = nrow(boundary_points),
    predictor_min = training_range[1],
    predictor_max = training_range[2],
    window_half_width = half_width,
    mahalanobis_probability = 0.98,
    density_probability = 0.05,
    loess_span = 2
  )
  
  list(
    predictor = predictor,
    model = boundary_model,
    training_range = training_range,
    observations = df_clean,
    boundary_points = boundary_points,
    boundary_line = boundary_line,
    summary = model_summary,
    plot = boundary_plot
  )
}

# ------------------------------------------------------------
# 5. Generate the Clay boundary first
# ------------------------------------------------------------

clay_result <- fit_boundary(
  data = df,
  predictor = "Clay",
  x_label = "Clay (%)"
)

print(clay_result$plot)
print(clay_result$summary)

# ------------------------------------------------------------
# 6. Generate the EC boundary next
# ------------------------------------------------------------

# Add your dataset's EC measurement unit to the axis label.
ec_result <- fit_boundary(
  data = df,
  predictor = "EC",
  x_label = "EC"
)

print(ec_result$plot)
print(ec_result$summary)

# ------------------------------------------------------------
# 7. Function: save one model and its outputs
# ------------------------------------------------------------

output_directory <- "outputs"

dir.create(
  output_directory,
  showWarnings = FALSE,
  recursive = TRUE
)

save_boundary <- function(result, output_directory) {
  
  prefix <- result$predictor
  
  write.csv(
    result$observations,
    file.path(
      output_directory,
      paste0(prefix, "_observations.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    result$boundary_points,
    file.path(
      output_directory,
      paste0(prefix, "_boundary_points.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    result$boundary_line,
    file.path(
      output_directory,
      paste0(prefix, "_boundary_predictions.csv")
    ),
    row.names = FALSE
  )
  
  saveRDS(
    list(
      rainfall_region = selected_region,
      predictor = result$predictor,
      model = result$model,
      training_range = result$training_range,
      settings = list(
        mahalanobis_probability = 0.98,
        density_probability = 0.05,
        density_grid_size = 100,
        n_windows = 15,
        window_half_width = result$summary$window_half_width,
        loess_span = 2,
        loess_degree = 2
      )
    ),
    file.path(
      output_directory,
      paste0(prefix, "_model.rds")
    )
  )
  
  ggsave(
    filename = file.path(
      output_directory,
      paste0(prefix, "_boundary_plot.pdf")
    ),
    plot = result$plot,
    width = 6,
    height = 4.5,
    units = "in"
  )
}

# ------------------------------------------------------------
# 8. Save Clay and EC outputs separately
# ------------------------------------------------------------

save_boundary(clay_result, output_directory)
save_boundary(ec_result, output_directory)

# Save a combined model summary.
model_summary <- bind_rows(
  clay_result$summary,
  ec_result$summary
)

write.csv(
  model_summary,
  file.path(output_directory, "model_summary.csv"),
  row.names = FALSE
)

# Record R and package versions.
writeLines(
  capture.output(sessionInfo()),
  file.path(output_directory, "session_info.txt")
)
