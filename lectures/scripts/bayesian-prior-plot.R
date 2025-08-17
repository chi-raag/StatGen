library(ggplot2)
library(dplyr)
library(gganimate)
library(viridis)

p <- seq(0, 1, length.out = 1000)

# Prior (Beta distribution from previous experiment)
prior <- dbeta(p, 3 + 1, 17 + 1)

# Likelihood (Beta distribution from current data)
likelihood <- dbeta(p, 11 + 1, 16 + 1)

# Posterior (Beta distribution combining both)
posterior <- dbeta(p, 14 + 1, 33 + 1)

# Combine into a data frame
plot_data <- data.frame(
  p = rep(p, 3),
  density = c(prior, likelihood, posterior),
  type = rep(c("Prior", "Likelihood", "Posterior"), each = length(p))
)

# Add flat prior distribution
flat_prior <- dbeta(p, 1, 1)
flat_posterior <- dbeta(p, 12, 17)

plot_data <- rbind(
  plot_data,
  data.frame(
    p = rep(p, 3),
    density = c(flat_prior, likelihood, flat_posterior),
    type = rep(c("Prior", "Likelihood", "Posterior"), each = length(p))
  )
)

# Add a column to distinguish between informative and flat scenarios
plot_data$scenario <- c(rep("Informative", 3000), rep("Flat", 3000))

# Create plot
# Plot for informative prior
p1 <- ggplot(
  plot_data %>% filter(scenario == "Informative"),
  aes(x = p, y = density, color = type)
) +
  geom_line(size = 2) +
  geom_vline(xintercept = 11 / 27, linetype = "dashed", color = "red") +
  labs(
    x = "Allele Frequency (p)",
    y = "Density",
    title = "Bayesian Update with Informative Prior"
  ) +
  scale_color_viridis_d() +
  theme_minimal(base_size = 15)

ggsave(
  plot = p1,
  filename = "presentations/session-01/_output/informative_prior_plot.png",
  width = 7,
  height = 4,
  units = "in",
  dpi = 150
)

# Plot for flat prior with different line types for overlapping curves
p2 <- ggplot(
  plot_data %>% filter(scenario == "Flat"),
  aes(x = p, y = density, color = type, linetype = type)
) +
  geom_line(size = 2) +
  scale_linetype_manual(values = c("Prior" = "dotted", "Likelihood" = "solid", "Posterior" = "dashed")) +
  geom_vline(xintercept = 11 / 27, linetype = "dotdash", color = "red") +
  labs(
    x = "Allele Frequency (p)",
    y = "Density",
    title = "Bayesian Update with Flat Prior"
  ) +
  theme_minimal(base_size = 15) +
  scale_color_viridis_d() +
  annotate("text",
    x = 0.7, y = max(flat_posterior) * 0.8,
    label = "Likelihood ≈ Posterior", color = "gray30"
  )

ggsave(
  plot = p2,
  filename = "presentations/session-01/_output/flat_prior_plot.png",
  width = 7,
  height = 4,
  units = "in",
  dpi = 150
)
# Create sequential Bayesian updating animation with visible observations
create_sequential_animation <- function() {
  # Simulate observed data: 11 successes, 16 failures
  observations <- c(rep(1, 11), rep(0, 16))
  observations <- sample(observations, replace = FALSE) # 1 = success, 0 = failure
  n_obs <- length(observations)

  # Starting with informative prior: Beta(4, 18)
  alpha_start <- 4
  beta_start <- 18

  # Create data for each observation
  all_data <- data.frame()
  all_obs_data <- data.frame()

  for (i in 0:n_obs) {
    if (i == 0) {
      # Initial prior
      alpha_current <- alpha_start
      beta_current <- beta_start
      obs_text <- "Starting Prior"
    } else {
      # Update based on observations up to point i
      successes <- sum(observations[1:i])
      failures <- i - successes
      alpha_current <- alpha_start + successes
      beta_current <- beta_start + failures

      last_obs <- if (observations[i] == 1) "Success" else "Failure"
      obs_text <- sprintf("After %d obs (last: %s)", i, last_obs)
    }

    # Calculate density for current posterior
    current_density <- dbeta(p, alpha_current, beta_current)

    frame_data <- data.frame(
      p = p,
      density = current_density,
      observation = i,
      alpha = alpha_current,
      beta = beta_current,
      frame_label = obs_text,
      successes = if (i == 0) 0 else sum(observations[1:i]),
      failures = if (i == 0) 0 else i - sum(observations[1:i])
    )

    all_data <- rbind(all_data, frame_data)

    # Create observation data for this frame
    if (i > 0) {
      obs_frame <- data.frame(
        obs_num = 1:i,
        obs_value = observations[1:i],
        obs_type = ifelse(observations[1:i] == 1, "Success", "Failure"),
        x_pos = seq(0.02, 0.98, length.out = max(i, 2))[1:i],
        y_pos = rep(max(current_density) * 1.05, i),
        observation = i
      )
      all_obs_data <- rbind(all_obs_data, obs_frame)
    }
  }

  # Only include observation data if we have observations
  if (nrow(all_obs_data) > 0) {
    p_anim <- ggplot(all_data, aes(x = p, y = density)) +
      geom_line(color = "#2E8B57") +
      geom_area(fill = "#2E8B57", alpha = 0.3) +
      geom_point(
        data = all_obs_data,
        aes(x = x_pos, y = y_pos, color = obs_type, shape = obs_type),
        size = 10, inherit.aes = FALSE
      ) +
      scale_color_manual(
        values = c("Success" = "blue", "Failure" = "forestgreen"),
        name = "Observations"
      ) +
      scale_shape_manual(
        values = c("Success" = 19, "Failure" = 4),
        name = "Observations"
      ) +
      labs(
        x = "Allele Frequency (p)",
        y = "Density",
        title = "Sequential Bayesian Updating"
      ) +
      theme_minimal(base_size = 30) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "top"
      ) +
      xlim(0, 1) +
      ylim(0, max(all_data$density) * 1.15) +
      transition_states(observation, transition_length = 1, state_length = 2) +
      ease_aes("cubic-in-out")
  } else {
    p_anim <- ggplot(all_data, aes(x = p, y = density)) +
      geom_line(color = "#2E8B57", size = 1.5) +
      geom_area(fill = "#2E8B57", alpha = 0.3) +
      labs(
        x = "Allele Frequency (p)",
        y = "Density",
        title = "Sequential Bayesian Updating: {closest_state$frame_label[1]}",
        subtitle = "Alpha = {closest_state$alpha[1]}, Beta = {closest_state$beta[1]}"
      ) +
      theme_minimal(base_size = 24) +
      theme(
        panel.grid.minor = element_blank()
      ) +
      xlim(0, 1) +
      ylim(0, max(all_data$density) * 1.15) +
      transition_states(observation, transition_length = 1, state_length = 5) +
      ease_aes("cubic-in-out")
  }

  return(p_anim)
}

# Create the sequential animation
sequential_anim <- create_sequential_animation()

# Render sequential animation
sequential_gif <- animate(sequential_anim,
  width = 1600, height = 900,
  fps = 10, duration = 30,
  renderer = gifski_renderer("sequential_bayesian_update.gif")
)
