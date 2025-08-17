library(ggplot2)
library(dplyr)

n <- 27
x <- 11
# Create a sequence of possible p values
p_seq <- seq(0, 1, length.out = 200)
# Compute the binomial likelihood for each p
likelihood <- dbinom(x, n, p_seq)
# Plot the likelihood function
ggplot(data.frame(p = p_seq, likelihood = likelihood), aes(x = p, y = likelihood)) +
    geom_line(size = 1) +
    labs(
        x = "Probability of Success (p)",
        y = "Likelihood",
        title = "Binomial Likelihood Function"
    ) +
    geom_vline(xintercept = x / n, color = "red", linetype = "dashed") +
    theme_light(base_size = 18) +
    theme(panel.grid = element_blank())

ggsave("presentations/session-01/_output/likelihood_plot.png",
    width = 8,
    height = 6,
    dpi = 150,
    unit = "in"
)
