library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Define the range for x values
x_values <- seq(-5, 5, length.out = 10000)  # Use 10000 for smooth lines

# Calculate densities for each shape parameter
density1 <- dpgnorm(x_values, p = 1)
density2 <- dpgnorm(x_values, p = 2)
density10 <- dpgnorm(x_values, p = 10)

# Combine the density results into a data frame
df <- data.frame(
  x = rep(x_values, 3),
  density = c(density1, density2, density10),
  shape_parameter = factor(rep(c("1", "2", "10"), each = length(x_values)),
                           levels = c("1", "2", "10"))
)

# Use RColorBrewer palette
color_palette <- brewer.pal(3, "Set1")

# Plot with all density functions on the same plot
p <- ggplot(df, aes(x = x, y = density, color = shape_parameter)) +
  geom_line() +  # Increase line thickness for better visibility
  scale_color_manual(values = color_palette) +  # Apply custom color palette
  theme_bw() +
  labs(x = "Value", y = "Density", title = "Density Functions of the p-generalized Normal Distribution for Different Shape Parameters", color = "shape parameter") +
  theme(
    legend.position = "bottom",  # Move legend to bottom  # Remove legend title
    legend.title = element_text(size = 12),  # Set legend title size
    plot.title = element_text(size = 16),
    legend.text = element_text(size = 12),  # Adjust legend text size
    panel.border = element_blank(),  # Remove panel border
    axis.line = element_blank(),  # Remove axis lines
    axis.ticks = element_blank()  # Remove axis ticks
  )

# Print the plot
print(p)
