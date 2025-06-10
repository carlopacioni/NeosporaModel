library(ggplot2)
library(dplyr)

# Colors
farm_colors <- c(
  "W" = "#1b9e77",
  "NC" = "#d95f02",
  "G" = "#7570b3",
  "PH" = "#e7298a",
  "NM" = "#66a61e"
)

# Dams data 
data_dam <- data.frame(
  Farm = c("W", "NC", "NC", "NC", "NC", "G", "PH", "PH", "NM", "NM"),
  SERO = c("S+", "S-", "S-", "S+", "S+", "S+", "S-", "S+", "S+", "S+"),
  PCR = c("P-", "P-", "P+", "P-", "P+", "P-", "P-", "P-", "P-", "P+"),
  Count = c(3, 78, 14, 1, 1, 21, 10, 78, 23, 1)
)

# Calves data 
data_calf <- data.frame(
  Farm = c("W", "G", "G", "PH", "NM", "NM"),
  SERO = c("S+", "S-", "S+", "S+", "S-", "S+"),
  PCR = c("P-", "P-", "P-", "P-", "P-", "P-"),
  Count = c(1, 9, 17, 39, 7, 9)
)

# Plot for Dams
plot_dam <- ggplot(data_dam, aes(x = Farm, y = Count, fill = Farm)) +
  geom_bar(stat = "identity") +
  facet_grid(PCR ~ SERO) +
  theme_minimal() +
  labs(
    title = "PCR vs SERO Results by Farm (Dams)",
    x = "Farm",
    y = "Number of Samples"
  ) +
  scale_fill_manual(values = farm_colors) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

# Plot for Calves
plot_calf <- ggplot(data_calf, aes(x = Farm, y = Count, fill = Farm)) +
  geom_bar(stat = "identity") +
  facet_grid(PCR ~ SERO) +
  theme_minimal() +
  labs(
    title = "PCR vs SERO Results by Farm (Calves)",
    x = "Farm",
    y = "Number of Samples"
  ) +
  scale_fill_manual(values = farm_colors) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

# Print
print(plot_dam)
print(plot_calf)
