
# ====== Randomly distributed matrices ======

## Library import ====

library(dplyr)
library(gridExtra)
library(ggplot2)
library(reshape2)

library(vegan)
library(permute)
library(lattice)

## Parameters ####
nrows <- 200
ncols <- 20
area <- "both"
distr_type <- "powerlaw"
alpha <- 0.6


## Random matrices ####

if(distr_type == "uniform") {
  if(area == "agents") {
    probs <- runif(nrows, 0, 1)
    prob_vector <- rep(probs, each = ncols)
  } else if(area == "items") {
    probs <- runif(ncols, 0, 1)
    prob_vector <- rep(probs, times = nrows)
  } else {  # both
    prob_vector <- runif(nrows * ncols, 0, 1)
  }

} else if(distr_type == "powerlaw") {
  if(area == "agents") {
    u <- runif(nrows)
    probs <- u^(1/(alpha + 1))
    prob_vector <- rep(probs, each = ncols)
  } else if(area == "ietms") {
    u <- runif(ncols)
    probs <- u^(1/(alpha + 1))
    prob_vector <- rep(probs, times = nrows)
  } else {  # both
    u <- runif(nrows * ncols)
    prob_vector <- u^(1/(alpha + 1))
  }

} else if(distr_type == "gaussian") {
  if(area == "agents") {
    x <- rnorm(nrows, 0.5, 0.2)
  } else if(area == "items") {
    x <- rnorm(ncols, 0.5, 0.2)
  } else {  # both
    x <- rnorm(nrows * ncols, 0.5, 0.2)
  }
  prob_vector <- pmin(pmax(x, 0), 1)  # Clamp to [0,1]
}


## Matrix Generation ####
binary_matrix <- matrix(rbinom(nrows * ncols, size = 1, prob = prob_vector),
                        nrow = nrows, ncol = ncols,
                        dimnames = list(paste0("agent_", 1:num_agents),
                                        paste0("item_", 1:num_items)))

print("Final probabilities:")
print(prob_vector)
print("Binary matrix:")
print(binary_matrix)

item_prevalence <- colSums(binary_matrix)
print(item_prevalence)
inventory_size <- rowSums(binary_matrix)
print(inventory_size)

## Plotting Matrix ####

# Sorting function implementation
sort_matrix <- function(M, iter) {
  for(i in 1:iter) {
    M <- M[, order(colSums(M), decreasing = TRUE)]
    M <- M[order(rowSums(M)), ]
  }
  return(M)
}


# Matrix plotting function
plot_matrix <- function(matrix_data, title) {
  
  df <- reshape2::melt(matrix_data)
  colnames(df) <- c("Agent", "CulturalItem", "Value")
  df$Agent <- factor(df$Agent, levels = rownames(matrix_data))
  df$CulturalItem <- factor(df$CulturalItem, levels = colnames(matrix_data))
  
  ggplot(df, aes(x = CulturalItem, y = Agent, fill = factor(Value))) +
    geom_tile() +
    scale_fill_manual(values = c("0" = "cornsilk", "1" = "olivedrab")) +
    theme_minimal() +
    labs(title = title, x = "Items", y = "Agents") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
}

# Generate and plot matrices
ordered_matrix <- sort_matrix(binary_matrix, 50)
p1 <- plot_matrix(binary_matrix, "Original Matrix")
p2 <- plot_matrix(ordered_matrix, "Sorted Matrix")
grid.arrange(p1, p2, ncol = 2)



## Plot rarity, prevalence, frequency ####

### 1. Item Prevalence Distribution ----
prevalence_plot <- ggplot(data.frame(Prevalence = item_prevalence), 
                          aes(x = Prevalence)) +
  geom_histogram(fill = "olivedrab", color = "olivedrab",
                 alpha = 0.9, bins = 20) +
  labs(title = "Item Prevalence Distribution",
       x = "Number of agentss possessing item",
       y = "Count") +
  theme_minimal()
prevalence_plot


### 1. bis Item prevalence ----
# Create a data frame with item IDs and their prevalence
item_data <- data.frame(
  Item = factor(1:length(item_prevalence)),  # Convert to factor for discrete axis
  Prevalence = item_prevalence)

prevalence_barplot <- ggplot(item_data, 
                             aes(x = reorder(Item, -Prevalence), 
                                 y = Prevalence)) +
  geom_bar(stat = "identity", fill = "olivedrab", 
           alpha = 0.9, width = 0.8) +
  # Add text labels at the top of bars
  geom_text(aes(label = Prevalence), vjust = -0.5, 
            color = "olivedrab4", size = 3) +
  labs(title = "Item Prevalence Distribution",
       x = "Items",
       y = "Number of Agents Possessing Item") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))  # Add space for top labels
prevalence_barplot


### 2. Inventory size Distribution ----
inventory_plot <- ggplot(data.frame(Inventory = inventory_size), 
                         aes(x = Inventory)) +
  geom_histogram(fill = "olivedrab", color = "olivedrab", 
                 alpha = 0.9, bins = 14) +
  labs(title = "Inventory size Distribution",
       x = "Number of items possessed by agents",
       y = "Count") +
  theme_minimal()
inventory_plot


### 3. Average Rarity vs Inventory Size ----
# Calculate item rarity (inverse of prevalence normalized)
item_rarity <- 1 - (item_prevalence/max(item_prevalence))

# Calculate agent-level metrics
agent_stats <- data.frame(
  Agent = rownames(binary_matrix),
  InventorySize = rowSums(binary_matrix),
  AvgRarity = apply(binary_matrix, 1, function(x) {
    mean(item_rarity[which(x == 1)])
  }))

# Remove agents with empty inventories
rarity_plot <- ggplot(agent_stats, aes(x = InventorySize, y = AvgRarity)) +
  geom_point(color = "olivedrab", alpha = 0.7) +
  geom_smooth(method = "loess", color = "olivedrab", se = FALSE) +
  labs(title = "Average Item Rarity vs Inventory Size",
       x = "Number of items possessed (Inventory Size)",
       y = "Average Rarity Score") +
  theme_minimal()
rarity_plot


### 3. bis Prevalence vs Average Inventory Size ----
# Calculate average inventory size for agents possessing each item
avg_inventory <- apply(binary_matrix, 2, function(item_col) {
  mean(rowSums(binary_matrix)[item_col == 1])
})

# Create data frame for plotting
item_stats <- data.frame(
  Prevalence = item_prevalence,
  AvgInventory = avg_inventory
)

# Create the scatter plot
prevalence_inventory_plot <- ggplot(item_stats, aes(x = Prevalence, y = AvgInventory)) +
  geom_point(color = "olivedrab", alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", color = "gray50") +
  labs(title = "Item Prevalence vs Average Inventory Size",
       x = "Item prevalence :\nnumber of agents possessing the item",
       y = "Average inventory size\nof agents possessing the item") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

prevalence_inventory_plot



### Saving plots ----
# Create plot directory with timestamp
plot_dir <- paste0("plots_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
dir.create(plot_dir)

# Save combined plot
combined_plot <- grid.arrange(p1, p2, ncol = 2)
ggsave(file.path(plot_dir, "00_matrices.png"), 
       combined_plot, width = 12, height = 8, dpi = 300)
ggsave(file.path(plot_dir, "03_item_prevalence_bars.png"), 
       prevalence_barplot, width = 8, height = 6)
ggsave(file.path(plot_dir, "04_item_prevalence.png"), 
       prevalence_plot, width = 8, height = 6)
ggsave(file.path(plot_dir, "05_inventory_plot.png"), 
       inventory_plot, width = 8, height = 6)
ggsave(file.path(plot_dir, "06_rarity_vs_inventory.png"), 
       rarity_plot, width = 8, height = 6)
ggsave(file.path(plot_dir, "07_prevalence_vs_inventory.png"), 
       prevalence_inventory_plot, width = 8, height = 6)


## Nestedness tests ####

### NODF ----
# r00

out_nodf_r00 <- oecosimu(binary_matrix, nestednodf, "r00", alt = "greater")
print(out_nodf_r00)
densityplot(permustats(out_nodf_r00), main = "NODF Permutations - r00",
            layout = c(1, 3), aspect = "fill", as.table = TRUE)

png(file.path(plot_dir, "01_NODF_r00.png"), width = 400, height = 800)
print(densityplot(permustats(out_nodf_r00), main = "NODF Permutations - r00",
                  layout = c(1, 3), aspect = "fill", as.table = TRUE))
dev.off()


# r0
out_nodf_r0 <- oecosimu(binary_matrix, nestednodf, "r0", alt = "greater")
print(out_nodf_r0)
densityplot(permustats(out_nodf_r0), main = "NODF Permutations - r0",
            layout = c(1, 3), aspect = "fill", as.table = TRUE)

png(file.path(plot_dir, paste0("01_NODF_r0.png")), width = 400, height = 800)
print(densityplot(permustats(out_nodf_r0), main = "NODF Permutations - r0",
                  layout = c(1, 3), aspect = "fill", as.table = TRUE))
dev.off()

# r1
out_nodf_r1 <- oecosimu(binary_matrix, nestednodf, "r1", alt = "greater")
print(out_nodf_r1)
densityplot(permustats(out_nodf_r1), main = "NODF Permutations - r1",
            layout = c(1, 3), aspect = "fill", as.table = TRUE)

png(file.path(plot_dir, paste0("01_NODF_r1.png")), width = 400, height = 800)
print(densityplot(permustats(out_nodf_r1), main = "NODF Permutations - r1",
                  layout = c(1, 3), aspect = "fill", as.table = TRUE))
dev.off()

# c0
out_nodf_c0 <- oecosimu(binary_matrix, nestednodf, "c0", alt = "greater")
print(out_nodf_c0)
densityplot(permustats(out_nodf_c0), main = "NODF Permutations - c0",
            layout = c(1, 3), aspect = "fill", as.table = TRUE)

png(file.path(plot_dir, paste0("01_NODF_c0.png")), width = 400, height = 800)
print(densityplot(permustats(out_nodf_c0), main = "NODF Permutations - c0",
                  layout = c(1, 3), aspect = "fill", as.table = TRUE))
dev.off()

# swap
# sequential algorithm for binary matrices that changes 
#the matrix structure, but does not influence marginal sums 
#(Gotelli & Entsminger 2003). This inspects 2 × 2 submatrices so
#long that a swap can be done.
out_nodf_swap <- oecosimu(binary_matrix, nestednodf, "swap", alt = "greater")
print(out_nodf_swap)
densityplot(permustats(out_nodf_swap), main = "NODF Permutations - swap",
            layout = c(1, 3), aspect = "fill", as.table = TRUE)

png(file.path(plot_dir, paste0("01_NODF_swap.png")), width = 400, height = 800)
print(densityplot(permustats(out_nodf_swap), main = "NODF Permutations - swap",
                  layout = c(1, 3), aspect = "fill", as.table = TRUE))
dev.off()

# curveball
out_nodf_cb <- oecosimu(binary_matrix, nestednodf, "curveball", alt = "greater")
print(out_nodf_cb)
densityplot(permustats(out_nodf_cb), main = "NODF Permutations - curveball",
            layout = c(1, 3), aspect = "fill", as.table = TRUE)

png(file.path(plot_dir, paste0("01_NODF_cb.png")), width = 400, height = 800)
print(densityplot(permustats(out_nodf_cb), main = "NODF Permutations - cb",
                  layout = c(1, 3), aspect = "fill", as.table = TRUE))
dev.off()


### Temp ----
# r00
out_temp_r00 <- oecosimu(binary_matrix, nestedtemp, "r00", alt = "less")
print(out_temp_r00)
densityplot(permustats(out_temp_r00), main = "TEMP Permutations - r00")

png(file.path(plot_dir, paste0("02_Temp_r00.png")), width = 800, height = 600)
print(densityplot(permustats(out_temp_r00), main = "TEMP Permutations - r00"))
dev.off()

# r0
out_temp_r0 <- oecosimu(binary_matrix, nestedtemp, "r0", alt = "less")
print(out_temp_r0)
densityplot(permustats(out_temp_r0), main = "TEMP Permutations - r0")

png(file.path(plot_dir, paste0("02_Temp_r0.png")), width = 800, height = 600)
print(densityplot(permustats(out_temp_r0), main = "TEMP Permutations - r0"))
dev.off()

# r1
out_temp_r1 <- oecosimu(binary_matrix, nestedtemp, "r1", alt = "less")
print(out_temp_r1)
densityplot(permustats(out_temp_r1), main = "TEMP Permutations - r1")

png(file.path(plot_dir, paste0("02_Temp_r1.png")), width = 800, height = 600)
print(densityplot(permustats(out_temp_r1), main = "TEMP Permutations - r1"))
dev.off()

# c0
out_temp_c0 <- oecosimu(binary_matrix, nestedtemp, "c0", alt = "less")
print(out_temp_c0)
densityplot(permustats(out_temp_c0), main = "TEMP Permutations - c0")

png(file.path(plot_dir, paste0("02_Temp_c0.png")), width = 800, height = 600)
print(densityplot(permustats(out_temp_c0), main = "TEMP Permutations - c0"))
dev.off()

# swap
out_temp_swap <- oecosimu(binary_matrix, nestedtemp, "swap", alt = "less")
print(out_temp_swap)
densityplot(permustats(out_temp_swap), main = "TEMP Permutations - swap")

png(file.path(plot_dir, paste0("02_Temp_swap.png")), width = 800, height = 600)
print(densityplot(permustats(out_temp_swap), main = "TEMP Permutations - swap"))
dev.off()

# curveball
out_temp_cb <- oecosimu(binary_matrix, nestedtemp, "curveball", alt = "less")
print(out_temp_cb)
densityplot(permustats(out_temp_cb), main = "TEMP Permutations - curveball")

png(file.path(plot_dir, paste0("02_Temp_cb.png")), width = 800, height = 600)
print(densityplot(permustats(out_temp_cb), main = "TEMP Permutations - curveball"))
dev.off()
