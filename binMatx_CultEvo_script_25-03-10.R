
# ====== Cultural Evolution ======

# tests :
appeal_scores <- c(-0.2, 0.6, 1.1, 0.01, 0.99)
items <- c(1, 2, 3, 4, 5)
# After clamping:
appeal_scores <- pmax(pmin(appeal_scores, 0.99), 0.01)
print(appeal_scores)  # [0.01, 0.60, 0.99, 0.01, 0.99]

appeal_scores <- 1 - appeal_scores
appeal_scores
appeal_scores <- appeal_scores/sum(appeal_scores)
appeal_scores
ba <- sample(items, 1, prob = appeal_scores)
ba

## Library import ====

library(dplyr)
library(gridExtra)
library(ggplot2)
library(reshape2)

library(vegan)
library(permute)
library(lattice)


## Set parameters ====

num_agents <- 200
num_items <- 20
num_generations <- 50

transmission_rate <- 0.8 # between 0 and 1
mutation_rate <- 0.01

set_quota <- TRUE
quota_mean <- 5

set_appeal <- FALSE
appeal_mean <- 0.4  # between 0 and 1

burn_in <- 0


## Cultural evolution simulation ====

### Initialize binary matrix ----
binary_matrix <- matrix(0, nrow = num_agents, ncol = num_items,
                        dimnames = list(paste0("agent_", 1:num_agents),
                                        paste0("item_", 1:num_items)))

### Assign quotas if enabled ----
if(set_quota) {
  std_dev_q <- quota_mean/2
  quotas <- rnorm(num_agents, mean = quota_mean, sd = std_dev_q)
  quotas <- pmin(pmax(round(quotas), 1), num_items)
  names(quotas) <- paste0("agent_", 1:num_agents)
  print(quotas)
}

### Initialize appeal scores if enabled ----
if(set_appeal) {
  # Generate scores with normal distribution and clamp to (0,1)
  std_dev_a <- appeal_mean/2
  appeal_scores <- rnorm(num_items, mean = appeal_mean, sd = std_dev_a)
  appeal_scores <- pmax(pmin(appeal_scores, 0.99), 0.01)
  names(appeal_scores) <- paste0("item_", 1:num_items)
  print(appeal_scores)
}

### Burn-in period (mutations only) ----
for(gen in 1:burn_in) {
  
  # Mutation process
  for(agent in 1:num_agents) {
    # 0 to 1 mutation
    if(runif(1) < mutation_rate) {
      available_0 <- which(binary_matrix[agent, ] == 0)
      if(length(available_0) > 0) {
        binary_matrix[agent, sample(available_0, 1)] <- 1
      }
    }
    # 1 to 0 mutation
    if(runif(1) < mutation_rate) {
      available_1 <- which(binary_matrix[agent, ] == 1)
      if(length(available_1) > 0) {
        binary_matrix[agent, sample(available_1, 1)] <- 0
      }
    }
  }
}

### Main simulation loop ----
for(gen in (burn_in + 1):num_generations) {
  
  # Apply mutations
  for(agent in 1:num_agents) {
    # 0 to 1 mutation
    if(runif(1) < mutation_rate) {
      available_0 <- which(binary_matrix[agent, ] == 0)
      if(length(available_0) > 0) {
        binary_matrix[agent, sample(available_0, 1)] <- 1
      }
    }
    # 1 to 0 mutation
    if(runif(1) < mutation_rate) {
      available_1 <- which(binary_matrix[agent, ] == 1)
      if(length(available_1) > 0) {
        binary_matrix[agent, sample(available_1, 1)] <- 0
      }
    }
  }
  
  ### Pairwise interactions ----
  num_pairs <- ceiling(transmission_rate * num_agents^2)
  shuffled_agents <- sample(rownames(binary_matrix))
  
  for(pair in 1:num_pairs) {
    if(2*pair > length(shuffled_agents)) break
    
    agent_A <- shuffled_agents[2*pair - 1]
    agent_B <- shuffled_agents[2*pair]
    
    # Skip invalid pairs
    if(identical(agent_A, agent_B)) next #same agent
    if(all(binary_matrix[agent_A, ] == 0) && all(binary_matrix[agent_B, ] == 0)) next #empty collections
    if(identical(binary_matrix[agent_A, ], binary_matrix[agent_B, ])) next #same collection profiles
    
    # Find unique items
    items_A <- which(binary_matrix[agent_A, ] == 1 & binary_matrix[agent_B, ] == 0)
    items_B <- which(binary_matrix[agent_B, ] == 1 & binary_matrix[agent_A, ] == 0)
    
    # Determine transmission direction
    if(length(items_A) > 0 && length(items_B) == 0) {
      giver <- agent_A
      receiver <- agent_B
      items_to_transmit <- items_A
    } else if(length(items_B) > 0 && length(items_A) == 0) {
      giver <- agent_B
      receiver <- agent_A
      items_to_transmit <- items_B
    } else if(length(items_A) > 0 && length(items_B) > 0) {
      if(runif(1) < 0.5) {
        giver <- agent_A
        receiver <- agent_B
        items_to_transmit <- items_A
      } else {
        giver <- agent_B
        receiver <- agent_A
        items_to_transmit <- items_B
      }
    } else {
      next
    }

    ### Perform transmission ----
    if (length(items_to_transmit) == 1) {
      transmitted_item <- items_to_transmit[0]
    }
    # Take appeal into account if enabled
    else if(set_appeal) {
      # Get appeal weights
      item_weights <- appeal_scores[names(items_to_transmit)]
      # Normalize weights so that it sums to 1
      item_weights <- item_weights/sum(item_weights)
      transmitted_item <- sample(items_to_transmit, 1, prob = item_weights)
    } else {
      transmitted_item <- sample(items_to_transmit, 1)
    }
    binary_matrix[receiver, transmitted_item] <- 1
    
    ### Enforce quotas ----
    if(set_quota) {
      for(agent in c(agent_A, agent_B)) {
        current_items <- which(binary_matrix[agent, ] == 1)
        quota_limit <- quotas[[agent]]
        
        while(length(current_items) > quota_limit) {
          if(length(current_items) == 1) {  # Directly drop the only item
            dropped <- current_items[1]
          } else if(set_appeal) {
            drop_weights <- 1 - appeal_scores[names(current_items)]
            drop_weights <- drop_weights/sum(drop_weights)
            dropped <- sample(current_items, 1, prob = drop_weights)
          } else {
            dropped <- sample(current_items, 1)
          }
          binary_matrix[agent, dropped] <- 0
          current_items <- which(binary_matrix[agent, ] == 1)
        }
      }
    }
  }
}


## Calculate final results ####
item_prevalence <- colSums(binary_matrix)
binary_matrix_df <- as.data.frame(binary_matrix)
colnames(binary_matrix_df) <- paste0("Item_", 1:num_items)
rownames(binary_matrix_df) <- paste0("Agent_", 1:num_agents)


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
  colnames(df) <- c("Agent", "CulturalItem", "Value")  # Renamed to avoid conflict
  # Convert to factors with proper ordering
  df$Agent <- factor(df$Agent, levels = rownames(matrix_data))  # Invert levels for y-axis
  df$CulturalItem <- factor(df$CulturalItem, levels = colnames(matrix_data))
  
  # Plot with inverted y-axis
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


### Saving plots ----
# Create plot directory with timestamp
plot_dir <- paste0("plots_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
dir.create(plot_dir)

# Save combined plot
combined_plot <- grid.arrange(p1, p2, ncol = 2)
ggsave(file.path(plot_dir, "00_matrices.png"), 
       combined_plot, width = 12, height = 8, dpi = 300)


## Nestedness tests ####

### NODF ----
# r00
out_nodf_r00 <- oecosimu(binary_matrix, nestednodf, "r00", alt = "greater")
print(out_nodf_r00)
densityplot(permustats(out_nodf_r00), main = "NODF Permutations - r00")

png(file.path(plot_dir, paste0("01_NODF_r00.png")), width = 800, height = 600)
print(densityplot(permustats(out_nodf_r00), main = "NODF Permutations - r00"))
dev.off()

# r0
out_nodf_r0 <- oecosimu(binary_matrix, nestednodf, "r0", alt = "greater")
print(out_nodf_r0)
densityplot(permustats(out_nodf_r0), main = "NODF Permutations - r0")

png(file.path(plot_dir, paste0("01_NODF_r0.png")), width = 800, height = 600)
print(densityplot(permustats(out_nodf_r0), main = "NODF Permutations - r0"))
dev.off()

# r1
out_nodf_r1 <- oecosimu(binary_matrix, nestednodf, "r1", alt = "greater")
print(out_nodf_r1)
densityplot(permustats(out_nodf_r1), main = "NODF Permutations - r1")

png(file.path(plot_dir, paste0("01_NODF_r1.png")), width = 800, height = 600)
print(densityplot(permustats(out_nodf_r1), main = "NODF Permutations - r1"))
dev.off()

# c0
out_nodf_c0 <- oecosimu(binary_matrix, nestednodf, "c0", alt = "greater")
print(out_nodf_c0)
densityplot(permustats(out_nodf_c0), main = "NODF Permutations - c0")

png(file.path(plot_dir, paste0("01_NODF_c0.png")), width = 800, height = 600)
print(densityplot(permustats(out_nodf_c0), main = "NODF Permutations - c0"))
dev.off()

# swap
out_nodf_swap <- oecosimu(binary_matrix, nestednodf, "swap", alt = "greater")
print(out_nodf_swap)
densityplot(permustats(out_nodf_swap), main = "NODF Permutations - swap")

png(file.path(plot_dir, paste0("01_NODF_swap.png")), width = 800, height = 600)
print(densityplot(permustats(out_nodf_swap), main = "NODF Permutations - swap"))
dev.off()

# curveball
out_nodf_cb <- oecosimu(binary_matrix, nestednodf, "curveball", alt = "greater")
print(out_nodf_cb)
densityplot(permustats(out_nodf_cb), main = "NODF Permutations - curveball")

png(file.path(plot_dir, paste0("01_NODF_cb.png")), width = 800, height = 600)
print(densityplot(permustats(out_nodf_cb), main = "NODF Permutations - cb"))
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
