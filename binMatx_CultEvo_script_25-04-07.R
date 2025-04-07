
# ====== Cultural Evolution ======


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
num_generations <- 500
burn_in <- 0 # between 0 and num_generations

transmission_rate <- 0.05 # between 0 and 1
mutation_rate <- 0.001 # between 0 and 1, has to be very low

set_quota <- TRUE
quota_mean <- num_items/2 # between 0 and num_items

set_appeal <- FALSE # appeal = selection
gateway_threshold <- 0.5 # between 0 and 1 : strength of gateway grows with the threshold


## Cultural evolution simulation ====

### Initialize binary matrix ----
binary_matrix <- matrix(0, nrow = num_agents, ncol = num_items,
                        dimnames = list(paste0("agent_", 1:num_agents),
                                        paste0("item_", 1:num_items)))

binary_matrices <- list()

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
  appeal_scores <- rnorm(num_items, mean = 0.5, sd = 0.2)
  appeal_scores <- pmax(pmin(appeal_scores, 0.99), 0.01)
  names(appeal_scores) <- paste0("item_", 1:num_items)
}


### Burn-in period (mutations only) ----
for(gen in 1:burn_in) {
  
  # Mutation process
  for(agent in 1:num_agents) {
    for(item in 1:num_items) {
      if (runif(1) < mutation_rate){
        if (binary_matrix[agent, item]==0){
          
          # If item_1 then no gateway
          if (item == 1) {
            binary_matrix[agent, item] <- 1
          # Mutation with gateway
          } else {
            prev_sum <- sum(binary_matrix[agent, 1:(item-1)]) + 1
            num_prev <- item - 1
            # Apply gateway threshold
            if(prev_sum >= (gateway_threshold * num_prev)) {
              binary_matrix[agent, item] <- 1
            }
          }
        } else {
          binary_matrix[agent, item] <- 0
        } } } } }


### Main simulation loop ----
for(gen in (burn_in + 1):num_generations) {
  
  # Mutation process
  for(agent in 1:num_agents) {
    for(item in 1:num_items) {
      if (runif(1) < mutation_rate){
        if (binary_matrix[agent, item]==0){
          
          # If item_1 then no gateway
          if (item == 1) {
            binary_matrix[agent, item] <- 1
            # Mutation with gateway
          } else {
            prev_sum <- sum(binary_matrix[agent, 1:(item-1)]) + 1
            num_prev <- item - 1
            # Apply gateway threshold
            if(prev_sum >= (gateway_threshold * num_prev)) {
              binary_matrix[agent, item] <- 1
            }
          }
        } else {
          binary_matrix[agent, item] <- 0
        } } } }
  
  ### Pairwise interactions ----
  num_pairs <- ceiling(transmission_rate * num_agents^2)
  shuffled_agents <- sample(rownames(binary_matrix))
  
  for(pair in 1:floor(length(shuffled_agents)/2)) {
    agent_A <- shuffled_agents[2*pair - 1]
    agent_B <- shuffled_agents[2*pair]

    # Skip invalid pairs
    if (any(is.na(c(agent_A, agent_B)))) next
       if(identical(agent_A, agent_B)) {
         next
       }
       if (all(binary_matrix[agent_A, ] == 0) && all(binary_matrix[agent_B, ] == 0)) {
         next
       }
       if (identical(binary_matrix[agent_A, ], binary_matrix[agent_B, ])) {
         next
       }
    
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
    # if there is only one item 
    if (length(items_to_transmit) == 1) {
      transmitted_item <- items_to_transmit[1]
      
    # Take appeal into account if enabled
    } else if(set_appeal) {
      # Get appeal weights
      item_weights <- appeal_scores[names(items_to_transmit)]
      # Normalize weights so that it sums to 1
      item_weights <- item_weights/sum(item_weights)
      transmitted_item <- sample(items_to_transmit, 1, prob = item_weights)
    } else {
      transmitted_item <- sample(items_to_transmit, 1)
    }

    # Gateway check before transmission
    # If item_1 then no gateway
    if (item == 1) {
      binary_matrix[agent, item] <- 1
    # Otherwise apply gateway
    } else {
      item_index <- which(colnames(binary_matrix) == paste0("item_", transmitted_item))
      previous_sum <- sum(binary_matrix[receiver, 1:(item_index-1)]) + 1
      num_prev_items <- item_index - 1
      
      if(previous_sum >= (gateway_threshold * num_prev_items)) {
        binary_matrix[receiver, transmitted_item] <- 1
        }
    }
    
    
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
        } } } }

  # Store the current matrix in the list
  binary_matrices[[gen]] <- binary_matrix
}


## Calculate final results ####
item_prevalence <- colSums(binary_matrix)
inventory_size <- rowSums(binary_matrix)
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


### 1. bis. Item prevalence ----
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
                 alpha = 0.9, bins = 20) +
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

# Rarity plot
rarity_plot <- ggplot(agent_stats, aes(x = InventorySize, y = AvgRarity)) +
  geom_point(color = "olivedrab", alpha = 0.7) +
  geom_smooth(method = "loess", color = "olivedrab", se = FALSE) +
  labs(title = "Average Item Rarity vs Inventory Size",
       x = "Number of items possessed (Inventory Size)",
       y = "Average Rarity Score") +
  theme_minimal()
rarity_plot


### 3.bis Prevalence vs Average Inventory Size ----
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



### Save plots ----
# Define your base directory
base_dir <- "C:/Users/lucil/OneDrive/Documents/CPES 2/Stage/my codes"
# Create plot directory with timestamp
plot_dir <- file.path(base_dir, paste0("plots_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
dir.create(plot_dir)



#Save all above plots
combined_plot <- grid.arrange(p1, p2, ncol = 2)
ggsave(file.path(plot_dir, "00_matrices.png"), 
       combined_plot, width = 12, height = 8, dpi = 300)
ggsave(file.path(plot_dir, "01_item_prevalence_bars.png"), 
       prevalence_barplot, width = 8, height = 6)
ggsave(file.path(plot_dir, "02_item_prevalence.png"), 
       prevalence_plot, width = 8, height = 6)
ggsave(file.path(plot_dir, "03_inventory_plot.png"), 
       inventory_plot, width = 8, height = 6)
ggsave(file.path(plot_dir, "04_rarity_vs_inventory.png"), 
       rarity_plot, width = 8, height = 6)
ggsave(file.path(plot_dir, "05_prevalence_vs_inventory.png"), 
       prevalence_inventory_plot, width = 8, height = 6)


## Nestedness tests ####

### NODF ----
# r00
out_nodf_r00 <- oecosimu(binary_matrix, nestednodf, "r00", alt = "greater")
print(out_nodf_r00)
densityplot(permustats(out_nodf_r00), main = "NODF Permutations - r00",
            layout = c(1, 3), aspect = "fill", as.table = TRUE)

png(file.path(plot_dir, paste0("06_NODF_r00.png")), width = 400, height = 800)
print(densityplot(permustats(out_nodf_r00), main = "NODF Permutations - r00",
                  layout = c(1, 3), aspect = "fill", as.table = TRUE))
dev.off()

# r0
out_nodf_r0 <- oecosimu(binary_matrix, nestednodf, "r0", alt = "greater")
print(out_nodf_r0)
densityplot(permustats(out_nodf_r0), main = "NODF Permutations - r0",
            layout = c(1, 3), aspect = "fill", as.table = TRUE)

png(file.path(plot_dir, paste0("06_NODF_r0.png")), width = 400, height = 800)
print(densityplot(permustats(out_nodf_r0), main = "NODF Permutations - r0",
                  layout = c(1, 3), aspect = "fill", as.table = TRUE))
dev.off()

# r1
out_nodf_r1 <- oecosimu(binary_matrix, nestednodf, "r1", alt = "greater")
print(out_nodf_r1)
densityplot(permustats(out_nodf_r1), main = "NODF Permutations - r1",
            layout = c(1, 3), aspect = "fill", as.table = TRUE)

png(file.path(plot_dir, paste0("06_NODF_r1.png")), width = 400, height = 800)
print(densityplot(permustats(out_nodf_r1), main = "NODF Permutations - r1",
                  layout = c(1, 3), aspect = "fill", as.table = TRUE))
dev.off()

# c0
out_nodf_c0 <- oecosimu(binary_matrix, nestednodf, "c0", alt = "greater")
print(out_nodf_c0)
densityplot(permustats(out_nodf_c0), main = "NODF Permutations - c0",
            layout = c(1, 3), aspect = "fill", as.table = TRUE)

png(file.path(plot_dir, paste0("06_NODF_c0.png")), width = 400, height = 800)
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

png(file.path(plot_dir, paste0("06_NODF_swap.png")), width = 400, height = 800)
print(densityplot(permustats(out_nodf_swap), main = "NODF Permutations - swap",
                  layout = c(1, 3), aspect = "fill", as.table = TRUE))
dev.off()

# curveball
out_nodf_cb <- oecosimu(binary_matrix, nestednodf, "curveball", alt = "greater")
print(out_nodf_cb)
densityplot(permustats(out_nodf_cb), main = "NODF Permutations - curveball",
            layout = c(1, 3), aspect = "fill", as.table = TRUE)

png(file.path(plot_dir, paste0("06_NODF_cb.png")), width = 400, height = 800)
print(densityplot(permustats(out_nodf_cb), main = "NODF Permutations - cb",
                  layout = c(1, 3), aspect = "fill", as.table = TRUE))
dev.off()


### Temp ----
# r00
out_temp_r00 <- oecosimu(binary_matrix, nestedtemp, "r00", alt = "less")
print(out_temp_r00)
densityplot(permustats(out_temp_r00), main = "TEMP Permutations - r00")

png(file.path(plot_dir, paste0("07_Temp_r00.png")), width = 800, height = 600)
print(densityplot(permustats(out_temp_r00), main = "TEMP Permutations - r00"))
dev.off()

# r0
out_temp_r0 <- oecosimu(binary_matrix, nestedtemp, "r0", alt = "less")
print(out_temp_r0)
densityplot(permustats(out_temp_r0), main = "TEMP Permutations - r0")

png(file.path(plot_dir, paste0("07_Temp_r0.png")), width = 800, height = 600)
print(densityplot(permustats(out_temp_r0), main = "TEMP Permutations - r0"))
dev.off()

# r1
out_temp_r1 <- oecosimu(binary_matrix, nestedtemp, "r1", alt = "less")
print(out_temp_r1)
densityplot(permustats(out_temp_r1), main = "TEMP Permutations - r1")

png(file.path(plot_dir, paste0("07_Temp_r1.png")), width = 800, height = 600)
print(densityplot(permustats(out_temp_r1), main = "TEMP Permutations - r1"))
dev.off()

# c0
out_temp_c0 <- oecosimu(binary_matrix, nestedtemp, "c0", alt = "less")
print(out_temp_c0)
densityplot(permustats(out_temp_c0), main = "TEMP Permutations - c0")

png(file.path(plot_dir, paste0("07_Temp_c0.png")), width = 800, height = 600)
print(densityplot(permustats(out_temp_c0), main = "TEMP Permutations - c0"))
dev.off()

# swap
out_temp_swap <- oecosimu(binary_matrix, nestedtemp, "swap", alt = "less")
print(out_temp_swap)
densityplot(permustats(out_temp_swap), main = "TEMP Permutations - swap")

png(file.path(plot_dir, paste0("07_Temp_swap.png")), width = 800, height = 600)
print(densityplot(permustats(out_temp_swap), main = "TEMP Permutations - swap"))
dev.off()

# curveball
out_temp_cb <- oecosimu(binary_matrix, nestedtemp, "curveball", alt = "less")
print(out_temp_cb)
densityplot(permustats(out_temp_cb), main = "TEMP Permutations - curveball")

png(file.path(plot_dir, paste0("07_Temp_cb.png")), width = 800, height = 600)
print(densityplot(permustats(out_temp_cb), main = "TEMP Permutations - curveball"))
dev.off()
