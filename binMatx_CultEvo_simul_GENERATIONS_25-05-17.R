# ====== Cultural Evolution ======

## 1. Library import ====
library(tidyverse)
library(vegan)
library(permute)
library(lattice)
library(foreach)
library(doParallel)
library(parallel)


## 2. Cultural evolution simulation function ====
run_simulation <- function(
    num_agents, num_items, num_generations, 
    interaction_rate, mutation_rate, 
    set_quota, quota_ratio, 
    set_appeal, appeal_mean, appeal_std_dev,
    gateway_threshold, burn_in=0) 
  {
  
  ### A. Initialize binary matrix ----
  binary_matrix <- matrix(0, nrow = num_agents, ncol = num_items,
    dimnames = list(paste0("agent_", 1:num_agents), paste0("item_", 1:num_items)))

  # Assign quotas if enabled
  quotas <- NULL
  if (set_quota) {
    quota_mean <- quota_ratio*num_items
    std_dev_q <- quota_mean / 2
    
    quotas <- rnorm(num_agents, mean = quota_mean, sd = std_dev_q)
    quotas <- pmin(pmax(round(quotas), 1), num_items)
    names(quotas) <- paste0("agent_", 1:num_agents)
  }
  
  # Initialize appeal scores if enabled
  appeal_scores <- NULL
  if (set_appeal) {
    appeal_scores <- rnorm(num_items, appeal_mean, appeal_std_dev)
    appeal_scores <- pmax(pmin(appeal_scores, 0.99), 0.01)
    names(appeal_scores) <- paste0("item_", 1:num_items)
  }
  
  ### B. Burn-in period (mutations only) ----
  for (gen in 1:burn_in) {
    for (agent in 1:num_agents) {
      for (item in 1:num_items) {
        if (runif(1) < mutation_rate) {
          if (binary_matrix[agent, item] == 0) {
            if (item == 1) {
              binary_matrix[agent, item] <- 1
            } else {
              prev_sum <- sum(binary_matrix[agent, 1:(item - 1)]) + 1
              num_prev <- item - 1
              if (prev_sum >= (gateway_threshold * num_prev)) {
                binary_matrix[agent, item] <- 1
              }
            }
          } else {
            binary_matrix[agent, item] <- 0
          }
        }
      }
    }
  }
  
  ### C. Main simulation loop ----
  for (gen in (burn_in + 1):num_generations) {
    # Mutation process
    for (agent in 1:num_agents) {
      for (item in 1:num_items) {
        if (runif(1) < mutation_rate) {
          if (binary_matrix[agent, item] == 0) {
            if (item == 1) {
              binary_matrix[agent, item] <- 1
            } else {
              prev_sum <- sum(binary_matrix[agent, 1:(item - 1)]) + 1
              num_prev <- item - 1
              if (prev_sum >= (gateway_threshold * num_prev)) {
                binary_matrix[agent, item] <- 1
              }
            }
          } else {
            binary_matrix[agent, item] <- 0
          }
        }
      }
    }
    
    ### D. Pairwise interactions ----
    num_pairs <- ceiling(interaction_rate * num_agents^2)
    shuffled_agents <- sample(rownames(binary_matrix))
    
    for (pair in 1:floor(length(shuffled_agents) / 2)) {
      agent_A <- shuffled_agents[2 * pair - 1]
      agent_B <- shuffled_agents[2 * pair]
      
      if (any(is.na(c(agent_A, agent_B))) || 
          identical(agent_A, agent_B) || 
          all(binary_matrix[agent_A, ] == 0) && all(binary_matrix[agent_B, ] == 0) || 
          identical(binary_matrix[agent_A, ], binary_matrix[agent_B, ])) next
      
      items_A <- which(binary_matrix[agent_A, ] == 1 & binary_matrix[agent_B, ] == 0)
      items_B <- which(binary_matrix[agent_B, ] == 1 & binary_matrix[agent_A, ] == 0)
      
      # Determine interaction direction
      if (length(items_A) > 0 && length(items_B) == 0) {
        giver <- agent_A
        receiver <- agent_B
        items_to_transmit <- items_A
      } else if (length(items_B) > 0 && length(items_A) == 0) {
        giver <- agent_B
        receiver <- agent_A
        items_to_transmit <- items_B
      } else if (length(items_A) > 0 && length(items_B) > 0) {
        if (runif(1) < 0.5) {
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
      
      # Perform interaction
      if (length(items_to_transmit) == 1) {
        transmitted_item <- items_to_transmit[1]
      } else if (set_appeal) {
        item_weights <- appeal_scores[names(items_to_transmit)]
        item_weights <- item_weights / sum(item_weights)
        transmitted_item <- sample(items_to_transmit, 1, prob = item_weights)
      } else {
        transmitted_item <- sample(items_to_transmit, 1)
      }
      
      # Gateway check before interaction
      if (transmitted_item == 1) {
        binary_matrix[receiver, transmitted_item] <- 1
      } else {
        item_index <- transmitted_item
        prev_sum <- sum(binary_matrix[receiver, 1:(item_index - 1)]) + 1
        num_prev_items <- item_index - 1
        if (prev_sum >= (gateway_threshold * num_prev_items)) {
          binary_matrix[receiver, transmitted_item] <- 1
        }
      }
      
      # Enforce quotas
      if (set_quota) {
        for (agent in c(agent_A, agent_B)) {
          current_items <- which(binary_matrix[agent, ] == 1)
          quota_limit <- quotas[[agent]]
          while (length(current_items) > quota_limit) {
            if (length(current_items) == 1) {
              dropped <- current_items[1]
            } else if (set_appeal) {
              drop_weights <- 1 - appeal_scores[names(current_items)]
              drop_weights <- drop_weights / sum(drop_weights)
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
  
  ### E. Calculate final results ----
  item_prevalence <- colSums(binary_matrix)
  inventory_size <- rowSums(binary_matrix)
  
  return(list(
    binary_matrix = binary_matrix,
    item_prevalence = item_prevalence,
    inventory_size = inventory_size
  ))
}

### F. Sort function ----
sort_matrix <- function(M, iter) {
  for(i in 1:iter) {
    M <- M[, order(colSums(M), decreasing = TRUE)]
    M <- M[order(rowSums(M)), ]
  }
  return(M)
}


## 3. Optimized simulation loop ====

### A. Parameters lists ----
# Agents and items as list of pairs
matrix_size <- list(c(500, 200), c(50, 20))
num_gen <- c(25, 50, 100, 150, 200, 250, 300)
interact_rate <- c(0.001, 0.005, 0.01, 0.05, 0.1)
mutat_rate <- c(0.001, 0.005, 0.01, 0.05, 0.1)
baselines <- c("r00", "r0", "r1", "c0", "curveball")

### B. Parallelization Setup ----
cores <- detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

### C. Pre-calculate Parameter Combinations ----
param_grid <- expand.grid(
  matrix_size = matrix_size,
  num_gen = num_gen,
  interact_rate = interact_rate,
  mutat_rate = mutat_rate,
  baselines = baselines,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE  # Ajout crucial
)

### D. Vectorized Result Collection ----
results_list <- foreach(idx = 1:nrow(param_grid), .packages = c("vegan", "permute", "lattice"), .combine = rbind) %dopar% {
  p <- param_grid[idx, ]
  a <- p$matrix_size[[1]]
  
  # Conversion explicite en caractère
  current_baseline <- as.character(p$baselines)
  
  ### E. Optimized Simulation Call ----
  res <- run_simulation(
    num_agents = a[1],
    num_items = a[2],
    num_generations = p$num_gen,
    interaction_rate = p$interact_rate,
    mutation_rate = p$mutat_rate,
    set_quota = TRUE,
    quota_ratio = 0.5,
    set_appeal = FALSE,
    appeal_mean = NA,
    appeal_std_dev = 0.2,
    gateway_threshold = 0
  )
  
  ### F. Fast Matrix Processing ----
  binary_matrix <- res$binary_matrix
  row_sums <- rowSums(binary_matrix)  # Cache row sums
  
  # Calculate nestedness
  out_nodf <- oecosimu(binary_matrix, nestednodf, method = current_baseline, alt = "greater", nsimul = 100)# Reduced permutations
  p_val <- out_nodf$oecosimu$pval[3]
  
  # Calculate regression coefficient
  avg_inventory <- colMeans(matrix(row_sums, nrow = nrow(binary_matrix), ncol = ncol(binary_matrix)) * binary_matrix)
  coeff <- coef(lm(avg_inventory ~ colSums(binary_matrix)))[2]
  
  ### G. Return Compact Result ----
  data.frame(
    num_agents = a[1],
    num_items = a[2],
    num_generations = p$num_gen,
    interaction_rate = p$interact_rate,
    mutation_rate = p$mutat_rate,
    set_quota = TRUE,
    quota_ratio = 0.5,
    set_appeal = FALSE,
    appeal_mean = NA_real_,
    appeal_std_dev = 0.2,
    gateway_threshold = 0,
    burn_in = 0,
    metric = 'NODF',
    null_model = p$baselines,
    p_value = p_val,
    coeff_plot = coeff,
    stringsAsFactors = FALSE
  )
}

### H. Cleanup and Save ----
stopCluster(cl)
write.csv(results_list, "nestedness_results.csv", row.names = FALSE)
