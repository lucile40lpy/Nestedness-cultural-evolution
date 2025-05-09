# ====== Cultural Evolution ======

## 1. Library import ====
library(dplyr)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(progress)

library(vegan)
library(permute)
library(lattice)


## 3. Cultural evolution simulation function ====
run_simulation <- function(
    num_agents, num_items, num_generations, 
    interaction_rate, mutation_rate, 
    set_quota, quota_ratio, 
    set_appeal, appeal_mean, appeal_std_dev,
    gateway_threshold, burn_in=0) 
  {
  
  ### Initialize binary matrix ----
  binary_matrix <- matrix(0, nrow = num_agents, ncol = num_items,
    dimnames = list(paste0("agent_", 1:num_agents), paste0("item_", 1:num_items)))

  ### Assign quotas if enabled ----
  quotas <- NULL
  if (set_quota) {
    quota_mean <- quota_ratio*num_items
    std_dev_q <- quota_mean / 2
    
    quotas <- rnorm(num_agents, mean = quota_mean, sd = std_dev_q)
    quotas <- pmin(pmax(round(quotas), 1), num_items)
    names(quotas) <- paste0("agent_", 1:num_agents)
  }
  
  ### Initialize appeal scores if enabled ----
  appeal_scores <- NULL
  if (set_appeal) {
    appeal_scores <- rnorm(num_items, appeal_mean, appeal_std_dev)
    appeal_scores <- pmax(pmin(appeal_scores, 0.99), 0.01)
    names(appeal_scores) <- paste0("item_", 1:num_items)
  }
  
  ### Burn-in period (mutations only) ----
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
  
  ### Main simulation loop ----
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
    
    ### Pairwise interactions ----
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
      
      ### Perform interaction ----
      if (length(items_to_transmit) == 1) {
        transmitted_item <- items_to_transmit[1]
      } else if (set_appeal) {
        item_weights <- appeal_scores[names(items_to_transmit)]
        item_weights <- item_weights / sum(item_weights)
        transmitted_item <- sample(items_to_transmit, 1, prob = item_weights)
      } else {
        transmitted_item <- sample(items_to_transmit, 1)
      }
      
      # Gateway check before interaction (corrected from original code)
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
      
      ### Enforce quotas ----
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
  
  ### Calculate final results ----
  item_prevalence <- colSums(binary_matrix)
  inventory_size <- rowSums(binary_matrix)
  
  return(list(
    binary_matrix = binary_matrix,
    item_prevalence = item_prevalence,
    inventory_size = inventory_size
  ))
}

## Sort function ----
sort_matrix <- function(M, iter) {
  for(i in 1:iter) {
    M <- M[, order(colSums(M), decreasing = TRUE)]
    M <- M[order(rowSums(M)), ]
  }
  return(M)
}

## 4. Simulation loop ====


### Parameters lists ----
# Agents and items as list of pairs
matrix_size <- list(c(500, 200), c(50, 20))

num_gen <- c(25, 50, 100, 150, 200, 250, 300, 500)
interact_rate <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5)
mutat_rate <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5)


### Load existing CSV file ----
df_results <- read.csv("simulation_results_empty.csv", sep = ';', stringsAsFactors = FALSE)


### Progress bar ----
total_iter <- length(matrix_size) * length(interact_rate) * 
  length(mutat_rate) * length(num_gen)

pb <- progress_bar$new(
  format = "[:bar] :percent | ETA: :eta | Elapsed: :elapsedfull",
  total = total_iter,
  clear = FALSE,
  width = 100
)


print(start)
### Simulation loops ----
metric = 'NODF'
for(a in matrix_size) {
  for(i in interact_rate) {
    for(m in mutat_rate) {
      for(g in num_gen) {
        
        #### Run simulation ----
        res <- run_simulation(
          num_agents = a[1],
          num_items = a[2],
          num_generations = g,
          interaction_rate = i,
          mutation_rate = m,
          set_quota = TRUE,
          quota_ratio = 0.5,
          set_appeal = FALSE,
          appeal_mean = NA,
          appeal_std_dev = 0.2,
          gateway_threshold = 0
        )
        # Get results
        matrix = res$binary_matrix
        item_prevalence = res$item_prevalence
        inventory_size = res$inventory_size
        sorted_matrix = sort_matrix(matrix, 20)*
        
        
        #### Calculate nestedness p-value ----
        if (metric == 'NODF'){
          out_nodf <- oecosimu(matrix, nestednodf, 'r00', alt = "greater")
          p_val = out_nodf$oecosimu$pval[3]
        } else {
          out_temp <- oecosimu(binary_matrix, nestedtemp, "r00", alt = "less")
          p_val = out_temp$oecosimu$pval
        }
        
        #### Calculate plot coeff ----
        # Calculate average inventory size for agents possessing each item
        avg_inventory <- apply(binary_matrix, 2, function(item_col) {
          mean(rowSums(binary_matrix)[item_col == 1])
        })
        # Create data frame for plotting
        item_stats <- data.frame(
          Prevalence = item_prevalence,
          AvgInventory = avg_inventory
        )
        # Fit linear model
        lm_model <- lm(AvgInventory ~ Prevalence, data = item_stats)
        coeff <- coefficients(lm_model)[2]
        
        #### Create new row ----
        new_row <- list(
            num_agents = a[1],
            num_items = a[2],
            num_generations = g,
            interaction_rate = i,
            mutation_rate = m,
            set_quota = TRUE,
            quota_ratio = 0.5,
            set_appeal = FALSE,
            appeal_mean = NA_real_,
            appeal_std_dev = 0.2,
            gateway_threshold = 0,
            burn_in = 0,
            metric = 'NODF',
            null_model = 'r00',
            p_value = p_val,
            matrix = matrix,
            sorted_matrix = sorted_matrix,
            coeff_plot = coeff
            )
        
        # Append to results
        df_results <- rbind(df_results, as.data.frame(new_row))
      }
    }
  }
}

### Save Results ====
View(df_results)
write.csv(results_df, "nestedness_results.csv", row.names = FALSE)


