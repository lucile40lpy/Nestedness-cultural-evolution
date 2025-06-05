# ======== Cultural Evolution ========

## ==== 1. Library import ====
library(tidyverse)
library(progress)

library(vegan)
library(permute)
library(lattice)


## ==== 2. Cultural evolution simulation model ====
run_simulation <- function(
    num_agents, num_items, num_generations, 
    interaction_rate, mutation_rate, 
    set_quota, quota_ratio, 
    set_appeal, appeal_mean, appeal_std_dev,
    gateway_threshold, burn_in=0) 
  {
  
  ### ---- A. Initialize binary matrix ----
  binary_matrix <- matrix(0, nrow = num_agents, ncol = num_items,
    dimnames = list(paste0("agent_", 1:num_agents), paste0("item_", 1:num_items)))
  
  ### ---- B. Assign quotas ----
  quotas <- NULL
  if (set_quota) {
    quota_mean <- quota_ratio*num_items
    std_dev_q <- quota_mean / 2
    
    quotas <- rnorm(num_agents, mean = quota_mean, sd = std_dev_q)
    quotas <- pmin(pmax(round(quotas), 1), num_items)
    names(quotas) <- paste0("agent_", 1:num_agents)
  }
  
  ### ---- C. Initialize appeal scores ----
  appeal_scores <- NULL
  if (set_appeal) {
    appeal_scores <- rnorm(num_items, appeal_mean, appeal_std_dev)
    appeal_scores <- pmax(pmin(appeal_scores, 0.99), 0.01)
    names(appeal_scores) <- paste0("item_", 1:num_items)
  }
  
  ### ---- D. Burn-in period (mutations only) ----
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
  
  ### ---- E. Main simulation loop ----
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
    
    ### ---- F. Pairwise interactions ----
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
      
      ### ---- G. Perform interaction ----
      if (length(items_to_transmit) == 1) {
        transmitted_item <- items_to_transmit[1]
      } else if (set_appeal) {
        item_weights <- appeal_scores[names(items_to_transmit)]
        item_weights <- item_weights / sum(item_weights)
        transmitted_item <- sample(items_to_transmit, 1, prob = item_weights)
      } else {
        transmitted_item <- sample(items_to_transmit, 1)
      }
      
      ### ---- H. Gateway check ----
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
      
      ### ---- I. Enforce quotas ----
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
  
  ### ---- J. Final results ----
  return(list(
    binary_matrix = binary_matrix))
}


## ==== 3. Functions ====

### ---- A. Sort function ----
sort_matrix <- function(M, iter) {
  for(i in 1:iter) {
    M <- M[, order(colSums(M), decreasing = TRUE)]
    M <- M[order(rowSums(M)), ]
  }
  return(M)
}

### ---- B. Compute correlation ----
compute_cor_coef <- function(mat) {
  avg_inventory <- apply(mat, 2, function(item_col) {
    mean(rowSums(mat)[item_col == 1])})
  
  item_stats <- data.frame(
    Prevalence = colSums(mat),
    AvgInventory = avg_inventory
  )
  return(cor(item_stats$Prevalence, item_stats$AvgInventory))
}

### ---- C. Clean matrix ----
clean_matrix <- function(mat) {
  # Remove empty rows (agents with no items)
  non_empty_rows <- rowSums(mat) > 0
  cleaned <- mat[non_empty_rows, , drop = FALSE]
  
  # Remove empty columns (items with no adoptions)
  non_empty_cols <- colSums(cleaned) > 0
  cleaned <- cleaned[, non_empty_cols, drop = FALSE]
  
  # Return cleaned matrix and dimensions
  return(list(
    matrix = cleaned,
    num_agents = sum(non_empty_rows),
    num_items = sum(non_empty_cols)
  ))
}

### ---- D. Compute nodf p-value ----
options(mc.cores = max(2, parallel::detectCores() - 2))

compute_p_val_nodf <- function(mat, b) {
  tryCatch({
    if(nrow(mat) > 1 && ncol(mat) > 1) {
      out_nodf <- oecosimu(comm = mat, nestfun = nestednodf, 
                           method = b, alt = "greater", nsimul = 200,
                           batchsize = 50, parallel = TRUE)
      out_nodf$oecosimu$pval[3]
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
}

### ---- E. Compute temp p-value ----
compute_p_val_temp <- function(mat, b) {
  tryCatch({
    if(nrow(mat) > 1 && ncol(mat) > 1) {
      out_temp <- oecosimu(comm = mat, nestfun = nestedtemp, 
                           method = b, alt = "less", nsimul = 200,
                           batchsize = 50, parallel = TRUE)
      out_temp$oecosimu$pval
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
}


## ==== 4. Simulation loop ====

### ---- A. Parameters lists ----
matrix_size <- list(c(500, 200), c(50, 20))

interact_rate <- c(0.001, 0.005, 0.01, 0.05, 0.1)
mutat_rate <- c(0.001, 0.005, 0.01)

quota_ratio <- c(0.1, 0.5)
gateway_threshold <- c(0, 0.1, 0.5, 0.9, 1)
set_appeal <- c(FALSE, TRUE)
appeal_mean <- list(NA, c(0.1, 0.5, 0.9))  # NA when set_appeal=FALSE

baselines <- c('r0', 'r1', 'curveball')

### ---- B. Initialize empty results dataframe ----
df_results <- data.frame(
  num_agents = integer(),
  num_items = integer(),
  num_generations = integer(),
  interaction_rate = numeric(),
  mutation_rate = numeric(),
  set_quota = logical(),
  quota_ratio = numeric(),
  set_appeal = logical(),
  appeal_mean = numeric(),
  appeal_std_dev = numeric(),
  gateway_threshold = numeric(),
  burn_in = integer(),
  num_agents_final = integer(),
  num_items_final = integer(),
  coef_cor = numeric(),
  p_value_temp_r0 = numeric(),
  p_value_nodf_r0 = numeric(),
  p_value_temp_r1 = numeric(),
  p_value_nodf_r1 = numeric(),
  p_value_temp_curveball = numeric(),
  p_value_nodf_curveball = numeric(),
  stringsAsFactors = FALSE
)

# Ensure matrices directory exists
if (!dir.exists("matrices_final")) {
  dir.create("matrices_final")
}


### ---- C. Progress bar ----
total_iter <- length(matrix_size) * length(interact_rate) * length(mutat_rate) * 
  length(quota_ratio) * length(gateway_threshold) * length(set_appeal) * 
  length(baselines) * 3  # 3 appeal means when set_appeal=TRUE

pb <- progress_bar$new(
  format = "[:bar] :percent | ETA: :eta | Elapsed: :elapsedfull",
  total = total_iter,
  clear = FALSE,
  width = 100
)

# OR
#counter <- 0

### ---- D. Simulation loops ----
# general parameters
for (a in matrix_size) {
  for (i in interact_rate) {
    for (m in mutat_rate) {
      
      # parameters of interest
      for (quota in quota_ratio) {
        for (gateway in gateway_threshold) {
          for (appeal_flag in set_appeal) {
            appeal_means <- if(appeal_flag) c(0.1, 0.5, 0.9) else NA
            for (am in appeal_means) {
              
              ### ---- E. Run model ----
              res <- run_simulation(
                num_agents = a[1],
                num_items = a[2],
                num_generations = 50,
                interaction_rate = 0.05,
                mutation_rate = 0.01,
                set_quota = TRUE,
                quota_ratio = quota,
                set_appeal = appeal_flag,
                appeal_mean = if(appeal_flag) am else NA,
                appeal_std_dev = 0.2,
                gateway_threshold = gateway
              )
              
              ### ---- F. Sort and clean matrices ----
              final_clean <- clean_matrix(res$binary_matrix)
              binary_matrix_clean <- sort_matrix(final_clean$matrix, 10)
              
              num_agents_f <- final_clean$num_agents
              num_items_f <- final_clean$num_items
              
              ### ---- G. Coefficient of correlation ----
              cor_coef <- compute_cor_coef(binary_matrix_clean)
              
              ### ---- H. Create new row ----
              new_row <- data.frame(
                num_agents = a[1],
                num_items = a[2],
                num_generations = 50,
                interaction_rate = i,
                mutation_rate = m,
                quota_ratio = quota,
                gateway_threshold = gateway,
                set_quota = TRUE,
                set_appeal = appeal_flag,
                appeal_mean = ifelse(appeal_flag, am, NA),
                appeal_std_dev = 0.2,
                num_agents_final = num_agents_f,
                num_items_final = num_items_f,
                coef_cor = cor_coef, 
                stringsAsFactors = FALSE
              )
              
              ### ---- I. Nestedness test ----
              for (b in baselines) {
                # NODF
                p_val_nodf <- compute_p_val_nodf(binary_matrix_clean, b)
                new_row[[paste0("p_value_nodf_", b)]] <- p_val_nodf
                # Temp
                p_val_temp <- compute_p_val_temp(binary_matrix_clean, b)
                new_row[[paste0("p_value_temp_", b)]] <- p_val_temp
                
                # Progress bar
                pb$tick()
                # OR
                #counter <- counter + 1
                #print(counter)
              }
              
              ### ---- J. Append to results and export matrices ----
              df_results <- rbind(df_results, new_row)
              
              matrix_id <- paste0("matrices_final/",a[1],"_",a[2],"_",i,"_",m,"_",quota,"_",gateway,"_",appeal_flag,"_",am, ".csv")
              write.csv(binary_matrix_clean, matrix_id)
            }
          }
        }
      }
    }
  }
}

### ---- K. Save Results ----
write.csv2(df_results, "nestedness_results_fixed_gens.csv", row.names = FALSE)

