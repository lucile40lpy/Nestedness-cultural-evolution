# ====== Cultural Evolution ======

## 1. Library import ====
library(permute)
library(lattice)
library(vegan)
library(progress)
library(parallel)


## 2. Cultural evolution simulation function ====
cultural_evolution <- function(
    num_agents, num_items, num_generations, 
    interaction_rate, mutation_rate, 
    set_quota, quota_ratio, 
    set_appeal, appeal_mean, appeal_std_dev,
    gateway_threshold, burn_in=0) 
  {
  
  # Initialize matrix storage list
  matrix_list <- list()
  
  ### A. Initialize binary matrix ----
  binary_matrix <- matrix(0, nrow = num_agents, ncol = num_items,
    dimnames = list(paste0("agent_", 1:num_agents), 
                    paste0("item_", 1:num_items)))

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
      
      # Store matrix every 10 generations
      if (gen %% 10 == 0) {
        matrix_list[[as.character(gen)]] <- binary_matrix
      }
    }
  }
  
  ### E. Return final results ----
  return(list(
    binary_matrix = binary_matrix,
    matrix_list = matrix_list
  ))
}

## 3. Functions ====

### A. Sort matrix ----
sort_matrix <- function(M, iter) {
  for(i in 1:iter) {
    M <- M[, order(colSums(M), decreasing = TRUE)]
    M <- M[order(rowSums(M)), ]
  }
  return(M)
}

### B. Compute correlation ----
compute_cor_coef <- function(mat) {
  avg_inventory <- apply(mat, 2, function(item_col) {
    mean(rowSums(mat)[item_col == 1])})
  
  item_stats <- data.frame(
    Prevalence = colSums(mat),
    AvgInventory = avg_inventory
  )
  return(cor(item_stats$Prevalence, item_stats$AvgInventory))
}

### C. Clean matrix ----
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

### D. Compute nodf p-value ----
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

## 4. Optimized simulation loop ====

### A. Parameters lists ----
matrix_size <- list(c(500, 200), c(50, 20))
num_gen <- c(30, 50, 100)
interact_rate <- c(0.001, 0.005, 0.05, 0.1)
mutat_rate <- c(0.001)
baselines <- c("r0", "r1", "curveball")


### B. Initialize empty results dataframe ----

# Define generation labels for columns
gen_labels <- paste0(seq(10, 100, 10))

# Create empty dataframe with dynamic columns
df_cols <- list(
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
  metric = character(),
  null_model = character(),
  p_value_final = numeric(),
  coeff_cor_final = numeric(),
  num_agents_final = integer(),
  num_items_final = integer()
)

# Add generation-specific columns
for (gen in gen_labels) {
  df_cols[[paste0("p_value_", gen)]] <- numeric()
  df_cols[[paste0("coeff_cor_", gen)]] <- numeric()
  df_cols[[paste0("num_agents_", gen)]] <- integer()
  df_cols[[paste0("num_items_", gen)]] <- integer()
}

df_results <- do.call(data.frame, c(df_cols, stringsAsFactors = FALSE))


### C. Progress bar ----
total_iter <- length(matrix_size) * length(interact_rate) * 
  length(mutat_rate) * length(num_gen) * length(baselines)

pb <- progress_bar$new(
 format = "[:bar] :percent | ETA: :eta | Elapsed: :elapsedfull",
 total = total_iter,
 clear = FALSE,
 width = 100
 )

### D. Simulation loops ----
metric = 'NODF'
for(a in matrix_size) {
  for(i in interact_rate) {
    for(m in mutat_rate) {
      for(g in num_gen) {
        
        ### E. Run simulation ----
        res <- cultural_evolution(
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
        
        ### F. Sort and clean matrices ----
        # Clean final matrix
        final_clean <- clean_matrix(res$binary_matrix)
        binary_matrix_clean <- sort_matrix(final_clean$matrix, 10)
        
        num_agents_final <- final_clean$num_agents
        num_items_final <- final_clean$num_items
        
        # Clean stored matrices
        stored_matrices <- list()
        stored_dims <- list()
        stored_generations <- names(res$matrix_list)
        
        for(gen_str in stored_generations) {
          clean_res <- clean_matrix(res$matrix_list[[gen_str]])
          stored_matrices[[gen_str]] <- sort_matrix(clean_res$matrix, 10)
          stored_dims[[gen_str]] <- list(
            agents = clean_res$num_agents,
            items = clean_res$num_items)
        }
        
        ### G. Compute coefficient of correlation ----
        # For final matrix
        cor_coef_final <- compute_cor_coef(binary_matrix_clean)
        
        # For stored matrices
        stored_cor_coefs <- list()
        for (gen_str in stored_generations) {
          stored_cor_coefs[[gen_str]] <- compute_cor_coef(stored_matrices[[gen_str]])
        }
        
        ### H. Compute nestedness ----
        for(b in baselines) {
          # For final matrix
          p_val_final <- compute_p_val_nodf(binary_matrix_clean, b)
          
          ### I. Create new row with all columns ----
          new_row <- data.frame(
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
            null_model = b,
            p_value_final = p_val_final,
            coeff_cor_final = cor_coef_final,
            num_agents_final = num_agents_final,
            num_items_final = num_items_final,
            stringsAsFactors = FALSE
          )
          
          ### J. Add generation-specific values ----
          for(gen in gen_labels) {
            # Initialize with NA
            new_row[[paste0("p_value_", gen)]] <- NA_real_
            new_row[[paste0("coeff_cor_", gen)]] <- NA_real_
            new_row[[paste0("num_agents_", gen)]] <- NA_integer_
            new_row[[paste0("num_items_", gen)]] <- NA_integer_
            
            # Fill if exists
            if(gen %in% stored_generations) {
              # Compute nestedness for this generation
              gen_mat <- stored_matrices[[gen]]
              p_val_gen <- compute_p_val_nodf(gen_mat, b)
              
              # Add values to row
              new_row[[paste0("p_value_", gen)]] <- p_val_gen
              new_row[[paste0("coeff_cor_", gen)]] <- stored_cor_coefs[[gen]]
              new_row[[paste0("num_agents_", gen)]] <- stored_dims[[gen]]$agents
              new_row[[paste0("num_items_", gen)]] <- stored_dims[[gen]]$items
            }
          }
          
          ### K. Append to results and export matrices ----
          df_results <- rbind(df_results, new_row)
          
          # Final matrix
          matrix_id <- paste0("matrices/",a[1],"_",a[2],"_",i,"_",m,"_",g, ".csv")
          write.csv(binary_matrix_clean, matrix_id)
          
          # Matrices for each stored generation
          for (gen_str in stored_generations) {
            gen_matrix <- stored_matrices[[gen_str]]
            matrix_id <- paste0("matrices/",a[1],"_",a[2],"_",i,"_",m,"_",g,"_",gen_str, ".csv")
            write.csv(gen_matrix, matrix_id)
          }
          
          pb$tick()
        }
      }
    }
  }
}

### L. Save results ----
write.csv(df_results, "nestedness_results.csv", row.names = FALSE)