#' fastLDSC: fast (matrix-based) LD Score Regression
#'
#' @import data.table
#' @importFrom stats qnorm dnorm var cov na.omit
#' @importFrom Rcpp sourceCpp evalCpp
#' @useDynLib fastLDSC, .registration = TRUE
#'
#' @param trait_vec_1 Vector of sumstat file paths.
#' @param trait_vec_2 Optional vector of sumstat file paths. If NULL, symmetric mode.
#' @param trait_names_1 Optional display names for trait_vec_1.
#' @param trait_names_2 Optional display names for trait_vec_2.
#' @param overlap_mat Optional overlap proportion matrix (rows=trait1, cols=trait2).
#' @param ld LD reference directory (by_chr=TRUE) or file path (by_chr=FALSE).
#' @param wld Weight directory or file path.
#' @param ld_prefix Prefix for per-chromosome ldscore files (e.g., "eur_w_ld_chr").
#' @param w_prefix Prefix for per-chromosome weight files. Defaults to ld_prefix when NULL.
#' @param by_chr Logical. If TRUE, read per-chromosome files 1:chr.
#' @param sep_weights Logical. If TRUE, weights are separate files.
#' @param chr Max chromosome index.
#' @param n.blocks Jackknife blocks.
#' @param stand,select,chisq.max LDSC parameters.
#'
#' @return A list containing the following matrices (rows correspond to \code{trait_vec_1}, columns to \code{trait_vec_2}):
#' \item{S}{Genetic Covariance Matrix. Diagonal elements represent SNP heritability, off-diagonal elements represent genetic covariance.}
#' \item{I}{Intercept Matrix. Represents the LD Score regression intercept (confounding + sample overlap).}
#' \item{RhoE}{Environmental Covariance Matrix. Derived from the intercept and sample overlap.}
#' \item{SE_S}{Standard Error matrix for S.}
#' \item{SE_I}{Standard Error matrix for I.}
#' \item{SE_RhoE}{Standard Error matrix for RhoE.}
#' \item{trait_names_1}{Vector of trait names for rows.}
#' \item{trait_names_2}{Vector of trait names for columns.}
#'
#' @export
fast_ldsc <- function(trait_vec_1, trait_vec_2 = NULL, 
                      trait_names_1 = NULL, trait_names_2 = NULL,
                      sample.prev_1 = NULL, population.prev_1 = NULL,
                      sample.prev_2 = NULL, population.prev_2 = NULL,
                      overlap_mat = NULL,
                      ld, wld, 
                      ld_prefix = "", w_prefix = NULL, 
                      by_chr = TRUE, sep_weights = FALSE, chr = 22,
                      n.blocks = 200, stand = FALSE, select = FALSE, chisq.max = NA) {
  
  # ============================================================================
  # 0. Initialization & Validation
  # ============================================================================
  
  # Handle defaults for symmetric mode
  if (is.null(trait_vec_2)) {
    trait_vec_2 <- trait_vec_1
    sample.prev_2 <- sample.prev_1
    population.prev_2 <- population.prev_1
    if (is.null(trait_names_2)) trait_names_2 <- trait_names_1
  }
  
  # Handle names
  if (is.null(trait_names_1)) trait_names_1 <- basename(trait_vec_1)
  if (is.null(trait_names_2)) trait_names_2 <- basename(trait_vec_2)
  
  if (length(trait_names_1) != length(trait_vec_1)) stop("Length mismatch: trait_names_1 vs trait_vec_1")
  if (length(trait_names_2) != length(trait_vec_2)) stop("Length mismatch: trait_names_2 vs trait_vec_2")
  
  # Prevalences
  if (is.null(sample.prev_1)) sample.prev_1 <- rep(NA, length(trait_vec_1))
  if (is.null(population.prev_1)) population.prev_1 <- rep(NA, length(trait_vec_1))
  if (is.null(sample.prev_2)) sample.prev_2 <- rep(NA, length(trait_vec_2))
  if (is.null(population.prev_2)) population.prev_2 <- rep(NA, length(trait_vec_2))
  
  is_symmetric <- identical(trait_vec_1, trait_vec_2)
  if (is_symmetric) {
    message(">>> [Start] LDSC (Symmetric Mode: Calculating Upper Triangle Only)")
  } else {
    message(">>> [Start] LDSC (Rectangular Mode: Calculating Full Grid)")
  }
  
  # Mapping Unique Traits
  unique_traits <- unique(c(trait_vec_1, trait_vec_2))
  n.unique.traits <- length(unique_traits)
  
  df_map <- rbind(
    data.frame(trait = trait_vec_1, samp = sample.prev_1, pop = population.prev_1, stringsAsFactors = FALSE),
    data.frame(trait = trait_vec_2, samp = sample.prev_2, pop = population.prev_2, stringsAsFactors = FALSE)
  )
  df_map <- df_map[!duplicated(df_map$trait), ]
  map_idx <- match(unique_traits, df_map$trait)
  unique_samp_prev <- df_map$samp[map_idx]
  unique_pop_prev  <- df_map$pop[map_idx]
  
  map_vec1 <- match(trait_vec_1, unique_traits)
  map_vec2 <- match(trait_vec_2, unique_traits)
  
  # Initialize Result Matrices
  n_rows <- length(trait_vec_1)
  n_cols <- length(trait_vec_2)
  
  S_mat    <- matrix(NA, n_rows, n_cols)    
  I_mat    <- matrix(NA, n_rows, n_cols)        
  RhoE_mat <- matrix(NA, n_rows, n_cols)
  SE_S     <- matrix(NA, n_rows, n_cols)
  SE_I     <- matrix(NA, n_rows, n_cols)
  SE_RhoE  <- matrix(NA, n_rows, n_cols)
  
  dimnames(S_mat) <- dimnames(I_mat) <- dimnames(RhoE_mat) <- list(trait_names_1, trait_names_2)
  dimnames(SE_S) <- dimnames(SE_I) <- dimnames(SE_RhoE) <- list(trait_names_1, trait_names_2)
  
  # ============================================================================
  # Step 1: Read Reference Panel
  # ============================================================================
  message(">>> [Step 1] Loading Reference Panel...")
  t_step1 <- Sys.time()
  
  # Define file index based on selection
  if (isFALSE(select)) file_idx <- 1:chr
  else if (select == "ODD") file_idx <- seq(1, chr, 2)
  else if (select == "EVEN") file_idx <- seq(2, chr, 2)
  else file_idx <- select
  
  if (by_chr) {
    # --- Mode A: Split Files ---
    x <- rbindlist(lapply(file_idx, function(i) {
      f <- file.path(ld, paste0(ld_prefix, i, ".l2.ldscore.gz"))
      if (!file.exists(f)) stop(paste("LD file missing:", f))
      fread(f, showProgress = FALSE)
    }))
    
    if (sep_weights) {
      w <- rbindlist(lapply(file_idx, function(i) {
        f <- file.path(wld, paste0(w_prefix, i, ".l2.ldscore.gz"))
        if (!file.exists(f)) stop(paste("Weight file missing:", f))
        fread(f, showProgress = FALSE)
      }))
    } else { w <- copy(x) }
    
    # Read M files
    m_list <- rbindlist(lapply(file_idx, function(i) {
      f <- file.path(ld, paste0(ld_prefix, i, ".l2.M_5_50"))
      if (!file.exists(f)) stop(paste("M file missing:", f))
      fread(f, header = FALSE, showProgress = FALSE)
    }))
    m_scalar <- sum(m_list[[1]])
    
  } else {
    # --- Mode B: Single File ---
    if (!file.exists(ld)) stop(paste("LD file missing:", ld))
    x <- fread(ld, showProgress = FALSE)
    
    if (sep_weights) {
      if (!file.exists(wld)) stop(paste("Weight file missing:", wld))
      w <- fread(wld, showProgress = FALSE)
    } else { w <- copy(x) }
    
    # Infer M file
    m_file <- sub("\\.ldscore\\.gz$", ".M_5_50", ld)
    if (!file.exists(m_file)) m_file <- sub("\\.ldscore$", ".M_5_50", sub("\\.gz$", "", ld))
    
    if (file.exists(m_file)) {
      m_dt <- fread(m_file, header = FALSE, showProgress = FALSE)
      m_scalar <- sum(m_dt[[1]])
    } else {
      stop(sprintf("Error: M file not found. Expected: %s", m_file))
    }
  }
  
  setkey(x, SNP)
  if ("L2" %in% names(w)) setnames(w, "L2", "wLD")
  setkey(w, SNP)
  
  message(sprintf("    Step 1 completed in %.2fs", difftime(Sys.time(), t_step1, units = "secs")))
  
  # ============================================================================
  # Step 2: Common SNPs
  # ============================================================================
  message(">>> [Step 2] Finding Global Common SNPs...")
  t_step2 <- Sys.time()
  
  common_snps <- x$SNP
  
  for (i in 1:n.unique.traits) {
    # Read with column check to prevent "Object 'Z' not found"
    dt <- fread(unique_traits[i], showProgress = FALSE)
    
    req_cols <- c("SNP", "N", "Z", "A1")
    if (!all(req_cols %in% names(dt))) {
      stop(sprintf("Trait %d (%s) missing required columns. Found: %s", 
                   i, basename(unique_traits[i]), paste(names(dt), collapse=", ")))
    }
    
    dt <- dt[, ..req_cols] # Subset safely
    dt <- na.omit(dt)
    
    # Safe filtering using dt$ syntax or data.table scoping
    current_max_chisq <- if (is.na(chisq.max)) max(0.001 * max(dt$N), 80) else chisq.max
    dt <- dt[dt$Z^2 <= current_max_chisq]
    
    common_snps <- intersect(common_snps, dt$SNP)
    if (length(common_snps) == 0) stop(sprintf("Error: No common SNPs left after processing trait %d!", i))
  }
  
  if (length(common_snps) < 10000) stop("Error: Too few common SNPs (<10k)!")
  message(sprintf("    Final Analysis SNP count: %d", length(common_snps)))
  message(sprintf("    Step 2 completed in %.2fs", difftime(Sys.time(), t_step2, units = "secs")))
  
  # ============================================================================
  # Step 3: Matrix Construction
  # ============================================================================
  message(">>> [Step 3] Loading Data & Aligning...")
  t_step3 <- Sys.time()
  
  x_subset <- x[J(common_snps), nomatch = NULL]
  w_subset <- w[J(common_snps), nomatch = NULL]
  
  setorder(x_subset, CHR, BP)
  setorder(w_subset, CHR, BP)
  
  final_snps <- x_subset$SNP
  L2_vec <- x_subset$L2
  wLD_vec <- w_subset$wLD
  
  # Clean up large objects
  rm(x, w, x_subset, w_subset); gc()
  
  Z_mat <- matrix(NA, nrow = length(final_snps), ncol = n.unique.traits)
  N_mat <- matrix(NA, nrow = length(final_snps), ncol = n.unique.traits)
  Chi2_mat <- matrix(NA, nrow = length(final_snps), ncol = n.unique.traits)
  
  base_A1 <- NULL
  
  for (i in 1:n.unique.traits) {
    dt <- fread(unique_traits[i], select = c("SNP", "N", "Z", "A1"), showProgress = FALSE)
    dt <- na.omit(dt)
    
    # Recalculate chi-sq threshold logic for consistency
    current_max_chisq <- if (is.na(chisq.max)) max(0.001 * max(dt$N), 80) else chisq.max
    dt <- dt[dt$Z^2 <= current_max_chisq]
    
    setkey(dt, SNP)
    dt_aligned <- dt[J(final_snps), nomatch = NULL]
    if (nrow(dt_aligned) != length(final_snps)) dt_aligned <- unique(dt_aligned, by = "SNP")
    
    if (i == 1) {
      base_A1 <- dt_aligned$A1
      align_factor <- 1
    } else {
      align_factor <- ifelse(dt_aligned$A1 == base_A1, 1, -1)
    }
    
    Z_mat[, i] <- dt_aligned$Z * align_factor
    N_mat[, i] <- dt_aligned$N
    Chi2_mat[, i] <- dt_aligned$Z^2
  }
  rm(base_A1); gc()
  message(sprintf("    Step 3 completed in %.2fs", difftime(Sys.time(), t_step3, units = "secs")))
  
  # ============================================================================
  # Step 4: Weight Calculation
  # ============================================================================
  message(">>> [Step 4] Pre-calculating Weights...")
  t_step4 <- Sys.time()
  
  PreCalc <- vector("list", n.unique.traits)
  Mean_N_vec <- colMeans(N_mat)
  Liab_Factors <- rep(NA, n.unique.traits)
  
  for (i in 1:n.unique.traits) {
    cur_Chi2 <- Chi2_mat[, i]
    cur_N <- N_mat[, i]
    
    mean_chi <- mean(cur_Chi2)
    mean_ld_n <- mean(L2_vec * cur_N)
    
    tot_agg <- (m_scalar * (mean_chi - 1)) / mean_ld_n
    tot_agg <- min(max(tot_agg, 0), 1)
    
    c_val <- tot_agg * cur_N / m_scalar
    ld_c <- pmax(L2_vec, 1)
    w_ld_c <- pmax(wLD_vec, 1)
    
    het_w <- 1 / (2 * (1 + (c_val * ld_c))^2)
    oc_w <- 1 / w_ld_c
    
    w_final <- het_w * oc_w
    init_w <- sqrt(w_final)
    weights <- init_w / sum(init_w)
    
    X_weighted <- cbind(L2_vec * weights, weights)
    PreCalc[[i]] <- list(X = X_weighted, init_w = init_w, weights = weights)
    
    # Liab Factor
    K <- unique_pop_prev[i]; P <- unique_samp_prev[i]
    if (!is.na(K) && !is.na(P)) {
      thd <- qnorm(1 - K)
      zv <- dnorm(thd)
      Liab_Factors[i] <- (K^2 * (1 - K)^2) / (P * (1 - P) * zv^2)
    }
  }
  message(sprintf("    Step 4 completed in %.2fs", difftime(Sys.time(), t_step4, units = "secs")))
  
  # ============================================================================
  # Step 5: Optimized Pairwise Analysis
  # ============================================================================
  message(">>> [Step 5] Running Pairwise Analysis...")
  t_step5 <- Sys.time()
  
  # 5.1 Define Tasks
  tasks <- list()
  idx_counter <- 1
  
  for (r in 1:n_rows) {
    # If symmetric, only do upper triangle (c >= r), else full grid
    c_start <- if (is_symmetric) r else 1
    
    for (c in c_start:n_cols) {
      tasks[[idx_counter]] <- list(
        z1 = map_vec1[r], z2 = map_vec2[c],
        row = r, col = c
      )
      idx_counter <- idx_counter + 1
    }
  }
  message(sprintf("    Total number of pairs to analyze: %d", length(tasks)))
  
  # 5.2 Worker Function (Refactored for DRY)
  run_pair_worker <- function(task) {
    j <- task$z1
    k <- task$z2
    r_idx <- task$row
    c_idx <- task$col
    
    cur_overlap <- if (!is.null(overlap_mat)) overlap_mat[r_idx, c_idx] else 0.0
    
    # Setup Data based on pair type
    if (j == k) {
      # Heritability (Diagonal)
      X_mat <- PreCalc[[j]]$X
      y_vec <- Chi2_mat[, j] * PreCalc[[j]]$weights
      N_bar <- Mean_N_vec[j]
      eff_overlap <- 1 # Self-overlap is 100%
    } else {
      # Covariance (Off-diagonal)
      w_comb_num <- PreCalc[[j]]$init_w + PreCalc[[k]]$init_w
      w_cov <- w_comb_num / sum(w_comb_num)
      X_mat <- cbind(L2_vec * w_cov, w_cov)
      
      ZZ <- Z_mat[, j] * Z_mat[, k]
      y_vec <- ZZ * w_cov
      N_bar <- sqrt(Mean_N_vec[j] * Mean_N_vec[k])
      eff_overlap <- cur_overlap
    }
    
    scale_factor <- m_scalar / N_bar
    
    # Call C++ Engine
    res <- run_jackknife_engine(X_mat, as.matrix(y_vec), n.blocks, eff_overlap, scale_factor)
    
    # Parse Results
    beta_hat <- res$reg[1]
    alpha_hat <- res$reg[2]
    
    C <- stats::cov(res$pseudo) / n.blocks
    var_beta <- C[1, 1]
    var_alpha <- C[2, 2]
    
    S_hat <- beta_hat * scale_factor
    I_hat <- alpha_hat
    var_S <- (scale_factor^2) * var_beta
    var_I <- var_alpha
    
    rho_e_est <- res$rho_e
    var_rho_e <- stats::var(as.vector(res$pseudo_rho_e)) / n.blocks
    
    return(list(
      S = S_hat, I = I_hat, 
      SE_S = sqrt(var_S), SE_I = sqrt(var_I),
      RhoE = rho_e_est, SE_RhoE = sqrt(var_rho_e),
      row = r_idx, col = c_idx
    ))
  }
  
  # 5.3 Execute
  results_list <- lapply(tasks, run_pair_worker)
  
  # 5.4 Reconstruct Matrices
  for (res in results_list) {
    if (is.null(res) || inherits(res, "try-error")) next
    r <- res$row; c <- res$col
    
    # Fill (r, c)
    S_mat[r, c]    <- res$S
    I_mat[r, c]    <- res$I
    RhoE_mat[r, c] <- res$RhoE
    SE_S[r, c]     <- res$SE_S
    SE_I[r, c]     <- res$SE_I
    SE_RhoE[r, c]  <- res$SE_RhoE
    
    # Fill (c, r) if symmetric
    if (is_symmetric && r != c) {
      S_mat[c, r]    <- res$S
      I_mat[c, r]    <- res$I
      RhoE_mat[c, r] <- res$RhoE
      SE_S[c, r]     <- res$SE_S
      SE_I[c, r]     <- res$SE_I
      SE_RhoE[c, r]  <- res$SE_RhoE
    }
  }
  
  return(list(
    S = S_mat, I = I_mat, RhoE = RhoE_mat,
    SE_S = SE_S, SE_I = SE_I, SE_RhoE = SE_RhoE,
    trait_names_1 = trait_names_1,
    trait_names_2 = trait_names_2
  ))
}