############################################################
# Reproduction script (current package API version, v4 robust ID fix)
# Remotely Close Associations:
# Openness to Experience and Semantic Memory Structure
#
# Data directory: D:/R/R_course/fuxian
# Main input files:
#   1. FINAL fluency.csv
#   2. FINAL open.csv
#
# Main output directory:
#   D:/R/R_course/fuxian/outputs
#
# Why this version exists:
# The authors' public code used older function names such as
# semnetcleaner(), autoDeStr(), cosine(), partboot(), and partboot.test().
# In current SemNetCleaner / SemNeT versions, several functions were renamed,
# moved, or replaced. This script uses the current API:
#   textcleaner(), similarity(), bootSemNeT(), test.bootSemNeT(), etc.
############################################################

############################
# 0. Basic settings
############################

rm(list = ls())

base_dir <- "D:/R/R_course/fuxian"
setwd(base_dir)

out_dir <- file.path(base_dir, "outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

set.seed(1234)

# For debugging, set FALSE first. For formal reproduction, set TRUE.
RUN_BOOTSTRAP <- FALSE

# For debugging, set FALSE first. Turn TRUE after the main pipeline runs.
RUN_SENSITIVITY <- FALSE

# The paper's bootstrap used 1000. For testing, temporarily use 50 or 100.
N_BOOT <- 1000

# The authors' code used 4 cores. This will be capped after packages load
# so the script does not request more cores than the local machine has.
N_CORES <- 4

# Mimic the authors' code where possible: the 50% condition used weighted = TRUE.
MATCH_SOURCE_BOOTSTRAP <- TRUE

############################
# 1. Install and load packages
############################

install_from_cran_or_runiverse <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(
      pkg,
      dependencies = TRUE,
      repos = c(
        "https://alexchristensen.r-universe.dev",
        "https://cloud.r-project.org"
      )
    )
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

install_from_cran_or_runiverse("psych")
install_from_cran_or_runiverse("SemNetCleaner")
install_from_cran_or_runiverse("SemNeT")

# NetworkToolbox is not strictly required for the current API route, but it is
# loaded for compatibility with the authors' original TMFG-based workflow.
if (!requireNamespace("NetworkToolbox", quietly = TRUE)) {
  install.packages("NetworkToolbox", dependencies = TRUE)
}
suppressPackageStartupMessages(library(NetworkToolbox))

needed_current_functions <- c(
  "textcleaner",
  "finalize", "equate", "similarity", "TMFG", "semnetmeas",
  "bootSemNeT", "test.bootSemNeT", "convert2cytoscape"
)

missing_current_functions <- needed_current_functions[
  !vapply(needed_current_functions, exists, logical(1), mode = "function")
]

if (length(missing_current_functions) > 0) {
  stop(
    "The following current-API functions were not found: ",
    paste(missing_current_functions, collapse = ", "),
    "\nPlease restart R and run this script from the first line. If the problem remains, run:\n",
    "install.packages('SemNetCleaner', repos = c('https://alexchristensen.r-universe.dev','https://cloud.r-project.org'))\n",
    "install.packages('SemNeT', repos = c('https://alexchristensen.r-universe.dev','https://cloud.r-project.org'))\n"
  )
}

############################
# 2. Utility functions
############################

as_numeric_binary_df <- function(x) {
  x <- as.data.frame(x, check.names = FALSE)
  x[] <- lapply(x, function(z) as.numeric(as.character(z)))
  x[is.na(x)] <- 0
  x
}

add_binary_columns_if_missing <- function(dat, cols) {
  for (cc in cols) {
    if (!(cc %in% colnames(dat))) dat[[cc]] <- 0
  }
  dat
}

get_tmfg_adjacency <- function(sim_mat) {
  sim_mat <- as.matrix(sim_mat)
  if (nrow(sim_mat) != ncol(sim_mat)) {
    stop("TMFG input is not square: ", nrow(sim_mat), " x ", ncol(sim_mat))
  }
  out <- SemNeT::TMFG(sim_mat)
  if (is.list(out) && "A" %in% names(out)) return(out$A)
  as.matrix(out)
}

# Current SemNeT uses similarity(data, method = "cosine").
# The authors' old code used cosine(data, addConstant = .01). The helper below
# uses the current cosine similarity and then adds the same small constant.
cosine_with_constant <- function(x, addConstant = 0.01) {
  x <- as_numeric_binary_df(x)
  if (nrow(x) < 2 || ncol(x) < 4) {
    stop("Not enough rows/columns to compute a stable cosine/TMFG network: ",
         nrow(x), " rows x ", ncol(x), " columns")
  }
  sim <- SemNeT::similarity(as.matrix(x), method = "cosine")
  sim <- as.matrix(sim)
  sim[!is.finite(sim)] <- 0
  if (!is.null(addConstant) && addConstant != 0) {
    sim <- sim + addConstant
    sim[sim > 1] <- 1
  }
  diag(sim) <- 1
  sim
}

semnet_measures <- function(A) {
  A <- as.matrix(A)
  fml <- names(formals(SemNeT::semnetmeas))
  if ("swm" %in% fml) {
    return(SemNeT::semnetmeas(A, swm = "rand"))
  }
  SemNeT::semnetmeas(A, weighted = FALSE)
}

extract_equated <- function(eq, i) {
  # equate() returns a list. Names vary across versions and depend on object names.
  if (!is.list(eq) || length(eq) < i || is.null(eq[[i]])) {
    stop("equate() did not return element ", i, ". Returned names: ", paste(names(eq), collapse = ", "))
  }
  as_numeric_binary_df(eq[[i]])
}

safe_alpha <- function(dat, scale_name) {
  dat <- as.data.frame(dat)
  dat[] <- lapply(dat, function(z) suppressWarnings(as.numeric(as.character(z))))
  if (any(vapply(dat, function(z) all(is.na(z)), logical(1)))) {
    warning("Skipping alpha for ", scale_name, ": at least one item column is entirely NA after numeric conversion.")
    return(data.frame(scale = scale_name, raw_alpha = NA_real_, std_alpha = NA_real_, stringsAsFactors = FALSE))
  }
  out <- suppressWarnings(psych::alpha(dat))
  data.frame(
    scale = scale_name,
    raw_alpha = unname(out$total[["raw_alpha"]]),
    std_alpha = unname(out$total[["std.alpha"]]),
    stringsAsFactors = FALSE
  )
}

############################
# 3. Read data
############################

fluency_file <- file.path(base_dir, "FINAL fluency.csv")
open_file    <- file.path(base_dir, "FINAL open.csv")

if (!file.exists(fluency_file)) stop("Cannot find file: ", fluency_file)
if (!file.exists(open_file))    stop("Cannot find file: ", open_file)

raw <- read.csv(
  fluency_file,
  as.is = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

latent <- read.csv(
  open_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(raw) < 2) {
  stop("FINAL fluency.csv has fewer than 2 columns. Check delimiter and file format.")
}
if (nrow(raw) == 0 || nrow(latent) == 0) {
  stop("One of the input files has 0 rows. Check file paths and CSV format.")
}

# The authors' public code used raw[, 1] as participant ID. Some CSV files
# have a blank or nonstandard first column name, so indexing by column name can
# return NULL. Use position-based indexing for robustness.
raw_id <- raw[[1]]
if (length(raw_id) != nrow(raw) || is.null(raw_id)) {
  warning("The first column in FINAL fluency.csv could not be used as an ID. Using row numbers as IDs.")
  raw_id <- seq_len(nrow(raw))
}
raw_id <- as.character(raw_id)
if (anyDuplicated(raw_id)) {
  warning("Participant IDs in FINAL fluency.csv are duplicated. Row order will be used; row names will be made unique for textcleaner.")
}

latent_id_col <- if ("id" %in% names(latent)) "id" else names(latent)[1]
latent_id <- latent[[latent_id_col]]

if (!("no_int" %in% names(latent))) {
  stop("FINAL open.csv does not contain no_int. The authors used latent$no_int for grouping.")
}

latent[["no_int"]] <- suppressWarnings(as.numeric(as.character(latent[["no_int"]])))
if (anyNA(latent[["no_int"]])) {
  stop("FINAL open.csv no_int contains NA or non-numeric values after conversion. Fix no_int before grouping.")
}

if (nrow(raw) != nrow(latent)) {
  stop("FINAL fluency.csv and FINAL open.csv have different row counts.")
}

# Reorder latent file only when the ID columns are valid and clearly match.
# If not, keep the original row order, which is how the authors' code combined files.
if (!is.null(latent_id) && length(latent_id) == length(raw_id) && !all(raw_id == as.character(latent_id))) {
  latent_id_chr <- as.character(latent_id)
  if (!anyDuplicated(raw_id) && !anyDuplicated(latent_id_chr) && setequal(raw_id, latent_id_chr)) {
    message("ID order differs. Reordering FINAL open.csv to match FINAL fluency.csv.")
    latent <- latent[match(raw_id, latent_id_chr), ]
  } else {
    warning("ID columns do not match cleanly or contain duplicates. Continuing by row order, consistent with the authors' cbind workflow.")
  }
}

message("Data loaded:")
message("  FINAL fluency.csv: ", nrow(raw), " rows x ", ncol(raw), " columns")
message("  FINAL open.csv:    ", nrow(latent), " rows x ", ncol(latent), " columns")

############################
# 4. Reliability checks for openness scales
############################

alpha_results <- list()
if (ncol(latent) >= 41) {
  alpha_results <- rbind(
    safe_alpha(latent[, 10:19], "BFAS Openness"),
    safe_alpha(latent[, 20:29], "BFAS Intellect"),
    safe_alpha(latent[, 30:41], "NEO Openness")
  )
  write.csv(alpha_results, file.path(out_dir, "00_alpha_results.csv"), row.names = FALSE)
}

############################
# 5. Clean verbal fluency data with current SemNetCleaner API
############################

# Current textcleaner() expects participant IDs as row names rather than as a data column.
# The fluency response columns are all columns except the first ID column.
fluency_responses <- raw[, -1, drop = FALSE]
if (ncol(fluency_responses) == 0) {
  stop("fluency_responses has 0 columns. FINAL fluency.csv was probably read with the wrong delimiter.")
}
fluency_responses[] <- lapply(fluency_responses, function(z) {
  z <- trimws(as.character(z))
  z[z == ""] <- "99"
  z
})
rownames(fluency_responses) <- make.unique(raw_id)

run_textcleaner <- function(dictionary_value) {
  args <- list(
    data = fluency_responses,
    miss = 99,
    partBY = "row",
    dictionary = dictionary_value,
    spelling = "US",
    keepStrings = FALSE,
    allowPunctuations = "-",
    allowNumbers = FALSE,
    lowercase = TRUE
  )
  fml <- names(formals(SemNetCleaner::textcleaner))
  if ("type" %in% fml) args$type <- "fluency"
  args <- args[names(args) %in% fml]
  do.call(SemNetCleaner::textcleaner, args)
}

# The animal dictionary is appropriate for the animal verbal fluency task.
# If dictionary loading fails in a local installation, the script falls back to NULL.
cleaned <- tryCatch(
  run_textcleaner("animals"),
  error = function(e) {
    message("textcleaner(..., dictionary = 'animals') failed. Retrying with dictionary = NULL.")
    message("Original error: ", conditionMessage(e))
    run_textcleaner(NULL)
  }
)

resp <- cleaned$responses$clean

cat("names(cleaned$responses):\n")
print(names(cleaned$responses))

if ("binary" %in% names(cleaned$responses)) {
  bin <- cleaned$responses$binary
} else if ("corrected" %in% names(cleaned$responses)) {
  bin <- SemNetCleaner::resp2bin(cleaned$responses$corrected)
} else {
  bin <- SemNetCleaner::resp2bin(resp)
}

cat("dim(resp): "); print(dim(resp))
cat("dim(bin): "); print(dim(bin))

if (is.null(bin) || nrow(bin) == 0 || ncol(bin) == 0) {
  stop("Could not obtain a valid binary matrix from textcleaner output.")
}
if (nrow(bin) != nrow(raw)) {
  stop("Binary matrix and FINAL fluency.csv have different row counts: ", nrow(bin), " vs ", nrow(raw))
}
rownames(bin) <- rownames(fluency_responses)

saveRDS(cleaned, file.path(out_dir, "01_cleaned_textcleaner_object.rds"))
write.csv(as.data.frame(resp), file.path(out_dir, "01_cleaned_responses.csv"), row.names = FALSE)
write.csv(as.data.frame(bin),  file.path(out_dir, "01_binary_textcleaner.csv"), row.names = FALSE)

con <- as_numeric_binary_df(bin)

############################
# 6. Manual corrections from the authors' public code
############################

# 6.1 Split catefrog into cat + frog, then remove catefrog.
if ("catefrog" %in% colnames(con)) {
  catefrog_rows <- which(con[["catefrog"]] != 0)
  con <- add_binary_columns_if_missing(con, c("cat", "frog"))
  con[catefrog_rows, c("cat", "frog")] <- 1
  con <- con[, colnames(con) != "catefrog", drop = FALSE]
  message("Handled catefrog. Rows affected: ", paste(catefrog_rows, collapse = ", "))
} else {
  message("No catefrog column found. Skipping catefrog correction.")
}

# 6.2 Manual completion for row 386.
animals_386 <- c(
  "cat", "dog", "fish", "elephant", "tiger", "zebra",
  "monkey", "giraffe", "lion", "dolphin", "chicken", "squirrel"
)
if (nrow(con) >= 386) {
  con <- add_binary_columns_if_missing(con, animals_386)
  con[386, animals_386] <- 1
}

# 6.3 Manual completion for row 499.
animals_499 <- c(
  "dog", "cat", "mouse", "moose", "horse",
  "lion", "tiger", "bear", "deer", "pig", "cow"
)
if (nrow(con) >= 499) {
  con <- add_binary_columns_if_missing(con, animals_499)
  con[499, animals_499] <- 1
}

zero_cases <- which(rowSums(con) == 0)
write.csv(
  data.frame(zero_case_row = zero_cases),
  file.path(out_dir, "03_zero_response_cases_after_manual_fix.csv"),
  row.names = FALSE
)
if (length(zero_cases) > 0) {
  warning("After manual corrections, zero-response cases remain: ", paste(zero_cases, collapse = ", "))
}

finalFull <- con
saveRDS(finalFull, file.path(out_dir, "03_final_cleaning_file_currentAPI.rds"))
write.csv(finalFull, file.path(out_dir, "03_final_cleaning_file_currentAPI.csv"), row.names = FALSE)

############################
# 7. Combine latent variable with binary responses and split groups
############################

if (nrow(finalFull) != nrow(latent)) {
  stop("Cleaned binary matrix and FINAL open.csv have different row counts: ",
       nrow(finalFull), " vs ", nrow(latent),
       ". Check whether textcleaner removed or duplicated any rows.")
}

comb <- data.frame(
  latent = latent[["no_int"]],
  id = raw_id,
  finalFull,
  check.names = FALSE
)
if (anyNA(comb$latent) || any(!is.finite(comb$latent))) {
  stop("Grouping variable latent/no_int contains NA or non-finite values.")
}

write.csv(comb, file.path(out_dir, "04_combined_latent_and_binary_responses.csv"), row.names = FALSE)

n_total <- nrow(comb)
if (n_total %% 2 != 0) stop("Sample size is odd; cannot split into equal high/low groups.")
if (n_total != 516) warning("Current n is not 516. The script will still split into equal halves. n = ", n_total)

n_half <- n_total / 2
comb_sorted <- comb[order(comb$latent), ]

low_group  <- comb_sorted[1:n_half, ]
high_group <- comb_sorted[(n_half + 1):n_total, ]

deLow  <- low_group[,  -c(1, 2), drop = FALSE]
deHigh <- high_group[, -c(1, 2), drop = FALSE]

write.csv(low_group,  file.path(out_dir, "05_low_group_full.csv"), row.names = FALSE)
write.csv(high_group, file.path(out_dir, "05_high_group_full.csv"), row.names = FALSE)
write.csv(deLow,      file.path(out_dir, "05_low_group_binary_only.csv"), row.names = FALSE)
write.csv(deHigh,     file.path(out_dir, "05_high_group_binary_only.csv"), row.names = FALSE)

message("Groups created:")
message("  Low openness group:  ", nrow(low_group))
message("  High openness group: ", nrow(high_group))

############################
# 8. Behavioral analyses
############################

sumAll <- rowSums(comb[, -c(1, 2), drop = FALSE])
cor_total_response <- cor.test(sumAll, comb$latent)

sumLow  <- rowSums(deLow)
sumHigh <- rowSums(deHigh)
ttest_total_response <- t.test(sumHigh, sumLow, var.equal = TRUE)

pooled_sd <- sqrt(
  ((length(sumHigh) - 1) * var(sumHigh) + (length(sumLow) - 1) * var(sumLow)) /
    (length(sumHigh) + length(sumLow) - 2)
)
cohen_d_high_minus_low <- (mean(sumHigh) - mean(sumLow)) / pooled_sd

onlyH <- deHigh[, colSums(deHigh) >= 1, drop = FALSE]
onlyL <- deLow[,  colSums(deLow)  >= 1, drop = FALSE]

uniH <- colnames(onlyH)
uniL <- colnames(onlyL)
uniT <- unique(c(uniH, uniL))

chitest <- data.frame(
  high = as.integer(uniT %in% uniH),
  low  = as.integer(uniT %in% uniL)
)

mcnemar_unique <- mcnemar.test(chitest$high, chitest$low)
mcnemar_phi <- sqrt(unname(mcnemar_unique$statistic) / length(uniT))

behavior_summary <- data.frame(
  statistic = c(
    "N_total", "N_low", "N_high",
    "unique_total", "unique_high", "unique_low",
    "unique_high_only", "unique_low_only",
    "mean_total_responses_high", "mean_total_responses_low",
    "cor_total_responses_with_latent_r", "cor_total_responses_with_latent_p",
    "t_total_responses", "t_total_responses_df", "t_total_responses_p",
    "cohen_d_high_minus_low",
    "mcnemar_chisq", "mcnemar_p", "mcnemar_phi"
  ),
  value = c(
    n_total, nrow(low_group), nrow(high_group),
    length(uniT), length(uniH), length(uniL),
    length(setdiff(uniH, uniL)), length(setdiff(uniL, uniH)),
    mean(sumHigh), mean(sumLow),
    unname(cor_total_response$estimate), cor_total_response$p.value,
    unname(ttest_total_response$statistic), unname(ttest_total_response$parameter),
    ttest_total_response$p.value,
    cohen_d_high_minus_low,
    unname(mcnemar_unique$statistic), mcnemar_unique$p.value, mcnemar_phi
  )
)

write.csv(behavior_summary, file.path(out_dir, "06_behavior_summary.csv"), row.names = FALSE)
write.csv(
  data.frame(response = uniT, high = chitest$high, low = chitest$low),
  file.path(out_dir, "06_unique_response_presence_by_group.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    cor_total_response = cor_total_response,
    ttest_total_response = ttest_total_response,
    mcnemar_unique = mcnemar_unique,
    behavior_summary = behavior_summary
  ),
  file.path(out_dir, "06_behavior_tests.rds")
)

############################
# 9. Semantic network analysis
############################

finLow  <- SemNeT::finalize(deLow, minCase = 2)
finHigh <- SemNeT::finalize(deHigh, minCase = 2)

eq <- SemNeT::equate(finLow, finHigh)

low_net_input  <- extract_equated(eq, 1)
high_net_input <- extract_equated(eq, 2)

saveRDS(eq, file.path(out_dir, "07_equated_network_input_object.rds"))
write.csv(low_net_input,  file.path(out_dir, "07_low_equated_network_input.csv"), row.names = FALSE)
write.csv(high_net_input, file.path(out_dir, "07_high_equated_network_input.csv"), row.names = FALSE)

cosLow  <- cosine_with_constant(low_net_input,  addConstant = .01)
cosHigh <- cosine_with_constant(high_net_input, addConstant = .01)

write.csv(cosLow,  file.path(out_dir, "08_cosine_low.csv"))
write.csv(cosHigh, file.path(out_dir, "08_cosine_high.csv"))

netLow  <- get_tmfg_adjacency(cosLow)
netHigh <- get_tmfg_adjacency(cosHigh)

write.csv(netLow,  file.path(out_dir, "09_TMFG_network_low.csv"))
write.csv(netHigh, file.path(out_dir, "09_TMFG_network_high.csv"))

measLow  <- semnet_measures(netLow)
measHigh <- semnet_measures(netHigh)

saveRDS(list(low = measLow, high = measHigh), file.path(out_dir, "10_semantic_network_measures.rds"))
capture.output(list(low = measLow, high = measHigh), file = file.path(out_dir, "10_semantic_network_measures.txt"))

############################
# 10. Sensitivity checks from the authors' code
############################

if (RUN_SENSITIVITY) {
  sensitivity_results <- list()
  
  lowOne  <- deLow[,  colSums(deLow)  >= 1, drop = FALSE]
  highOne <- deHigh[, colSums(deHigh) >= 1, drop = FALSE]
  
  sensitivity_results$noEqLow <- semnet_measures(
    get_tmfg_adjacency(cosine_with_constant(lowOne, addConstant = 0))
  )
  
  sensitivity_results$noEqHigh <- semnet_measures(
    get_tmfg_adjacency(cosine_with_constant(highOne, addConstant = 0))
  )
  
  eq2 <- SemNeT::equate(lowOne, highOne)
  eq2Low  <- extract_equated(eq2, 1)
  eq2High <- extract_equated(eq2, 2)
  
  sensitivity_results$eqLow <- semnet_measures(
    get_tmfg_adjacency(cosine_with_constant(eq2Low, addConstant = .01))
  )
  
  sensitivity_results$eqHigh <- semnet_measures(
    get_tmfg_adjacency(cosine_with_constant(eq2High, addConstant = .01))
  )
  
  saveRDS(sensitivity_results, file.path(out_dir, "11_sensitivity_checks.rds"))
  capture.output(sensitivity_results, file = file.path(out_dir, "11_sensitivity_checks.txt"))
}

############################
# 11. Partial network bootstrap using current SemNeT API
############################

if (RUN_BOOTSTRAP) {
  
  run_bootSemNeT <- function(prop) {
    weighted_arg <- FALSE
    if (MATCH_SOURCE_BOOTSTRAP && isTRUE(all.equal(prop, 0.50))) weighted_arg <- TRUE
    
    SemNeT::bootSemNeT(
      high_net_input,
      low_net_input,
      prop = prop,
      iter = N_BOOT,
      sim = "cosine",
      cores = N_CORES,
      method = "TMFG",
      type = "node",
      weighted = weighted_arg
    )
  }
  
  nodedrop <- list(
    fifty   = run_bootSemNeT(0.50),
    sixty   = run_bootSemNeT(0.60),
    seventy = run_bootSemNeT(0.70),
    eighty  = run_bootSemNeT(0.80),
    ninety  = run_bootSemNeT(0.90)
  )
  
  saveRDS(nodedrop, file.path(out_dir, "12_nodedrop_bootSemNeT_object.rds"))
  
  boot_tests <- list(
    ninety  = SemNeT::test.bootSemNeT(nodedrop$ninety),
    eighty  = SemNeT::test.bootSemNeT(nodedrop$eighty),
    seventy = SemNeT::test.bootSemNeT(nodedrop$seventy),
    sixty   = SemNeT::test.bootSemNeT(nodedrop$sixty),
    fifty   = SemNeT::test.bootSemNeT(nodedrop$fifty)
  )
  
  saveRDS(boot_tests, file.path(out_dir, "13_bootSemNeT_tests.rds"))
  capture.output(boot_tests, file = file.path(out_dir, "13_bootSemNeT_tests.txt"))
  
  # Plot export is wrapped in try() because plot.bootSemNeT output format can vary by version.
  plot_obj <- try(
    plot(
      nodedrop$fifty, nodedrop$sixty, nodedrop$seventy,
      nodedrop$eighty, nodedrop$ninety,
      groups = c("High", "Low")
    ),
    silent = TRUE
  )
  
  saveRDS(plot_obj, file.path(out_dir, "14_bootSemNeT_plot_object.rds"))
}

############################
# 12. Cytoscape input files
############################

highCyto <- SemNeT::convert2cytoscape(netHigh)
lowCyto  <- SemNeT::convert2cytoscape(netLow)

write.csv(highCyto, file.path(out_dir, "15_high_open_cyto.csv"), row.names = FALSE)
write.csv(lowCyto,  file.path(out_dir, "15_low_open_cyto.csv"), row.names = FALSE)

############################
# 13. Session information
############################

capture.output(sessionInfo(), file = file.path(out_dir, "99_session_info.txt"))

message("Reproduction script finished. Outputs saved to: ", out_dir)
message("Check these files first:")
message("  06_behavior_summary.csv")
message("  10_semantic_network_measures.txt")
message("  13_bootSemNeT_tests.txt, if RUN_BOOTSTRAP = TRUE")
message("  15_high_open_cyto.csv / 15_low_open_cyto.csv")
