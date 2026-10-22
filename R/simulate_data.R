# simulate_data.R
# Creates a SIMULATED lipidomics dataset (no real patient data) shaped like
# a quantitative lipidomics table: samples x lipid species, in nmol/mL.
# Run once from the project folder:  source("R/simulate_data.R")

set.seed(42)

n_per_group <- 150
groups <- rep(c("Control", "Cancer"), each = n_per_group)
n <- length(groups)
sex <- sample(c("F", "M"), n, replace = TRUE)

# Lipid species, grouped by lipid class
lipids <- data.frame(
  lipid = c(paste0("CE ",  c("16:0", "18:1", "18:2", "20:4", "22:6")),
            paste0("Cer ", c("34:1", "36:1", "40:1", "41:1", "42:1", "42:2")),
            paste0("LPC ", c("16:0", "18:0", "18:1", "18:2", "20:4")),
            paste0("PC ",  c("32:0", "34:1", "34:2", "36:2", "36:3", "36:4", "38:4", "38:6")),
            paste0("SM ",  c("32:1", "34:1", "36:1", "39:1", "40:1", "41:1", "42:1", "42:2")),
            paste0("TG ",  c("50:1", "52:2", "52:3", "54:3", "54:4"))),
  stringsAsFactors = FALSE
)
lipids$class <- sub(" .*", "", lipids$lipid)

# Typical concentration level per class (nmol/mL, made up but plausible order)
class_level <- c(CE = 800, Cer = 1.5, LPC = 40, PC = 250, SM = 60, TG = 150)

# Simulated "cancer effect" on log scale: very long-chain SM and Cer lower,
# some LPC lower, small random effects elsewhere
effect <- setNames(rnorm(nrow(lipids), 0, 0.04), lipids$lipid)
effect[c("SM 39:1", "SM 40:1", "SM 41:1", "SM 42:1", "Cer 41:1", "Cer 42:1")] <- -0.35
effect[c("LPC 18:2", "CE 16:0", "PC 36:2", "PC 36:3")] <- -0.20
effect[c("TG 52:2", "TG 54:3")] <- 0.10

# Correlated noise: a shared factor within each lipid class + individual noise
class_factor <- sapply(unique(lipids$class), function(cl) rnorm(n, 0, 0.25))
mat <- sapply(seq_len(nrow(lipids)), function(j) {
  cl <- lipids$class[j]
  base <- log(class_level[cl] * runif(1, 0.3, 2))
  sex_eff <- ifelse(sex == "M", rnorm(1, 0, 0.08), 0)
  log_conc <- base + class_factor[, cl] + sex_eff +
    ifelse(groups == "Cancer", effect[j], 0) + rnorm(n, 0, 0.30)
  round(exp(log_conc), 4)
})
colnames(mat) <- lipids$lipid

wide <- data.frame(sample_id = sprintf("S%03d", seq_len(n)),
                   group = groups, sex = sex, mat,
                   check.names = FALSE)

dir.create("data", showWarnings = FALSE)
write.csv(wide, "data/demo_lipids.csv", row.names = FALSE)
message("Wrote data/demo_lipids.csv: ", n, " samples x ", nrow(lipids), " lipids")
