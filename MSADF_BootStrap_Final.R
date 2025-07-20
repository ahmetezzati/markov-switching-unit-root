# نصب و بارگذاری بسته‌ها
if (!require("MSwM")) install.packages("MSwM")
library(MSwM)

# --- داده‌ها ---
y <- read.table("clipboard", header = FALSE)
y <- ts(y)
dy <- diff(y)
dy_lag1 <- c(NA, dy[-length(dy)])
dy <- dy[-1]
dy_lag1 <- dy_lag1[-1]

df_sim <- data.frame(dy = dy, dy_lag1 = dy_lag1)
model_sim <- lm(dy ~ dy_lag1, data = df_sim)
sw <- rep(TRUE, length(coef(model_sim)) + 1)
msm_sim <- msmFit(model_sim, k = 3, sw = sw)

coefs <- msm_sim@Coef
intercepts <- coefs[, "(Intercept)"]
phis <- coefs[, "dy_lag1"]
sds <- sqrt(msm_sim@std)
P <- msm_sim@transMat

# --- تابع شبیه‌سازی برای 3 رژیم ---
simulate_y_MS <- function(n, intercepts, phis, sds, P, y0 = y[1]) {
  dy_sim <- numeric(n)
  regimes <- numeric(n)
  
  regimes[1] <- sample(1:3, 1)
  dy_sim[1] <- rnorm(1, mean = intercepts[regimes[1]], sd = sds[regimes[1]])
  
  for (t in 2:n) {
    prev_regime <- regimes[t - 1]
    regimes[t] <- sample(1:3, 1, prob = P[prev_regime, ])
    r <- regimes[t]
    dy_sim[t] <- intercepts[r] + phis[r] * dy_sim[t - 1] + rnorm(1, 0, sds[r])
  }
  
  y_sim <- cumsum(c(y0, dy_sim))
  return(y_sim)
}

# --- تولید سری‌های مصنوعی ---
set.seed(123)
simulated_y_list <- lapply(1:100, function(i) simulate_y_MS(length(dy), intercepts, phis, sds, P))

# --- استخراج t-value ---
extract_t_value <- function(summary_output, regime_num) {
  regime_start <- grep(paste0("Regime ", regime_num), summary_output)[1]
  next_regime <- grep(paste0("Regime ", regime_num + 1), summary_output)
  
  if (length(next_regime) == 0) {
    regime_lines <- summary_output[regime_start:length(summary_output)]
  } else {
    regime_lines <- summary_output[regime_start:(next_regime[1] - 1)]
  }
  
  line <- grep("y_lag1\\(S\\)", regime_lines, value = TRUE)
  if (length(line) == 0) return(NA)
  fields <- strsplit(trimws(line), "\\s+")[[1]]
  if (length(fields) < 4) return(NA)
  
  return(as.numeric(fields[4]))
}

# --- سری واقعی ---
y_real <- y
y_lag1_real <- y_real[-length(y_real)]
dy_real <- diff(y_real)
dy_lag1_real <- c(NA, dy_real[-length(dy_real)])

y_lag1_real <- y_lag1_real[-1]
dy_real <- dy_real[-1]
dy_lag1_real <- dy_lag1_real[-1]

df_real <- data.frame(dy = dy_real, y_lag1 = y_lag1_real, dy_lag1 = dy_lag1_real)
model_real <- lm(dy ~ y_lag1 + dy_lag1, data = df_real)
sw <- rep(TRUE, length(coef(model_real)) + 1)
msm_real <- msmFit(model_real, k = 3, sw = sw)

summary_real <- capture.output(summary(msm_real))
t_real <- sapply(1:3, function(r) extract_t_value(summary_real, r))

# --- آزمون روی سری‌های مصنوعی ---
t_sim_all <- matrix(NA, nrow = length(simulated_y_list), ncol = 3)

for (i in seq_along(simulated_y_list)) {
  y_sim <- simulated_y_list[[i]][-1]
  y_lag1_sim <- y_sim[-length(y_sim)]
  dy_sim_clean <- diff(y_sim)
  dy_lag1_sim <- c(NA, dy_sim_clean[-length(dy_sim_clean)])
  
  y_lag1_sim <- y_lag1_sim[-1]
  dy_sim_clean <- dy_sim_clean[-1]
  dy_lag1_sim <- dy_lag1_sim[-1]
  
  df_sim <- data.frame(dy = dy_sim_clean, y_lag1 = y_lag1_sim, dy_lag1 = dy_lag1_sim)
  model_sim <- lm(dy ~ y_lag1 + dy_lag1, data = df_sim)
  
  msm_sim <- try(msmFit(model_sim, k = 3, sw = sw), silent = TRUE)
  if (inherits(msm_sim, "try-error")) next
  
  summary_sim <- capture.output(summary(msm_sim))
  t_sim_all[i, ] <- sapply(1:3, function(r) extract_t_value(summary_sim, r))
}
