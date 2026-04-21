# 加载必要的包
library(mvtnorm)
library(tidyverse)
library(nhanesA)
suppressWarnings(suppressMessages(library(Bessel)))
library(glmmTMB)
library(copula)
library(statmod)

# # 获取并整理数据
# demo_raw <- nhanes('DEMO_J')
# biopro_raw <- nhanes('BIOPRO_J')
# bmx_raw <- nhanes('BMX_J')
# ghb_raw <- nhanes('GHB_J')
# alq_raw <- nhanes('ALQ_J')
# bpx_raw <- nhanes('BPX_J')
# 
# demo_sub <- demo_raw %>% select(SEQN, RIDAGEYR, RIAGENDR)
# biopro_sub <- biopro_raw %>% select(SEQN, LBXSATSI, LBXSASSI)
# bmx_sub <- bmx_raw %>% select(SEQN, BMXBMI)
# ghb_sub <- ghb_raw %>% select(SEQN, LBXGH)
# alq_sub <- alq_raw %>% select(SEQN, ALQ130)
# bpx_sub <- bpx_raw %>% select(SEQN, BPXSY1)
# 
# # 按受访者编号(SEQN)进行内连接合并
# df_merged <- demo_sub %>%
#   inner_join(biopro_sub, by = "SEQN") %>%
#   inner_join(bmx_sub, by = "SEQN") %>% 
#   inner_join(ghb_sub, by = "SEQN") %>% 
#   inner_join(alq_sub, by = "SEQN") %>% 
#   inner_join(bpx_sub, by = "SEQN")
# 
# df_clean <- df_merged %>%
#   # 剔除包含 NA 缺失值的样本
#   na.omit() %>%
#   # 转换为易于理解的列名
#   rename(
#     Age = RIDAGEYR,
#     Gender = RIAGENDR,
#     BMI = BMXBMI,
#     ALT = LBXSATSI,
#     AST = LBXSASSI,
#     HbA1c = LBXGH,
#     Alcohol_Cups = ALQ130,
#     SysBP = BPXSY1
#   ) %>%
#   # 将性别转换为 0(女) 和 1(男) 的虚拟变量
#   mutate(Gender = ifelse(Gender == "Male", 1, 0))
# 
# # 去掉编号列并转换单位
# df <- df_clean[, -1]
# df$ALT <- df$ALT / 10
# df$AST <- df$AST / 10
# 
# # 因变量画图
# par(mfrow = c(2, 1))
# hist(df$ALT, breaks = 30)
# hist(df$AST, breaks = 30)
# 
# 
# # 数据集划分训练集与测试集
# df$ID <- 1:nrow(df)
# 
# saveRDS(df, file = "NHANES.rds")

df <- readRDS("NHANES.rds")
n <- nrow(df)

set.seed(2026)
train_idx <- sample(1:n, size = 0.8 * n)

train_data <- df[train_idx, ]
test_data  <- df[-train_idx, ]

n_train <- nrow(train_data)
n_test <- nrow(test_data)
d <- 2

# 提取测试集的响应变量真实值
Y_test <- as.matrix(test_data[, c("ALT", "AST")])

# 建立用于存储结果的数据框
results_train <- data.frame(
  Model = character(), 
  Train_LogLik = numeric(), 
  Train_AIC = numeric(), 
  Train_BIC = numeric(), 
  stringsAsFactors = FALSE
)

results_test <- data.frame(
  Model = character(),
  Test_LogLik = numeric(),
  Test_AIC = numeric(),
  Test_BIC = numeric(),
  MAE_ALT = numeric(),
  MAE_AST = numeric(),
  stringsAsFactors = FALSE
)

# 回归公式 (复用)
formula_eq <- ~ Age + Gender + BMI + HbA1c + SysBP + Alcohol_Cups

# ==========================================
# 2. 模型 1: 独立 Gamma 均值回归
# ==========================================
# 在训练集上拟合
fit_g_alt <- glm(update(ALT ~ ., formula_eq), family = Gamma(link = "log"), data = train_data)
fit_g_ast <- glm(update(AST ~ ., formula_eq), family = Gamma(link = "log"), data = train_data)

# 训练集效果
ll.gamma <- as.numeric(logLik(fit_g_alt) + logLik(fit_g_ast))
p.gamma <- length(coef(fit_g_alt)) * d + d
AIC.gamma <- -2 * ll.gamma + 2 * p.gamma
BIC.gamma <- -2 * ll.gamma + log(n_train) * p.gamma

# 在测试集上预测条件均值 mu
pred_g_alt <- predict(fit_g_alt, newdata = test_data, type = "response")
pred_g_ast <- predict(fit_g_ast, newdata = test_data, type = "response")

# 计算 MAE
mae_g_alt <- mean(abs(Y_test[, "ALT"] - pred_g_alt))
mae_g_ast <- mean(abs(Y_test[, "AST"] - pred_g_ast))

# 计算测试集对数似然值
# 注意：需要提取训练集的离散参数 (dispersion) 来计算测试集密度
shape_alt <- 1 / summary(fit_g_alt)$dispersion
shape_ast <- 1 / summary(fit_g_ast)$dispersion

ll_g_alt <- sum(dgamma(Y_test[, "ALT"], shape = shape_alt, rate = shape_alt / pred_g_alt, log = TRUE))
ll_g_ast <- sum(dgamma(Y_test[, "AST"], shape = shape_ast, rate = shape_ast / pred_g_ast, log = TRUE))

test_ll_gamma <- ll_g_alt + ll_g_ast

# ==========================================
# 3. 模型 2: 独立逆高斯 (Inverse Gaussian) 均值回归
# ==========================================
# 在训练集上拟合
fit_ig_alt <- glm(update(ALT ~ ., formula_eq), family = inverse.gaussian(link = "log"), data = train_data)
fit_ig_ast <- glm(update(AST ~ ., formula_eq), family = inverse.gaussian(link = "log"), data = train_data)

# 训练集效果
ll.ig <- as.numeric(logLik(fit_ig_alt) + logLik(fit_ig_ast))
p.ig <- length(coef(fit_ig_alt)) * d + d
AIC.ig <- -2 * ll.ig + 2 * p.ig
BIC.ig <- -2 * ll.ig + log(n_train) * p.ig

# 在测试集上预测条件均值 mu
pred_ig_alt <- predict(fit_ig_alt, newdata = test_data, type = "response")
pred_ig_ast <- predict(fit_ig_ast, newdata = test_data, type = "response")

# 计算 MAE
mae_ig_alt <- mean(abs(Y_test[, "ALT"] - pred_ig_alt))
mae_ig_ast <- mean(abs(Y_test[, "AST"] - pred_ig_ast))

# 计算测试集对数似然值
# IG 对数密度公式: -0.5 * log(2 * pi * phi * y^3) - (y - mu)^2 / (2 * phi * mu^2 * y)
calc_ig_ll <- function(y, mu, phi) {
  sum(-0.5 * log(2 * pi * phi * y^3) - (y - mu)^2 / (2 * phi * mu^2 * y))
}

phi_alt <- summary(fit_ig_alt)$dispersion
phi_ast <- summary(fit_ig_ast)$dispersion

ll_ig_alt <- calc_ig_ll(Y_test[, "ALT"], pred_ig_alt, phi_alt)
ll_ig_ast <- calc_ig_ll(Y_test[, "AST"], pred_ig_ast, phi_ast)

test_ll_ig <- ll_ig_alt + ll_ig_ast

# ==========================================

# GeGa均值回归模型
IGaGa.MLE <- function(alpha, Beta, Var.tau, x, z, ep = 1e-8) {
  # 收敛指标与迭代次数
  index <- 0
  k <- 1
  
  # 样本量n
  n <- nrow(x)
  d <- ncol(x)
  
  # 确定τ分布的参数λ
  lambda <- 1 / Var.tau + 2
  
  log_f <- function(x, z, alpha, Beta, lambda) {
    d <- ncol(Beta)
    alpha.plus <- d * alpha
    f <- function(data) {
      x <- data[1:d]
      z <- data[-(1:d)]
      mu <- exp(z %*% Beta)
      term1 <- gamma(alpha.plus + lambda) * (lambda - 1)^lambda / gamma(lambda)
      term2 <- prod(alpha^alpha * x^(alpha - 1) / gamma(alpha) / mu^alpha)
      term3 <- (sum(alpha * x / mu) + lambda - 1)^(-alpha.plus - lambda)
      term1 * term2 * term3
    }
    pdf <- apply(cbind(x, z), 1, f)
    return(sum(log(pdf)))
  }
  
  # US算法求c+log(s)-ψ(s)=0的零点
  US <- function(s, c) {
    index <- 0
    k <- 0
    
    s_new <- s
    
    while (index == 0) {
      s <- s_new
      q <- c + log(s) - digamma(s) + pi^2 * s / 6 - 1 / s
      s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
      
      if (abs(s_new - s) < 1e-5) {
        index <- 1
        break
      }
      k <- k + 1
    }
    return(s_new)
  }
  # US算法求c+1/(s-1)+log(s)-ψ(s)=0的零点
  US.IGa <- function(s, c) {
    index <- 0
    k <- 0
    
    s_new <- s
    
    while (index == 0) {
      s <- s_new
      q <- c + log(s - 1) + 1 / (s - 1) - digamma(s) + pi^2 * s / 6 - 1 / s
      s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
      
      if (abs(s_new - s) < 1e-5) {
        index <- 1
        break
      }
      k <- k + 1
    }
    return(s_new)
  }
  
  # GIG的矩
  GIG.m <- function(x, z, alpha, Beta, lambda) {
    d <- ncol(x)
    alpha.plus <- d * alpha
    compute.b <- function(data) {
      x <- data[1:d]
      z <- data[-(1:d)]
      mu <- exp(z %*% Beta)
      sum(alpha * x / mu) + lambda - 1
    }
    a <- alpha.plus + lambda
    b <- apply(cbind(x, z), 1, compute.b)
    
    result1 <- a / b
    result2 <- log(b) - digamma(a)
    
    
    return(cbind(result1, result2))
  }
  
  alpha_new <- alpha
  Beta_new <- Beta
  lambda_new <- lambda
  loglikeli_new <- log_f(x, z, alpha, Beta, lambda)
  
  while (k <= 1000) {
    alpha <- alpha_new
    Beta <- Beta_new
    lambda <- lambda_new
    mu <- exp(z %*% Beta)
    loglikeli <- loglikeli_new
    
    gig.vals <- GIG.m(x, z, alpha, Beta, lambda)
    b1 <- gig.vals[, 1]
    d1 <- gig.vals[, 2]
    Wb1 <- diag(b1)
    Bb1 <- mean(b1)
    Bd1 <- mean(d1)
    # 更新αj
    c <- 1 + mean(log(x) - log(mu) - apply(x, 2, function(s) b1 * s) / mu) - Bd1
    alpha_new <- US(alpha, c)
    cat("Alpha:", alpha_new, "\n")
    # 更新βj
    for (i in 1:d) {
      r <- log(mu[, i]) + x[, i] / mu[, i] - 1 / b1
      Beta_new[, i] <- solve(t(z) %*% Wb1 %*% z) %*% t(z) %*% Wb1 %*% r
    }
    cat("Beta:", Beta_new, "\n")
    # 更新λ
    c <- 1 - Bb1 - Bd1
    lambda_new <- US.IGa(lambda, c)
    cat("Var(tau):", 1 / (lambda_new - 2), "\n")
    # 更新似然函数
    loglikeli_new <- log_f(x, z, alpha_new, Beta_new, lambda_new)
    cat("Loglikelihood:", loglikeli_new, "\n")
    
    if (abs((loglikeli_new - loglikeli) / loglikeli) < ep) {
      index <- 1
      break
    }
    cat("Iteration number:", k, "\n")
    k <- k + 1
  }
  
  result <- list(Alpha = alpha_new, 
                 Beta = Beta_new, 
                 Lambda = lambda_new, 
                 Loglikelihood = loglikeli_new,
                 number = k, 
                 index = index)
  return(result)
}

GIGGa.MLE <- function(alpha, Beta, Var.tau, x, z, ep = 1e-8) {
  # 收敛指标与迭代次数
  index <- 0
  k <- 1
  
  # 样本量n
  n <- nrow(x)
  d <- ncol(x)
  
  # GIGGa
  GIG.lambda <- function(Var.tau) {
    f <- function(lambda) {
      besselK(lambda, 2) * besselK(lambda, 0) / (besselK(lambda, 1))^2 - 1 - Var.tau
    }
    return(uniroot(f, c(1, 5), extendInt = "yes")$root)
  }
  a <- function(lambda) lambda * besselK(lambda, 1) / besselK(lambda, 0)
  b <- function(lambda) lambda * besselK(lambda, 0) / besselK(lambda, 1)
  
  log_f <- function(x, z, alpha, Beta, lambda) {
    d <- ncol(Beta)
    alpha.plus <- d * alpha
    f <- function(data) {
      x <- data[1:d]
      z <- data[-(1:d)]
      mu <- exp(z %*% Beta)
      term1 <- prod(alpha^alpha * x^(alpha - 1) / gamma(alpha) / mu^alpha)
      term2 <- besselK(sqrt(a(lambda) * (b(lambda) + 2 * sum(alpha * x / mu))), alpha.plus) / besselK(lambda, 0)
      term3 <- (a(lambda) / (b(lambda) + 2 * sum(alpha * x / mu)))^(alpha.plus / 2)
      term1 * term2 * term3
    }
    pdf <- apply(cbind(x, z), 1, f)
    return(sum(log(pdf)))
  }
  
  # 确定τ分布的参数λ
  lambda <- GIG.lambda(Var.tau)
  
  # US算法求c+log(s)-ψ(s)=0的零点
  US <- function(s, c) {
    index <- 0
    k <- 0
    
    s_new <- s
    
    while (index == 0) {
      s <- s_new
      q <- c + log(s) - digamma(s) + pi^2 * s / 6 - 1 / s
      s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
      
      if (abs(s_new - s) < 1e-5) {
        index <- 1
        break
      }
      k <- k + 1
    }
    return(s_new)
  }
  
  # 牛顿法求解λ的迭代
  solve.lambda <- function(lam.ini, a, b) {
    opt <- function(lambda) {
      2 * log(besselK(lambda, 0)) + lambda * (a * besselK(lambda, 1) / besselK(lambda, 0) + b * besselK(lambda, 0) / besselK(lambda, 1))
    }
    
    fit <- optim(
      par = lam.ini,
      fn = opt,
      method = "L-BFGS-B",
      lower = 1e-6
    )
    
    return(fit$par)
  }
  
  # GIG的矩
  GIG.m <- function(x, z, alpha, Beta, lambda) {
    d <- ncol(x)
    alpha.plus <- d * alpha
    
    compute.b <- function(data) {
      x <- data[1:d]
      z <- data[-(1:d)]
      mu <- exp(z %*% Beta)
      b(lambda) + 2 * sum(alpha * x / mu)
    }
    a.value <- a(lambda)
    b.value <- apply(cbind(x, z), 1, compute.b)
    p <- -alpha.plus
    
    result1 <- BesselK(sqrt(a.value * b.value), p + 1, expo = TRUE) / BesselK(sqrt(a.value * b.value), p, expo = TRUE) * (b.value / a.value)^0.5
    result2 <- BesselK(sqrt(a.value * b.value), p - 1, expo = TRUE) / BesselK(sqrt(a.value * b.value), p, expo = TRUE) * (b.value / a.value)^(-0.5)
    ep <- 1e-6
    result3 <- (log(BesselK(sqrt(a.value * b.value), p + ep, expo = FALSE)) - log(BesselK(sqrt(a.value * b.value), p, expo = FALSE))) / ep + 0.5 * log(b.value / a.value)
    
    return(cbind(result1, result2, result3))
  }
  
  alpha_new <- alpha
  Beta_new <- Beta
  lambda_new <- lambda
  loglikeli_new <- log_f(x, z, alpha, Beta, lambda)
  
  while (k <= 1000) {
    alpha <- alpha_new
    Beta <- Beta_new
    lambda <- lambda_new
    mu <- exp(z %*% Beta)
    loglikeli <- loglikeli_new
    
    gig.vals <- GIG.m(x, z, alpha, Beta, lambda)
    a1 <- gig.vals[, 1]
    b1 <- gig.vals[, 2]
    d1 <- gig.vals[, 3]
    Wb1 <- diag(b1)
    Ba1 <- mean(a1)
    Bb1 <- mean(b1)
    Bd1 <- mean(d1)
    # 更新αj
    c <- 1 + mean(log(x) - log(mu) - apply(x, 2, function(s) b1 * s) / mu) - Bd1
    alpha_new <- US(alpha, c)
    cat("Alpha:", alpha_new, "\n")
    # 更新βj
    for (i in 1:d) {
      r <- log(mu[, i]) + x[, i] / mu[, i] - 1 / b1
      Beta_new[, i] <- solve(t(z) %*% Wb1 %*% z) %*% t(z) %*% Wb1 %*% r
    }
    cat("Beta:", Beta_new, "\n")
    # 更新λ
    lambda_new <- solve.lambda(lambda, Ba1, Bb1)
    cat("Var(tau):", besselK(lambda, 2) * besselK(lambda, 0) / (besselK(lambda, 1))^2 - 1, "\n")
    # 更新似然函数
    loglikeli_new <- log_f(x, z, alpha_new, Beta_new, lambda_new)
    cat("Loglikelihood:", loglikeli_new, "\n")
    
    if (abs((loglikeli_new - loglikeli) / loglikeli) < ep) {
      index <- 1
      break
    }
    cat("Iteration number:", k, "\n")
    k <- k + 1
  }
  
  result <- list(Alpha = alpha_new, 
                 Beta = Beta_new, 
                 Lambda = lambda_new, 
                 Loglikelihood = loglikeli_new,
                 number = k, 
                 index = index)
  return(result)
}

IGauGa.MLE <- function(alpha, Beta, Var.tau, y, X, ep = 1e-8) {
  index <- 0
  k <- 1
  
  n <- nrow(y)
  d <- ncol(y)
  
  lambda <- 1 / Var.tau
  
  # 定义对数似然函数
  log_f3 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    f <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      term1 <- exp(lambda) * sqrt(2 * lambda / pi)
      term2 <- prod(alpha^alpha * y^(alpha - 1) / gamma(alpha) / mu^alpha)
      term3 <- (lambda / (lambda + 2 * sum(alpha * y / mu)))^(alpha.plus / 2 + 1 / 4)
      term4 <- besselK(sqrt(lambda * (lambda + 2 * sum(alpha * y / mu))), alpha.plus + 1 / 2, expo = FALSE)
      term1 * term2 * term3 * term4
    }
    pdf <- apply(rbind(t(y), X), 2, f)
    return(sum(log(pdf)))
  }
  
  # US算法求c+log(s)-ψ(s)=0的零点
  US <- function(s, c) {
    index <- 0
    k <- 0
    
    s_new <- s
    
    while (index == 0) {
      s <- s_new
      q <- c + log(s) - digamma(s) + pi^2 * s / 6 - 1 / s
      s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
      
      if (abs(s_new - s) < 1e-5) {
        index <- 1
        break
      }
      k <- k + 1
    }
    return(s_new)
  }
  
  # GIG的矩
  f3 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    compute.b <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      lambda + 2 * sum(alpha * y / mu)
    }
    a <- lambda
    b <- apply(rbind(t(y), X), 2, compute.b)
    p <- -alpha.plus - 1 / 2
    
    result1 <- BesselK(sqrt(a * b), p + 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^0.5
    result2 <- BesselK(sqrt(a * b), p - 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^(-0.5)
    ep <- 1e-6
    result3 <- (log(BesselK(sqrt(a * b), p + ep, expo = FALSE)) - log(BesselK(sqrt(a * b), p, expo = FALSE))) / ep + 0.5 * log(b / a)
    
    return(cbind(result1, result2, result3))
  }
  
  alpha_new <- alpha
  Beta_new <- Beta
  lambda_new <- lambda
  loglikeli_new <- log_f3(y, X, alpha, Beta, lambda)
  
  while (k <= 1000) {
    alpha <- alpha_new
    Beta <- Beta_new
    lambda <- lambda_new
    mu <- exp(t(X) %*% Beta)
    loglikeli <- loglikeli_new
    
    Moment <- f3(y, X, alpha, Beta, lambda)
    a3 <- Moment[, 1]
    b3 <- Moment[, 2]
    d3 <- Moment[, 3]
    Ba3 <- mean(a3)
    Bb3 <- mean(b3)
    Wb3 <- diag(b3)
    Bd3 <- mean(d3)
    # 更新αj
    c <- 1 + mean(log(y) - log(mu) - apply(y, 2, function(x) b3 * x) / mu) - Bd3
    alpha_new <- US(alpha, c)
    cat("Alpha:", alpha_new, "\n")
    # 更新βj
    for (j in 1:d) {
      r <- log(mu[, j]) + y[, j] / mu[, j] - 1 / b3
      Beta_new[, j] <- solve(X %*% Wb3 %*% t(X)) %*% X %*% Wb3 %*% r
    }
    cat("Beta:", Beta_new, "\n")
    # 更新λ
    lambda_new <- 1 / (Ba3 + Bb3 - 2)
    cat("Lambda:", lambda_new, "\n")
    # 更新似然函数
    loglikeli_new <- log_f3(y, X, alpha_new, Beta_new, lambda_new)
    cat("Loglikelihood:", loglikeli_new, "\n")
    
    if (abs((loglikeli_new - loglikeli) / loglikeli) < ep) {
      index <- 1
      break
    }
    cat("Iteration number:", k, "\n")
    k <- k + 1
  }
  
  result <- list(Alpha = alpha_new, 
                 Beta = Beta_new, 
                 Lambda = lambda_new, 
                 Loglikelihood = loglikeli_new,
                 number = k, 
                 index = index)
  return(result)
}

RIGGa.MLE <- function(alpha, Beta, Var.tau, y, X, ep = 1e-8) {
  index <- 0
  k <- 1
  
  n <- nrow(y)
  d <- ncol(y)
  
  lambda <- (3 + sqrt(1 + 4 * Var.tau)) / (4 - 2 * Var.tau)
  
  # 定义对数似然函数
  log_f4 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    f <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      term1 <- exp(1 / (lambda - 1)) * sqrt(2 * lambda / (pi * (lambda - 1)))
      term2 <- prod(alpha^alpha * y^(alpha - 1) / gamma(alpha) / mu^alpha)
      term3 <- ((lambda / (lambda - 1)) / ((lambda * (lambda - 1))^(-1) + 2 * sum(alpha * y / mu)))^(-alpha.plus / 2 + 1 / 4)
      term4 <- lambda / (lambda - 1) * (1 / (lambda * (lambda - 1)) + 2 * sum(alpha * y / mu))
      term5 <- besselK(sqrt(term4), -alpha.plus + 1 / 2, expo = FALSE)
      term1 * term2 / term3 * term5
    }
    pdf <- apply(rbind(t(y), X), 2, f)
    return(sum(log(pdf)))
  }
  
  # US算法求c+log(s)-ψ(s)=0的零点
  US <- function(s, c) {
    index <- 0
    k <- 0
    
    s_new <- s
    
    while (index == 0) {
      s <- s_new
      q <- c + log(s) - digamma(s) + pi^2 * s / 6 - 1 / s
      s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
      
      if (abs(s_new - s) < 1e-5) {
        index <- 1
        break
      }
      k <- k + 1
    }
    return(s_new)
  }
  
  # GIG的矩
  f4 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    compute.b <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      1 / (lambda * (lambda - 1)) + 2 * sum(alpha * y / mu)
    }
    a <- lambda / (lambda - 1)
    b <- apply(rbind(t(y), X), 2, compute.b)
    p <- -alpha.plus + 1 / 2
    
    result1 <- BesselK(sqrt(a * b), p + 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^0.5
    result2 <- BesselK(sqrt(a * b), p - 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^(-0.5)
    ep <- 1e-6
    result3 <- (log(BesselK(sqrt(a * b), p + ep, expo = FALSE)) - log(BesselK(sqrt(a * b), p, expo = FALSE))) / ep + 0.5 * log(b / a)
    
    return(cbind(result1, result2, result3))
  }
  
  alpha_new <- alpha
  Beta_new <- Beta
  lambda_new <- lambda
  loglikeli_new <- log_f4(y, X, alpha, Beta, lambda)
  
  while (k <= 1000) {
    alpha <- alpha_new
    Beta <- Beta_new
    lambda <- lambda_new
    mu <- exp(t(X) %*% Beta)
    loglikeli <- loglikeli_new
    
    Moment <- f4(y, X, alpha, Beta, lambda)
    a4 <- Moment[, 1]
    b4 <- Moment[, 2]
    d4 <- Moment[, 3]
    Ba4 <- mean(a4)
    Bb4 <- mean(b4)
    Wb4 <- diag(b4)
    Bd4 <- mean(d4)
    # 更新αj
    c <- 1 + mean(log(y) - log(mu) - apply(y, 2, function(x) b4 * x) / mu) - Bd4
    alpha_new <- US(alpha, c)
    cat("Alpha:", alpha_new, "\n")
    # 更新βj
    for (j in 1:d) {
      r <- log(mu[, j]) + y[, j] / mu[, j] - 1 / b4
      Beta_new[, j] <- solve(X %*% Wb4 %*% t(X)) %*% X %*% Wb4 %*% r
    }
    cat("Beta:", Beta_new, "\n")
    # 更新λ
    lambda_new <- (1 + 2 * Bb4 + sqrt((1 + 2 * Bb4)^2 - 4 * Bb4 * (3 - Ba4))) / (2 * (3 - Ba4))
    cat("Lambda:", lambda_new, "\n")
    # 更新似然函数
    loglikeli_new <- log_f4(y, X, alpha_new, Beta_new, lambda_new)
    cat("Loglikelihood:", loglikeli_new, "\n")
    
    if (abs((loglikeli_new - loglikeli) / loglikeli) < ep) {
      index <- 1
      break
    }
    cat("Iteration number:", k, "\n")
    k <- k + 1
  }
  
  result <- list(Alpha = alpha_new, 
                 Beta = Beta_new, 
                 Lambda = lambda_new, 
                 Loglikelihood = loglikeli_new,
                 number = k, 
                 index = index)
  return(result)
}

LNGa.MLE <- function(alpha, Beta, Var.tau, y, X, Type = "LN", ep = 1e-6) {
  index <- 0
  k <- 1
  
  n <- nrow(y)
  d <- ncol(y)
  
  if (Type == "Ga") {
    lambda <- 1 / Var.tau
  } else if (Type == "IGa") {
    lambda <- 1 / Var.tau + 2
  } else if (Type == "IGau") {
    lambda <- 1 / Var.tau
  } else if (Type == "RIG") {
    lambda <- (3 + sqrt(1 + 4 * Var.tau)) / (4 - 2 * Var.tau)
  } else if (Type == "LN") {
    lambda <- log(Var.tau + 1)
  }
  
  # 定义对数似然函数
  log_f1 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    f <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      term1 <- 2 * lambda^((lambda + alpha.plus) / 2) / gamma(lambda)
      term2 <- prod(alpha^alpha * y^(alpha - 1) / gamma(alpha) / mu^alpha)
      term3 <- sum(alpha * y / mu)^((lambda - alpha.plus) / 2)
      term4 <- besselK(sqrt(4 * lambda * sum(alpha * y / mu)), lambda - alpha.plus, expo = FALSE)
      term1 * term2 * term3 * term4
    }
    pdf <- apply(rbind(t(y), X), 2, f)
    return(sum(log(pdf)))
  }
  
  log_f2 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    f <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      term1 <- gamma(lambda + alpha.plus) * (lambda - 1)^lambda / gamma(lambda)
      term2 <- prod(alpha^alpha * y^(alpha - 1) / gamma(alpha) / mu^alpha)
      term3 <- (sum(alpha * y / mu) + lambda - 1)^(-lambda - alpha.plus)
      term1 * term2 * term3
    }
    pdf <- apply(rbind(t(y), X), 2, f)
    return(sum(log(pdf)))
  }
  
  log_f3 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    f <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      term1 <- exp(lambda) * sqrt(2 * lambda / pi)
      term2 <- prod(alpha^alpha * y^(alpha - 1) / gamma(alpha) / mu^alpha)
      term3 <- (lambda / (lambda + 2 * sum(alpha * y / mu)))^(alpha.plus / 2 + 1 / 4)
      term4 <- besselK(sqrt(lambda * (lambda + 2 * sum(alpha * y / mu))), alpha.plus + 1 / 2, expo = FALSE)
      term1 * term2 * term3 * term4
    }
    pdf <- apply(rbind(t(y), X), 2, f)
    return(sum(log(pdf)))
  }
  
  log_f4 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    f <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      term1 <- exp(1 / (lambda - 1)) * sqrt(2 * lambda / (pi * (lambda - 1)))
      term2 <- prod(alpha^alpha * y^(alpha - 1) / gamma(alpha) / mu^alpha)
      term3 <- ((lambda / (lambda - 1)) / ((lambda * (lambda - 1))^(-1) + 2 * sum(alpha * y / mu)))^(-alpha.plus / 2 + 1 / 4)
      term4 <- lambda / (lambda - 1) * (1 / (lambda * (lambda - 1)) + 2 * sum(alpha * y / mu))
      term5 <- besselK(sqrt(term4), -alpha.plus + 1 / 2, expo = FALSE)
      term1 * term2 / term3 * term5
    }
    pdf <- apply(rbind(t(y), X), 2, f)
    return(sum(log(pdf)))
  }
  
  log_f5 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    f <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      h5 <- function(tau) {
        tau^(-alpha.plus - 1) * exp(-sum(alpha * y / mu) / tau - (log(tau) + lambda / 2)^2 / (2 * lambda))
      }
      term1 <- 1 / sqrt(2 * pi * lambda)
      term2 <- prod(alpha^alpha * y^(alpha - 1) / gamma(alpha) / mu^alpha)
      term3 <- integrate(h5, 0, Inf)$value
      term1 * term2 * term3
    }
    pdf <- apply(rbind(t(y), X), 2, f)
    return(sum(log(pdf)))
  }
  
  # US算法求c+log(s)-ψ(s)=0的零点
  US <- function(s, c) {
    index <- 0
    k <- 0
    
    s_new <- s
    
    while (index == 0) {
      s <- s_new
      q <- c + log(s) - digamma(s) + pi^2 * s / 6 - 1 / s
      s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
      
      if (abs(s_new - s) < 1e-5) {
        index <- 1
        break
      }
      k <- k + 1
    }
    return(s_new)
  }
  # US算法求c+1/(s-1)+log(s)-ψ(s)=0的零点
  US.IGa <- function(s, c) {
    index <- 0
    k <- 0
    
    s_new <- s
    
    while (index == 0) {
      s <- s_new
      q <- c + log(s - 1) + 1 / (s - 1) - digamma(s) + pi^2 * s / 6 - 1 / s
      s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
      
      if (abs(s_new - s) < 1e-5) {
        index <- 1
        break
      }
      k <- k + 1
    }
    return(s_new)
  }
  
  # GIG的矩
  f1 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    compute.b <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      2 * sum(alpha * y / mu)
    }
    a <- 2 * lambda
    b <- apply(rbind(t(y), X), 2, compute.b)
    p <- lambda - alpha.plus
    
    result1 <- BesselK(sqrt(a * b), p + 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^0.5
    result2 <- BesselK(sqrt(a * b), p - 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^(-0.5)
    ep <- 1e-6
    result3 <- (log(BesselK(sqrt(a * b), p + ep, expo = FALSE)) - log(BesselK(sqrt(a * b), p, expo = FALSE))) / ep + 0.5 * log(b / a)
    
    return(cbind(result1, result2, result3))
  }
  
  f2 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    compute.b <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      sum(alpha * y / mu) + lambda - 1
    }
    a <- lambda + alpha.plus
    b <- apply(rbind(t(y), X), 2, compute.b)
    
    result2 <- a / b
    result3 <- log(b) - digamma(a)
    
    return(cbind(result2, result3))
  }
  
  f3 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    compute.b <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      lambda + 2 * sum(alpha * y / mu)
    }
    a <- lambda
    b <- apply(rbind(t(y), X), 2, compute.b)
    p <- -alpha.plus - 1 / 2
    
    result1 <- BesselK(sqrt(a * b), p + 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^0.5
    result2 <- BesselK(sqrt(a * b), p - 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^(-0.5)
    ep <- 1e-6
    result3 <- (log(BesselK(sqrt(a * b), p + ep, expo = FALSE)) - log(BesselK(sqrt(a * b), p, expo = FALSE))) / ep + 0.5 * log(b / a)
    
    return(cbind(result1, result2, result3))
  }
  
  f4 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    compute.b <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      1 / (lambda * (lambda - 1)) + 2 * sum(alpha * y / mu)
    }
    a <- lambda / (lambda - 1)
    b <- apply(rbind(t(y), X), 2, compute.b)
    p <- -alpha.plus + 1 / 2
    
    result1 <- BesselK(sqrt(a * b), p + 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^0.5
    result2 <- BesselK(sqrt(a * b), p - 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^(-0.5)
    ep <- 1e-6
    result3 <- (log(BesselK(sqrt(a * b), p + ep, expo = FALSE)) - log(BesselK(sqrt(a * b), p, expo = FALSE))) / ep + 0.5 * log(b / a)
    
    return(cbind(result1, result2, result3))
  }
  
  f5 <- function(y, X, alpha, Beta, lambda) {
    d <- ncol(y)
    alpha.plus <- d * alpha
    compute.moment <- function(data) {
      y <- data[1:d]
      X <- data[-(1:d)]
      mu <- exp(t(X) %*% Beta)
      h5 <- function(tau) {
        tau^(-alpha.plus - 1) * exp(-sum(alpha * y / mu) / tau - (log(tau) + lambda / 2)^2 / (2 * lambda))
      }
      int <- integrate(h5, 0, Inf)$value
      b5 <- function(tau) {
        tau^(-alpha.plus - 2) * exp(-sum(alpha * y / mu) / tau - (log(tau) + lambda / 2)^2 / (2 * lambda))
      }
      d5 <- function(tau) {
        log(tau) * tau^(-alpha.plus - 1) * exp(-sum(alpha * y / mu) / tau - (log(tau) + lambda / 2)^2 / (2 * lambda))
      }
      s5 <- function(tau) {
        (log(tau))^2 * tau^(-alpha.plus - 1) * exp(-sum(alpha * y / mu) / tau - (log(tau) + lambda / 2)^2 / (2 * lambda))
      }
      result1 <- integrate(b5, 0, Inf)$value / int
      result2 <- integrate(d5, 0, Inf)$value / int
      result3 <- integrate(s5, 0, Inf)$value / int
      return(c(result1, result2, result3))
    }
    result <- apply(rbind(t(y), X), 2, compute.moment) %>% t()
    
    return(result)
  }
  
  alpha_new <- alpha
  Beta_new <- Beta
  lambda_new <- lambda
  if (Type == "Ga") {
    loglikeli_new <- log_f1(y, X, alpha, Beta, lambda)
  } else if (Type == "IGa") {
    loglikeli_new <- log_f2(y, X, alpha, Beta, lambda)
  } else if (Type == "IGau") {
    loglikeli_new <- log_f3(y, X, alpha, Beta, lambda)
  } else if (Type == "RIG") {
    loglikeli_new <- log_f4(y, X, alpha, Beta, lambda)
  } else if (Type == "LN") {
    loglikeli_new <- log_f5(y, X, alpha, Beta, lambda)
  }
  
  while (k <= 1000) {
    alpha <- alpha_new
    Beta <- Beta_new
    lambda <- lambda_new
    mu <- exp(t(X) %*% Beta)
    loglikeli <- loglikeli_new
    
    if (Type == "Ga") {
      Moment <- f1(y, X, alpha, Beta, lambda)
      a1 <- Moment[, 1]
      b1 <- Moment[, 2]
      d1 <- Moment[, 3]
      Ba1 <- mean(a1)
      Wb1 <- diag(b1)
      Bd1 <- mean(d1)
      # 更新αj
      c <- 1 + mean(log(y) - log(mu) - apply(y, 2, function(x) b1 * x) / mu) - Bd1
      alpha_new <- US(alpha, c)
      cat("Alpha:", alpha_new, "\n")
      # 更新βj
      for (j in 1:d) {
        r <- log(mu[, j]) + y[, j] / mu[, j] - 1 / b1
        Beta_new[, j] <- solve(X %*% Wb1 %*% t(X)) %*% X %*% Wb1 %*% r
      }
      cat("Beta:", Beta_new, "\n")
      # 更新λ
      c <- 1 + Bd1 - Ba1
      lambda_new <- US(lambda, c)
      cat("Lambda:", lambda_new, "\n")
      # 更新似然函数
      loglikeli_new <- log_f1(y, X, alpha_new, Beta_new, lambda_new)
      cat("Loglikelihood:", loglikeli_new, "\n")
    } else if (Type == "IGa") {
      Moment <- f2(y, X, alpha, Beta, lambda)
      b2 <- Moment[, 1]
      d2 <- Moment[, 2]
      Wb2 <- diag(b2)
      Bb2 <- mean(b2)
      Bd2 <- mean(d2)
      # 更新αj
      c <- 1 + mean(log(y) - log(mu) - apply(y, 2, function(x) b2 * x) / mu) - Bd2
      alpha_new <- US(alpha, c)
      cat("Alpha:", alpha_new, "\n")
      # 更新βj
      for (j in 1:d) {
        r <- log(mu[, j]) + y[, j] / mu[, j] - 1 / b2
        Beta_new[, j] <- solve(X %*% Wb2 %*% t(X)) %*% X %*% Wb2 %*% r
      }
      cat("Beta:", Beta_new, "\n")
      # 更新λ
      c <- 1 - Bb2 - Bd2
      lambda_new <- US.IGa(lambda, c)
      cat("Lambda:", lambda_new, "\n")
      # 更新似然函数
      loglikeli_new <- log_f2(y, X, alpha_new, Beta_new, lambda_new)
      cat("Loglikelihood:", loglikeli_new, "\n")
    } else if (Type == "IGau") {
      Moment <- f3(y, X, alpha, Beta, lambda)
      a3 <- Moment[, 1]
      b3 <- Moment[, 2]
      d3 <- Moment[, 3]
      Ba3 <- mean(a3)
      Bb3 <- mean(b3)
      Wb3 <- diag(b3)
      Bd3 <- mean(d3)
      # 更新αj
      c <- 1 + mean(log(y) - log(mu) - apply(y, 2, function(x) b3 * x) / mu) - Bd3
      alpha_new <- US(alpha, c)
      cat("Alpha:", alpha_new, "\n")
      # 更新βj
      for (j in 1:d) {
        r <- log(mu[, j]) + y[, j] / mu[, j] - 1 / b3
        Beta_new[, j] <- solve(X %*% Wb3 %*% t(X)) %*% X %*% Wb3 %*% r
      }
      cat("Beta:", Beta_new, "\n")
      # 更新λ
      lambda_new <- 1 / (Ba3 + Bb3 - 2)
      cat("Lambda:", lambda_new, "\n")
      # 更新似然函数
      loglikeli_new <- log_f3(y, X, alpha_new, Beta_new, lambda_new)
      cat("Loglikelihood:", loglikeli_new, "\n")
    } else if (Type == "RIG") {
      Moment <- f4(y, X, alpha, Beta, lambda)
      a4 <- Moment[, 1]
      b4 <- Moment[, 2]
      d4 <- Moment[, 3]
      Ba4 <- mean(a4)
      Bb4 <- mean(b4)
      Wb4 <- diag(b4)
      Bd4 <- mean(d4)
      # 更新αj
      c <- 1 + mean(log(y) - log(mu) - apply(y, 2, function(x) b4 * x) / mu) - Bd4
      alpha_new <- US(alpha, c)
      cat("Alpha:", alpha_new, "\n")
      # 更新βj
      for (j in 1:d) {
        r <- log(mu[, j]) + y[, j] / mu[, j] - 1 / b4
        Beta_new[, j] <- solve(X %*% Wb4 %*% t(X)) %*% X %*% Wb4 %*% r
      }
      cat("Beta:", Beta_new, "\n")
      # 更新λ
      lambda_new <- (1 + 2 * Bb4 + sqrt((1 + 2 * Bb4)^2 - 4 * Bb4 * (3 - Ba4))) / (2 * (3 - Ba4))
      cat("Lambda:", lambda_new, "\n")
      # 更新似然函数
      loglikeli_new <- log_f4(y, X, alpha_new, Beta_new, lambda_new)
      cat("Loglikelihood:", loglikeli_new, "\n")
    } else if (Type == "LN") {
      Moment <- f5(y, X, alpha, Beta, lambda)
      b5 <- Moment[, 1]
      d5 <- Moment[, 2]
      s5 <- Moment[, 3]
      Bb5 <- mean(b5)
      Bd5 <- mean(d5)
      Wb5 <- diag(b5)
      Bs5 <- mean(s5)
      # 更新αj
      c <- 1 + mean(log(y) - log(mu) - apply(y, 2, function(x) b5 * x) / mu) - Bd5
      alpha_new <- US(alpha, c)
      cat("Alpha:", alpha_new, "\n")
      # 更新βj
      for (j in 1:d) {
        r <- log(mu[, j]) + y[, j] / mu[, j] - 1 / b5
        Beta_new[, j] <- solve(X %*% Wb5 %*% t(X)) %*% X %*% Wb5 %*% r
      }
      cat("Beta:", Beta_new, "\n")
      # 更新λ
      lambda_new <- 2 * (sqrt(1 + Bs5) - 1)
      cat("Lambda:", lambda_new, "\n")
      # 更新似然函数
      loglikeli_new <- log_f5(y, X, alpha_new, Beta_new, lambda_new)
      cat("Loglikelihood:", loglikeli_new, "\n")
    }
    
    if (abs((loglikeli_new - loglikeli) / loglikeli) < ep) {
      index <- 1
      break
    }
    cat("Iteration number:", k, "\n")
    k <- k + 1
  }
  
  result <- list(Alpha = alpha_new, 
                 Beta = Beta_new, 
                 Lambda = lambda_new, 
                 Loglikelihood = loglikeli_new,
                 number = k, 
                 index = index)
  return(result)
}

X_train <- as.matrix(train_data %>% select(ALT, AST))
Z_train <- as.matrix(cbind(Intercept = 1, train_data %>% select(Age, Gender, BMI, HbA1c, SysBP, Alcohol_Cups)))

X_test <- as.matrix(test_data %>% select(ALT, AST))
Z_test <- as.matrix(cbind(Intercept = 1, test_data %>% select(Age, Gender, BMI, HbA1c, SysBP, Alcohol_Cups)))

n_train <- nrow(X_train)
d <- ncol(X_train)
p <- ncol(Z_train)
p.gega <- p * d + 2

alpha <- 1
Var.tau <- 1
Beta <- matrix(0, nrow = p, ncol = d)

IGaGa.fit <- IGaGa.MLE(alpha, Beta, Var.tau, X_train, Z_train)
GIGGa.fit <- GIGGa.MLE(alpha, Beta, Var.tau, X_train, Z_train)
LNGa.fit <- LNGa.MLE(alpha, Beta, Var.tau, X_train, t(Z_train))
IGauGa.fit <- IGauGa.MLE(alpha, Beta, Var.tau, X_train, t(Z_train))
RIGGa.fit <- RIGGa.MLE(alpha, Beta, Var.tau, X_train, t(Z_train))

# 训练集效果
IGaGa.AIC <- p.gega * 2 - 2 * IGaGa.fit$Loglikelihood
GIGGa.AIC <- p.gega * 2 - 2 * GIGGa.fit$Loglikelihood
LNGa.AIC <- p.gega * 2 - 2 * LNGa.fit$Loglikelihood
IGauGa.AIC <- p.gega * 2 - 2 * IGauGa.fit$Loglikelihood
RIGGa.AIC <- p.gega * 2 - 2 * RIGGa.fit$Loglikelihood

IGaGa.BIC <- p.gega * log(n_train) - 2 * IGaGa.fit$Loglikelihood
GIGGa.BIC <- p.gega * log(n_train) - 2 * GIGGa.fit$Loglikelihood
LNGa.BIC <- p.gega * log(n_train) - 2 * LNGa.fit$Loglikelihood
IGauGa.BIC <- p.gega * log(n_train) - 2 * IGauGa.fit$Loglikelihood
RIGGa.BIC <- p.gega * log(n_train) - 2 * RIGGa.fit$Loglikelihood

results_train[1, ] <- list("IGaGa", IGaGa.fit$Loglikelihood, IGaGa.AIC, IGaGa.BIC)
results_train[2, ] <- list("GIGGa", GIGGa.fit$Loglikelihood, GIGGa.AIC, GIGGa.BIC)
results_train[3, ] <- list("LNGa", LNGa.fit$Loglikelihood, LNGa.AIC, LNGa.BIC)
results_train[4, ] <- list("IGauGa", IGauGa.fit$Loglikelihood, IGauGa.AIC, IGauGa.BIC)
results_train[5, ] <- list("RIGGa", RIGGa.fit$Loglikelihood, RIGGa.AIC, RIGGa.BIC)
results_train[6, ] <- list("Independent Gamma", ll.gamma, AIC.gamma, BIC.gamma)
results_train[7, ] <- list("Independent IG", ll.ig, AIC.ig, BIC.ig)

# 测试集效果
IGaGa.loglikeli <- function(x, z, alpha, Beta, lambda) {
  d <- ncol(Beta)
  alpha.plus <- d * alpha
  f <- function(data) {
    x <- data[1:d]
    z <- data[-(1:d)]
    mu <- exp(z %*% Beta)
    term1 <- gamma(alpha.plus + lambda) * (lambda - 1)^lambda / gamma(lambda)
    term2 <- prod(alpha^alpha * x^(alpha - 1) / gamma(alpha) / mu^alpha)
    term3 <- (sum(alpha * x / mu) + lambda - 1)^(-alpha.plus - lambda)
    term1 * term2 * term3
  }
  pdf <- apply(cbind(x, z), 1, f)
  return(sum(log(pdf)))
}

GIGGa.loglikeli <- function(x, z, alpha, Beta, lambda) {
  a <- function(lambda) lambda * besselK(lambda, 1) / besselK(lambda, 0)
  b <- function(lambda) lambda * besselK(lambda, 0) / besselK(lambda, 1)
  
  d <- ncol(Beta)
  alpha.plus <- d * alpha
  f <- function(data) {
    x <- data[1:d]
    z <- data[-(1:d)]
    mu <- exp(z %*% Beta)
    term1 <- prod(alpha^alpha * x^(alpha - 1) / gamma(alpha) / mu^alpha)
    term2 <- besselK(sqrt(a(lambda) * (b(lambda) + 2 * sum(alpha * x / mu))), alpha.plus) / besselK(lambda, 0)
    term3 <- (a(lambda) / (b(lambda) + 2 * sum(alpha * x / mu)))^(alpha.plus / 2)
    term1 * term2 * term3
  }
  pdf <- apply(cbind(x, z), 1, f)
  return(sum(log(pdf)))
}

LNGa.loglikeli <- function(y, X, alpha, Beta, lambda) {
  d <- ncol(y)
  alpha.plus <- d * alpha
  f <- function(data) {
    y <- data[1:d]
    X <- data[-(1:d)]
    mu <- exp(t(X) %*% Beta)
    h5 <- function(tau) {
      tau^(-alpha.plus - 1) * exp(-sum(alpha * y / mu) / tau - (log(tau) + lambda / 2)^2 / (2 * lambda))
    }
    term1 <- 1 / sqrt(2 * pi * lambda)
    term2 <- prod(alpha^alpha * y^(alpha - 1) / gamma(alpha) / mu^alpha)
    term3 <- integrate(h5, 0, Inf)$value
    term1 * term2 * term3
  }
  pdf <- apply(rbind(t(y), X), 2, f)
  return(sum(log(pdf)))
}

IGauGa.loglikeli <- function(y, X, alpha, Beta, lambda) {
  d <- ncol(y)
  alpha.plus <- d * alpha
  f <- function(data) {
    y <- data[1:d]
    X <- data[-(1:d)]
    mu <- exp(t(X) %*% Beta)
    term1 <- exp(lambda) * sqrt(2 * lambda / pi)
    term2 <- prod(alpha^alpha * y^(alpha - 1) / gamma(alpha) / mu^alpha)
    term3 <- (lambda / (lambda + 2 * sum(alpha * y / mu)))^(alpha.plus / 2 + 1 / 4)
    term4 <- besselK(sqrt(lambda * (lambda + 2 * sum(alpha * y / mu))), alpha.plus + 1 / 2, expo = FALSE)
    term1 * term2 * term3 * term4
  }
  pdf <- apply(rbind(t(y), X), 2, f)
  return(sum(log(pdf)))
}

RIGGa.loglikeli <- function(y, X, alpha, Beta, lambda) {
  d <- ncol(y)
  alpha.plus <- d * alpha
  f <- function(data) {
    y <- data[1:d]
    X <- data[-(1:d)]
    mu <- exp(t(X) %*% Beta)
    term1 <- exp(1 / (lambda - 1)) * sqrt(2 * lambda / (pi * (lambda - 1)))
    term2 <- prod(alpha^alpha * y^(alpha - 1) / gamma(alpha) / mu^alpha)
    term3 <- ((lambda / (lambda - 1)) / ((lambda * (lambda - 1))^(-1) + 2 * sum(alpha * y / mu)))^(-alpha.plus / 2 + 1 / 4)
    term4 <- lambda / (lambda - 1) * (1 / (lambda * (lambda - 1)) + 2 * sum(alpha * y / mu))
    term5 <- besselK(sqrt(term4), -alpha.plus + 1 / 2, expo = FALSE)
    term1 * term2 / term3 * term5
  }
  pdf <- apply(rbind(t(y), X), 2, f)
  return(sum(log(pdf)))
}

IGaGa.test <- IGaGa.loglikeli(X_test, Z_test, IGaGa.fit$Alpha, IGaGa.fit$Beta, IGaGa.fit$Lambda)
IGaGa.MAE <- colMeans(abs(exp(Z_test %*% IGaGa.fit$Beta) - X_test))

GIGGa.test <- GIGGa.loglikeli(X_test, Z_test, GIGGa.fit$Alpha, GIGGa.fit$Beta, GIGGa.fit$Lambda)
GIGGa.MAE <- colMeans(abs(exp(Z_test %*% GIGGa.fit$Beta) - X_test))

IGauGa.test <- IGauGa.loglikeli(X_test, t(Z_test), IGauGa.fit$Alpha, IGauGa.fit$Beta, IGauGa.fit$Lambda)
IGauGa.MAE <- colMeans(abs(exp(Z_test %*% IGauGa.fit$Beta) - X_test))

RIGGa.test <- RIGGa.loglikeli(X_test, t(Z_test), RIGGa.fit$Alpha, RIGGa.fit$Beta, RIGGa.fit$Lambda)
RIGGa.MAE <- colMeans(abs(exp(Z_test %*% RIGGa.fit$Beta) - X_test))

LNGa.test <- LNGa.loglikeli(X_test, t(Z_test), LNGa.fit$Alpha, LNGa.fit$Beta, LNGa.fit$Lambda)
LNGa.MAE <- colMeans(abs(exp(Z_test %*% LNGa.fit$Beta) - X_test))

results_test[1, ] <- list("IGaGa", IGaGa.test, 
                          -2 * IGaGa.test + 2 * p.gega, 
                          -2 * IGaGa.test + log(n_test) * p.gega,
                          IGaGa.MAE[1], IGaGa.MAE[2])
results_test[2, ] <- list("GIGGa", GIGGa.test, 
                          -2 * GIGGa.test + 2 * p.gega, 
                          -2 * GIGGa.test + log(n_test) * p.gega,
                          GIGGa.MAE[1], GIGGa.MAE[2])
results_test[3, ] <- list("LNGa", LNGa.test, 
                          -2 * LNGa.test + 2 * p.gega, 
                          -2 * LNGa.test + log(n_test) * p.gega,
                          LNGa.MAE[1], LNGa.MAE[2])
results_test[4, ] <- list("IGauGa", IGauGa.test, 
                          -2 * IGauGa.test + 2 * p.gega, 
                          -2 * IGauGa.test + log(n_test) * p.gega,
                          IGauGa.MAE[1], IGauGa.MAE[2])
results_test[5, ] <- list("RIGGa", RIGGa.test, 
                          -2 * RIGGa.test + 2 * p.gega, 
                          -2 * RIGGa.test + log(n_test) * p.gega,
                          RIGGa.MAE[1], RIGGa.MAE[2])

results_test[6, ] <- list("Independent Gamma", test_ll_gamma, 
                          -2 * test_ll_gamma + 2 * p.gamma, 
                          -2 * test_ll_gamma + log(n_test) * p.gamma,
                          mae_g_alt, mae_g_ast)

results_test[7, ] <- list("Independent IG", test_ll_ig, 
                          -2 * test_ll_ig + 2 * p.ig, 
                          -2 * test_ll_ig + log(n_test) * p.ig,
                          mae_ig_alt, mae_ig_ast)


# 模型 2: 基于 Copula 的多元广义线性模型
# ==========================================
cat("\n--- 拟合 Copula 多元 GLM (Gaussian Copula) ---\n")

# 1. 拟合边缘独立的 Gamma GLM
formula_eq <- ~ Age + Gender + BMI + HbA1c + SysBP + Alcohol_Cups
fit_m_alt <- glm(update(ALT ~ ., formula_eq), family = Gamma(link = "log"), data = train_data)
fit_m_ast <- glm(update(AST ~ ., formula_eq), family = Gamma(link = "log"), data = train_data)

# 2. 提取 Gamma 形状参数并计算 PIT (概率积分变换)，得到 U 分布 [0,1]
shape_alt <- 1 / summary(fit_m_alt)$dispersion
shape_ast <- 1 / summary(fit_m_ast)$dispersion

u_alt <- pgamma(train_data$ALT, shape=shape_alt, rate=shape_alt/fitted(fit_m_alt))
u_ast <- pgamma(train_data$AST, shape=shape_ast, rate=shape_ast/fitted(fit_m_ast))

# 极值截断(防止对数密度计算时出现 Inf)
U_train <- cbind(u_alt, u_ast)
U_train <- pmax(pmin(U_train, 1-1e-7), 1e-7)

# 3. 拟合高斯 Copula (提取依赖结构)
norm_cop <- normalCopula(dim = d, dispstr = "un") # un 表示非结构化相关矩阵
fit_cop <- fitCopula(norm_cop, pobs(U_train), method = "mpl")

# 4. 训练集联合对数似然 (Joint LogLik = 边缘 LL 之和 + Copula LL)
ll_marginals_train <- as.numeric(logLik(fit_m_alt) + logLik(fit_m_ast))
ll_cop_train <- sum(dCopula(U_train, fit_cop@copula, log = TRUE))
ll_joint_train <- ll_marginals_train + ll_cop_train

p_copula <- length(coef(fit_m_alt))*d + d + d # 边缘参数 + 边缘散布 + 3个Copula相关系数
aic_copula_train <- -2 * ll_joint_train + 2 * p_copula
bic_copula_train <- -2 * ll_joint_train + log(n_train) * p_copula

# 5. 测试集预测与指标计算
# MAE 与独立 GLM 完全相同 (因为均值方程独立)
pred_cop_alt <- predict(fit_m_alt, newdata = test_data, type = "response")
pred_cop_ast <- predict(fit_m_ast, newdata = test_data, type = "response")

mae_cop_alt <- mean(abs(test_data$ALT - pred_cop_alt))
mae_cop_ast <- mean(abs(test_data$AST - pred_cop_ast))

# 计算测试集似然
# 计算测试集的 U 分布
u_alt_test <- pgamma(test_data$ALT, shape=shape_alt, rate=shape_alt/pred_cop_alt)
u_ast_test <- pgamma(test_data$AST, shape=shape_ast, rate=shape_ast/pred_cop_ast)

U_test <- pmax(pmin(cbind(u_alt_test, u_ast_test), 1-1e-7), 1e-7)

ll_marginals_test <- sum(dgamma(test_data$ALT, shape=shape_alt, rate=shape_alt/pred_cop_alt, log=TRUE)) +
  sum(dgamma(test_data$AST, shape=shape_ast, rate=shape_ast/pred_cop_ast, log=TRUE))
ll_cop_test <- sum(dCopula(U_test, fit_cop@copula, log = TRUE))
test_ll_copula <- ll_marginals_test + ll_cop_test

test_aic_copula <- -2 * test_ll_copula + 2 * p_copula
test_bic_copula <- -2 * test_ll_copula + log(n_test) * p_copula

cat(sprintf("Copula - Train LL: %.2f | Train AIC: %.2f | Train BIC: %.2f\n", ll_joint_train, aic_copula_train, bic_copula_train))
cat(sprintf("Copula - Test LL: %.2f | Test AIC: %.2f | Test BIC: %.2f\n", test_ll_copula, test_aic_copula, test_bic_copula))

results_train[8, ] <- list("Copula.gamma", ll_joint_train, aic_copula_train, bic_copula_train)
results_test[8, ] <- list("Copula.gamma", test_ll_copula, test_aic_copula, test_bic_copula, 
                           mae_cop_alt, mae_cop_ast)

# 如果尚未安装 statmod，请先取消注释并安装
# install.packages("statmod")

cat("\n--- 拟合 Copula 多元 GLM (边际分布: 逆高斯 Inverse Gaussian) ---\n")

# ==========================================
# 1. 在训练集上拟合边缘独立的 Inverse Gaussian GLM
# ==========================================
formula_eq <- ~ Age + Gender + BMI + HbA1c + SysBP + Alcohol_Cups

fit_ig_alt <- glm(update(ALT ~ ., formula_eq), family = inverse.gaussian(link = "log"), data = train_data)
fit_ig_ast <- glm(update(AST ~ ., formula_eq), family = inverse.gaussian(link = "log"), data = train_data)


# 提取逆高斯分布的散布参数 (phi)
phi_alt <- summary(fit_ig_alt)$dispersion
phi_ast <- summary(fit_ig_ast)$dispersion

# ==========================================
# 2. 计算概率积分变换 (PIT)，生成 U 分布
# 使用 statmod::pinvgauss 稳健计算逆高斯 CDF
# ==========================================
u_ig_alt <- pinvgauss(train_data$ALT, mean = fitted(fit_ig_alt), dispersion = phi_alt)
u_ig_ast <- pinvgauss(train_data$AST, mean = fitted(fit_ig_ast), dispersion = phi_ast)

# 极值截断，防止丢入 Copula 算 log-density 时出现无穷大
U_ig_train <- pmax(pmin(cbind(u_ig_alt, u_ig_ast), 1-1e-7), 1e-7)

# ==========================================
# 3. 拟合高斯 Copula (提取依赖结构)
# ==========================================
norm_cop_ig <- normalCopula(dim = d, dispstr = "un") # 非结构化相关矩阵
fit_cop_ig <- fitCopula(norm_cop_ig, pobs(U_ig_train), method = "mpl")

# ==========================================
# 4. 训练集联合对数似然、AIC 与 BIC
# ==========================================
# 使用 statmod::dinvgauss 计算边缘对数似然
ll_ig_marg_train_alt <- sum(dinvgauss(train_data$ALT, mean = fitted(fit_ig_alt), dispersion = phi_alt, log = TRUE))
ll_ig_marg_train_ast <- sum(dinvgauss(train_data$AST, mean = fitted(fit_ig_ast), dispersion = phi_ast, log = TRUE))
ll_ig_marginals_train <- ll_ig_marg_train_alt + ll_ig_marg_train_ast

# Copula 结构对数似然
ll_ig_cop_train <- sum(dCopula(U_ig_train, fit_cop_ig@copula, log = TRUE))

# 训练集联合总似然
ll_ig_joint_train <- ll_ig_marginals_train + ll_ig_cop_train

# 参数个数 = 3*(回归系数) + 3*(散布参数 phi) + 3*(Copula相关系数)
p_ig_copula <- length(coef(fit_ig_alt)) * d + d + d 
aic_ig_copula_train <- -2 * ll_ig_joint_train + 2 * p_ig_copula
bic_ig_copula_train <- -2 * ll_ig_joint_train + log(n_train) * p_ig_copula

cat(sprintf("IG-Copula - Train LL: %.2f | Train AIC: %.2f\n", ll_ig_joint_train, aic_ig_copula_train))

# ==========================================
# 5. 测试集预测与指标计算
# ==========================================
# MAE: 预测点依然是均值 (由于 link="log"，预测结果与普通独立 IG 完全一致)
pred_ig_cop_alt <- predict(fit_ig_alt, newdata = test_data, type = "response")
pred_ig_cop_ast <- predict(fit_ig_ast, newdata = test_data, type = "response")

mae_ig_cop_alt <- mean(abs(test_data$ALT - pred_ig_cop_alt))
mae_ig_cop_ast <- mean(abs(test_data$AST - pred_ig_cop_ast))

# 计算测试集似然
# 计算测试集的 U 分布
u_ig_test_alt <- pinvgauss(test_data$ALT, mean = pred_ig_cop_alt, dispersion = phi_alt)
u_ig_test_ast <- pinvgauss(test_data$AST, mean = pred_ig_cop_ast, dispersion = phi_ast)
U_ig_test <- pmax(pmin(cbind(u_ig_test_alt, u_ig_test_ast), 1-1e-7), 1e-7)

# 测试集边缘似然和 Copula 似然
ll_ig_marg_test_alt <- sum(dinvgauss(test_data$ALT, mean = pred_ig_cop_alt, dispersion = phi_alt, log = TRUE))
ll_ig_marg_test_ast <- sum(dinvgauss(test_data$AST, mean = pred_ig_cop_ast, dispersion = phi_ast, log = TRUE))
ll_ig_marginals_test <- ll_ig_marg_test_alt + ll_ig_marg_test_ast

ll_ig_cop_test <- sum(dCopula(U_ig_test, fit_cop_ig@copula, log = TRUE))

# 测试集联合总似然
test_ll_ig_copula <- ll_ig_marginals_test + ll_ig_cop_test
test_aic_ig_copula <- -2 * test_ll_ig_copula + 2 * p_ig_copula
test_bic_ig_copula <- -2 * test_ll_ig_copula + log(n_test) * p_ig_copula

results_train[9, ] <- list("Copula.ig", ll_ig_joint_train, aic_ig_copula_train, bic_ig_copula_train)
results_test[9, ] <- list("Copula.ig", test_ll_ig_copula, test_aic_ig_copula, test_bic_ig_copula, 
                           mae_ig_cop_alt, mae_ig_cop_ast)


# # ==========================================
# # 4. 模型 3: 多元对数正态回归模型
# # ==========================================
# # 构建训练集和测试集的对数响应矩阵
# Y_train_log <- log(as.matrix(train_data[, c("ALT", "AST")]))
# Y_test_log  <- log(Y_test)
# 
# # 在训练集拟合多元线性回归
# fit_mln <- lm(Y_train_log ~ Age + Gender + BMI + HbA1c + SysBP + Alcohol_Cups, data = train_data)
# 
# # 提取训练集的协方差矩阵 Sigma
# Sigma_hat <- crossprod(residuals(fit_mln)) / nrow(train_data)
# 
# # 训练集效果
# ll.mln <- - (n_train * d / 2) * log(2 * pi) - (n_train / 2) * log(det(Sigma_hat)) - (n_train * d / 2) - sum(Y_train_log) 
# p.mln <- nrow(coef(fit_mln)) * d + d * (d + 1) / 2
# AIC.mln <- -2 * ll.mln + 2 * p.mln
# BIC.mln <- -2 * ll.mln + log(n_train) * p.mln
# 
# # 在测试集上预测对数均值
# pred_mln_log <- predict(fit_mln, newdata = test_data)
# 
# # 计算 MAE (需要转换回原始尺度)
# # 对数正态期望 E(X) = exp(mu + diag(Sigma)/2)
# pred_mln_orig <- exp(pred_mln_log + matrix(rep(diag(Sigma_hat)/2, n_test), nrow=n_test, byrow=TRUE))
# 
# mae_mln_alt <- mean(abs(Y_test[, "ALT"] - pred_mln_orig[, 1]))
# mae_mln_ast <- mean(abs(Y_test[, "AST"] - pred_mln_orig[, 2]))
# 
# # 计算测试集对数似然值
# # 计算多元正态密度 + 雅可比调整 (Jacobian Adjustment: -sum(log(y)))
# ll_mln_log_scale <- 0
# for (i in 1:n_test) {
#   ll_mln_log_scale <- ll_mln_log_scale + dmvnorm(Y_test_log[i, ], mean = pred_mln_log[i, ], sigma = Sigma_hat, log = TRUE)
# }
# test_ll_mln <- ll_mln_log_scale - sum(Y_test_log)
# 
# results_train[10, ] <- list("Multivariate Lognormal", ll.mln, AIC.mln, BIC.mln)
# 
# results_test[10, ] <- list("Multivariate Lognormal", test_ll_mln, 
#                           -2 * test_ll_mln + 2 * p.mln, 
#                           -2 * test_ll_mln + log(n_test) * p.mln,
#                           mae_mln_alt, mae_mln_ast)
# 
# # 模型 1: 多元广义线性混合模型 (Multivariate GLMM)
# # ==========================================
# cat("\n--- 拟合多元 GLMM (Shared Random Effect) ---\n")
# 
# # 1. 宽表转长表 (GLMM 必需)
# train_long <- pivot_longer(train_data, cols=c("ALT", "AST"), 
#                            names_to="Response", values_to="Value")
# test_long  <- pivot_longer(test_data, cols=c("ALT", "AST"), 
#                            names_to="Response", values_to="Value")
# 
# # 2. 拟合模型
# # - 均值方程：Response 及其与协变量的交互项，确保每个酶有自己的 \beta
# # - 随机效应：(1 | ID) 表示同一个患者拥有一个共享的随机截距（模仿你的 \tau）
# # - 散布参数：dispformula = ~ Response 允许三个酶拥有不同的 Gamma 形状参数 \alpha
# fit_glmm <- glmmTMB(Value ~ Response + Response:(Age + Gender + BMI + HbA1c + SysBP + Alcohol_Cups) + (1 | ID),
#                     dispformula = ~ Response,
#                     family = Gamma(link = "log"), 
#                     data = train_long, 
#                     control = glmmTMBControl(
#                       optCtrl = list(
#                         iter.max = 1000,    # 增加最大迭代次数（默认300）
#                         eval.max = 2000,    # 增加函数评估次数
#                         rel.tol = 1e-6,     # 提高收敛精度（默认1e-10）
#                         x.tol = 1e-8        # 参数变化容差
#                       )
#                     ))
# 
# # 3. 提取训练集指标
# ll_glmm_train <- logLik(fit_glmm)
# p_glmm <- attr(ll_glmm_train, "df") # 提取参数个数
# aic_glmm_train <- AIC(fit_glmm)
# bic_glmm_train <- BIC(fit_glmm)
# 
# # 4. 测试集预测与指标计算
# # 预测均值 (注意：为了严格对比，这里仅预测固定效应的期望，因为新患者的随机效应不可知)
# pred_glmm_test <- predict(fit_glmm, newdata = test_long, type = "response", re.form = NA, allow.new.levels = TRUE)
# test_long$Pred <- pred_glmm_test
# 
# # 计算测试集 MAE
# mae_glmm_alt <- mean(abs(test_long$Value[test_long$Response=="ALT"] - test_long$Pred[test_long$Response=="ALT"]))
# mae_glmm_ast <- mean(abs(test_long$Value[test_long$Response=="AST"] - test_long$Pred[test_long$Response=="AST"]))
# 
# # 计算测试集对数似然 (Test LogLik)
# # 在 GLMM 中，精确的测试集似然需要对新个体的随机效应进行积分，但为了与你的 Ge-Ga 平级对比，
# # 我们可以直接计算新数据在模型下的边缘似然。使用模型的 dispersion 估值。
# cat("\n--- 开始执行 GLMM 边缘似然的蒙特卡洛积分计算 ---\n")
# 
# # 1. 提取随机效应的标准差 (sigma_u)
# # 从 glmmTMB 模型对象中提取患者共享随机截距的标准差
# sigma_u <- attr(VarCorr(fit_glmm)$cond$ID, "stddev")
# 
# # 2. 计算测试集上的固定效应部分和散布参数
# # type="link" 返回 Z*beta (对数尺度上的均值预测)
# test_long$eta_fixed <- predict(fit_glmm, newdata = test_long, type = "link", re.form = NA, allow.new.levels = TRUE)
# # type="disp" 返回特定于各维度(ALT/AST/ALP)的散布参数 phi
# test_long$phi <- predict(fit_glmm, newdata = test_long, type = "disp", allow.new.levels = TRUE)
# test_long$shape <- 1 / test_long$phi
# 
# # 获取测试集中所有的独立患者 ID
# unique_ids <- unique(test_long$ID)
# 
# # 蒙特卡洛抽样次数 (10000 次已经足够达到非常高的精度)
# K <- 10000
# 
# # 定义稳健的 Log-Sum-Exp 函数以防止数值下溢 (Underflow)
# # 计算 log(mean(exp(x))) 的安全版本
# log_mean_exp <- function(x) {
#   max_x <- max(x)
#   # 提取最大值以确保 exp() 内的数值不会极度负向导致全为 0
#   max_x + log(mean(exp(x - max_x)))
# }
# 
# # 3. 循环计算每个患者的边缘似然
# ll_glmm_test <- 0
# 
# # 开启进度条 (因为硬算需要一点时间)
# pb <- txtProgressBar(min = 0, max = length(unique_ids), style = 3)
# 
# for (i in seq_along(unique_ids)) {
#   current_id <- unique_ids[i]
#   
#   # 提取该患者在 3 个维度上的真实值、固定预测子、形状参数
#   patient_data <- test_long[test_long$ID == current_id, ]
#   
#   # 从 N(0, sigma_u^2) 中随机抽取 K 个潜在的随机效应 u_i
#   u_sim <- rnorm(K, mean = 0, sd = sigma_u)
#   
#   # 建立一个长度为 K 的向量，储存该患者在这 K 种宇宙下的联合对数似然
#   log_L_k <- numeric(K)
#   
#   # 针对该患者的每个生化指标 (ALT, AST, ALP) 进行累加
#   for (j in 1:nrow(patient_data)) {
#     y_val <- patient_data$Value[j]
#     shape_val <- patient_data$shape[j]
#     
#     # 结合固定效应与模拟的随机效应，反解出均值 mu
#     # mu_k = exp(Z*beta + u)
#     mu_k <- exp(patient_data$eta_fixed[j] + u_sim)
#     
#     # 计算条件对数似然
#     # 注意 rate = shape / mu
#     ll_cond <- dgamma(y_val, shape = shape_val, rate = shape_val / mu_k, log = TRUE)
#     
#     # 各维度条件独立，故将其相加
#     log_L_k <- log_L_k + ll_cond
#   }
#   
#   # 针对这 K 个模拟的联合对数似然，执行蒙特卡洛积分求均值
#   # L_i = 1/K * sum(L_i_k)  =>  log(L_i) = log_mean_exp(log_L_k)
#   patient_marginal_ll <- log_mean_exp(log_L_k)
#   
#   # 累加到总测试集似然中
#   ll_glmm_test <- ll_glmm_test + patient_marginal_ll
#   
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# 
# # 4. 重新计算 GLMM 的真实测试集指标
# aic_glmm_test <- -2 * ll_glmm_test + 2 * p_glmm
# bic_glmm_test <- -2 * ll_glmm_test + log(n_test) * p_glmm
# 
# results_train[11, ] <- list("GLMM", ll_glmm_train, aic_glmm_train, bic_glmm_train)
# results_test[11, ] <- list("GLMM", ll_glmm_test, aic_glmm_test, bic_glmm_test, 
#                           mae_glmm_alt, mae_glmm_ast)


# 5. 输出最终结果
# ==========================================
print("测试集评估结果 (Test Set Performance):")
print(results_train)
print(results_test)

results_train[, 2:4] %>% round(digits = 4)
results_test[, 2:6] %>% round(digits = 4)

