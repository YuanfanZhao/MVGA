# 导入需要的package
suppressMessages(suppressWarnings(library(tidyverse)))
suppressWarnings(suppressMessages(library(SuppDists)))
suppressWarnings(suppressMessages(library(Bessel)))
suppressWarnings(suppressMessages(library(GIGrvg)))
suppressMessages(suppressWarnings(library(copula)))
suppressWarnings(suppressMessages(library(fitdistrplus)))
suppressWarnings(suppressMessages(library(actuar)))

# CR多元Gamma分布的EM算法实现
# ---------------------------------------------------------

#' 牛顿法求解 alpha: C - digamma(alpha) = 0
#' @param C 目标常数
#' @param init_val alpha的初始值
#' @param tol 容差
solve_alpha_newton <- function(C, init_val = 1, tol = 1e-6) {
  alpha <- init_val
  for (iter in 1:100) {
    # 目标函数 f(alpha) = digamma(alpha) - C = 0
    f_val <- digamma(alpha) - C
    # 导数 f'(alpha) = trigamma(alpha)
    f_deriv <- trigamma(alpha)
    
    alpha_new <- alpha - f_val / f_deriv
    
    # 防止alpha变为非正数
    if (alpha_new <= 0) {
      alpha_new <- alpha / 2 
    }
    
    if (abs(alpha_new - alpha) < tol) {
      return(alpha_new)
    }
    alpha <- alpha_new
  }
  return(alpha)
}

#' E步：计算单个样本的条件期望和对数似然
#' 基于密度 f(s|x_i, theta) ∝ s^(alpha0-1) * prod((x_ij - s)^(alphaj-1)) * exp((d-1)*beta*s)
e_step <- function(x_i, alpha0, alpha, beta) {
  d <- length(x_i)
  m_i <- min(x_i) # 积分上限 m_i = min(x_i1, ..., x_id)
  
  # 定义未归一化的对数后验密度函数，使用对数防止数值溢出
  log_f <- function(s) {
    # 避免边界产生log(0)的错误
    s <- pmax(pmin(s, m_i - 1e-10), 1e-10)
    term1 <- (alpha0 - 1) * log(s)
    term2 <- sum((alpha - 1) * log(x_i - s))
    term3 <- (d - 1) * beta * s
    return(term1 + term2 + term3)
  }
  log_f_vec <- Vectorize(log_f)
  
  # 寻找最大值以进行数值稳定化 (Log-Sum-Exp 技巧连续版)
  opt <- optimize(log_f_vec, interval = c(0, m_i), maximum = TRUE)
  max_log_f <- opt$objective
  
  # 稳定的未归一化密度函数
  f_unnorm <- function(s) {
    exp(log_f_vec(s) - max_log_f)
  }
  
  # 计算积分项
  # 注意：在实际复杂数据中，integrate函数可能会遇到细微的计算瓶颈，
  # 若有报错可替换为高斯积分或网格积分。
  I_norm <- integrate(f_unnorm, lower = 0, upper = m_i, stop.on.error = FALSE)$value
  I_s <- integrate(function(s) s * f_unnorm(s), lower = 0, upper = m_i, stop.on.error = FALSE)$value
  I_logs <- integrate(function(s) log(s) * f_unnorm(s), lower = 0, upper = m_i, stop.on.error = FALSE)$value
  
  I_logx_s <- numeric(d)
  for(j in 1:d) {
    I_logx_s[j] <- integrate(function(s) log(x_i[j] - s) * f_unnorm(s), 
                             lower = 0, upper = m_i, stop.on.error = FALSE)$value
  }
  
  # 计算条件期望 E(Y_i0 | x_i, theta) 等
  E_s <- I_s / I_norm
  E_logs <- I_logs / I_norm
  E_logx_s <- I_logx_s / I_norm
  
  # 计算该样本的边缘对数似然贡献 (Observed Log-Likelihood)
  # L(x_i) = int f_com(x_i, s) ds 
  log_C <- (alpha0 + sum(alpha)) * log(beta) - lgamma(alpha0) - sum(lgamma(alpha)) - beta * sum(x_i)
  log_L_i <- log_C + max_log_f + log(I_norm)
  
  return(list(E_s = E_s, E_logs = E_logs, E_logx_s = E_logx_s, log_L_i = log_L_i))
}

#' 主函数：多元Gamma分布的EM算法
#' @param X 数据矩阵 (n行 d列)
#' @param max_iter 最大迭代次数
#' @param tol 收敛容差
cr_gamma_em <- function(X, max_iter = 1000, tol = 1e-3) {
  n <- nrow(X)
  d <- ncol(X)
  
  # 1. 粗略初始化参数
  means <- colMeans(X)
  vars <- apply(X, 2, var)
  beta_init <- mean(means / vars)
  alpha_init <- means * beta_init
  alpha0_init <- min(alpha_init) / 3 
  alpha_init <- alpha_init - alpha0_init
  
  theta <- list(alpha0 = alpha0_init, alpha = alpha_init, beta = beta_init)
  log_lik_history <- numeric()
  
  cat("开始EM算法迭代...\n")
  
  for (iter in 1:max_iter) {
    E_s <- numeric(n)
    E_logs <- numeric(n)
    E_logx_s <- matrix(0, nrow = n, ncol = d)
    log_L_total <- 0
    
    # ---------------- E步 ----------------
    for (i in 1:n) {
      res <- e_step(X[i, ], theta$alpha0, theta$alpha, theta$beta)
      E_s[i] <- res$E_s
      E_logs[i] <- res$E_logs
      E_logx_s[i, ] <- res$E_logx_s
      log_L_total <- log_L_total + res$log_L_i
    }
    
    log_lik_history <- c(log_lik_history, log_L_total)
    
    # 打印当前迭代信息
    cat(sprintf("Iter %d: Log-Likelihood = %.4f\n", iter, log_L_total))
    
    # 检查收敛
    if (iter > 1 && abs(log_L_total - log_lik_history[iter - 1]) / abs(log_lik_history[iter - 1]) < tol) {
      cat("算法已收敛。\n")
      break
    }
    
    # ---------------- M步 ----------------
    # 基于完全数据对数似然的偏导数等于0的方程进行更新
    
    # 更新 alpha0
    C0 <- log(theta$beta) + mean(E_logs)
    theta$alpha0 <- solve_alpha_newton(C0, init_val = theta$alpha0)
    
    # 更新 alpha_j (j = 1...d)
    for (j in 1:d) {
      Cj <- log(theta$beta) + mean(E_logx_s[, j])
      theta$alpha[j] <- solve_alpha_newton(Cj, init_val = theta$alpha[j])
    }
    
    # 更新 beta
    sum_alpha_all <- theta$alpha0 + sum(theta$alpha)
    sum_X <- sum(X)
    sum_E_s <- sum(E_s)
    # 根据文档偏导数推导出的 beta 解析解
    theta$beta <- (n * sum_alpha_all) / (sum_X - (d - 1) * sum_E_s)
  }
  
  return(list(
    MLE_Estimates = theta,
    Final_Log_Likelihood = log_L_total,
    Iterations = iter,
    History = log_lik_history
  ))
}

# MVGA拟合函数
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

LNGa.MLE <- function(alpha, Beta, Var.tau, y, X, Type = "LN", ep = 1e-8) {
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

# 多元对数正态分布
fit_multivariate_lognormal <- function(data) {
  # 转换为矩阵，并检查数据
  data <- as.matrix(data)
  if (any(data <= 0)) {
    stop("所有数据必须为正数，因为对数正态分布要求 x > 0")
  }
  
  n <- nrow(data)
  d <- ncol(data)
  
  # 取对数
  log_data <- log(data)
  
  # 极大似然估计：多元正态的MLE
  mu_hat <- colMeans(log_data)
  # 注意：cov() 返回无偏估计（除以 n-1），MLE应除以 n
  Sigma_hat <- cov(log_data) * (n - 1) / n
  
  # 计算对数似然
  # 1. 多元正态部分的对数似然
  # 使用Cholesky分解提高数值稳定性
  chol_Sigma <- chol(Sigma_hat)
  log_det_Sigma <- 2 * sum(log(diag(chol_Sigma)))
  
  # 计算马氏距离
  # 对每个观测，计算 (logx - mu)' Sigma^{-1} (logx - mu)
  # 利用前向/后向替换求解
  centered <- t(t(log_data) - mu_hat)  # 中心化矩阵
  # 求解 Sigma^{-1} * centered'，更高效的方法是使用backsolve
  # 令 L = chol(Sigma)，则 Sigma = L' L，Sigma^{-1} = L^{-1} L'^{-1}
  # 马氏距离 = (L'^{-1} (x-mu))' (L'^{-1} (x-mu))
  # 即先解 L' v = (x-mu)，然后距离 = sum(v^2)
  v <- backsolve(chol_Sigma, t(centered), upper.tri = FALSE, transpose = TRUE)
  # v 是 d x n 矩阵，每列对应一个观测的 v
  mahalanobis_dist <- colSums(v^2)
  
  ll_norm <- - (n * d / 2) * log(2 * pi) - (n / 2) * log_det_Sigma - 0.5 * sum(mahalanobis_dist)
  
  # 2. 雅可比项：每个观测变换的贡献是 -sum(log(x_i))，对所有观测求和
  jacobian_term <- -sum(log_data)
  
  loglik <- ll_norm + jacobian_term
  
  # 返回结果
  list(
    mu = mu_hat,
    Sigma = Sigma_hat,
    loglik = loglik,
    n = n,
    d = d
  )
}

# Copula拟合
fit_positive_copula <- function(data, marginal_dists = rep("gamma", ncol(data)), copula_type = "gaussian") {
  drig <- function(x, mu, lambda) {
    term1 <- sqrt(lambda/(2 * pi * x))
    term2 <- exp(-lambda * (1 - mu * x)^2 / (2 * mu^2 * x))
    term1 * term2
  }
  
  prig <- function(q, mu, lambda) {
    1 - pinvgauss(1/q, mean = mu, shape = lambda)
  }
  
  qrig <- function(p, mu, lambda) {
    1 / qinvgauss(1 - p, mean = mu, shape = lambda)
  }
  
  d <- ncol(data)
  n <- nrow(data)
  results <- list()
  
  # 1. 拟合边缘分布
  marginal_fits <- vector("list", d)
  u_matrix <- matrix(0, n, d)
  marginal_loglik <- 0
  
  for (i in 1:d) {
    x <- data[, i]
    dist <- marginal_dists[i]
    
    # 根据分布类型设置正确的初始参数名称
    start_params <- switch(
      dist,
      "gamma" = list(shape = mean(x)^2/var(x), rate = mean(x)/var(x)),
      "invgamma" = list(shape = 2 + mean(x)^2/var(x), scale = mean(x)*(1 + mean(x)^2/var(x))),
      "invgauss" = list(mean = mean(x), shape = mean(x)^3/var(x)),
      "rig" = list(mu = 1/mean(x), lambda = mean(x)^3/var(x)),
      "lnorm" = list(meanlog = mean(log(x)), sdlog = sd(log(x))),
      stop("不支持的分布类型: ", dist)
    )
    
    if (dist == "rig") {
      fit <- tryCatch({
        fitdistrplus::fitdist(x, "rig", start = start_params, method = "mle", lower = c(0, 0))
      }, error = function(e) {
        NULL
      })
      if (!is.null(fit)) {
        # 转换为均匀分布
        u_matrix[, i] <- prig(x, fit$estimate["mu"], fit$estimate["lambda"])
        marginal_loglik <- marginal_loglik + fit$loglik
      } else {
        u_matrix[, i] <- pobs(x)
      }
    } else {
      fit <- fitdistrplus::fitdist(x, dist, method = "mle", start = start_params)
      
      # 转换为均匀分布
      cdf_fun <- switch(
        dist,
        "gamma" = pgamma,
        "invgamma" = pinvgamma,
        "invgauss" = pinvgauss,
        "lnorm" = plnorm,
        stop("不支持的分布类型")
      )
      
      u_matrix[, i] <- do.call(cdf_fun, c(list(x), as.list(fit$estimate)))
      marginal_loglik <- marginal_loglik + fit$loglik
      
    }
    marginal_fits[[i]] <- fit
  }
  
  # 2. 创建Copula对象
  copula_obj <- switch(
    copula_type,
    "gaussian" = normalCopula(rep(0.1, d*(d-1)/2), dim = d, dispstr = "un"),
    "t" = tCopula(rep(0.1, d*(d-1)/2), dim = d, dispstr = "un", df = 5, df.fixed = FALSE),
    "clayton" = claytonCopula(2, dim = d),
    "gumbel" = gumbelCopula(2, dim = d),
    "frank" = frankCopula(2, dim = d),
    stop("不支持的Copula类型: ", copula_type)
  )
  
  # 3. 拟合Copula参数
  fit_copula <- fitCopula(
    copula_obj,
    u_matrix,
    method = "ml"
  )
  
  # 4. 计算总对数似然
  copula_loglik <- fit_copula@loglik
  total_loglik <- marginal_loglik + copula_loglik
  
  # 5. 准备结果输出
  results$marginal_fits <- marginal_fits
  results$copula_fit <- fit_copula
  # results$u_matrix <- u_matrix
  results$total_loglik <- total_loglik
  # results$copula_loglik <- copula_loglik
  results$AIC <- -2 * total_loglik + 2 * (length(fit_copula@estimate) + sum(sapply(marginal_fits, function(f) length(f$estimate))))
  results$BIC <- -2 * total_loglik + log(n) * (length(fit_copula@estimate) + sum(sapply(marginal_fits, function(f) length(f$estimate))))
  class(results) <- "positive_copula_fit"
  return(results)
}

# 拟合多元对数t分布
t.MLE <- function(X) {
  fit <- QRM::fit.mst(X)
  return(list(
    mu = fit$mu,
    Sigma = fit$Sigma,
    nv = fit$df, 
    loglik = fit$ll.max - sum(X)
  ))
}

# 准备数据
data <- readxl::read_xlsx("数据.xlsx")
data <- data[, 1:6]
data <- data %>% filter(TOTSLF23 != 0)
data <- data %>% filter(TOTMCR23 != 0)
# data <- data %>% filter(TOTMCD23 != 0)
data <- data %>% filter(TOTPRV23 != 0)
data <- data / 1000

x <- data[, c(3, 4, 6)] %>% as.matrix()
n <- nrow(x)
d <- ncol(x)

# 数据统计
Mean <- colMeans(x)
Corr <- cor(x)
Cov <- cov(x)

alpha <- 1
Var.tau <- 1
z <- rep(1, n)
Beta <- matrix(rep(0, d), ncol = d)

IGaGa.fit <- IGaGa.MLE(alpha, Beta, Var.tau, x, z)
GIGGa.fit <- GIGGa.MLE(alpha, Beta, Var.tau, x, z)
LNGa.fit <- LNGa.MLE(alpha, Beta, Var.tau, x, t(z))
IGauGa.fit <- IGauGa.MLE(alpha, Beta, Var.tau, x, t(z))
RIGGa.fit <- RIGGa.MLE(alpha, Beta, Var.tau, x, t(z))

para <- 2 + d

IGaGa.AIC <- para * 2 - 2 * IGaGa.fit$Loglikelihood
GIGGa.AIC <- para * 2 - 2 * GIGGa.fit$Loglikelihood
LNGa.AIC <- para * 2 - 2 * LNGa.fit$Loglikelihood
IGauGa.AIC <- para * 2 - 2 * IGauGa.fit$Loglikelihood
RIGGa.AIC <- para * 2 - 2 * RIGGa.fit$Loglikelihood

IGaGa.BIC <- para * log(n) - 2 * IGaGa.fit$Loglikelihood
GIGGa.BIC <- para * log(n) - 2 * GIGGa.fit$Loglikelihood
LNGa.BIC <- para * log(n) - 2 * LNGa.fit$Loglikelihood
IGauGa.BIC <- para * log(n) - 2 * IGauGa.fit$Loglikelihood
RIGGa.BIC <- para * log(n) - 2 * RIGGa.fit$Loglikelihood

Lognormal.fit <- fit_multivariate_lognormal(x)
Lognormal.AIC <- (d + sum(1:d)) * 2 - 2 * Lognormal.fit$loglik
Lognormal.BIC <- (d + sum(1:d)) * log(n) - 2 * Lognormal.fit$loglik

logx <- log(x)
t.fit <- t.MLE(logx)
t.AIC <- (d + sum(1:d) + 1) * 2 - 2 * t.fit$loglik
t.BIC <- (d + sum(1:d) + 1) * log(n) - 2 * t.fit$loglik

Copula1.fit <- fit_positive_copula(data = x, marginal_dists = c("gamma", "gamma", "gamma"), copula_type = "gaussian")
Copula2.fit <- fit_positive_copula(data = x, marginal_dists = c("invgamma", "invgamma", "invgamma"), copula_type = "gaussian")
Copula3.fit <- fit_positive_copula(data = x, marginal_dists = c("invgauss", "invgauss", "invgauss"), copula_type = "gaussian")
Copula4.fit <- fit_positive_copula(data = x, marginal_dists = c("lnorm", "lnorm", "lnorm"), copula_type = "gaussian")

Copula5.fit <- fit_positive_copula(data = x, marginal_dists = c("gamma", "gamma", "gamma"), copula_type = "t")
Copula6.fit <- fit_positive_copula(data = x, marginal_dists = c("invgamma", "invgamma", "invgamma"), copula_type = "t")
Copula7.fit <- fit_positive_copula(data = x, marginal_dists = c("invgauss", "invgauss", "invgauss"), copula_type = "t")
Copula8.fit <- fit_positive_copula(data = x, marginal_dists = c("lnorm", "lnorm", "lnorm"), copula_type = "t")

Copula9.fit <- fit_positive_copula(data = x, marginal_dists = c("gamma", "gamma", "gamma"), copula_type = "clayton")
Copula10.fit <- fit_positive_copula(data = x, marginal_dists = c("invgamma", "invgamma", "invgamma"), copula_type = "clayton")
Copula11.fit <- fit_positive_copula(data = x, marginal_dists = c("invgauss", "invgauss", "invgauss"), copula_type = "clayton")
Copula12.fit <- fit_positive_copula(data = x, marginal_dists = c("lnorm", "lnorm", "lnorm"), copula_type = "clayton")

Copula13.fit <- fit_positive_copula(data = x, marginal_dists = c("gamma", "gamma", "gamma"), copula_type = "gumbel")
Copula14.fit <- fit_positive_copula(data = x, marginal_dists = c("invgamma", "invgamma", "invgamma"), copula_type = "gumbel")
Copula15.fit <- fit_positive_copula(data = x, marginal_dists = c("invgauss", "invgauss", "invgauss"), copula_type = "gumbel")
Copula16.fit <- fit_positive_copula(data = x, marginal_dists = c("lnorm", "lnorm", "lnorm"), copula_type = "gumbel")

Copula17.fit <- fit_positive_copula(data = x, marginal_dists = c("gamma", "gamma", "gamma"), copula_type = "frank")
Copula18.fit <- fit_positive_copula(data = x, marginal_dists = c("invgamma", "invgamma", "invgamma"), copula_type = "frank")
Copula19.fit <- fit_positive_copula(data = x, marginal_dists = c("invgauss", "invgauss", "invgauss"), copula_type = "frank")
Copula20.fit <- fit_positive_copula(data = x, marginal_dists = c("lnorm", "lnorm", "lnorm"), copula_type = "frank")

Copula.loglikelihood <- c(Copula1.fit$total_loglik, Copula2.fit$total_loglik, Copula3.fit$total_loglik, Copula4.fit$total_loglik, 
                          Copula5.fit$total_loglik, Copula6.fit$total_loglik, Copula7.fit$total_loglik, Copula8.fit$total_loglik, 
                          Copula9.fit$total_loglik, Copula10.fit$total_loglik, Copula11.fit$total_loglik, Copula12.fit$total_loglik, 
                          Copula13.fit$total_loglik, Copula14.fit$total_loglik, Copula15.fit$total_loglik, Copula16.fit$total_loglik, 
                          Copula17.fit$total_loglik, Copula18.fit$total_loglik, Copula19.fit$total_loglik, Copula20.fit$total_loglik)

which.max(Copula.loglikelihood)

CR.fit <- cr_gamma_em(x, tol = 1e-5)
CR.AIC <- 2 * 5 - 2 * CR.fit$Final_Log_Likelihood
CR.BIC <- 5 * log(n) - 2 * CR.fit$Final_Log_Likelihood

Result <- data.frame(Distribution = c("IGaGa", "GIGGa", "LNGa", "IGauGa", "RIGGa", "Lognormal", "CR", "logt", "Copula"), 
                     Loglikelihood = c(IGaGa.fit$Loglikelihood, GIGGa.fit$Loglikelihood, LNGa.fit$Loglikelihood, IGauGa.fit$Loglikelihood, RIGGa.fit$Loglikelihood, Lognormal.fit$loglik, CR.fit$Final_Log_Likelihood, Copula8.fit$total_loglik, t.fit$loglik) %>% round(digits = 4), 
                     AIC = c(IGaGa.AIC, GIGGa.AIC, LNGa.AIC, IGauGa.AIC, RIGGa.AIC, Lognormal.AIC, CR.AIC, Copula8.fit$AIC, t.AIC) %>% round(digits = 4),
                     BIC = c(IGaGa.BIC, GIGGa.BIC, LNGa.BIC, IGauGa.BIC, RIGGa.BIC, Lognormal.BIC, CR.BIC, Copula8.fit$BIC, t.BIC) %>% round(digits = 4))
Result %>% knitr::kable()

# Mean
abs(exp(IGaGa.fit$Beta) - colMeans(x)) %>% sum() %>% round(digits = 4)
abs(exp(GIGGa.fit$Beta) - colMeans(x)) %>% sum() %>% round(digits = 4)
abs(exp(LNGa.fit$Beta) - colMeans(x)) %>% sum() %>% round(digits = 4)
abs(exp(IGauGa.fit$Beta) - colMeans(x)) %>% sum() %>% round(digits = 4)
abs(exp(RIGGa.fit$Beta) - colMeans(x)) %>% sum() %>% round(digits = 4)
exp(Lognormal.fit$mu + 0.5 * sqrt(diag(Lognormal.fit$Sigma))) %>% round(digits = 4)
((CR.fit$MLE_Estimates$alpha0 + CR.fit$MLE_Estimates$alpha) / CR.fit$MLE_Estimates$beta) %>% round(digits = 4)

# MLE and CI
# GIG-Ga的参数估计

# GIG反解λ=====================================================================
GIG.lambda <- function(Var.tau) {
  f <- function(lambda) {
    besselK(lambda, 2) * besselK(lambda, 0) / (besselK(lambda, 1))^2 - 1 - Var.tau
  }
  return(uniroot(f, c(1, 5), extendInt = "yes")$root)
}
a <- function(lambda) lambda * besselK(lambda, 1) / besselK(lambda, 0)
b <- function(lambda) lambda * besselK(lambda, 0) / besselK(lambda, 1)

# 生成IGaGa随机数===============================================================
rGIGGa <- function(n, z, alpha, Beta, Var.tau) {
  # 确定τ分布的参数λ并生成τ分布的随机数
  lambda <- GIG.lambda(Var.tau)
  tau <- rgig(n, 0, b(lambda), a(lambda))
  
  # 计算Gamma分布的速率参数
  d <- ncol(Beta)
  rate.para <- alpha / apply(exp(z %*% Beta), 2, function(x) x * tau)
  
  # 生成样本x
  x <- matrix(0, n, d)
  for (i in 1:d) {
    x[, i] <- sapply(rate.para[, i], function(s) rgamma(1, alpha, s))
  }
  return(x)
}

# 定义对数似然函数
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

# NEM算法估计MLE================================================================
GIGGa.MLE <- function(alpha, Beta, Var.tau, x, z, ep = 1e-8) {
  # 收敛指标与迭代次数
  index <- 0
  k <- 1
  
  # 样本量n
  n <- nrow(x)
  d <- ncol(x)
  
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

# =========================================================
# 参数打包：alpha, Beta(p×d), lambda -> theta 向量
# =========================================================
pack_theta <- function(alpha, Beta, lambda) {
  c(alpha, as.vector(Beta), lambda)
}

# =========================================================
# 参数解包：theta -> alpha, Beta, lambda
# =========================================================
unpack_theta <- function(theta, p, d) {
  alpha  <- theta[1]
  Beta   <- matrix(theta[2:(1 + p * d)], nrow = p, ncol = d)
  lambda <- theta[length(theta)]
  list(alpha = alpha, Beta = Beta, lambda = lambda)
}

# =========================================================
# 标准参数向量形式的 log-likelihood
# =========================================================
loglik_theta <- function(theta, x, z, p, d) {
  par <- unpack_theta(theta, p, d)
  log_f(
    x      = x,
    z      = z,
    alpha  = par$alpha,
    Beta   = par$Beta,
    lambda = par$lambda
  )
}

# =========================================================
# 每个观测的 score contribution
# 输出: n × dim(theta) 矩阵
# =========================================================
score_contributions_theta <- function(theta, x, z, p, d) {
  n <- nrow(x)
  dim_theta <- length(theta)
  
  score_mat <- matrix(0, nrow = n, ncol = dim_theta)
  
  for (i in 1:n) {
    grad_i <- numDeriv::grad(
      func = function(th) {
        par <- unpack_theta(th, p, d)
        
        # 单个观测 log-density
        xi <- matrix(x[i, ], nrow = 1)
        zi <- matrix(z[i, ], nrow = 1)
        
        log_f(
          x      = xi,
          z      = zi,
          alpha  = par$alpha,
          Beta   = par$Beta,
          lambda = par$lambda
        )
      },
      theta,
      method = "Richardson",
      method.args = list(eps = 1e-10, d = 1e-4, r = 6)
    )
    
    score_mat[i, ] <- grad_i
  }
  
  colnames(score_mat) <- c(
    "alpha",
    paste0("beta_", 1:(p * d)),
    "lambda"
  )
  
  return(score_mat)
}

# =========================================================
# Hessian matrix of log-likelihood
# =========================================================
hessian_theta <- function(theta, x, z, p, d) {
  numDeriv::hessian(
    func = function(th) loglik_theta(th, x, z, p, d),
    theta,
    method = "Richardson",
    method.args = list(eps = 1e-10, d = 1e-4, r = 6)
  )
}

# =========================================================
# Sandwich covariance matrix
# =========================================================
sandwich_covariance <- function(theta_hat, x, z, p, d) {
  # score contributions
  S <- score_contributions_theta(theta_hat, x, z, p, d)
  
  # B matrix
  B <- crossprod(S)
  
  # Hessian
  H <- -hessian_theta(theta_hat, x, z, p, d)
  
  # Sandwich
  Sigma <- solve(H) %*% B %*% solve(H)
  
  return(Sigma)
}

# =========================================================
# Asymptotic confidence intervals
# =========================================================
asymptotic_CI <- function(theta_hat, Sigma, level = 0.95) {
  z_alpha <- qnorm(1 - (1 - level) / 2)
  se <- sqrt(diag(Sigma))
  
  lower <- theta_hat - z_alpha * se
  upper <- theta_hat + z_alpha * se
  
  list(
    lower = lower,
    upper = upper,
    se    = se
  )
}

delta_CI_Vartau <- function(lambda_hat, var_lambda, level = 0.95) {
  z_alpha <- qnorm(1 - (1 - level) / 2)
  
  g_hat <- besselK(lambda_hat, 2) * besselK(lambda_hat, 0) / (besselK(lambda_hat, 1))^2 - 1
  g_prime <- -(besselK(lambda_hat, 1) + besselK(lambda_hat, 3)) * besselK(lambda_hat, 0) / 2 / (besselK(lambda_hat, 1))^2 - besselK(lambda_hat, 2) / besselK(lambda_hat, 1) + besselK(lambda_hat, 2) * besselK(lambda_hat, 0) * (besselK(lambda_hat, 0) + besselK(lambda_hat, 2)) / (besselK(lambda_hat, 1))^3
  
  se_g <- abs(g_prime) * sqrt(var_lambda)
  
  lower <- g_hat - z_alpha * se_g
  upper <- g_hat + z_alpha * se_g
  
  c(lower = lower, upper = upper)
}

# Import the data
data <- readxl::read_xlsx("数据.xlsx")
data <- data[, 1:6]
data <- data %>% filter(TOTSLF23 != 0)
data <- data %>% filter(TOTMCR23 != 0)
data <- data %>% filter(TOTPRV23 != 0)
data <- data / 1000

x <- data[, c(3, 4, 6)] %>% as.matrix()
n <- nrow(x)
d <- ncol(x)

alpha <- 1
Var.tau <- 1
z <- rep(1, n) %>% as.matrix()
Beta <- matrix(rep(0, d), ncol = d)

Res <- GIGGa.MLE(alpha, Beta, Var.tau, x, z)


theta_hat <- pack_theta(Res$Alpha, Res$Beta, Res$Lambda)

# MLE
result_hat <- c(Res$Alpha, exp(Res$Beta), besselK(Res$Lambda, 2) * besselK(Res$Lambda, 0) / (besselK(Res$Lambda, 1))^2 - 1)
result_hat %>% round(digits = 4)

# 渐近CI
p <- ncol(z)
Sigma_hat <- sandwich_covariance(theta_hat, x, z, p, d)

CI <- asymptotic_CI(theta_hat, Sigma_hat)

# μ的置信区间
mu_hat <- exp(Res$Beta)
se_g <- mu_hat * sqrt(diag(Sigma_hat[2:4, 2:4]))
mu.lower <- mu_hat - qnorm(1 - (1 - 0.95) / 2) * se_g
mu.upper <- mu_hat + qnorm(1 - (1 - 0.95) / 2) * se_g

# Var(τ)的置信区间
idx_lambda <- length(theta_hat)
CI_Vartau <- delta_CI_Vartau(
  lambda_hat = theta_hat[idx_lambda],
  var_lambda = Sigma_hat[idx_lambda, idx_lambda]
)

lower.bound <- CI$lower
upper.bound <- CI$upper
lower.bound[2:4] <- mu.lower
upper.bound[2:4] <- mu.upper
lower.bound[idx_lambda] <- CI_Vartau[1]
upper.bound[idx_lambda] <- CI_Vartau[2]

lower.bound %>% round(digits = 4)
upper.bound %>% round(digits = 4)

# Bootstrap CI
mean.std.CI <- function(lasample) {
  G <- dim(lasample)[1]
  lasort <- apply(lasample, 2, sort)
  indexx <- floor(c(0.025 * G, 0.975 * G))
  laL <- (lasort[indexx[1], ] + lasort[indexx[1] + 1, ]) / 2
  laU <- (lasort[indexx[2], ] + lasort[indexx[2] + 1, ]) / 2
  results <- c(laL, laU)
  return(results)
}

Boot.GIGGa <- function(n, alpha, Beta, Var.tau) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(SuppDists)))
  suppressWarnings(suppressMessages(library(Bessel)))
  suppressWarnings(suppressMessages(library(GIGrvg)))
  suppressWarnings(suppressMessages(library(numDeriv)))
  suppressWarnings(suppressMessages(library(MASS)))
  rGIGGa <- function(n, z, alpha, Beta, Var.tau) {
    # 确定τ分布的参数λ并生成τ分布的随机数
    lambda <- GIG.lambda(Var.tau)
    tau <- rgig(n, 0, b(lambda), a(lambda))
    
    # 计算Gamma分布的速率参数
    d <- ncol(Beta)
    rate.para <- alpha / apply(exp(z %*% Beta), 2, function(x) x * tau)
    
    # 生成样本x
    x <- matrix(0, n, d)
    for (i in 1:d) {
      x[, i] <- sapply(rate.para[, i], function(s) rgamma(1, alpha, s))
    }
    return(x)
  }
  
  g <- 1
  G <- 6
  
  boot.estimation <- list(Alpha.boot = numeric(length = G), 
                          Mu.boot = matrix(nrow = G, ncol = 3), 
                          Var.tau.boot = numeric(length = G))
  
  while (g <= G) {
    boot.sample <- rGIGGa(n, z, alpha, Beta, Var.tau)
    Res <- tryCatch({
      GIGGa.MLE(alpha, Beta, Var.tau, boot.sample, z)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    if (!all(is.na(Res))) {
      boot.estimation[[1]][g] <- Res$Alpha
      boot.estimation[[2]][g, ] <- exp(Res$Beta)
      boot.estimation[[3]][g] <- besselK(Res$Lambda, 2) * besselK(Res$Lambda, 0) / (besselK(Res$Lambda, 1))^2 - 1
      g <- g + 1
    }
  }
  
  return(boot.estimation)
}

num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Boot.GIGGa(n, Res$Alpha, Res$Beta, besselK(Res$Lambda, 2) * besselK(Res$Lambda, 0) / (besselK(Res$Lambda, 1))^2 - 1)
stopCluster(cl)

Result.alpha <- result_part[[1]][[1]]
Result.mu <- result_part[[1]][[2]]
Result.var.tau <- result_part[[1]][[3]]

for (i in 2:num_cores) {
  Result.alpha <- c(Result.alpha, result_part[[i]][[1]])
  Result.mu <- rbind(Result.mu, result_part[[i]][[2]])
  Result.var.tau <- c(Result.var.tau, result_part[[i]][[3]])
}

mean.std.CI(Result.alpha %>% as.matrix()) %>% round(digits = 4)
mean.std.CI(Result.mu %>% as.matrix()) %>% round(digits = 4)
mean.std.CI(Result.var.tau %>% as.matrix()) %>% round(digits = 4)


# IGauGa的参数估计=================================================================================================================================
log_f <- function(y, X, alpha, Beta, lambda) {
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
  pdf <- apply(cbind(y, X), 1, f)
  return(sum(log(pdf)))
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

# =========================================================
# 参数打包：alpha, Beta(p×d), lambda -> theta 向量
# =========================================================
pack_theta <- function(alpha, Beta, lambda) {
  c(alpha, as.vector(Beta), lambda)
}

# =========================================================
# 参数解包：theta -> alpha, Beta, lambda
# =========================================================
unpack_theta <- function(theta, p, d) {
  alpha  <- theta[1]
  Beta   <- matrix(theta[2:(1 + p * d)], nrow = p, ncol = d)
  lambda <- theta[length(theta)]
  list(alpha = alpha, Beta = Beta, lambda = lambda)
}

# =========================================================
# 标准参数向量形式的 log-likelihood
# =========================================================
loglik_theta <- function(theta, x, z, p, d) {
  par <- unpack_theta(theta, p, d)
  log_f(
    y      = x,
    X      = z,
    alpha  = par$alpha,
    Beta   = par$Beta,
    lambda = par$lambda
  )
}

# =========================================================
# 每个观测的 score contribution
# 输出: n × dim(theta) 矩阵
# =========================================================
score_contributions_theta <- function(theta, x, z, p, d) {
  n <- nrow(x)
  dim_theta <- length(theta)
  
  score_mat <- matrix(0, nrow = n, ncol = dim_theta)
  
  for (i in 1:n) {
    grad_i <- numDeriv::grad(
      func = function(th) {
        par <- unpack_theta(th, p, d)
        
        # 单个观测 log-density
        xi <- matrix(x[i, ], nrow = 1)
        zi <- matrix(z[i, ], nrow = 1)
        
        log_f(
          y      = xi,
          X      = zi,
          alpha  = par$alpha,
          Beta   = par$Beta,
          lambda = par$lambda
        )
      },
      theta,
      method = "Richardson",
      method.args = list(eps = 1e-10, d = 1e-4, r = 6)
    )
    
    score_mat[i, ] <- grad_i
  }
  
  colnames(score_mat) <- c(
    "alpha",
    paste0("beta_", 1:(p * d)),
    "lambda"
  )
  
  return(score_mat)
}

# =========================================================
# Hessian matrix of log-likelihood
# =========================================================
hessian_theta <- function(theta, x, z, p, d) {
  numDeriv::hessian(
    func = function(th) loglik_theta(th, x, z, p, d),
    theta,
    method = "Richardson",
    method.args = list(eps = 1e-10, d = 1e-4, r = 6)
  )
}

# =========================================================
# Sandwich covariance matrix
# =========================================================
sandwich_covariance <- function(theta_hat, x, z, p, d) {
  # score contributions
  S <- score_contributions_theta(theta_hat, x, z, p, d)
  
  # B matrix
  B <- crossprod(S)
  
  # Hessian
  H <- -hessian_theta(theta_hat, x, z, p, d)
  
  # Sandwich
  Sigma <- solve(H) %*% B %*% solve(H)
  
  return(Sigma)
}

# =========================================================
# Asymptotic confidence intervals
# =========================================================
asymptotic_CI <- function(theta_hat, Sigma, level = 0.95) {
  z_alpha <- qnorm(1 - (1 - level) / 2)
  se <- sqrt(diag(Sigma))
  
  lower <- theta_hat - z_alpha * se
  upper <- theta_hat + z_alpha * se
  
  list(
    lower = lower,
    upper = upper,
    se    = se
  )
}

delta_CI_Vartau <- function(lambda_hat, var_lambda, level = 0.95) {
  z_alpha <- qnorm(1 - (1 - level) / 2)
  
  g_hat <- 1 / lambda_hat
  
  se_g <- sqrt(var_lambda) / lambda_hat^2
  
  lower <- g_hat - z_alpha * se_g
  upper <- g_hat + z_alpha * se_g
  
  c(lower = lower, upper = upper)
}

# Import the data
data <- readxl::read_xlsx("数据.xlsx")
data <- data[, 1:6]
data <- data %>% filter(TOTSLF23 != 0)
data <- data %>% filter(TOTMCR23 != 0)
data <- data %>% filter(TOTPRV23 != 0)
data <- data / 1000

x <- data[, c(3, 4, 6)] %>% as.matrix()
n <- nrow(x)
d <- ncol(x)

alpha <- 1
Var.tau <- 1
z <- rep(1, n) %>% as.matrix()
Beta <- matrix(rep(0, d), ncol = d)

Res <- IGauGa.MLE(alpha, Beta, Var.tau, x, t(z))


theta_hat <- pack_theta(Res$Alpha, Res$Beta, Res$Lambda)

# MLE
result_hat <- c(Res$Alpha, exp(Res$Beta), 1 / Res$Lambda)
result_hat %>% round(digits = 4)

# 渐近CI
p <- ncol(z)
Sigma_hat <- sandwich_covariance(theta_hat, x, z, p, d)

CI <- asymptotic_CI(theta_hat, Sigma_hat)

# μ的置信区间
mu_hat <- exp(Res$Beta)
se_g <- mu_hat * sqrt(diag(Sigma_hat[2:4, 2:4]))
mu.lower <- mu_hat - qnorm(1 - (1 - 0.95) / 2) * se_g
mu.upper <- mu_hat + qnorm(1 - (1 - 0.95) / 2) * se_g

# Var(τ)的置信区间
idx_lambda <- length(theta_hat)
CI_Vartau <- delta_CI_Vartau(
  lambda_hat = theta_hat[idx_lambda],
  var_lambda = Sigma_hat[idx_lambda, idx_lambda]
)

lower.bound <- CI$lower
upper.bound <- CI$upper
lower.bound[2:4] <- mu.lower
upper.bound[2:4] <- mu.upper
lower.bound[idx_lambda] <- CI_Vartau[1]
upper.bound[idx_lambda] <- CI_Vartau[2]

lower.bound %>% round(digits = 4)
upper.bound %>% round(digits = 4)

# Bootstrap CI
mean.std.CI <- function(lasample) {
  G <- dim(lasample)[1]
  lasort <- apply(lasample, 2, sort)
  indexx <- floor(c(0.025 * G, 0.975 * G))
  laL <- (lasort[indexx[1], ] + lasort[indexx[1] + 1, ]) / 2
  laU <- (lasort[indexx[2], ] + lasort[indexx[2] + 1, ]) / 2
  results <- c(laL, laU)
  return(results)
}

Boot.IGauGa <- function(n, alpha, Beta, Var.tau) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(SuppDists)))
  suppressWarnings(suppressMessages(library(Bessel)))
  suppressWarnings(suppressMessages(library(GIGrvg)))
  suppressWarnings(suppressMessages(library(numDeriv)))
  suppressWarnings(suppressMessages(library(MASS)))
  rIGauGa <- function(n, z, alpha, Beta, Var.tau) {
    # 确定τ分布的参数λ并生成τ分布的随机数
    lambda <- 1 / Var.tau
    tau <- rinvGauss(n, 1, lambda)
    
    # 计算Gamma分布的速率参数
    d <- ncol(Beta)
    rate.para <- alpha / apply(exp(z %*% Beta), 2, function(x) x * tau)
    
    # 生成样本x
    x <- matrix(0, n, d)
    for (i in 1:d) {
      x[, i] <- sapply(rate.para[, i], function(s) rgamma(1, alpha, s))
    }
    return(x)
  }
  
  g <- 1
  G <- 6
  
  boot.estimation <- list(Alpha.boot = numeric(length = G), 
                          Mu.boot = matrix(nrow = G, ncol = 3), 
                          Var.tau.boot = numeric(length = G))
  
  while (g <= G) {
    boot.sample <- rIGauGa(n, z, alpha, Beta, Var.tau)
    Res <- tryCatch({
      IGauGa.MLE(alpha, Beta, Var.tau, boot.sample, t(z))
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    if (!all(is.na(Res))) {
      boot.estimation[[1]][g] <- Res$Alpha
      boot.estimation[[2]][g, ] <- exp(Res$Beta)
      boot.estimation[[3]][g] <- 1 / Res$Lambda
      g <- g + 1
    }
  }
  
  return(boot.estimation)
}

num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Boot.IGauGa(n, Res$Alpha, Res$Beta, 1 / Res$Lambda)
stopCluster(cl)

Result.alpha <- result_part[[1]][[1]]
Result.mu <- result_part[[1]][[2]]
Result.var.tau <- result_part[[1]][[3]]

for (i in 2:num_cores) {
  Result.alpha <- c(Result.alpha, result_part[[i]][[1]])
  Result.mu <- rbind(Result.mu, result_part[[i]][[2]])
  Result.var.tau <- c(Result.var.tau, result_part[[i]][[3]])
}

mean.std.CI(Result.alpha %>% as.matrix()) %>% round(digits = 4)
mean.std.CI(Result.mu %>% as.matrix()) %>% round(digits = 4)
mean.std.CI(Result.var.tau %>% as.matrix()) %>% round(digits = 4)

# RIGGa的参数估计============================================================================================
log_f <- function(y, X, alpha, Beta, lambda) {
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
  pdf <- apply(cbind(y, X), 1, f)
  return(sum(log(pdf)))
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

# =========================================================
# 参数打包：alpha, Beta(p×d), lambda -> theta 向量
# =========================================================
pack_theta <- function(alpha, Beta, lambda) {
  c(alpha, as.vector(Beta), lambda)
}

# =========================================================
# 参数解包：theta -> alpha, Beta, lambda
# =========================================================
unpack_theta <- function(theta, p, d) {
  alpha  <- theta[1]
  Beta   <- matrix(theta[2:(1 + p * d)], nrow = p, ncol = d)
  lambda <- theta[length(theta)]
  list(alpha = alpha, Beta = Beta, lambda = lambda)
}

# =========================================================
# 标准参数向量形式的 log-likelihood
# =========================================================
loglik_theta <- function(theta, x, z, p, d) {
  par <- unpack_theta(theta, p, d)
  log_f(
    y      = x,
    X      = z,
    alpha  = par$alpha,
    Beta   = par$Beta,
    lambda = par$lambda
  )
}

# =========================================================
# 每个观测的 score contribution
# 输出: n × dim(theta) 矩阵
# =========================================================
score_contributions_theta <- function(theta, x, z, p, d) {
  n <- nrow(x)
  dim_theta <- length(theta)
  
  score_mat <- matrix(0, nrow = n, ncol = dim_theta)
  
  for (i in 1:n) {
    grad_i <- numDeriv::grad(
      func = function(th) {
        par <- unpack_theta(th, p, d)
        
        # 单个观测 log-density
        xi <- matrix(x[i, ], nrow = 1)
        zi <- matrix(z[i, ], nrow = 1)
        
        log_f(
          y      = xi,
          X      = zi,
          alpha  = par$alpha,
          Beta   = par$Beta,
          lambda = par$lambda
        )
      },
      theta,
      method = "Richardson",
      method.args = list(eps = 1e-10, d = 1e-4, r = 6)
    )
    
    score_mat[i, ] <- grad_i
  }
  
  colnames(score_mat) <- c(
    "alpha",
    paste0("beta_", 1:(p * d)),
    "lambda"
  )
  
  return(score_mat)
}

# =========================================================
# Hessian matrix of log-likelihood
# =========================================================
hessian_theta <- function(theta, x, z, p, d) {
  numDeriv::hessian(
    func = function(th) loglik_theta(th, x, z, p, d),
    theta,
    method = "Richardson",
    method.args = list(eps = 1e-10, d = 1e-4, r = 6)
  )
}

# =========================================================
# Sandwich covariance matrix
# =========================================================
sandwich_covariance <- function(theta_hat, x, z, p, d) {
  # score contributions
  S <- score_contributions_theta(theta_hat, x, z, p, d)
  
  # B matrix
  B <- crossprod(S)
  
  # Hessian
  H <- -hessian_theta(theta_hat, x, z, p, d)
  
  # Sandwich
  Sigma <- solve(H) %*% B %*% solve(H)
  
  return(Sigma)
}

# =========================================================
# Asymptotic confidence intervals
# =========================================================
asymptotic_CI <- function(theta_hat, Sigma, level = 0.95) {
  z_alpha <- qnorm(1 - (1 - level) / 2)
  se <- sqrt(diag(Sigma))
  
  lower <- theta_hat - z_alpha * se
  upper <- theta_hat + z_alpha * se
  
  list(
    lower = lower,
    upper = upper,
    se    = se
  )
}

delta_CI_Vartau <- function(lambda_hat, var_lambda, level = 0.95) {
  z_alpha <- qnorm(1 - (1 - level) / 2)
  
  g_hat <- (2 * lambda_hat - 1) * (lambda_hat - 1) / lambda_hat^2
  
  se_g <- sqrt(var_lambda) * abs((3 * lambda_hat - 2) / lambda_hat^3)
  
  lower <- g_hat - z_alpha * se_g
  upper <- g_hat + z_alpha * se_g
  
  c(lower = lower, upper = upper)
}

# Import the data
data <- readxl::read_xlsx("数据.xlsx")
data <- data[, 1:6]
data <- data %>% filter(TOTSLF23 != 0)
data <- data %>% filter(TOTMCR23 != 0)
data <- data %>% filter(TOTPRV23 != 0)
data <- data / 1000

x <- data[, c(3, 4, 6)] %>% as.matrix()
n <- nrow(x)
d <- ncol(x)

alpha <- 1
Var.tau <- 1
z <- rep(1, n) %>% as.matrix()
Beta <- matrix(rep(0, d), ncol = d)

Res <- RIGGa.MLE(alpha, Beta, Var.tau, x, t(z))


theta_hat <- pack_theta(Res$Alpha, Res$Beta, Res$Lambda)

# MLE
result_hat <- c(Res$Alpha, exp(Res$Beta), (2 * Res$Lambda - 1) * (Res$Lambda - 1) / Res$Lambda^2)
result_hat %>% round(digits = 4)

# 渐近CI
p <- ncol(z)
Sigma_hat <- sandwich_covariance(theta_hat, x, z, p, d)

CI <- asymptotic_CI(theta_hat, Sigma_hat)

# μ的置信区间
mu_hat <- exp(Res$Beta)
se_g <- mu_hat * sqrt(diag(Sigma_hat[2:4, 2:4]))
mu.lower <- mu_hat - qnorm(1 - (1 - 0.95) / 2) * se_g
mu.upper <- mu_hat + qnorm(1 - (1 - 0.95) / 2) * se_g

# Var(τ)的置信区间
idx_lambda <- length(theta_hat)
CI_Vartau <- delta_CI_Vartau(
  lambda_hat = theta_hat[idx_lambda],
  var_lambda = Sigma_hat[idx_lambda, idx_lambda]
)

lower.bound <- CI$lower
upper.bound <- CI$upper
lower.bound[2:4] <- mu.lower
upper.bound[2:4] <- mu.upper
lower.bound[idx_lambda] <- CI_Vartau[1]
upper.bound[idx_lambda] <- CI_Vartau[2]

lower.bound %>% round(digits = 4)
upper.bound %>% round(digits = 4)

# Bootstrap CI
mean.std.CI <- function(lasample) {
  G <- dim(lasample)[1]
  lasort <- apply(lasample, 2, sort)
  indexx <- floor(c(0.025 * G, 0.975 * G))
  laL <- (lasort[indexx[1], ] + lasort[indexx[1] + 1, ]) / 2
  laU <- (lasort[indexx[2], ] + lasort[indexx[2] + 1, ]) / 2
  results <- c(laL, laU)
  return(results)
}

Boot.RIGGa <- function(n, alpha, Beta, Var.tau) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(SuppDists)))
  suppressWarnings(suppressMessages(library(Bessel)))
  suppressWarnings(suppressMessages(library(GIGrvg)))
  suppressWarnings(suppressMessages(library(numDeriv)))
  suppressWarnings(suppressMessages(library(MASS)))
  rRIGGa <- function(n, z, alpha, Beta, Var.tau) {
    # 确定τ分布的参数λ并生成τ分布的随机数
    lambda <- (3 + sqrt(1 + 4 * Var.tau)) / (4 - 2 * Var.tau)
    tau <- 1 / rinvGauss(n, lambda, lambda / (lambda - 1))
    
    # 计算Gamma分布的速率参数
    d <- ncol(Beta)
    rate.para <- alpha / apply(exp(z %*% Beta), 2, function(x) x * tau)
    
    # 生成样本x
    x <- matrix(0, n, d)
    for (i in 1:d) {
      x[, i] <- sapply(rate.para[, i], function(s) rgamma(1, alpha, s))
    }
    return(x)
  }
  
  g <- 1
  G <- 6
  
  boot.estimation <- list(Alpha.boot = numeric(length = G), 
                          Mu.boot = matrix(nrow = G, ncol = 3), 
                          Var.tau.boot = numeric(length = G))
  
  while (g <= G) {
    boot.sample <- rRIGGa(n, z, alpha, Beta, Var.tau)
    Res <- tryCatch({
      RIGGa.MLE(alpha, Beta, Var.tau, boot.sample, t(z))
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    if (!all(is.na(Res))) {
      boot.estimation[[1]][g] <- Res$Alpha
      boot.estimation[[2]][g, ] <- exp(Res$Beta)
      boot.estimation[[3]][g] <- (2 * Res$Lambda - 1) * (Res$Lambda - 1) / Res$Lambda^2
      g <- g + 1
    }
  }
  
  return(boot.estimation)
}

num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Boot.RIGGa(n, Res$Alpha, Res$Beta, (2 * Res$Lambda - 1) * (Res$Lambda - 1) / Res$Lambda^2)
stopCluster(cl)

Result.alpha <- result_part[[1]][[1]]
Result.mu <- result_part[[1]][[2]]
Result.var.tau <- result_part[[1]][[3]]

for (i in 2:num_cores) {
  Result.alpha <- c(Result.alpha, result_part[[i]][[1]])
  Result.mu <- rbind(Result.mu, result_part[[i]][[2]])
  Result.var.tau <- c(Result.var.tau, result_part[[i]][[3]])
}

mean.std.CI(Result.alpha %>% as.matrix()) %>% round(digits = 4)
mean.std.CI(Result.mu %>% as.matrix()) %>% round(digits = 4)
mean.std.CI(Result.var.tau %>% as.matrix()) %>% round(digits = 4)

# LNGa的参数估计===========================================================================================

log_f <- function(y, X, alpha, Beta, lambda) {
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
  pdf <- apply(cbind(y, X), 1, f)
  return(sum(log(pdf)))
}

LNGa.MLE <- function(alpha, Beta, Var.tau, y, X, Type = "LN", ep = 1e-8) {
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

# =========================================================
# 参数打包：alpha, Beta(p×d), lambda -> theta 向量
# =========================================================
pack_theta <- function(alpha, Beta, lambda) {
  c(alpha, as.vector(Beta), lambda)
}

# =========================================================
# 参数解包：theta -> alpha, Beta, lambda
# =========================================================
unpack_theta <- function(theta, p, d) {
  alpha  <- theta[1]
  Beta   <- matrix(theta[2:(1 + p * d)], nrow = p, ncol = d)
  lambda <- theta[length(theta)]
  list(alpha = alpha, Beta = Beta, lambda = lambda)
}

# =========================================================
# 标准参数向量形式的 log-likelihood
# =========================================================
loglik_theta <- function(theta, x, z, p, d) {
  par <- unpack_theta(theta, p, d)
  log_f(
    y      = x,
    X      = z,
    alpha  = par$alpha,
    Beta   = par$Beta,
    lambda = par$lambda
  )
}

# =========================================================
# 每个观测的 score contribution
# 输出: n × dim(theta) 矩阵
# =========================================================
score_contributions_theta <- function(theta, x, z, p, d) {
  n <- nrow(x)
  dim_theta <- length(theta)
  
  score_mat <- matrix(0, nrow = n, ncol = dim_theta)
  
  for (i in 1:n) {
    grad_i <- numDeriv::grad(
      func = function(th) {
        par <- unpack_theta(th, p, d)
        
        # 单个观测 log-density
        xi <- matrix(x[i, ], nrow = 1)
        zi <- matrix(z[i, ], nrow = 1)
        
        log_f(
          y      = xi,
          X      = zi,
          alpha  = par$alpha,
          Beta   = par$Beta,
          lambda = par$lambda
        )
      },
      theta,
      method = "Richardson",
      method.args = list(eps = 1e-10, d = 1e-4, r = 6)
    )
    
    score_mat[i, ] <- grad_i
  }
  
  colnames(score_mat) <- c(
    "alpha",
    paste0("beta_", 1:(p * d)),
    "lambda"
  )
  
  return(score_mat)
}

# =========================================================
# Hessian matrix of log-likelihood
# =========================================================
hessian_theta <- function(theta, x, z, p, d) {
  numDeriv::hessian(
    func = function(th) loglik_theta(th, x, z, p, d),
    theta,
    method = "Richardson",
    method.args = list(eps = 1e-10, d = 1e-4, r = 6)
  )
}

# =========================================================
# Sandwich covariance matrix
# =========================================================
sandwich_covariance <- function(theta_hat, x, z, p, d) {
  # score contributions
  S <- score_contributions_theta(theta_hat, x, z, p, d)
  
  # B matrix
  B <- crossprod(S)
  
  # Hessian
  H <- -hessian_theta(theta_hat, x, z, p, d)
  
  # Sandwich
  Sigma <- solve(H) %*% B %*% solve(H)
  
  return(Sigma)
}

# =========================================================
# Asymptotic confidence intervals
# =========================================================
asymptotic_CI <- function(theta_hat, Sigma, level = 0.95) {
  z_alpha <- qnorm(1 - (1 - level) / 2)
  se <- sqrt(diag(Sigma))
  
  lower <- theta_hat - z_alpha * se
  upper <- theta_hat + z_alpha * se
  
  list(
    lower = lower,
    upper = upper,
    se    = se
  )
}

delta_CI_Vartau <- function(lambda_hat, var_lambda, level = 0.95) {
  z_alpha <- qnorm(1 - (1 - level) / 2)
  
  g_hat <- exp(lambda_hat) - 1
  
  se_g <- sqrt(var_lambda) * exp(lambda_hat)
  
  lower <- g_hat - z_alpha * se_g
  upper <- g_hat + z_alpha * se_g
  
  c(lower = lower, upper = upper)
}

# Import the data
data <- readxl::read_xlsx("MEPS.xlsx")
data <- data[, 1:6]
data <- data %>% filter(TOTSLF23 != 0)
data <- data %>% filter(TOTMCR23 != 0)
data <- data %>% filter(TOTPRV23 != 0)
data <- data / 1000

x <- data[, c(3, 4, 6)] %>% as.matrix()
n <- nrow(x)
d <- ncol(x)

alpha <- 1
Var.tau <- 1
z <- rep(1, n) %>% as.matrix()
Beta <- matrix(rep(0, d), ncol = d)

Res <- LNGa.MLE(alpha, Beta, Var.tau, x, t(z))


theta_hat <- pack_theta(Res$Alpha, Res$Beta, Res$Lambda)

# MLE
result_hat <- c(Res$Alpha, exp(Res$Beta), exp(Res$Lambda) - 1)
result_hat %>% round(digits = 4)

# 渐近CI
p <- ncol(z)
Sigma_hat <- sandwich_covariance(theta_hat, x, z, p, d)

CI <- asymptotic_CI(theta_hat, Sigma_hat)

# μ的置信区间
mu_hat <- exp(Res$Beta)
se_g <- mu_hat * sqrt(diag(Sigma_hat[2:4, 2:4]))
mu.lower <- mu_hat - qnorm(1 - (1 - 0.95) / 2) * se_g
mu.upper <- mu_hat + qnorm(1 - (1 - 0.95) / 2) * se_g

# Var(τ)的置信区间
idx_lambda <- length(theta_hat)
CI_Vartau <- delta_CI_Vartau(
  lambda_hat = theta_hat[idx_lambda],
  var_lambda = Sigma_hat[idx_lambda, idx_lambda]
)

lower.bound <- CI$lower
upper.bound <- CI$upper
lower.bound[2:4] <- mu.lower
upper.bound[2:4] <- mu.upper
lower.bound[idx_lambda] <- CI_Vartau[1]
upper.bound[idx_lambda] <- CI_Vartau[2]

lower.bound %>% round(digits = 4)
upper.bound %>% round(digits = 4)

# Bootstrap CI
mean.std.CI <- function(lasample) {
  G <- dim(lasample)[1]
  lasort <- apply(lasample, 2, sort)
  indexx <- floor(c(0.025 * G, 0.975 * G))
  laL <- (lasort[indexx[1], ] + lasort[indexx[1] + 1, ]) / 2
  laU <- (lasort[indexx[2], ] + lasort[indexx[2] + 1, ]) / 2
  results <- c(laL, laU)
  return(results)
}

Boot.LNGa <- function(n, alpha, Beta, Var.tau) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(SuppDists)))
  suppressWarnings(suppressMessages(library(Bessel)))
  suppressWarnings(suppressMessages(library(GIGrvg)))
  suppressWarnings(suppressMessages(library(numDeriv)))
  suppressWarnings(suppressMessages(library(MASS)))
  rLNGa <- function(n, z, alpha, Beta, Var.tau) {
    # 确定τ分布的参数λ并生成τ分布的随机数
    lambda <- log(Var.tau + 1)
    tau <- rlnorm(n, -lambda / 2, lambda)
    
    # 计算Gamma分布的速率参数
    d <- ncol(Beta)
    rate.para <- alpha / apply(exp(z %*% Beta), 2, function(x) x * tau)
    
    # 生成样本x
    x <- matrix(0, n, d)
    for (i in 1:d) {
      x[, i] <- sapply(rate.para[, i], function(s) rgamma(1, alpha, s))
    }
    return(x)
  }
  
  g <- 1
  G <- 6
  
  boot.estimation <- list(Alpha.boot = numeric(length = G), 
                          Mu.boot = matrix(nrow = G, ncol = 3), 
                          Var.tau.boot = numeric(length = G))
  
  while (g <= G) {
    boot.sample <- rLNGa(n, z, alpha, Beta, Var.tau)
    Res <- tryCatch({
      LNGa.MLE(alpha, Beta, Var.tau, boot.sample, t(z))
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    if (!all(is.na(Res))) {
      boot.estimation[[1]][g] <- Res$Alpha
      boot.estimation[[2]][g, ] <- exp(Res$Beta)
      boot.estimation[[3]][g] <- exp(Res$Lambda) - 1
      g <- g + 1
    }
  }
  
  return(boot.estimation)
}

num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Boot.LNGa(n, Res$Alpha, Res$Beta, exp(Res$Lambda) - 1)
stopCluster(cl)

Result.alpha <- result_part[[1]][[1]]
Result.mu <- result_part[[1]][[2]]
Result.var.tau <- result_part[[1]][[3]]

for (i in 2:num_cores) {
  Result.alpha <- c(Result.alpha, result_part[[i]][[1]])
  Result.mu <- rbind(Result.mu, result_part[[i]][[2]])
  Result.var.tau <- c(Result.var.tau, result_part[[i]][[3]])
}

mean.std.CI(Result.alpha %>% as.matrix()) %>% round(digits = 4)
mean.std.CI(Result.mu %>% as.matrix()) %>% round(digits = 4)
mean.std.CI(Result.var.tau %>% as.matrix()) %>% round(digits = 4)
