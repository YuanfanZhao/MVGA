
# 开始计时
start_time <- Sys.time()

# 导入需要的Package
library(doParallel)
library(foreach)
library(tidyverse)
library(abind)

# ==============================================================================================================================

# 定义Simualtion需要的函数
Simulation <- function(number, n, alpha, Beta, Var.tau) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(SuppDists)))
  suppressWarnings(suppressMessages(library(Bessel)))
  suppressWarnings(suppressMessages(library(numDeriv)))
  suppressWarnings(suppressMessages(library(MASS)))
  
  # 生成IGaGa随机数===============================================================
  rIGaGa <- function(n, z, alpha, Beta, Var.tau) {
    # 确定τ分布的参数λ并生成τ分布的随机数
    lambda <- 1 / Var.tau + 2
    tau <- 1 / rgamma(n, lambda, lambda - 1)
    
    # 计算Gamma分布的速率参数
    d <- ncol(Beta)
    rate.para <- alpha / apply(exp(z %*% Beta), 2, function(x) x * tau)
    
    # 生成样本x
    x <- matrix(0, n, d)
    for (i in 1:d) {
      x[, i] <- sapply(rate.para[, i], function(b) rgamma(1, alpha, b))
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
      term1 <- gamma(alpha.plus + lambda) * (lambda - 1)^lambda / gamma(lambda)
      term2 <- prod(alpha^alpha * x^(alpha - 1) / gamma(alpha) / mu^alpha)
      term3 <- (sum(alpha * x / mu) + lambda - 1)^(-alpha.plus - lambda)
      term1 * term2 * term3
    }
    pdf <- apply(cbind(x, z), 1, f)
    return(sum(log(pdf)))
  }
  
  # NEM算法估计MLE================================================================
  IGaGa.MLE <- function(alpha, Beta, Var.tau, x, z, ep = 1e-8) {
    # 收敛指标与迭代次数
    index <- 0
    k <- 1
    
    # 样本量n
    n <- nrow(x)
    d <- ncol(x)
    
    # 确定τ分布的参数λ
    lambda <- 1 / Var.tau + 2
    
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
  
  # =========================================================
  # Delta method CI for Var(tau)
  # Var(tau) = 1 / (lambda - 2)
  # =========================================================
  delta_CI_Vartau <- function(lambda_hat, var_lambda, level = 0.95) {
    z_alpha <- qnorm(1 - (1 - level) / 2)
    
    g_hat <- 1 / (lambda_hat - 2)
    g_prime <- -1 / (lambda_hat - 2)^2
    
    se_g <- abs(g_prime) * sqrt(var_lambda)
    
    lower <- g_hat - z_alpha * se_g
    upper <- g_hat + z_alpha * se_g
    
    c(lower = lower, upper = upper)
  }
  
  
  
  # 开始Simulation
  k <- 1
  
  result <- list(Alpha = matrix(nrow = number, ncol = length(alpha)), 
                 Beta = array(dim = c(nrow(Beta), ncol(Beta), number)),
                 Var.tau = numeric(length = number), 
                 Loglikelihood = numeric(length = number),
                 CR = matrix(nrow = number, ncol = 8),
                 Width = matrix(nrow = number, ncol = 8))
  
  while(k <= number) {
    z1 <- rep(1, n)
    z2 <- rnorm(n, 0, 1)
    z3 <- rpois(n, 2)
    z <- cbind(z1, z2, z3)
    # z <- rep(1, n)
    x <- rIGaGa(n, z, alpha, Beta, Var.tau)
    
    Res <- tryCatch({
      IGaGa.MLE(alpha, Beta, Var.tau, x, z)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    
    
    if (!all(is.na(Res))) {
      if (Res$index == 1) {
        result[[1]][k, ] <- Res$Alpha
        result[[2]][, , k] <- Res$Beta
        result[[3]][k] <- 1 / (Res$Lambda - 2)
        result[[4]][k] <- Res$Loglikelihood
        
        theta_true <- pack_theta(alpha, Beta, Var.tau)
        theta_hat <- pack_theta(Res$Alpha, Res$Beta, Res$Lambda)
        
        p <- ncol(z)
        d <- ncol(x)
        
        Sigma_hat <- sandwich_covariance(theta_hat, x, z, p, d)
        
        CI <- asymptotic_CI(theta_hat, Sigma_hat)
        
        # Var(tau) 的区间
        idx_lambda <- length(theta_hat)
        CI_Vartau <- delta_CI_Vartau(
          lambda_hat = theta_hat[idx_lambda],
          var_lambda = Sigma_hat[idx_lambda, idx_lambda]
        )
        
        lower.bound <- CI$lower
        upper.bound <- CI$upper
        lower.bound[idx_lambda] <- CI_Vartau[1]
        upper.bound[idx_lambda] <- CI_Vartau[2]
        
        result[[5]][k, ] <- ifelse(lower.bound <= theta_true & upper.bound >= theta_true, 1, 0)
        result[[6]][k, ] <- upper.bound - lower.bound
        k <- k + 1
      }
    }
  }
  
  return(result)
}

# ==============================================================================================================================

# 设定参数
num_cores <- detectCores()
all_times <- num_cores * 105
n <- 900
alpha <- 2
Beta <- matrix(c(1, 2, -1, 0, 3, -2), ncol = 2)
# Beta <- matrix(c(log(2), log(3)), ncol = 2)
Var.tau <- 0.5

# z <- rep(1, n)
# Beta <- matrix(c(log(2), log(3)), ncol = 2)
# ==============================================================================================================================

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Simulation(number = ceiling(all_times / num_cores), n, alpha, Beta, Var.tau)
stopCluster(cl)

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)

# ==============================================================================================================================

# 汇总拟合结果并保存
Result_Alpha <- result_part[[1]]$Alpha
Result_Beta <- result_part[[1]]$Beta
Result_Var.tau <- result_part[[1]]$Var.tau
Result_Loglikelihood <- result_part[[1]]$Loglikelihood
Result_CR <- result_part[[1]]$CR
Result_Width <- result_part[[1]]$Width

for (i in 2:num_cores) {
  Result_Alpha <- rbind(Result_Alpha, result_part[[i]]$Alpha)
  Result_Beta <- abind(Result_Beta, result_part[[i]]$Beta, along = 3)
  Result_Var.tau <- c(Result_Var.tau, result_part[[i]]$Var.tau)
  Result_Loglikelihood <- c(Result_Loglikelihood, result_part[[i]]$Loglikelihood)
  Result_CR <- rbind(Result_CR, result_part[[i]]$CR)
  Result_Width <- rbind(Result_Width, result_part[[i]]$Width)
}

result <- list(Alpha = Result_Alpha, 
               Beta = Result_Beta, 
               Var.tau = Result_Var.tau, 
               Loglikelihood = Result_Loglikelihood, 
               CR = Result_CR,
               Width = Result_Width)

# ==============================================================================================================================

# 计算A-MLE和A-Moment
MLE.Alpha <- colMeans(result$Alpha)
MLE.Beta <- apply(result$Beta, c(1, 2), mean)
MLE.Var.tau <- mean(result$Var.tau)
MLE.Loglikelihood <- mean(result$Loglikelihood)
CI.CR <- colMeans(Result_CR)
CI.Width <- colMeans(Result_Width)

# ==============================================================================================================================

# 计算MSE
MLE.Alpha.MSE <- colSums((result$Alpha - alpha)^2 / all_times)
MLE.Var.tau.MSE <- sum((result$Var.tau - Var.tau)^2 / all_times)
b <- array(dim = dim(result$Beta))
for (i in 1:all_times) {
  b[, , i] <- result$Beta[, , i] - Beta
}
MLE.Beta.MSE <- apply(b^2, c(1, 2), mean)

MLE.Alpha %>% round(digits = 4)
MLE.Alpha.MSE %>% round(digits = 4)
MLE.Beta %>% round(digits = 4)
MLE.Beta.MSE %>% round(digits = 4)
MLE.Var.tau %>% round(digits = 4)
MLE.Var.tau.MSE %>% round(digits = 4)
CI.CR %>% round(digits = 4)
CI.Width %>% round(digits = 4)
# ==============================================================================================================================

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)
