
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
  suppressWarnings(suppressMessages(library(numDeriv)))
  suppressWarnings(suppressMessages(library(MASS)))
  
  # 生成IGaGa随机数===============================================================
  rGeGa <- function(n, X, alpha, Beta, Var.tau, Type = "LN") {
    if (Type == "Ga") {
      lambda <- 1 / Var.tau
      tau <- rgamma(n, lambda, lambda)
    } else if (Type == "IGa") {
      lambda <- 1 / Var.tau + 2
      tau <- 1 / rgamma(n, lambda, lambda - 1)
    } else if (Type == "IGau") {
      lambda <- 1 / Var.tau
      tau <- rinvGauss(n, 1, lambda)
    } else if (Type == "RIG") {
      lambda <- (3 + sqrt(1 + 4 * Var.tau)) / (4 - 2 * Var.tau)
      tau <- 1 / rinvGauss(n, lambda, lambda / (lambda - 1))
    } else if (Type == "LN") {
      lambda <- log(Var.tau + 1)
      tau <- rlnorm(n, -lambda / 2, lambda)
    }
    d <- ncol(Beta)
    beta <- alpha / apply(exp(t(X) %*% Beta), 2, function(x) x * tau)
    y <- matrix(0, n, d)
    for (i in 1:d) {
      y[, i] <- apply(beta, 1, function(x) rgamma(1, alpha, x[i]))
    }
    return(y)
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
  
  # 定义对数似然函数
  NEM <- function(alpha, Beta, Var.tau, y, X, Type = "LN", ep = 1e-5) {
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
    log_f5(
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
          zi <- matrix(z[, i], ncol = 1)
          
          log_f5(
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
  
  # =========================================================
  # Delta method CI for Var(tau)
  # Var(tau) = 1 / (lambda - 2)
  # =========================================================
  delta_CI_Vartau <- function(lambda_hat, var_lambda, level = 0.95) {
    z_alpha <- qnorm(1 - (1 - level) / 2)
    
    g_hat <- exp(lambda_hat) - 1
    g_prime <- exp(lambda_hat)
    
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
    x <- rGeGa(n, t(z), alpha, Beta, Var.tau)
    
    Res <- tryCatch({
      NEM(alpha, Beta, Var.tau, x, t(z))
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    
    
    if (!all(is.na(Res))) {
      if (Res$index == 1) {
        result[[1]][k, ] <- Res$Alpha
        result[[2]][, , k] <- Res$Beta
        result[[3]][k] <- exp(Res$Lambda) - 1
        result[[4]][k] <- Res$Loglikelihood
        
        theta_true <- pack_theta(alpha, Beta, Var.tau)
        theta_hat <- pack_theta(Res$Alpha, Res$Beta, Res$Lambda)
        
        p <- ncol(z)
        d <- ncol(x)
        
        Sigma_hat <- sandwich_covariance(theta_hat, x, t(z), p, d)
        
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
n <- 300
alpha <- 3
Beta <- matrix(c(1, 3, -2, -1, 4, 0), ncol = 2)
# Beta <- matrix(c(log(2), log(3)), ncol = 2)
Var.tau <- exp(1) - 1
# Simulation(1, n, alpha, Beta, Var.tau)
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
