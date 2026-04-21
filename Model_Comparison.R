# 开始计时
start_time <- Sys.time()

# 导入需要的Package
library(doParallel)
library(foreach)
library(tidyverse)
library(abind)

# ==============================================================================

# 定义Simualtion需要的函数
Simulation <- function(number, size, alpha, Beta, Var.tau, Ratio) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(SuppDists)))
  suppressWarnings(suppressMessages(library(Bessel)))
  suppressWarnings(suppressMessages(library(GIGrvg)))
  
  GIG.lambda <- function(Var.tau) {
    f <- function(lambda) {
      besselK(lambda, 2) * besselK(lambda, 0) / (besselK(lambda, 1))^2 - 1 - Var.tau
    }
    return(uniroot(f, c(1, 5), extendInt = "yes")$root)
  }
  a <- function(lambda) lambda * besselK(lambda, 1) / besselK(lambda, 0)
  b <- function(lambda) lambda * besselK(lambda, 0) / besselK(lambda, 1)
  
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
  rIGauGa <- function(n, X, alpha, Beta, Var.tau, Type) {
    lambda <- 1 / Var.tau
    tau <- rinvGauss(n, 1, lambda)
    
    d <- ncol(Beta)
    beta <- alpha / apply(exp(t(X) %*% Beta), 2, function(x) x * tau)
    y <- matrix(0, n, d)
    for (i in 1:d) {
      y[, i] <- apply(beta, 1, function(x) rgamma(1, alpha, x[i]))
    }
    return(y)
  }
  rRIGGa <- function(n, X, alpha, Beta, Var.tau, Type) {
    lambda <- (3 + sqrt(1 + 4 * Var.tau)) / (4 - 2 * Var.tau)
    tau <- 1 / rinvGauss(n, lambda, lambda / (lambda - 1))
    d <- ncol(Beta)
    beta <- alpha / apply(exp(t(X) %*% Beta), 2, function(x) x * tau)
    y <- matrix(0, n, d)
    for (i in 1:d) {
      y[, i] <- apply(beta, 1, function(x) rgamma(1, alpha, x[i]))
    }
    return(y)
  }
  rLNGa <- function(n, X, alpha, Beta, Var.tau, Type = "LN") {
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
    
    # 确定τ分布的参数λ
    lambda <- GIG.lambda(Var.tau)
    
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
  
  # NEM算法估计MLE================================================================
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
  
  # 开始Simulation
  k <- 1
  Loglikelihood <- matrix(nrow = number, ncol = 5)
  
  while(k <= number) {
    X <- rbind(rep(1, n))
    y1 <- rIGaGa(n, t(X), alpha, Beta, Var.tau)
    y2 <- rGIGGa(n, t(X), alpha, Beta, Var.tau)
    y4 <- rIGauGa(n, X, alpha, Beta, Var.tau)
    y5 <- rRIGGa(n, X, alpha, Beta, Var.tau)
    y3 <- rLNGa(n, X, alpha, Beta, Var.tau)
    z <- rmultinom(n, 1, Ratio)
    y <- y1 * z[1, ] + y2 * z[2, ] + y3 * z[3, ] + y4 * z[4, ] + y5 * z[5, ]
    para.num <- ncol(y) + 2
    
    Res1 <- tryCatch({
      IGaGa.MLE(alpha, Beta, Var.tau, y, t(X))
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    Res2 <- tryCatch({
      GIGGa.MLE(alpha, Beta, Var.tau, y, t(X))
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    Res4 <- tryCatch({
      IGauGa.MLE(alpha, Beta, Var.tau, y, X)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    Res5 <- tryCatch({
      RIGGa.MLE(alpha, Beta, Var.tau, y, X)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    Res3 <- tryCatch({
      LNGa.MLE(alpha, Beta, Var.tau, y, X)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    
    if (!all(is.na(Res1)) & !all(is.na(Res2)) & !all(is.na(Res3)) & !all(is.na(Res4)) & !all(is.na(Res5))) {
      if (Res1$index == 1 & Res3$index == 1 & Res4$index == 1 & Res5$index == 1) {
        Loglikelihood[k, 1] <- -2 * Res1$Loglikelihood + 2 * para.num
        Loglikelihood[k, 2] <- -2 * Res2$Loglikelihood + 2 * para.num
        Loglikelihood[k, 3] <- -2 * Res3$Loglikelihood + 2 * para.num
        Loglikelihood[k, 4] <- -2 * Res4$Loglikelihood + 2 * para.num
        Loglikelihood[k, 5] <- -2 * Res5$Loglikelihood + 2 * para.num
        k <- k + 1
      }
    }
  }
  
  return(Loglikelihood)
}

# ==============================================================================================================================

# 设定参数
num_cores <- detectCores()
all_times <- num_cores * 11
n <- 500
alpha <- 2
Beta <- matrix(c(log(2), log(1.5), log(0.8)), nrow = 1)
Var.tau <- 0.5
Ratio <- c(0, 0, 0.9, 0.1, 0)

# ==============================================================================================================================

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Simulation(number = ceiling(all_times / num_cores), n, alpha, Beta, Var.tau, Ratio)
stopCluster(cl)

# Simulation(number = 1, n, alpha, Beta, Var.tau, Ratio)
# ==============================================================================================================================

# 汇总拟合结果并保存
Result <- result_part[[1]]
for (i in 2:num_cores) {
  Result <- rbind(Result, result_part[[i]])
}
Rank <- apply(-Result, 1, rank) %>% t()

colMeans(Result, na.rm = TRUE) %>% round(digits = 1)
colMeans(Rank) %>% round(digits = 2)
# ==============================================================================================================================

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)
