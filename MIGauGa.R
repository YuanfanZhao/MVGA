# 生成MGeGa随机数===============================================================
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

# NEM算法估计MLE================================================================
IGauGa.MLE <- function(alpha, Beta, Var.tau, y, X, Type, ep = 1e-5) {
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