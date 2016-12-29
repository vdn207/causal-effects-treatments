library(MASS)

# Introducing a seed as suggested by the authors for reproducibility purposes. 
set.seed(18)

# Initializes nu by generating a pi via stick breaking and then drawing nu from a beta(pi, 1 - pi)
init_nu<-function(N, K, alpha){
  # Initial draw of pi's via stick-breaking
  v <- rbeta(n = K, shape1 = alpha, shape2 = 1)
  pi <- rep(0, K)
  pi[1] <- v[1]
  for (k in 2:K){
    pi[k] <- v[k] * pi[k - 1]
  }
  
  return(sapply(1:K, function(k) rbeta(N, pi[k], 1 - pi[k])))
}

# Initializes all the parameters by first drawing nu and an initial phi and then drawing all of the others
# conditional on these two.
initial_draw<-function(Y, X, N, D, K, alpha, a, b, sigmasq.A, sigmasq.n){  
  out<-list()
  
  out$nu <- init_nu(N, K, alpha)
  
  #Initialize phi based on nu
  init.phi <- apply(X, 2, function(x) coef(lm(x ~ -1 + out$nu)))
  init.phi[is.na(init.phi)] <- rnorm(sum(is.na(init.phi)), 0, 1)
  
  update <- update_A(X, N, D, K, sigmasq.A, sigmasq.n, init.phi, out$nu)
  out$phi <- update$phi
  out$big.Phi <- update$big.Phi
  
  update <- update_betatau(Y, K, N, a, b, out$nu)
  out$m <- update$m
  out$S <- update$S
  out$c <- update$c
  out$d <- update$d
  
  out$lambda <- update_pi(alpha, N, K, out$nu)
  
  return(out)
}

# Update phi and Phi
update_A<-function(X, N, D, K, sigmasq.A, sigmasq.n, phi, nu){  
  big.Phi <- lapply(apply(nu, 2, function(nuk) 1/sigmasq.A + sum(nuk)/sigmasq.n)^(-1), 
                    function(x) x*diag(D))
  for (k in 1:K){
    X.resid <- X - (nu%*%phi - nu[,k,drop=FALSE]%*%phi[k,,drop=FALSE])
    phi[k,] <- (1/sigmasq.n * colSums(nu[,k]*X.resid))*(1/sigmasq.A + sum(nu[,k,drop=FALSE])/sigmasq.n)^(-1)
  }
  
  out<-list()
  out$phi <- phi
  out$big.Phi <- big.Phi
  return(out)
}

# Update m, S, c, and d
update_betatau<-function(Y, K, N, a, b, nu){  
  ztz<-t(nu)%*%nu
  diag(ztz) <- apply(nu, 2, sum)
  
  out <- list()
  
  out$S <- ginv(ztz + diag(K))
  out$m <- out$S %*% t(nu) %*% Y
  
  out$c <- a + N/2
  out$d <- b + (t(Y)%*%Y - t(Y) %*% nu %*% out$S %*% t(nu) %*% Y)/2
  return(out)
}

# Update nu
update_Z<-function(nu, lambda, c, d, phi, big.Phi, m, S, sigmasq.n, Y, N, K, X){
  pi.component <- apply(lambda, 1, function(lam) digamma(lam[1]) - digamma(lam[2]))
  
  big.Phi.tr <- sapply(1:K, function(k) sum(diag(big.Phi[[k]])))
  X.component <- -1/(2*sigmasq.n) * t(apply(apply(phi, 1, function(ph) -2*ph%*%t(X) + sum(ph^2)), 1, 
                                            function(x) x + big.Phi.tr))
  Y.component <- as.numeric(-c/(2*d)) * t(apply(-2*Y%*%t(m), 1, function(y) y + d*diag(S)/(c-1) + m^2))
  
  for (k in 1:K){
    X.sumexcept <- -1/(2*sigmasq.n) * as.numeric(2*phi[k,]%*%t(nu%*%phi - nu[,k,drop=FALSE]%*%phi[k,,drop=FALSE]))
    Y.sumexcept <- -c/(2*d) * 2*m[k]*as.numeric(nu%*%m - nu[,k]*m[k])
    v <- pi.component[k] + X.component[,k] + X.sumexcept + Y.component[,k]  + Y.sumexcept
    nu[,k] <- (1 + exp(-v))^(-1)
    
  }  
  return(nu)
}

#Update lambda
update_pi<-function(alpha, N, K, nu){
  lambda <- matrix(rep(0, 2*K), nrow = K, ncol = 2)
  sumnu <- apply(nu, 2, sum)
  lambda[,1] <- alpha/K + sumnu
  lambda[,2] <- 1 + N - sumnu
  
  return (lambda)
}

# Find the absolute value of the biggest parameter change, and subtract the threhsold
convergence.check<-function(new.c, new.d, new.phi, new.big.Phi, new.lambda, new.m, new.S, new.nu, 
                            old.c, old.d, old.phi, old.big.Phi, old.lambda, old.m, old.S, old.nu,
                            threshold){
  return(max(abs(c(new.c - old.c, new.d - old.d, new.phi - old.phi, unlist(new.big.Phi) - unlist(old.big.Phi),
                   new.lambda - old.lambda, new.m - old.m, new.S - old.S, new.nu - old.nu))) - threshold)
  
}

# The main inference procedure.  This is the function used to estimate the sIBP
estimate_parameters<-function(Y, X, K, alpha, a, b, sigmasq.A, sigmasq.n, silent = FALSE){
  # Standardize Y and X
  Y <- (Y - mean(Y))/sd(Y)
  X.means <- apply(X, 2, mean)
  X.sd <- apply(X, 2, sd)
  X <- t(apply(X, 1, function(row) (row-X.means)/X.sd))
  
  N<-length(Y)
  D<-ncol(X)
  initial.params <- initial_draw(Y, X, N, D, K, alpha, a, b, sigmasq.A, sigmasq.n)
  
  # Initialize parameters
  c <- initial.params$c
  d <- initial.params$d
  phi <- initial.params$phi
  big.Phi <- initial.params$big.Phi
  lambda <- initial.params$lambda
  m <- initial.params$m
  S <- initial.params$S
  nu <- initial.params$nu
  
  iterations<-0
  converged<-FALSE
  
  # Update until convergence
  while (!converged){
    iterations<-iterations+1
    
    new.nu <- update_Z(nu, lambda, c, d, phi, big.Phi, m, S, sigmasq.n, Y, N, K, X)
    new.lambda <- update_pi(alpha, N, K, new.nu)
    update <- update_A(X, N, D, K, sigmasq.A, sigmasq.n, phi, new.nu)
    new.phi <- update$phi
    new.big.Phi <- update$big.Phi
    update <- update_betatau(Y, K, N, a, b, new.nu)
    new.c <- update$c
    new.d <- update$d
    new.m <- update$m
    new.S <- update$S
    
    threshold<-10^(-3)
    
    conv.dist <- convergence.check(new.c, new.d, new.phi, new.big.Phi, new.lambda, new.m, new.S, new.nu, 
                                   c, d, phi, big.Phi, lambda, m, S, nu, threshold)
    converged <- conv.dist < 0
    
    c <- new.c
    d <- new.d
    phi <- new.phi
    big.Phi <- new.big.Phi
    lambda <- new.lambda
    m <- new.m
    S <- new.S
    nu <- new.nu
    
    if (iterations %% 10 == 0){
      if (!silent){
      print(iterations)
      print(conv.dist)
      }
    }
  }
  
  return(list("nu"=nu, "m"=m, "S"=S, "lambda"=lambda, "phi"=phi, "big.Phi"=big.Phi, "c"=c, "d"=d))
}

# Code used for testing.  No longer in use.
sibp_llikelihood<-function(Y, X, K, alpha, a, b, sigmasq.A, sigmasq.n, params){
  pi <- params$lambda[,1]/(params$lambda[,1] + params$lambda[,2])
  tau <- params$c/params$d
  beta <- params$m
  A <- params$phi
  Z <- matrix(as.numeric(params$nu >= 0.5), nrow = nrow(params$nu))
  
  pi.component <- sum(dbeta(pi, alpha/K, 1, log = TRUE))
  A.component <- sum(apply(A, 1, function(a) sum(dnorm(a, mean = 0, sd = sqrt(sigmasq.A), log = TRUE))))
  Z.component <- sum(apply(Z, 1, function(z) sum(dbinom(z, size = 1, prob = pi, log = TRUE)))) 
  
  X.component <- sum(dnorm(as.matrix(X), mean = Z%*%A, sd = sqrt(sigmasq.n), log = TRUE))
  
  tau.component <- dgamma(tau, shape = a, rate = b, log = TRUE)
  beta.component <- sum(dnorm(beta, mean = 0, sd = 1/tau, log = TRUE))
  Y.component <- sum(dnorm(Y, mean = Z%*%beta, sd = 1/tau, log = TRUE))
  
  return(as.numeric(pi.component + Z.component + A.component + X.component + tau.component + 
                      beta.component + Y.component))
}

find_sigmasqn<-function(Y, X, K, alpha, a = 0.1, b=0.1, sigmasq.A=5, num.words = 15, iters = 1){
  paramseq<-seq(0.1, 0.8, by = 0.1)
  best.params <- list()
  best.sigmasq.n <- list()
  best.val <- -Inf
  for (sigmasq.n in paramseq){
    print(sigmasq.n)
    for (i in 1:iters){
      this.params <- estimate_parameters(Y, X, K, alpha, a, b, sigmasq.A, sigmasq.n)
      #this.val <- total_coherence(this.params, topic_coherence(this.params, X, num.words))
      this.val <- sibp_llikelihood(Y, X, K, alpha, a, b, sigmasq.A, sigmasq.n, this.params)
      if (this.val > best.val){
        best.val <- this.val
        best.params <- this.params
        best.sigmasq.n <- sigmasq.n
      }
    }
  }
  return(list("sigmasq.n" = best.sigmasq.n, "params" = best.params))
}

# Code used to automatically choose K.  Currently not in use.
harmonic<-function(N){
  return(sum(rep(1,N)/seq(1:N)))
}

# The CE metric reported in the paper.  Used to automatically evaluate the semantic coherence of topics.
topic_exclusivity<-function(params, X, num.words){
  top.words <- apply(params$phi, 1, order, decreasing = TRUE)[1:num.words,]
  C <- array()
  for (t in 1:ncol(params$nu)){
    C[t] <- 0
    docs.in.topic <- which(params$nu[,t] > 0.5)
    for (m in 2:num.words){
      for (l in 1:(m-1)){
        word1 <- top.words[m,t]
        word2 <- top.words[l,t]
        if (length(docs.in.topic) > 1){
          C[t] <- C[t] + cov(X[docs.in.topic,word1], X[docs.in.topic,word2])*length(docs.in.topic)
        }
        if (length(docs.in.topic) > 1 & length(docs.in.topic) < (nrow(X)-1)){
          C[t] <- C[t] - cov(X[-docs.in.topic,word1], X[-docs.in.topic,word2]) * 
            (nrow(X) - length(docs.in.topic))
        }
      }
    }
    C[t] <- C[t]/((num.words^2 - num.words)/2)
  }
  return(sum(C))
}

# Predicts the Z on the test set
predict_Z<-function(params, N, K, X, alpha, sigmasq.n){
  nu <- init_nu(N, K, alpha)
  phi <- params$phi
  big.Phi <- params$big.Phi
  lambda <- params$lambda
  
  converged <- FALSE
  while(!converged)
  {
    old.nu <- nu
    pi.component <- apply(lambda, 1, function(lam) digamma(lam[1]) - digamma(lam[2]))
    
    big.Phi.tr <- sapply(1:K, function(k) sum(diag(big.Phi[[k]])))
    X.component <- -1/(2*sigmasq.n) * t(apply(apply(phi, 1, function(ph) -2*ph%*%t(X) + sum(ph^2)), 1, 
                                              function(x) x + big.Phi.tr))
    
    for (k in 1:K){
      X.sumexcept <- -1/(2*sigmasq.n) * as.numeric(2*phi[k,]%*%t(nu%*%phi - nu[,k,drop=FALSE]%*%phi[k,,drop=FALSE]))
      v <- pi.component[k] + X.component[,k] + X.sumexcept
      nu[,k] <- (1 + exp(-v))^(-1)
      
    }  
    if (max(abs(nu - old.nu)) < 10^(-3)) converged <- TRUE
  }

  return(nu)
}
