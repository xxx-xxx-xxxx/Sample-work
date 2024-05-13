sim.logS <- function(S_0, mu, sigma, T_, N) {
  dt <- T_/N
  muT <- (mu - 0.5*sigma^2)

  logS_t <- vector("double", N)
  logS_t[1] <- as.numeric(log(S_0))
  
  BM <- rnorm(N,mean = 0, sd = sqrt(dt))
  
  for (i in 2:N) {
    logS_t[i] <- muT*dt + sigma*BM[i-1] + logS_t[i - 1]
  }
  return(logS_t)
}

parameters <- read.table("input.txt", sep=",")


#TEST, SET PARAMETERS
S_0 <- as.numeric(parameters[1])
B_0 <- as.numeric(parameters[2])
mu <- as.numeric(parameters[3])
sigma <- as.numeric(parameters[4])
r <- as.numeric(parameters[5])
T_ <- as.numeric(parameters[6])
p <- as.numeric(parameters[7])
N <- as.numeric(parameters[8])
M <- as.numeric(parameters[9])

dt <- T_/N
t <- 1:N

#B
f_1 <- function(t) (log(B_0) + r*dt*(t))
XV_diff <- vector("double", M)

for (j in 1:M){
  
  #S_t
  logS_t <- sim.logS(S_0,mu, sigma, T_, N)
  S <- exp(logS_t)
  #####################################################################################
  
  #v(S_t,t)
  v_t <- vector("double", N)
  v_t[1] <- exp(-r*(T_))*((S[1])^p)*exp(p*(r - 0.5*sigma^2)*(T_) + 0.5*(p^2)*sigma^2*(T_))
  
  for (i in 2:N) {
    v_t[i] <- exp(-r*(T_ - i*dt))*((S[i])^p)*exp(p*(r - 0.5*sigma^2)*(T_ - i*dt) + 0.5*(p^2)*sigma^2*(T_ - i*dt))
  }
  ###########################################################################################
  
  # Delta_t
  Delta_t <- vector("double", N)
  Delta_t[1] <- p*exp(-r*(T_))*((S[1])^(p - 1))*exp(p*(r - 0.5*sigma^2)*(T_) + 0.5*(p^2)*(sigma^2)*(T_))
  for (i in 2:N) {
    Delta_t[i] <- p*exp(-r*(T_ - i*dt))*((S[i])^(p - 1))*exp(p*(r - 0.5*sigma^2)*(T_ - i*dt) + 0.5*(p^2)*(sigma^2)*(T_ - i*dt))
  }
  ###########################################################################################
  
  # X_t
  X_t <- vector("double", N)
  X_t[1] <- v_t[1]
  for (i in 2:N) {
    X_t[i] <- as.numeric(Delta_t[i - 1])*(S[i] - S[i - 1]) + (X_t[i - 1] - Delta_t[i - 1]*S[i - 1])*(exp(f_1(i)) - exp(f_1(i - 1)))/exp(f_1(i - 1)) + X_t[i - 1]
  }
  
  XV_diff[j] <- X_t[N] - v_t[N]
}



png("partf.png")
hist(XV_diff, breaks=20)

dev.off()


