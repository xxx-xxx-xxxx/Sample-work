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

f_1 <- function(t) (log(B_0) + r*dt*(t))

parameters <- read.table("input.txt", sep=",")

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

logS_t <- sim.logS(S_0,mu, sigma, T_, N)
S <- exp(logS_t)

s_t <- data.frame(t,S)

png("partc.png")

plot(t, S, type="l", ylim = c(min(S), max(S)))
lines(t, exp(f_1(t)), type="l", ylim = c(1,max(exp(f_1(t)))), col = 2)
dev.off()






