library(usethis)
library(MASS)
set.seed(475646)

D <- 50
Intercept <- runif(D, 5, 9)
LFC <- c(rep(0, D*0.4), runif(D*0.6, 0.25, 2.5))
LipidConcentrationEff <- c(rep(0, D/2), runif(D/2, -3, 0))
B <- cbind(Intercept, LFC)

N <- 2
X <- rbind(rep(1, N), c(rep(0, N/2), rep(1, N/2)))

means <- B%*%X
A <- 2^apply(means, c(1,2), function(item) rnorm(1, item, 0))
mean(log2(colSums(A)[((N/2)+1):N]))-mean(log2(colSums(A)[1:(N/2)]))
mean(log2(colSums(A)[1:(N/2)]))

N <- 100
X <- rbind(rep(1, N), c(rep(0, N/2), rep(1, N/2)))
means <- B%*%X
A <- 2^apply(means, c(1,2), function(item) rnorm(1, item, 1))

R <- apply(A, 2, function(col) col/sum(col))
Y <- apply(R, 2, function(col) rmultinom(1, 5000, col))
data_sim_1 <- list(seq_counts=Y,
                   absolute_abundances=A,
                   lfc=LFC,
                   intercepts=Intercept)


setwd("..")
usethis::use_data(data_sim_1, overwrite=T)
