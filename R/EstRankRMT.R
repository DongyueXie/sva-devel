#'@title Estimate rank using RMT, taken from isva package
#'@param data.m matrix, features by samples
#'@param rmax maximum rank
#'@return rank
#'@export


EstRankRMT = function (data.m, rmax = 50){
  n=nrow(data.m)
  p=ncol(data.m)
  rmax = min(c(rmax,min(c(n,p))))
  
  M <- apply(data.m, 2, function(X) {
    (X - mean(X))/sqrt(var(X))
  })
  sigma2 <- var(as.vector(M))
  Q <- nrow(data.m)/ncol(data.m)
  ns <- ncol(data.m)
  lambdaMAX <- sigma2 * (1 + 1/Q + 2 * sqrt(1/Q))
  lambdaMIN <- sigma2 * (1 + 1/Q - 2 * sqrt(1/Q))
  delta <- lambdaMAX - lambdaMIN
  step = delta/ns
  # roundN <- 3
  # step <- round(delta/ns, roundN)
  # while (step == 0) {
  #   roundN <- roundN + 1
  #   step <- round(delta/ns, roundN)
  # }
  lambda.v <- seq(lambdaMIN, lambdaMAX, by = step)
  dens.v <- (Q/(2 * pi * sigma2)) * sqrt((lambdaMAX - lambda.v) * (lambda.v - lambdaMIN))/lambda.v
  thdens.o <- list(min = lambdaMIN, max = lambdaMAX, step = step,
                   lambda = lambda.v, dens = dens.v)
  #C <- 1/nrow(M) * t(M) %*% M
  #eigen.o <- eigen(C, symmetric = TRUE)
  eigen.o = svds(M/sqrt(nrow(M)),rmax)
  intdim <- sum(eigen.o$d^2 > thdens.o$max)
  return(intdim)
}





