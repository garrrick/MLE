mle_Newton_Raphson <- function(approx.weight, observation, censor, 
                               shape.x, scale.x, curetime.x, 
                               general.population.shape, general.population.scale, 
                               shape.initials, scale.initials, curetime.initials, 
                               smooth.parameter, tolerance = 10^(-7), stop.count = 100, unstable.stop.count = 10){
  # observation = obs; approx.weight = length(observation); censor = cen; shape.x = gxx; scale.x = pxx; curetime.x = cxx; general.population.shape = general.shape; general.population.scale = general.scale; shape.initials = initial.shape; scale.initials = initial.scale; curetime.initials = initial.curetime; smooth.parameter = ll.initial; tolerance = 10^(-3); stop.count = 50; unstable.stop.count = 10
  # x = c(shape.initials, scale.initials, curetime.initials)
  shape.length <- 1:length(shape.initials)
  scale.length <- (length(shape.initials)+1):(length(shape.initials)+length(scale.initials))
  Weibull.length <- c(shape.length, scale.length)
  curetime.length <- (length(shape.initials)+length(scale.initials)+1):(length(shape.initials)+length(scale.initials)+length(curetime.initials))
  library(matrixcalc)
  size <- length(observation)
  llh <- function(x){
    x1 <- x[1:length(shape.initials)]
    e.shape <- exp(crossprod(x1, t(shape.x)))
    x2 <- x[(length(shape.initials)+1):(length(shape.initials)+length(scale.initials))]
    e.scale <- exp(crossprod(x2, t(scale.x)))
    x3 <- x[(length(shape.initials)+length(scale.initials)+1):length(x)]
    e.ttc <- exp(crossprod(x3, t(curetime.x)))
    s.p1 <- exp(-(observation/e.scale)^e.shape)
    s.p1.ttc <- exp(-(e.ttc/e.scale)^e.shape)
    d.p1 <- (e.shape/e.scale)*((observation/e.scale)^(e.shape-1))*exp(-(observation/e.scale)^e.shape)
    h.p1 <- (e.shape/e.scale)*((observation/e.scale)^(e.shape-1))
    approx.sigma <- ssize^(-1/approx.weight)
    ind.approx <- (1+exp(-(e.ttc-observation)*approx.sigma^(-1)))^(-1)
    s.p <- pweibull(observation, shape = general.population.shape, scale = general.population.scale, lower.tail = FALSE)
    d.p <- dweibull(observation, shape = general.population.shape, scale = general.population.scale)
    h.p <- d.p/s.p
    penalty <- (2^(-1))*smooth.parameter*crossprod(x3)
    -(mean(censor*log(h.p + h.p1*ind.approx) +
             log((s.p1 - s.p1.ttc)*ind.approx + s.p1.ttc)) - penalty)
  }
  llh.indicator <- function(x){
    x1 <- x[1:length(shape.initials)]
    e.shape <- exp(crossprod(x1, t(shape.x)))
    x2 <- x[(length(shape.initials)+1):(length(shape.initials)+length(scale.initials))]
    e.scale <- exp(crossprod(x2, t(scale.x)))
    x3 <- x[(length(shape.initials)+length(scale.initials)+1):length(x)]
    e.ttc <- exp(crossprod(x3, t(curetime.x)))
    s.p1 <- exp(-(observation/e.scale)^e.shape)
    s.p1.ttc <- exp(-(e.ttc/e.scale)^e.shape)
    d.p1 <- (e.shape/e.scale)*((observation/e.scale)^(e.shape-1))*exp(-(observation/e.scale)^e.shape)
    h.p1 <- (e.shape/e.scale)*((observation/e.scale)^(e.shape-1))
    ind <- I(e.ttc - observation > 0)
    s.p <- pweibull(observation, shape = general.population.shape, scale = general.population.scale, lower.tail = FALSE)
    d.p <- dweibull(observation, shape = general.population.shape, scale = general.population.scale)
    h.p <- d.p/s.p
    penalty <- (2^(-1))*smooth.parameter*crossprod(x3)
    mean(censor*log(h.p + h.p1*ind) +
           log((s.p1 - s.p1.ttc)*ind + s.p1.ttc)) - penalty
  }
  hessian <- function(x){
    x1 <- matrix(x[1:length(shape.initials)], nc = 1) # length(shape.initials)*1
    e.shape <- exp(crossprod(t(shape.x), x1)) # sample.size*1
    x2 <- matrix(x[(length(shape.initials)+1):(length(shape.initials)+length(scale.initials))], nc = 1) # length(shape.initials)*1
    e.scale <- exp(crossprod(t(scale.x), x2)) # sample.size*1
    x3 <- matrix(x[(length(shape.initials)+length(scale.initials)+1):length(x)], nc = 1) # length(curetime.initials)*1
    e.ttc <- exp(crossprod(t(curetime.x), x3)) # sample.size*1

    d.shape <- matrix(e.shape, nc = length(shape.initials), nr = size, byrow = FALSE)*shape.x # sample.size*length(shape.initials)
    d.scale <- matrix(e.scale, nc = length(scale.initials), nr = size, byrow = FALSE)*scale.x # sample.size*length(scale.initials)
    d.curetime <- matrix(e.ttc, nc = length(curetime.initials), nr = size, byrow = FALSE)*curetime.x # sample.size*length(curetime.initials)

    s.p1 <- exp(-(observation/e.scale)^e.shape) # sample.size*1
    s.p1.matrix <- matrix(s.p1, nc = length(x), nr = size, byrow = FALSE)
    s.p1.ttc <- exp(-(e.ttc/e.scale)^e.shape) # sample.size*1
    s.p1.ttc.matrix <- matrix(s.p1.ttc, nc = length(x), nr = size, byrow = FALSE)
    d.p1 <- (e.shape/e.scale)*(observation/e.scale)^(e.shape-1)*exp(-(observation/e.scale)^e.shape) # sample.size*1
    d.p1.matrix <- matrix(d.p1, nc = length(x), nr = size, byrow = FALSE)
    d.p1.ttc <- (e.shape/e.scale)*(e.ttc/e.scale)^(e.shape-1)*exp(-(e.ttc/e.scale)^e.shape) # sample.size*1
    d.p1.ttc.matrix <- matrix(d.p1.ttc, nc = length(x), nr = size, byrow = FALSE)
    h.p1 <- (e.shape/e.scale)*(observation/e.scale)^(e.shape-1) # sample.size*1
    h.p1.matrix <- matrix(h.p1, nc = length(x), nr = size, byrow = FALSE)
    h.p1.ttc <- (e.shape/e.scale)*(e.ttc/e.scale)^(e.shape-1) # sample.size*1
    h.p1.ttc.matrix <- matrix(h.p1.ttc, nc = length(x), nr = size, byrow = FALSE)
    approx.sigma <- size^(-1/approx.weight)
    ind.approx <- (1+exp(-(e.ttc-observation)*approx.sigma^(-1)))^(-1)
    ind.approx.matrix <- matrix(ind.approx, nc = length(x), nr = size, byrow = FALSE)
    s.p <- matrix(pweibull(observation, shape = general.population.shape, scale = general.population.scale, lower.tail = FALSE), nc = 1, nr = size) # sample.size*1
    s.p.matrix <- matrix(s.p, nc = length(x), nr = size, byrow = FALSE)
    d.p <- matrix(dweibull(observation, shape = general.population.shape, scale = general.population.scale), nc = 1, nr = size) # sample.size*1
    d.p.matrix <- matrix(d.p, nc = length(x), nr = size, byrow = FALSE)
    h.p <- d.p/s.p # sample.size*1
    h.p.matrix <- matrix(h.p, nc = length(x), nr = size, byrow = FALSE)
    penalty <- (2^(-1))*smooth.parameter*crossprod(x3) # 1*1


    derivative.s.p1.1 <- matrix(-((observation/e.scale)^e.shape)*(log(observation)-log(e.scale))*s.p1, nc = length(shape.initials), nr = size, byrow = FALSE)*d.shape
    derivative.s.p1.2 <- matrix((observation/e.scale)*d.p1, nc = length(scale.initials), nr = size, byrow = FALSE)*d.scale
    derivative.s.p1 <- cbind(derivative.s.p1.1, 
                             derivative.s.p1.2) # sample.size*(length(shape.initials) + length(scale.initials))
    derivative.h.p1 <- cbind(matrix((1+e.shape*(log(observation)-log(e.scale)))*(e.scale^(-1))*(observation/e.scale)^(e.shape-1), nc = length(shape.initials), nr = size, byrow = FALSE)*d.shape, 
                             matrix(-(e.shape*(e.scale)^(-1))*h.p1, nc = length(scale.initials), nr = size, byrow = FALSE)*d.scale) # sample.size*(length(shape.initials) + length(scale.initials))
    derivative.ind.approx <- matrix(approx.sigma^(-1)*ind.approx*(1-ind.approx), nc = length(curetime.initials), nr = size, byrow = FALSE)*d.curetime # sample.size*length(curetime.initials)
    derivative.penalty <- matrix(smooth.parameter*x3, nc = length(curetime.initials), nr = size, byrow = TRUE) # sample.size*length(curetime.initials)
    derivative.s.p1.ttc1 <- matrix(-((e.ttc/e.scale)^e.shape)*(log(e.ttc)-log(e.scale))*s.p1.ttc, nc = length(shape.initials), nr = size, byrow = FALSE)*d.shape
    derivative.s.p1.ttc2 <- matrix((e.ttc/e.scale)*d.p1.ttc, nc = length(scale.initials), nr = size, byrow = FALSE)*d.scale
    derivative.s.p1.ttc3 <- matrix(-d.p1.ttc, nc = length(curetime.initials), nr = size, byrow = FALSE)*d.curetime
    derivative.s.p1.ttc <- cbind(derivative.s.p1.ttc1, 
                                 derivative.s.p1.ttc2,
                                 derivative.s.p1.ttc3) # sample.size*(length(shape.initials) + length(scale.initials) + length(curetime.initials))



    second.derivative.ind.approx.constant <- matrix(e.ttc*approx.sigma^(-1)*ind.approx*(1-ind.approx)*(1+approx.sigma^(-1)*e.ttc*(1-2*ind.approx)), nc = length(curetime.initials), nr = size, byrow = FALSE)
    second.derivative.h.p1.constant11 <- matrix(h.p1*((1+e.shape*(log(observation)-log(e.scale)))^2+e.shape*(log(observation)-log(e.scale))), nc = length(shape.initials), nr = size, byrow = FALSE)
    second.derivative.h.p1.constant22 <- matrix(e.shape^2*h.p1, nc = length(scale.initials), nr = size, byrow = FALSE)
    second.derivative.h.p1.constant12 <- matrix(-e.shape*h.p1*(2+e.shape*(log(observation)-log(e.scale))), nc = length(shape.initials), nr = size)
    second.derivative.s.p1.constant11 <- matrix(-d.p1*observation*(log(observation)-log(e.scale))*(1+e.shape*(log(observation)-log(e.scale))*(1-(observation/e.scale)^e.shape)), nc = length(shape.initials), nr = size, byrow = FALSE)
    second.derivative.s.p1.constant22 <- matrix(observation*d.p1*(observation*h.p1-e.shape), nc = length(scale.initials), nr = size, byrow = FALSE)
    second.derivative.s.p1.constant12 <- matrix(observation*d.p1*(1+e.shape*(log(observation)-log(e.scale))*(1-(observation/e.scale)^e.shape)), nc = length(shape.initials), nr = size, byrow = FALSE)
    second.derivative.s.p1.ttc.constant11 <- matrix(-d.p1.ttc*e.ttc*(log(e.ttc)-log(e.scale))*(1+e.shape*(log(e.ttc)-log(e.scale))*(1-(e.ttc/e.scale)^e.shape)), nc = length(shape.initials), nr = size, byrow = FALSE)
    second.derivative.s.p1.ttc.constant22 <- matrix(e.ttc*d.p1.ttc*(e.ttc*h.p1.ttc-e.shape), nc = length(scale.initials), nr = size, byrow = FALSE)
    second.derivative.s.p1.ttc.constant33 <- matrix(e.ttc*d.p1.ttc*(e.ttc*h.p1.ttc-e.shape), nc = length(curetime.initials), nr = size, byrow = FALSE)
    second.derivative.s.p1.ttc.constant12 <- matrix(e.ttc*d.p1.ttc*(1+e.shape*(log(e.ttc)-log(e.scale))*(1-(e.ttc/e.scale)^e.shape)), nc = length(shape.initials), nr = size, byrow = FALSE)
    second.derivative.s.p1.ttc.constant23 <- matrix(-e.ttc*d.p1.ttc*(e.ttc*h.p1.ttc-e.shape), nc = length(scale.initials), nr = size, byrow = FALSE)
    second.derivative.s.p1.ttc.constant13 <- matrix(-e.ttc*d.p1.ttc*(1+e.shape*(log(e.ttc)-log(e.scale))*(1-(e.ttc/e.scale)^e.shape)), nc = length(shape.initials), nr = size, byrow = FALSE)
    
    
    
    censor.matrix <- matrix(censor, nc = length(x), nr = size, byrow = FALSE) # sample.size*(length(shape.initials) + length(scale.initials) + length(curetime.initials))
    model.hazard <- h.p + h.p1*ind.approx
    model.survival <- (s.p1 - s.p1.ttc)*ind.approx + s.p1.ttc
    
    
    
    derivative.model.hazard <- cbind(
                                     matrix(ind.approx, nc = length(shape.initials) + length(scale.initials), nr = size, byrow = FALSE)*derivative.h.p1, 
                                     matrix(h.p1, nc = length(curetime.initials), nr = size, byrow = FALSE)*derivative.ind.approx
                                     ) # sample.size*(length(shape.initials) + length(scale.initials) + length(curetime.initials))
    derivative.model.survival <- cbind(
                                       matrix(ind.approx, nc = length(shape.initials) + length(scale.initials), nr = size, byrow = FALSE)*derivative.s.p1,
                                       matrix(s.p1 - s.p1.ttc, nc = length(curetime.initials), nr = size, byrow = FALSE)*derivative.ind.approx
                                       ) + 
                                 matrix(1-ind.approx, nc = length(x), nr = size, byrow = FALSE)*derivative.s.p1.ttc
    
    
    
    second.derivative.model.hazard.1 <- crossprod(matrix(-censor*model.hazard^(-2), nc = length(x), nr = size)*
                                                  derivative.model.hazard, derivative.model.hazard)
    
    
    second.derivative.model.survival.1 <- crossprod(matrix(-model.survival^(-2), nc = length(x), nr = size)*
                                                    derivative.model.survival, derivative.model.survival)
    
    
    second.derivative.model.hazard.2 <- matrix(NA, nc = length(x), nr = length(x))
    second.derivative.model.hazard.2[1:length(shape.initials),1:length(shape.initials)] <- 
      crossprod(matrix(censor*model.hazard^(-1)*ind.approx, nc = length(shape.initials), nr = size)*
                  second.derivative.h.p1.constant11*shape.x, shape.x)
    second.derivative.model.hazard.2[(length(shape.initials)+1):(length(shape.initials) + length(scale.initials)),(length(shape.initials)+1):(length(shape.initials) + length(scale.initials))] <- 
      crossprod(matrix(censor*model.hazard^(-1)*ind.approx, nc = length(scale.initials), nr = size)*
                  second.derivative.h.p1.constant22*scale.x, scale.x)
    second.derivative.model.hazard.2[(length(shape.initials) + length(scale.initials) + 1):length(x),(length(shape.initials) + length(scale.initials) + 1):length(x)] <- 
      crossprod(matrix(censor*model.hazard^(-1)*h.p1, nc = length(curetime.initials), nr = size)*
                  second.derivative.ind.approx.constant*curetime.x, curetime.x)
    second.derivative.model.hazard.2[1:(length(shape.initials) + length(scale.initials)),(length(shape.initials) + length(scale.initials) + 1):length(x)] <- 
      crossprod(matrix(censor*model.hazard^(-1), nc = length(shape.initials) + length(scale.initials), nr = size)*
                  derivative.h.p1, derivative.ind.approx)
    second.derivative.model.hazard.2[(length(shape.initials) + length(scale.initials) + 1):length(x), 1:(length(shape.initials) + length(scale.initials))] <- 
      t(crossprod(matrix(censor*model.hazard^(-1), nc = length(shape.initials) + length(scale.initials), nr = size)*
                  derivative.h.p1, derivative.ind.approx))
    second.derivative.model.hazard.2[1:length(shape.initials),(length(shape.initials) + 1):(length(shape.initials) + length(scale.initials))] <- 
      crossprod(matrix(censor*model.hazard^(-1)*ind.approx, nc = length(shape.initials), nr = size)*
                  second.derivative.h.p1.constant12*shape.x, scale.x)
    second.derivative.model.hazard.2[(length(shape.initials) + 1):(length(shape.initials) + length(scale.initials)),1:length(shape.initials)] <- 
      t(crossprod(matrix(censor*model.hazard^(-1)*ind.approx, nc = length(shape.initials), nr = size)*
                  second.derivative.h.p1.constant12*shape.x, scale.x))


    second.derivative.model.survival.2 <- matrix(NA, nc = length(x), nr = length(x))
    second.derivative.model.survival.2[1:length(shape.initials),1:length(shape.initials)] <- 
      crossprod(matrix(model.survival^(-1)*ind.approx, nc = length(shape.initials), nr = size)*
                (second.derivative.s.p1.constant11-second.derivative.s.p1.ttc.constant11)*shape.x, shape.x) + 
      crossprod(matrix(model.survival^(-1), nc = length(shape.initials), nr = size)*
                  second.derivative.s.p1.ttc.constant11*shape.x, shape.x)
    second.derivative.model.survival.2[(length(shape.initials)+1):(length(shape.initials) + length(scale.initials)),(length(shape.initials)+1):(length(shape.initials) + length(scale.initials))] <- 
      crossprod(matrix(model.survival^(-1)*ind.approx, nc = length(scale.initials), nr = size)*
                  (second.derivative.s.p1.constant22-second.derivative.s.p1.ttc.constant22)*scale.x, scale.x) + 
      crossprod(matrix(model.survival^(-1), nc = length(scale.initials), nr = size)*
                  second.derivative.s.p1.ttc.constant22*scale.x, scale.x)
    second.derivative.model.survival.2[(length(shape.initials) + length(scale.initials) + 1):length(x),(length(shape.initials) + length(scale.initials) + 1):length(x)] <-   
      crossprod(matrix(model.survival^(-1)*(1-ind.approx), nc = length(curetime.initials), nr = size)*
                  second.derivative.s.p1.ttc.constant33*curetime.x, curetime.x) + 
      crossprod(matrix(model.survival^(-1)*(s.p1-s.p1.ttc), nc = length(curetime.initials), nr = size)*
                  second.derivative.ind.approx.constant*curetime.x, curetime.x) -
      crossprod(matrix(model.survival^(-1), nc = length(curetime.initials), nr = size)*
                  derivative.s.p1.ttc3, derivative.ind.approx) -
      t(crossprod(matrix(model.survival^(-1), nc = length(curetime.initials), nr = size)*
                  derivative.s.p1.ttc3, derivative.ind.approx))
    second.derivative.model.survival.2[1:(length(shape.initials) + length(scale.initials)), (length(shape.initials) + length(scale.initials) + 1):length(x)] <- 
      crossprod(matrix(model.survival^(-1)*(1-ind.approx), nc = length(shape.initials) + length(scale.initials), nr = size)*
                  cbind(second.derivative.s.p1.ttc.constant13, second.derivative.s.p1.ttc.constant23)*cbind(shape.x, scale.x), curetime.x) + 
      crossprod(matrix(model.survival^(-1), nc = length(shape.initials) + length(scale.initials), nr = size)*
                  (derivative.s.p1 - cbind(derivative.s.p1.ttc1, derivative.s.p1.ttc2)), derivative.ind.approx)
    second.derivative.model.survival.2[(length(shape.initials) + length(scale.initials) + 1):length(x), 1:(length(shape.initials) + length(scale.initials))] <-
      t(crossprod(matrix(model.survival^(-1)*(1-ind.approx), nc = length(shape.initials) + length(scale.initials), nr = size)*
                  cbind(second.derivative.s.p1.ttc.constant13, second.derivative.s.p1.ttc.constant23)*cbind(shape.x, scale.x), curetime.x) + 
      crossprod(matrix(model.survival^(-1), nc = length(shape.initials) + length(scale.initials), nr = size)*
                  (derivative.s.p1 - cbind(derivative.s.p1.ttc1, derivative.s.p1.ttc2)), derivative.ind.approx))
    second.derivative.model.survival.2[1:length(shape.initials),(length(shape.initials) + 1):(length(shape.initials) + length(scale.initials))] <-
      crossprod(matrix(model.survival^(-1)*ind.approx, nc = length(shape.initials), nr = size)*
                  second.derivative.s.p1.constant12*shape.x, scale.x) +
      crossprod(matrix(model.survival^(-1)*(1-ind.approx), nc = length(shape.initials), nr = size)*
                  second.derivative.s.p1.ttc.constant12*shape.x, scale.x)
    second.derivative.model.survival.2[(length(shape.initials) + 1):(length(shape.initials) + length(scale.initials)),1:length(shape.initials)] <- 
      t(crossprod(matrix(model.survival^(-1)*ind.approx, nc = length(shape.initials), nr = size)*
                  second.derivative.s.p1.constant12*shape.x, scale.x) +
      crossprod(matrix(model.survival^(-1)*(1-ind.approx), nc = length(shape.initials), nr = size)*
                  second.derivative.s.p1.ttc.constant12*shape.x, scale.x))
    
      
    
    second.derivative.penalty <- matrix(0, nc = length(x), nr = length(x))
    second.derivative.penalty[((length(shape.initials)+length(scale.initials)+1):length(x)),((length(shape.initials)+length(scale.initials)+1):length(x))] <- smooth.parameter*diag(length(curetime.initials))
    
    
    
    (second.derivative.model.hazard.1 + 
      second.derivative.model.hazard.2 + 
      second.derivative.model.survival.1 + 
      second.derivative.model.survival.2)/size - 
      second.derivative.penalty
  }
  score <- function(x){
    x1 <- matrix(x[1:length(shape.initials)], nc = 1)
    e.shape <- exp(crossprod(t(shape.x), x1))
    x2 <- matrix(x[(length(shape.initials)+1):(length(shape.initials)+length(scale.initials))], nc = 1)
    e.scale <- exp(crossprod(t(scale.x), x2))
    x3 <- matrix(x[(length(shape.initials)+length(scale.initials)+1):length(x)], nc = 1)
    e.ttc <- exp(crossprod(t(curetime.x), x3))
    d.shape <- matrix(e.shape, nc = length(shape.initials), nr = size, byrow = FALSE)*shape.x
    d.scale <- matrix(e.scale, nc = length(scale.initials), nr = size, byrow = FALSE)*scale.x
    d.curetime <- matrix(e.ttc, nc = length(curetime.initials), nr = size, byrow = FALSE)*curetime.x
    s.p1 <- exp(-(observation/e.scale)^e.shape)
    s.p1.ttc <- exp(-(e.ttc/e.scale)^e.shape)
    d.p1 <- (e.shape/e.scale)*(observation/e.scale)^(e.shape-1)*exp(-(observation/e.scale)^e.shape)
    d.p1.ttc <- (e.shape/e.scale)*(e.ttc/e.scale)^(e.shape-1)*exp(-(e.ttc/e.scale)^e.shape)
    h.p1 <- (e.shape/e.scale)*(observation/e.scale)^(e.shape-1)
    h.p1.ttc <- (e.shape/e.scale)*(e.ttc/e.scale)^(e.shape-1)
    approx.sigma <- size^(-1/approx.weight)
    ind.approx <- (1+exp(-(e.ttc-observation)*approx.sigma^(-1)))^(-1)
    s.p <- matrix(pweibull(observation, shape = general.population.shape, scale = general.population.scale, lower.tail = FALSE), nc = 1, nr = size)
    d.p <- matrix(dweibull(observation, shape = general.population.shape, scale = general.population.scale), nc = 1, nr = size)
    h.p <- d.p/s.p
    penalty <- (2^(-1))*smooth.parameter*crossprod(x3)
    derivative.s.p1 <- cbind(matrix(-(observation/e.scale)^e.shape*log(observation/e.scale)*s.p1, nc = length(shape.initials), nr = size, byrow = FALSE)*d.shape, 
                             matrix((observation/e.scale)*d.p1, nc = length(scale.initials), nr = size, byrow = FALSE)*d.scale,
                             matrix(0, nc = length(curetime.initials), nr = size))
    derivative.h.p1 <- cbind(matrix((1+e.shape*(log(observation)-log(e.scale)))*(e.scale^(-1))*(observation/e.scale)^(e.shape-1), nc = length(shape.initials), nr = size, byrow = FALSE)*d.shape, 
                             matrix(-(e.shape*(e.scale)^(-1))*h.p1, nc = length(scale.initials), nr = size, byrow = FALSE)*d.scale,
                             matrix(0, nc = length(curetime.initials), nr = size))
    derivative.ind.approx <- cbind(matrix(0, nc = length(shape.initials), nr = size), 
                                   matrix(0, nc = length(scale.initials), nr = size),
                                   matrix(approx.sigma^(-1)*ind.approx*(1-ind.approx), nc = length(curetime.initials), nr = size, byrow = FALSE)*d.curetime)
    derivative.penalty <- cbind(matrix(0, nc = length(shape.initials), nr = size), 
                                matrix(0, nc = length(scale.initials), nr = size),
                                matrix(smooth.parameter*x3, nc = length(curetime.initials), nr = size, byrow = TRUE))
    derivative.s.p1.ttc <- cbind(matrix(-(e.ttc/e.scale)^e.shape*log(e.ttc/e.scale)*s.p1.ttc, nc = length(shape.initials), nr = size, byrow = FALSE)*d.shape, 
                                 matrix((e.ttc/e.scale)*d.p1.ttc, nc = length(scale.initials), nr = size, byrow = FALSE)*d.scale,
                                 matrix(-d.p1.ttc, nc = length(curetime.initials), nr = size, byrow = FALSE)*d.curetime)
    censor.matrix <- matrix(censor, nc = length(x), nr = size, byrow = FALSE)
    model.hazard <- matrix(h.p + h.p1*ind.approx, nc = length(x), nr = size, byrow = FALSE)
    model.survival <- matrix((s.p1 - s.p1.ttc)*ind.approx + s.p1.ttc, nc = length(x), nr = size, byrow = FALSE)
    matrix(colMeans(censor.matrix*model.hazard^(-1)*(matrix(h.p1, nc = length(x), nr = size)*derivative.ind.approx + 
                                                         matrix(ind.approx, nc = length(x), nr = size)*derivative.h.p1) + 
                      model.survival^(-1)*((derivative.s.p1 - derivative.s.p1.ttc)*matrix(ind.approx, nc = length(x), nr = size) + 
                                               (matrix(s.p1 - s.p1.ttc, nc = length(x), nr = size))*derivative.ind.approx + derivative.s.p1.ttc) - 
                      derivative.penalty), nc = 1, nr = length(x), byrow = FALSE)
  }
  Newton.Raphson <- function(x, rate) {
    x-rate*crossprod(solve(hessian(x)), score(x))
  }
  x0 <- x0.initials <- c(shape.initials, scale.initials, curetime.initials)
  unstable.count <- 0
  if (sum(is.na(x0.initials)) == 0 & 
      sum(is.nan(x0.initials)) == 0 & 
      sum(is.infinite(x0.initials)) == 0 & 
      !is.na(llh(x0.initials)) & 
      !is.nan(llh(x0.initials)) & 
      !is.infinite(llh(x0.initials)) &
      sum(is.na(hessian(x0.initials))) == 0 &
      sum(is.nan(hessian(x0.initials))) == 0 & 
      sum(is.infinite(hessian(x0.initials))) == 0) {
    if (!is.singular.matrix(hessian(x0.initials)) &
         is.positive.definite(round(-hessian(x0.initials), 4))) {
      x0 <- x0.initials
    } else {
      if (unstable.count == unstable.stop.count) {
        estimate <- x0.initials
        estimate_converge <- 'no'
        break
      } else {
        unstable.count <- unstable.count + 1
        x0[Weibull.length] <- 
          x0.initials[Weibull.length] + runif(length(x0.initials[Weibull.length]), -1, 1)/sqrt(crossprod(x0.initials[Weibull.length]+.001))
        x0[curetime.length] <- 
          x0.initials[curetime.length] + runif(length(x0.initials[curetime.length]), -1, -.5)/sqrt(crossprod(x0.initials[curetime.length]+.001))
      }
    }
  } else {
    if (unstable.count == unstable.stop.count) {
      estimate <- x0.initials
      estimate_converge <- 'no'
      break
    } else {
      unstable.count <- unstable.count + 1
      x0[Weibull.length] <- 
        x0.initials[Weibull.length] + runif(length(x0.initials[Weibull.length]), -1, 1)/sqrt(crossprod(x0.initials[Weibull.length]+.001))
      x0[curetime.length] <- 
        x0.initials[curetime.length] + runif(length(x0.initials[curetime.length]), -1, -.5)/sqrt(crossprod(x0.initials[curetime.length]+.001))
    }
  }
  like <- rr <- c()
  line_search <- function(r){
    logl <- llh(x0-r*crossprod(solve(hessian(x0)), score(x0)))
    if (is.na(logl) | is.nan(logl) | is.infinite(logl)) {
      logl <- llh(x0)
    }
    return(logl)
  }
  rate.points <- l <- c()
  seeds <- seq(-5,5, .5)
  for (ll in 1:length(seeds)) {
    rate.points[ll] <- optim(par = 1, fn = line_search, method = 'Brent', lower = -5+seeds[ll], upper = 5+seeds[ll])$par
    l[ll] <- line_search(rate.points[ll])
  }
  dis <- min(l)
  rr[1] <- rate <- rate.points[l == dis][1]
  if (is.na(rate) | is.nan(rate) | is.infinite(rate)) {
    rr[1] <- rate <- 1
  }
  count <- unstable.count <- 0
  repeat{
    like[count + 1] <- llh(x0)
    x1 <- Newton.Raphson(x = x0, rate = rate)
    if (sum(is.na(x1)) == 0 & 
        sum(is.nan(x1)) == 0 & 
        sum(is.infinite(x1)) == 0 & 
        !is.na(llh(x1)) & 
        !is.nan(llh(x1)) & 
        !is.infinite(llh(x1)) &
        sum(is.na(hessian(x1))) == 0 &
        sum(is.nan(hessian(x1))) == 0 & 
        sum(is.infinite(hessian(x1))) == 0) {
      if (!is.singular.matrix(hessian(x1)) &
           is.positive.definite(round(-hessian(x1), 4))) {
        if (sqrt(crossprod(score(x1))) <= tolerance) {
          count <- count + 1
          estimate <- x1
          estimate_converge <- 'yes'
          break
        } else if (count == stop.count) {
          estimate <- x1
          estimate_converge <- 'no'
          break
        } else {
          line_search <- function(r){
            logl <- llh(x1-r*crossprod(solve(hessian(x1)), score(x1)))
            if (is.na(logl) | is.nan(logl) | is.infinite(logl)) {
              logl <- llh(x1)
            }
            return(logl)
          }
          rate.points <- l <- c()
          seeds <- seq(-5,5, .5)
          for (ll in 1:length(seeds)) {
            rate.points[ll] <- optim(par = 1, fn = line_search, method = 'Brent', lower = -5+seeds[ll], upper = 5+seeds[ll])$par
            l[ll] <- line_search(rate.points[ll])
          }
          dis <- min(l)
          rr[count + 2] <- rate <- rate.points[l == dis][1]
          if (is.na(rate) | is.nan(rate) | is.infinite(rate)) {
            rr[count + 2] <- rate <- 1
          }
          count <- count + 1
          x0 <- x1
        }
      } else {
        if (unstable.count == unstable.stop.count) {
          estimate <- x0
          estimate_converge <- 'no'
          break
        } else {
          unstable.count <- unstable.count + 1
          count <- 0
          x0[Weibull.length] <- 
            x0.initials[Weibull.length] + runif(length(x0.initials[Weibull.length]), -1, 1)/sqrt(crossprod(x0.initials[Weibull.length]+.001))
          x0[curetime.length] <- 
            x0.initials[curetime.length] + runif(length(x0.initials[curetime.length]), -1, -.5)/sqrt(crossprod(x0.initials[curetime.length]+.001))
          line_search <- function(r){
            logl <- llh(x0-r*crossprod(solve(hessian(x0)), score(x0)))
            if (is.na(logl) | is.nan(logl) | is.infinite(logl)) {
              logl <- llh(x0)
            }
            return(logl)
          }
          rate.points <- l <- c()
          seeds <- seq(-5,5, .5)
          for (ll in 1:length(seeds)) {
            rate.points[ll] <- optim(par = 1, fn = line_search, method = 'Brent', lower = -5+seeds[ll], upper = 5+seeds[ll])$par
            l[ll] <- line_search(rate.points[ll])
          }
          dis <- min(l)
          rr[1] <- rate <- rate.points[l == dis][1]
          if (is.na(rate) | is.nan(rate) | is.infinite(rate)) {
            rr[1] <- rate <- 1
          }
        }
      }
    } else {
      if (unstable.count == unstable.stop.count) {
        estimate <- x0
        estimate_converge <- 'no'
        break
      } else {
        unstable.count <- unstable.count + 1
        count <- 0
        x0[Weibull.length] <- 
          x0.initials[Weibull.length] + runif(length(x0.initials[Weibull.length]), -1, 1)/sqrt(crossprod(x0.initials[Weibull.length]+.001))
        x0[curetime.length] <- 
          x0.initials[curetime.length] + runif(length(x0.initials[curetime.length]), -1, -.5)/sqrt(crossprod(x0.initials[curetime.length]+.001))
        line_search <- function(r){
          logl <- llh(x0-r*crossprod(solve(hessian(x0)), score(x0)))
          if (is.na(logl) | is.nan(logl) | is.infinite(logl)) {
            logl <- llh(x0)
          }
          return(logl)
        }
        rate.points <- l <- c()
        seeds <- seq(-1,1, .1)
        for (ll in 1:length(seeds)) {
          rate.points[ll] <- optim(par = 1, fn = line_search, method = 'Brent', lower = -5+seeds[ll], upper = 5+seeds[ll])$par
          l[ll] <- line_search(rate.points[ll])
        }
        dis <- min(l)
        rr[1] <- rate <- rate.points[l == dis][1]
        if (is.na(rate) | is.nan(rate) | is.infinite(rate)) {
          rr[1] <- rate <- 1
        }
      }
    }
  }
  est <- list(shape = estimate[shape.length], 
              scale = estimate[scale.length], 
              curetime = estimate[curetime.length], 
              Gradient = score(estimate),
              Observed.Score = score(estimate),
              Hessian = hessian(estimate),
              Observed.Information = -hessian(estimate),
              indicator.loglikelihood = llh.indicator(estimate),
              estimate_converge = estimate_converge)
  return(est)
}