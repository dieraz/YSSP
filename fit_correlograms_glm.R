### case 1: exponential
fit_correlogram<- function(z,r,auto=TRUE){
  
  z <- z[1:(length(z)/10)]
  rad <- r[1:(length(r)/10)]
  temp <- data.frame(cbind(rad,z))
  
  # plot data
  #plot(tempprev, xlab= 'Radii', ylab='Pairwise Density')
  we<-1/(1+rad)^3
  # fit non-linear model

  
  if (auto)
  {
    mod1 <- glm(z ~ rad, family = quasipoisson(link = "log"), data = temp, weights = we)
    c_mod1 <-coef(mod1)
    d_mod1 <- abs(1/c_mod1[2])
         
  }
  else
  {
    zcc <- 1 - z/max(z)
    temp <- data.frame(cbind(rad,zcc))
    mod1 <- glm(zcc ~ rad - 1, family = quasipoisson(link = "log"), data = temp, weights = we)
    c_mod1 <-coef(mod1)
    d_mod1 <- 1/c_mod1

  }
  
  
  return(d_mod1)
  
}

### case 2: exponential
fit_correlogram_glm1<- function(z,r,auto=TRUE){
  
  z <- z[1:(length(z)/10)]
  rad <- r[1:(length(r)/10)]
  temp <- data.frame(cbind(rad,z))
  we<-1/(1+rad)^3
  # fit non-linear model
  d_mod1 <- NA  
  
  if (sum(temp$z) | 0) {
    
  if (auto)
  {
    zac <- (z-min(z))/max(z)
    temp <- data.frame(cbind(rad,zac))
    mod1 <- glm(zac ~ rad - 1, family = quasipoisson(link = "log"), data = temp, weights = we)
    c_mod1 <-coef(mod1)
    d_mod1 <- abs(1/c_mod1)

  }
  else
  {
    zcc <- 1 - z/max(z)
    temp <- data.frame(cbind(rad,zcc))
    mod1 <- glm(zcc ~ rad - 1, family = quasipoisson(link = "log"), data = temp, weights = we)
    c_mod1 <-coef(mod1)
    d_mod1 <- abs(1/c_mod1)

  }
  
  }

  
  return(d_mod1)
}

### case 3: gaussian
fit_correlogram_glm2 <- function(z,r,auto=TRUE){
  
  z <- z[1:(length(z)/10)]
  rad <- r[1:(length(r)/10)]
  tempprev <- data.frame(cbind(rad,z))

  we<-1/(1+rad)^3

  
  if (auto)
  {
    zac <- (z-min(z))/max(z)
    rad <- rad^2
    temp <- data.frame(cbind(rad,zac))
    mod1 <- glm(zac ~ rad - 1, family = quasipoisson(link = "log"), data = temp, weights = we)
    c_mod1 <-coef(mod1)
    d_mod1 <- sqrt(2/c_mod1)

  }
  else
  {
    zcc <- 1 - z/max(z)
    rad <- rad^2
    temp <- data.frame(cbind(rad,zcc))
    mod1 <- glm(zcc ~ rad - 1, family = quasipoisson(link = "log"), data = temp, weights = we)
    c_mod1 <-coef(mod1)
    d_mod1 <- sqrt(2/c_mod1)

  }
  
  
  
  
}





