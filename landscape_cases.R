source('~/Desktop/IIASA/fit_correlograms_glm.R', echo=TRUE)
source('~/Desktop/IIASA/make_correlogram_20160731.R', echo=TRUE)
source('~/Desktop/IIASA/ex_prepare_inputs.R', echo=TRUE)

require(raster)
require(secr)
require(igraph)

tempmask <- make.mask(nx = 400, ny = 400, spacing = 10)
# p: controls fragmentation, A: expected proportion of habitat
prepreimage0 <- raster(randomHabitat(tempmask, p = 0.6, A = 1, minpatch = 10))
preimage0 <- as.matrix(prepreimage0)
preimage0[is.na(preimage0)] <- 2
image0 <- prepare_input(preimage0)
output0 <- input_2LCT(image0)

M01 <- output0$M1
M02 <- output0$M2
M0F <- output0$MF

dfim0 <- make_correlogram(M01,M02,M0F,return.vectors = TRUE)
z0 <- dfim0$z
r0 <- dfim0$r
L0 <- fit_correlogram_glm1(z0,r0,auto = FALSE) 

dfim1 <- make_correlogram(M01,M01,M0F,return.vectors = TRUE)
z1 <- dfim1$z
r1 <- dfim1$r
L1 <- fit_correlogram_glm1(z1,r1,auto = TRUE) # fit autocorrelograms (same LC)

dfim2 <- make_correlogram(M02,M02,M0F,return.vectors = TRUE)
z2 <- dfim2$z
r2 <- dfim2$r
L2 <- fit_correlogram_glm1(z2,r2,auto = TRUE) # fit autocorrelograms (same LC)

par(mfrow=c(2,2))
image(image0)
plot(r0,z0)
plot(r1,z1)
plot(r2,z2)
