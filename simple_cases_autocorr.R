source('~/Desktop/IIASA/fit_correlograms_glm.R', echo=TRUE)
source('~/Desktop/IIASA/make_correlogram_20160731.R', echo=TRUE)
source('~/Desktop/IIASA/ex_prepare_inputs.R', echo=TRUE)


preimage0<- array(1,c(12,12))
image0 <- prepare_input(preimage0)
dfim0 <- make_correlogram(image0,image0,image0,return.vectors = TRUE)
z0 <- dfim0$z
r0 <- dfim0$r
L0 <- fit_correlogram_glm1(z0,r0,auto = TRUE) # fit autocorrelograms (same LC)

preimage1<- array(1,c(25,25))
image1 <- prepare_input(preimage1)
dfim1 <- make_correlogram(image1,image1,image1,return.vectors = TRUE)
z1 <- dfim1$z
r1 <- dfim1$r
L1 <- fit_correlogram_glm1(z1,r1,auto = TRUE) # fit autocorrelograms (same LC)

preimage2<- array(1,c(50,50))
image2 <- prepare_input(preimage2)
dfim2 <- make_correlogram(image2,image2,image2,return.vectors = TRUE)
z2 <- dfim2$z
r2 <- dfim2$r
L2 <- fit_correlogram_glm1(z2,r2,auto = TRUE) # fit autocorrelograms (same LC)

preimage3<- array(1,c(75,75))
image3 <- prepare_input(preimage3)
dfim3 <- make_correlogram(image3,image3,image3,return.vectors = TRUE)
z3 <- dfim3$z
r3 <- dfim3$r
L3 <- fit_correlogram_glm1(z3,r3,auto = TRUE) # fit autocorrelograms (same LC)

preimage4<- array(1,c(100,100))
image4 <- prepare_input(preimage4)
dfim4 <- make_correlogram(image4,image4,image4,return.vectors = TRUE)
z4 <- dfim4$z
r4 <- dfim4$r
L4 <- fit_correlogram_glm1(z4,r4,auto = TRUE) # fit autocorrelograms (same LC)

preimage5<- array(1,c(200,200))
image5 <- prepare_input(preimage5)
dfim5 <- make_correlogram(image5,image5,image5,return.vectors = TRUE)
z5 <- dfim5$z
r5 <- dfim5$r
L5 <- fit_correlogram_glm1(z5,r5,auto = TRUE) # fit autocorrelograms (same LC)

par(mfrow=c(2,3))
image(image0)
image(image1)
image(image2)
plot(r0,z0)

plot(r1,z1)
plot(r2,z2)

par(mfrow=c(2,3))
image(image3)
image(image4)
image(image5)
plot(r3,z3)
plot(r4,z4)
plot(r5,z5)
