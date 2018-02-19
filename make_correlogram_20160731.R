### locations of n equally spaced points in [0,L] (aligned left,right, or middle)
loc0 <- function(n,L=1) (0:(n-1))*L/n
loc1 <- function(n,L=1) (1:n)*L/n
loc05 <- function(n,L=1) (1:n-0.5)*L/n

### bilinear interpolation
# returns values at (x,y) by inter- and extrapolating on matrix Z
# length(x)==length(y), 0<=x<=1, 0<=y<=1, dim(Z)==c(n,m)
bilinear_interpolation <- function(x,y,Z)
{
  stopifnot(length(x)==length(y))
  stopifnot(all(x>=0 & y>=0),all(x<=1 & y<=1))
  stopifnot(nrow(Z)>1 & ncol(Z)>1)

  m <- nrow(Z)
  n <- ncol(Z)

  i <- pmax(pmin(floor(m*x+0.5),m-1),1)
  j <- pmax(pmin(floor(n*y+0.5),n-1),1)

  zij <- diag(Z[i,j,drop=FALSE])
  zip1 <- diag(Z[i+1,j,drop=FALSE])
  zjp1 <- diag(Z[i,j+1,drop=FALSE])
  zijp1 <- diag(Z[i+1,j+1,drop=FALSE])

  t <- (x-loc05(m)[i])*m
  u <- (y-loc05(n)[j])*n

  return((1-t)*(1-u)*zij+(1-t)*u*zjp1+t*(1-u)*zip1+t*u*zijp1)
}

### returns resampled matrix, by default resized to nearest power of 2
resample <- function(A,
                     m=2^round(log2(nrow(A))),n=2^round(log2(ncol(A))),
                     intp=bilinear_interpolation)
{
  if (all(dim(A)==c(m,n)))
    # nothing to do
    return(A)
  else
    return(outer(loc05(m),loc05(n),intp,A))
}

shiftpeak <- function(A)
{
  m <- nrow(A)
  n <- ncol(A)
  rm <- 1:(m/2)
  rn <- 1:(n/2)
  return(rbind(cbind(A[rm+m/2,rn+n/2],A[rm+m/2,rn]),
               cbind(A[rm,rn+n/2],A[rm,rn])))
}

### correlogram
mkcf <- function(A,B=A,normalize=FALSE)
{
  stopifnot(all(dim(A)==dim(B)),is.logical(normalize))

  message("FFT'ing A ...")
  fft_A <- fft(A)

  if (identical(A,B))
    fft_B <- fft_A
  else
    fft_B <- fft(B)

  # 2D correlation function
  Z <- abs(fft(fft_A*Conj(fft_B),inverse=TRUE))
  if (normalize)
    Z <- Z/(length(A)*sqrt(length(Z)))
  return(Z)
}

### whole procedure
make_correlogram <- function(A,B=A,C=NULL,
                             px_width=1,px_height=px_width,
                             symmetrize.cross=TRUE,
                             keep.all=FALSE,
                             return.vectors=FALSE)
{
  # check arguments
  if (all(dim(A)!=2^round(log2(dim(A)))))
    stop("Please supply only maps with height and widths already powers of two")
  if (is.null(C))
    stop("Please always supply a municipality mask of the same size as the map(s).")
  stopifnot(all(dim(B)==dim(A)),all(dim(C)==dim(A)))

  # auto- or cross-correlation?
  do.auto <- identical(A,B)

  # remember map widths and heigths
  m <- nrow(A)
  n <- ncol(A)

  # will need this
  N <- matrix(1,m,n)

  if (do.auto)
    Z <- mkcf(A,C*A)/mkcf(N,C)
  else
  {
    if (symmetrize.cross)
    {
      A1 <- C*A
      B1 <- C*B
      Z <- (mkcf(A,B1)+mkcf(B,A1)-mkcf(A1,B1))/mkcf(N,C)
    }
    else
      Z <- mkcf(B,C*A)/mkcf(N,C)
  }

  message("Shifting peak to middle ...")
  rm <- 1:(m/2)
  rn <- 1:(n/2)
  Z <- rbind(cbind(Z[rm+m/2,rn+n/2],Z[rm+m/2,rn]),
             cbind(Z[rm,rn+n/2],Z[rm,rn]))

  # half height and width of map
  hh <- m*px_height/2
  hw <- n*px_width/2

  # corresponding radii
  R00 <- sqrt(outer(loc0(m/2,L=hh)^2,loc0(n/2,L=hw)^2,"+"))
  R01 <- sqrt(outer(loc0(m/2,L=hh)^2,loc1(n/2,L=hw)^2,"+"))
  R10 <- sqrt(outer(loc1(m/2,L=hh)^2,loc0(n/2,L=hw)^2,"+"))
  R11 <- sqrt(outer(loc1(m/2,L=hh)^2,loc1(n/2,L=hw)^2,"+"))

  # arrange so that radius matrix corresponds to shifted-peak Z matrix
  R <- rbind(cbind(R11[(m/2):1,(n/2):1],R10[(m/2):1,]),
             cbind(R01[,(n/2):1],R00))

  # for convenience, store distances along rows and columns also separately
  # the peak is now at 1+m/2,1+n/2
  rdist <- R[,1+n/2]
  rdist[n/2+rn] <- -rdist[n/2+rn]
  cdist <- R[1+m/2,]
  cdist[rm] <- -cdist[rm]

  stopifnot(all(dim(R)==dim(Z)))

  if (!keep.all)
  {
    # only keep half of original map dimensions
    m_ <- m/2
    n_ <- n/2
    # the peak is now at 1+m/2,1+n/2
    # so want to grab
    mr1 <- 1:(m_/2)+m_
    nr1 <- 1:(n_/2)+n_
    mr2 <- 1:(m_/2)+m_/2
    nr2 <- 1:(n_/2)+n_/2

    Z <- Z[c(mr2,mr1),c(nr2,nr1)]
    R <- R[c(mr2,mr1),c(nr2,nr1)]

    rdist <- rdist[c(mr2,mr1)]
    cdist <- cdist[c(nr2,nr1)]
  }

  if (return.vectors)
  {
    # sort vectors by radius
    s_ <- sort(as.vector(R),index.return=TRUE)
    return(data.frame(r=s_$x,z=as.vector(Z)[s_$ix]))
  }
  else
    return(list(Z=Z,R=R,rdist=rdist,cdist=cdist))
}
