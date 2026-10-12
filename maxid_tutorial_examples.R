# Spatial configuration helpers for maxid_regression_tutorial.html/.Rmd.
# First source maxid_regression_mcem.R and maxid_regression_fast_fit.R.
# This file defines functions; it does not simulate or fit automatically.
# Current observations use L=mu+c*log B+c*E, without population centering.
# elevatr is needed only when elevation is downloaded; supplied elevations
# allow the same spatial design and MCEM workflow to run without network access.

anisotropic_metric <- function(angle=30, axis_scales=c(1.8,0.8)) {
  if (length(axis_scales)!=2 || any(!is.finite(axis_scales)) || any(axis_scales<=0))
    stop("axis_scales must contain two positive values")
  a <- angle*pi/180
  R <- matrix(c(cos(a),sin(a),-sin(a),cos(a)),2,2)
  diag(1/axis_scales,2)%*%t(R)
}

make_spatial_design <- function(layout=c("grid","uniform","clustered"),
                                nx=25L,ny=20L,n=500L,nt=400L,
                                knot_nx=3L,knot_ny=3L,radii=6,
                                metric=diag(2),normalize=TRUE,
                                xlim=c(0,10),ylim=c(0,10),
                                locations=NULL,knots=NULL,elevation=NULL,seed=2026L) {
  layout <- match.arg(layout)
  if (length(xlim)!=2 || length(ylim)!=2 || any(!is.finite(c(xlim,ylim))) ||
      diff(xlim)<=0 || diff(ylim)<=0) stop("Invalid domain bounds")
  for (value in c(nx,ny,n,nt,knot_nx,knot_ny))
    if (!is.finite(value) || value<1 || value!=as.integer(value)) stop("Counts must be positive integers")
  if (nt<2) stop("At least two times are needed for an increasing GHG covariate")
  metric <- as.matrix(metric)
  if (!identical(dim(metric),c(2L,2L)) || any(!is.finite(metric)) || abs(det(metric))<1e-12)
    stop("metric must be a nonsingular 2 by 2 matrix")
  set.seed(seed)
  if (!is.null(locations)) xy <- as.matrix(locations)
  else if (layout=="grid") {
    xy <- as.matrix(expand.grid(x=seq(xlim[1],xlim[2],length.out=nx),
                               y=seq(ylim[1],ylim[2],length.out=ny)))
  } else if (layout=="uniform") {
    xy <- cbind(x=runif(n,xlim[1],xlim[2]),y=runif(n,ylim[1],ylim[2]))
  } else {
    # Four Gaussian clusters, rejecting points outside the domain.
    centers <- as.matrix(expand.grid(x=xlim[1]+c(.25,.75)*diff(xlim),
                                    y=ylim[1]+c(.25,.75)*diff(ylim)))
    xy <- matrix(numeric(),0,2)
    while (nrow(xy)<n) {
      m <- max(32L,n-nrow(xy)); group <- sample.int(4,m,replace=TRUE)
      proposal <- centers[group,,drop=FALSE]+
        cbind(rnorm(m,sd=.065*diff(xlim)),rnorm(m,sd=.065*diff(ylim)))
      keep <- proposal[,1]>=xlim[1] & proposal[,1]<=xlim[2] &
              proposal[,2]>=ylim[1] & proposal[,2]<=ylim[2]
      xy <- rbind(xy,proposal[keep,,drop=FALSE])
    }
    xy <- xy[seq_len(n),,drop=FALSE]
  }
  if (ncol(xy)!=2 || !nrow(xy) || any(!is.finite(xy))) stop("locations must be a finite two-column matrix")
  if (is.null(knots)) {
    knots <- as.matrix(expand.grid(
      x=seq(xlim[1]+.1*diff(xlim),xlim[2]-.1*diff(xlim),length.out=knot_nx),
      y=seq(ylim[1]+.1*diff(ylim),ylim[2]-.1*diff(ylim),length.out=knot_ny)))
  } else knots <- as.matrix(knots)
  if (ncol(knots)!=2 || !nrow(knots) || any(!is.finite(knots))) stop("Invalid knot coordinates")
  n <- nrow(xy); K <- nrow(knots)
  if (length(radii)==1) radii <- rep(radii,K)
  if (length(radii)!=K || any(!is.finite(radii)) || any(radii!=6))
    stop("Support radii are fixed at 6 for every knot")
  W <- vapply(seq_len(K),function(k) {
    delta <- sweep(xy,2,knots[k,],"-")%*%t(metric)
    r <- sqrt(rowSums(delta^2))/radii[k]
    pmax(1-r,0)^4*(1+4*r)
  },numeric(n))
  uncovered <- which(rowSums(W)<=0)
  unused <- which(colSums(W)<=0)
  if (length(uncovered)) stop(length(uncovered)," locations have no positive basis weight; add or reposition knots")
  if (length(unused)) stop("Unused knots: ",paste(unused,collapse=", "),"; remove or reposition them")
  if (normalize) W <- W/rowSums(W)
  if (is.null(elevation)) {
    if (!requireNamespace("elevatr",quietly=TRUE))
      stop("Install the elevatr package to download terrain, or supply elevation")
    pts  <- data.frame(x = xy[, 1]/2-115, y = xy[, 2]/2+30)       # lon, lat
    out  <- elevatr::get_elev_point(pts, prj = 4326, src = "aws", z = 6)
    elevation <- out$elevation                                     # meters, same row order as xy
  }
  if (length(elevation)!=n || any(!is.finite(elevation)) || sd(elevation)==0)
    stop("elevation must vary across locations")
  list(xy=xy,knots=knots,radii=radii,metric=metric,W=W,
       elev=as.numeric(scale(elevation)),ghg=as.numeric(scale(seq(0,1,length.out=nt))),
       n=n,nt=nt,layout=layout,xlim=xlim,ylim=ylim,normalized=normalize)
}

theta_pattern <- function(K) rep(c(0,0.0001,0,0.0006,0.005,0.0015,0,0.002,0.003),length.out=K)

simulate_spatial <- function(design=make_spatial_design(),theta=NULL,
                            beta=c(intercept=1,elevation=0.5,GHG=0.3),
                            cscale=0.12,alpha0=0.1,seed=43L) {
  if (!exists("simulate_model",mode="function"))
    stop("Source maxid_regression_mcem.R before calling simulate_spatial()")
  K <- ncol(design$W); n <- design$n; nt <- design$nt
  if (is.null(theta)) theta <- theta_pattern(K)
  if (length(theta)!=K || any(!is.finite(theta)) || any(theta<0) ||
      length(beta)!=3 || any(!is.finite(beta))) stop("Invalid theta/beta")
  if (!is.finite(cscale) || cscale<=0 || !is.finite(alpha0) || alpha0<=0)
    stop("Invalid noise scales")
  # Delegate to one simulator so tutorial examples and fitted likelihood agree.
  simulate_model(design,theta=theta,beta=beta,cscale=cscale,alpha0=alpha0,seed=seed)
}

plot_spatial_design <- function(d,main="Spatial design",supports=TRUE) {
  plot(d$xy,pch=16,cex=.55,col="#227C85",asp=1,
       xlim=d$xlim,ylim=d$ylim,xlab="x",ylab="y",main=main)
  if (supports) {
    a <- seq(0,2*pi,length.out=150)
    circle <- cbind(cos(a),sin(a))%*%t(solve(d$metric))
    for (k in seq_len(nrow(d$knots)))
      lines(sweep(d$radii[k]*circle,2,d$knots[k,],"+"),col="#ADB8C2")
  }
  points(d$knots,pch=4,col="#BA4D2F",cex=1.2,lwd=2)
}

example_spatial_designs <- function() {
  list(grid=make_spatial_design(),
       uniform=make_spatial_design("uniform"),
       clustered=make_spatial_design("clustered"),
       sparse_grid=make_spatial_design(nx=10,ny=10),
       dense_grid=make_spatial_design(nx=25,ny=20),
       dense_knots=make_spatial_design(knot_nx=4,knot_ny=4),
       anisotropic=make_spatial_design(metric=anisotropic_metric()))
}
