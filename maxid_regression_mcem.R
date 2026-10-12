# Regression with an uncentered max-id field and separate log-Laplace noise.
# Base R only. Run: Rscript maxid_regression_mcem.R
# B_it=sum_k w_ik^(1/alpha) Z_kt; S_it=B_it.
# E_it ~ Laplace(location=0, scale=alpha0); R_it=exp(E_it)*S_it.
# L_it = x_it'beta + c * log S_it + c * E_it.
# Equivalently L=mu+c*log B+d*U, d=c*alpha0, U~Laplace(0,1).
# Thus E[L_it] = x_it'beta + c*E_theta[log B_it].
# Simulation defaults alpha0=0.1, c=0.12; examples hold both scales fixed.
# Beta describes a regression location on the log response scale, not its marginal mean.
# alpha=1/2 is fixed; f(z)=exp(-1/(4z))/(2*sqrt(pi)*z^(3/2)).
# expPS density: f(z)*exp(-theta*z+sqrt(theta)); zeros are allowed.
# MCEM approximates marginal MLE; it is not exact finite-sample integration.
# Inference uses base R. Default terrain downloads require elevatr; supplied
# elevation values avoid that dependency and permit reproducible offline runs.

maxid_model_spec <- function() "uncentered_field_separate_Laplace_scale_v4"
wendland1 <- function(r) pmax(1-r,0)^4 * (1+4*r)
# W holds the spatial basis weights; A holds coefficients in the latent sum.
# Normalize W (if desired) BEFORE raising it to 1/alpha. Do not renormalize A.
power_weights <- function(W,alpha=0.5) {
  W <- as.matrix(W)
  if (length(alpha)!=1L || !is.finite(alpha) || alpha!=0.5)
    stop("This implementation fixes alpha at 1/2")
  if (!is.numeric(W) || !nrow(W) || !ncol(W) || any(!is.finite(W)) ||
      any(W<0) || any(rowSums(W)<=0)) stop("Invalid nonnegative basis weights")
  A <- W^(1/alpha)
  if (any(!is.finite(A)) || any(rowSums(A)<=0)) stop("Powered weights overflow or underflow")
  A
}
make_design <- function(n=500L, nt=400L, seed=2026L, elevation=NULL) {
  set.seed(seed)
  # Regular 25 x 20 lattice within a 10 x 10 spatial domain.
  if(n==500L) xy <- as.matrix(expand.grid(x=seq(0,10,length.out=25),
                                         y=seq(0,10,length.out=20)))
  else xy <- cbind(runif(n,0,10),runif(n,0,10))
  knots <- as.matrix(expand.grid(x=c(1,5,9),y=c(1,5,9)))
  W <- vapply(seq_len(9),function(k)
    wendland1(sqrt(rowSums(sweep(xy,2,knots[k,],'-')^2))/6),numeric(n))
  stopifnot(all(rowSums(W)>0))
  W <- W/rowSums(W) # fixed normalization, identically used in fitting
  if (is.null(elevation)) {
    if (!requireNamespace("elevatr",quietly=TRUE))
      stop("Install the elevatr package to download terrain, or supply elevation")
    pts <- data.frame(x=xy[,1]/2-115,y=xy[,2]/2+30)
    elevation <- elevatr::get_elev_point(pts,prj=4326,src="aws",z=6)$elevation
  }
  if (length(elevation)!=n || any(!is.finite(elevation)) || sd(elevation)==0)
    stop("elevation must be finite and vary across locations")
  elev <- as.numeric(scale(elevation))
  ghg <- as.numeric(scale(seq(0,1,length.out=nt)))
  list(xy=xy,knots=knots,radii=rep(6,9),W=W,elev=elev,ghg=ghg,n=n,nt=nt)
}
# Data-only scale initialization. The simulator's dat$cscale is never a fit initializer.
.maxid_initial_scale <- function(L,X,beta,alpha0,c_init=NULL) {
  if (!is.null(c_init)) {
    if (length(c_init)!=1L || !is.finite(c_init) || c_init<=0 ||
        !is.finite(c_init*alpha0)) stop("c_init must be a positive finite scalar")
    return(as.numeric(c_init))
  }
  residual <- as.vector(L)-as.numeric(X%*%beta)
  d <- median(abs(residual))/log(2)
  if (!is.finite(d) || d<=0) d <- mean(abs(residual))
  if (!is.finite(d) || d<=0) stop("No positive residual scale; supply nondegenerate data")
  d/alpha0
}
# Laplace-transform identity for m_B=E log B, B=sum A_k Z_k, A=W^(1/alpha).
# LT_B=exp(sum(sqrt(theta_k)-sqrt(theta_k+u*A_k))).
# Transform u=exp(v); deterministic trapezoidal quadrature.
mean_log_sum <- function(W,theta,ng=601L,alpha=0.5,vlim=c(-36,30)) {
  A <- power_weights(W,alpha)
  if (length(theta)!=ncol(A) || any(!is.finite(theta)) || any(theta<0) ||
      length(ng)!=1 || !is.finite(ng) || ng<3 || ng!=as.integer(ng) ||
      length(vlim)!=2 || any(!is.finite(vlim)) || vlim[1]>=vlim[2])
    stop("Invalid tilt or quadrature controls")
  v <- seq(vlim[1],vlim[2],length.out=ng); u <- exp(v)
  a <- matrix(0,nrow(A),ng)
  for(k in seq_len(ncol(A))) a <- a + sqrt(theta[k])-
    sqrt(theta[k]+outer(A[,k],u))
  f <- matrix(exp(-u),nrow(A),ng,byrow=TRUE)-exp(a)
  as.numeric((rowSums(f)-0.5*(f[,1]+f[,ng]))*(v[2]-v[1]))
}
# S=B: its population log mean does not depend on the noise scale alpha0.
# This is an optional diagnostic; neither simulation nor fitting subtracts it.
# Retain alpha0 as a compatibility argument; it is validated, not multiplied in.
mean_log_field <- function(W,theta,ng=601L,alpha=0.5,alpha0=0.1,vlim=c(-36,30)) {
  if (length(alpha0)!=1 || !is.finite(alpha0) || alpha0<=0) stop("Invalid alpha0 scale")
  mean_log_sum(W,theta,ng,alpha,vlim)
}
rps_half <- function(n) 1/(2*rnorm(n)^2)
# Exact rejection sampling, acceptance probability exp(-sqrt(theta)).
rexpPS_half <- function(n,theta) {
  if(theta==0) return(rps_half(n))
  out <- numeric(0)
  while(length(out)<n) {
    z <- rps_half(max(32L,ceiling((n-length(out))*exp(sqrt(theta))*1.1)))
    out <- c(out,z[log(runif(length(z))) < -theta*z])
  }
  out[seq_len(n)]
}
rlaplace <- function(n,b=1) {
  u <- runif(n)-0.5
  -b*sign(u)*log1p(-2*abs(u))
}
simulate_model <- function(d=make_design(),
 theta=c(0,0.0001,0,0.0006,0.005,0.0015,0,0.002,0.003),
 beta=c(intercept=1,elevation=0.5,GHG=0.3),cscale=0.12,alpha0=0.1,seed=43) {
  K <- ncol(d$W)
  if (length(theta)!=K || any(!is.finite(theta)) || any(theta<0) ||
      length(beta)!=3 || any(!is.finite(beta)) || length(cscale)!=1 ||
      !is.finite(cscale) || cscale<=0 || length(alpha0)!=1 ||
      !is.finite(alpha0) || alpha0<=0) stop("Invalid model parameters")
  set.seed(seed)
  Z <- t(vapply(theta,function(th) rexpPS_half(d$nt,th),numeric(d$nt)))
  A <- power_weights(d$W)
  mu <- outer(beta[1]+beta[2]*d$elev,rep(1,d$nt))+
    outer(rep(1,d$n),beta[3]*d$ghg)
  B <- A%*%Z; logB <- log(B); logS <- logB
  E <- matrix(rlaplace(d$n*d$nt,alpha0),d$n)
  spatial <- cscale*logS; noise <- cscale*E
  L <- mu+spatial+noise
  c(d,list(L=L,response=exp(L),Z=Z,theta=theta,beta=beta,
           alpha=0.5,cscale=cscale,alpha0=alpha0,effective_scale=cscale*alpha0,
           A=A,B=B,S=B,R=exp(logS+E),logB=logB,logS=logS,E=E,
           center=rep(0,d$n),center_base=rep(0,d$n),mu=mu,
           spatial=spatial,noise=noise,model_spec=maxid_model_spec()))
}
# Posterior sampler in q=log Z. Include the Jacobian: log f_Z(exp(q))+q.
# Coordinate Metropolis proposals update all independent times in parallel.
e_step <- function(dat,beta,theta,q,proposal,ndraw=100L,burn=150L,thin=4L) {
  W <- power_weights(dat$W); nt <- dat$nt; K <- ncol(W); n <- dat$n
  mu <- outer(beta[1]+beta[2]*dat$elev,rep(1,nt))+
    outer(rep(1,n),beta[3]*dat$ghg)
  target <- (dat$L-mu)/dat$cscale
  B <- W %*% exp(q)
  lpobs <- -colSums(abs(target-log(B)))/dat$alpha0
  logs <- array(NA_real_,c(n,nt,ndraw)); zsum <- numeric(K)
  acc <- tries <- numeric(K)
  total <- burn+ndraw*thin; draw <- 0L
  for(it in seq_len(total)) {
    for(k in seq_len(K)) {
      qp <- q[k,]+rnorm(nt,sd=proposal[k]); zp <- exp(qp); z <- exp(q[k,])
      Bp <- B+outer(W[,k],zp-z)
      valid <- apply(Bp,2,function(x) all(is.finite(x)&x>0))
      lpnew <- rep(-Inf,nt)
      lpnew[valid] <- -colSums(abs(target[,valid,drop=FALSE]-log(Bp[,valid,drop=FALSE])))/dat$alpha0
      # Constants in prior cancel; on q scale prior=-q/2-exp(-q)/4-theta*exp(q).
      ratio <- lpnew-lpobs-0.5*(qp-q[k,])-0.25*(exp(-qp)-exp(-q[k,]))-
        theta[k]*(zp-z)
      take <- !is.na(ratio)&log(runif(nt))<ratio
      q[k,take] <- qp[take]; B[,take] <- Bp[,take,drop=FALSE]; lpobs[take] <- lpnew[take]
      if(it<=burn) {acc[k]<-acc[k]+sum(take); tries[k]<-tries[k]+nt}
    }
    # Adapt only during burn-in; retained draws use fixed proposals.
    if(it<=burn && it%%50L==0L) {
      proposal <- pmin(1.5,pmax(0.005,proposal*exp(acc/tries-0.35)))
      acc[] <- tries[] <- 0
    }
    if(it>burn && (it-burn)%%thin==0L) {
      draw <- draw+1L; logs[,,draw] <- log(B); zsum <- zsum+rowSums(exp(q))
    }
  }
  # logs contain log B=log S; alpha0 is only the Laplace noise scale.
  list(q=q,proposal=proposal,logs=logs,meanZ=zsum/(ndraw*nt))
}
fit_mcem <- function(dat,iterations=35L,draws_start=60L,draws_end=400L,
                     seed=91L,verbose=TRUE,c_init=NULL) {
  if (!identical(dat$model_spec,maxid_model_spec()))
    stop("Regenerate data with the revised simulator; incompatible/missing model_spec")
  alpha <- dat[["alpha",exact=TRUE]]
  if (!is.null(alpha) && (length(alpha)!=1 || !is.finite(alpha) || alpha!=0.5))
    stop("This fitter fixes alpha at 1/2")
  if (length(dat$alpha0)!=1L || !is.finite(dat$alpha0) || dat$alpha0<=0)
    stop("alpha0 must be a positive fixed scalar")
  if (!is.null(dat$radii) && any(dat$radii!=6)) stop("Support radii must be fixed at 6")
  set.seed(seed); K <- ncol(dat$W)
  # Initialization uses no simulation truths. OLS is only a starting value.
  X <- cbind(1,rep(dat$elev,dat$nt),rep(dat$ghg,each=dat$n))
  beta <- as.numeric(qr.solve(X,as.vector(dat$L)))
  cscale <- .maxid_initial_scale(dat$L,X,beta,dat$alpha0,c_init)
  eta <- rep(sqrt(0.002),K) # eta=sqrt(theta): permits exact boundary zero
  q <- matrix(log(1),K,dat$nt); proposal <- rep(0.15,K)
  history <- matrix(NA_real_,iterations,3+K+4,
    dimnames=list(NULL,c('intercept','elevation','GHG',paste0('theta',1:K),
                        'cscale','effective_scale','Qgain','ndraw')))
  for(it in seq_len(iterations)) {
    nd <- round(draws_start*(draws_end/draws_start)^((it-1)/max(1,iterations-1)))
    dat$cscale <- cscale
    e <- e_step(dat,beta,eta^2,q,proposal,ndraw=nd,
                burn=if(it==1) 1000L else 150L)
    q <- e$q; proposal <- e$proposal
    # M-step uses expected complete-data log likelihood. f_PS term is fixed.
    # The uncentered observation location is independent of theta given Z.
    objective <- function(p) {
      b <- p[1:3]; et <- p[3+seq_len(K)]; th <- et^2
      dscale <- exp(p[4+K])*dat$alpha0
      if (!is.finite(dscale) || dscale<=0) return(1e100)
      mu <- outer(b[1]+b[2]*dat$elev,rep(1,dat$nt))+
        outer(rep(1,dat$n),b[3]*dat$ghg)
      target <- (dat$L-mu)/dscale
      if (any(!is.finite(target))) return(1e100)
      loss <- 0
      for(j in seq_len(nd)) loss <- loss+sum(abs(target-e$logs[,,j]/dat$alpha0))/nd
      # expPS normalization is exp(sqrt(theta)); do not omit it.
      log(dscale)+(loss+dat$nt*sum(th*e$meanZ-et))/(dat$n*dat$nt)
    }
    old <- c(beta,eta,log(cscale)); qold <- objective(old)
    opt <- optim(old,objective,method='L-BFGS-B',
      lower=c(rep(-Inf,3),rep(0,K),-Inf),upper=rep(Inf,4+K),
      control=list(maxit=100,factr=1e7,ndeps=rep(1e-5,4+K)))
    # Generalized EM: require improvement of THIS sampled Q function.
    if(is.finite(opt$value)&&opt$value<=qold) {
      beta<-opt$par[1:3];eta<-opt$par[3+seq_len(K)];cscale<-exp(opt$par[4+K])
    }
    history[it,] <- c(beta,eta^2,cscale,cscale*dat$alpha0,
                     qold-objective(c(beta,eta,log(cscale))),nd)
    if(verbose) cat(sprintf('MCEM %02d/%02d: beta=(%.3f, %.3f, %.3f), c=%.5f, Q gain=%.3g, draws=%d\n',
                           it,iterations,beta[1],beta[2],beta[3],cscale,history[it,'Qgain'],nd))
  }
  list(beta=setNames(beta,c('intercept','elevation','GHG')),theta=eta^2,
       cscale=cscale,history=history,q=q,proposal=proposal,
       model_spec=maxid_model_spec(),effective_scale=cscale*dat$alpha0,
       alpha=0.5,alpha0=dat$alpha0,
       method='Monte Carlo generalized EM with uncentered field and separate Laplace noise scale',
       note='Assess mixing and repeat with more draws; fixed iteration count is not a convergence certificate.')
}
# Repeat independent posterior chains at the final fit to assess Monte Carlo error.
posterior_check <- function(dat,fit,ndraw=500L,seeds=c(101,202)) {
  if (!identical(dat$model_spec,maxid_model_spec()) ||
      !identical(fit$model_spec,maxid_model_spec())) stop("Incompatible data or fitted model")
  if (!identical(dat$alpha0,fit$alpha0)) stop("alpha0 must match the fitted fixed value")
  dat$cscale <- fit$cscale
  lapply(seeds,function(seed) {
    set.seed(seed)
    e <- e_step(dat,fit$beta,fit$theta,fit$q,fit$proposal,ndraw,burn=500L,thin=8L)
    list(meanZ=e$meanZ,proposal=e$proposal)
  })
}
run_demo <- function(outdir='maxid_results',iterations=35L) {
  dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
  if (!exists("fit_mcem_fast",mode="function")) source("maxid_regression_fast_fit.R")
  dat <- simulate_model()
  fit <- fit_mcem_fast(dat,iterations=iterations,estimate_c=FALSE,c_fixed=0.12)
  tab <- data.frame(parameter=c(names(dat$beta),paste0('theta',seq_along(dat$theta)),'c'),
                    truth=c(dat$beta,dat$theta,dat$cscale),estimate=c(fit$beta,fit$theta,fit$cscale))
  print(tab,digits=5)
  # Quantify scale balance on the modeled (log) response scale.
  cat('SD regression location:',sd(as.vector(dat$mu)),
      '; SD spatial term:',sd(as.vector(dat$spatial)),
      '; SD noise:',sd(as.vector(dat$noise)), '\n')
  cat('No population centering is used in this model.\n')
  write.csv(tab,file.path(outdir,'estimates.csv'),row.names=FALSE)
  write.csv(fit$history,file.path(outdir,'mcem_history.csv'),row.names=FALSE)
  saveRDS(list(data=dat,fit=fit),file.path(outdir,'simulation_and_fit.rds'))
  pdf(file.path(outdir,'diagnostics.pdf'),width=10,height=7)
  par(mfrow=c(2,2)); plot(dat$xy,col=hcl.colors(50)[cut(dat$elev,50)],pch=16,
    main='Fixed elevation covariate',xlab='x',ylab='y'); points(dat$knots,pch=4,cex=2)
  plot(seq_along(dat$ghg),dat$ghg,type='l',lwd=2,
       main='Increasing GHG covariate',xlab='Time replicate',ylab='Standardized GHG')
  matplot(fit$history[,1:3],type='l',lty=1,xlab='MCEM iteration',ylab='beta',
          ylim=range(c(as.matrix(fit$history[,1:3]),dat$beta)))
  abline(h=dat$beta,lty=2,col=1:3)
  legend('topright',names(dat$beta),col=1:3,lty=1,bty='n',cex=.8)
  matplot(fit$history[,4:12],type='l',lty=1,xlab='MCEM iteration',ylab='theta',
          ylim=range(c(as.matrix(fit$history[,4:12]),dat$theta)))
  abline(h=dat$theta,lty=2,col=1:9)
  legend('topright',paste0('theta',1:9),col=1:9,lty=1,ncol=3,bty='n',cex=.7)
  dev.off(); invisible(list(data=dat,fit=fit,comparison=tab))
}
if(sys.nframe()==0L) run_demo()
