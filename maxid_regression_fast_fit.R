# Faster Monte Carlo ECM for the SAME uncentered model as maxid_regression_mcem.R.
# Base R only. Source the original file for make_design()/simulate_model().
# No automatic fitting or compilation happens when this file is sourced.
# fit_mcem_fast(dat) fits dat$L, not dat$response; alpha remains fixed at 1/2.
# Revised hierarchy: A=W^(1/alpha), B=A%*%Z, S=B,
# E~Laplace(0,scale=alpha0). Work with noise scale d=cscale*alpha0.
# Conditional location is mu+c*log B; no population centering is evaluated.
# Features:
#   * separate convex LAD, closed-form tilt, and optional positive-scale updates;
#   * equivalent location augmentation accelerates fixed-c fits;
#   * optional mean-log-B quadrature cache retained for diagnostics only;
#   * sorted Monte Carlo offsets for fast absolute-loss evaluation;
#   * joint log-latent proposals with a covariance for each temporal replicate;
#   * adaptation ONLY during warm-up, plus occasional coordinate moves.

.maxid_check <- function(dat) {
  required_spec <- "uncentered_field_separate_Laplace_scale_v4"
  if (!exists("maxid_model_spec",mode="function") ||
      !identical(maxid_model_spec(),required_spec))
    stop("Source the updated uncentered maxid_regression_mcem.R before this fitter")
  needed <- c("L", "W", "elev", "ghg", "alpha0", "model_spec")
  if (!all(needed %in% names(dat))) stop("dat must contain: ", paste(needed, collapse=", "))
  if (!identical(dat$model_spec,maxid_model_spec()))
    stop("Incompatible model_spec: regenerate data with the revised simulator")
  alpha <- dat[["alpha",exact=TRUE]] # $alpha would partially match alpha0 in R
  if (!is.null(alpha) && (length(alpha)!=1L ||
      !is.finite(alpha) || alpha!=0.5)) stop("This fitter fixes alpha at 1/2")
  L <- as.matrix(dat$L); W <- as.matrix(dat$W)
  if (!is.numeric(L) || !is.numeric(W) || any(!is.finite(L)) ||
      any(!is.finite(W)) || any(W < 0)) stop("L and nonnegative W must be finite numeric matrices")
  n <- nrow(L); nt <- ncol(L); K <- ncol(W)
  if (n < 1 || nt < 1 || K < 1 || nrow(W) != n ||
      any(rowSums(W) <= 0) || any(colSums(W) <= 0)) stop("Invalid/unused basis functions")
  if (length(dat$elev) != n || length(dat$ghg) != nt ||
      any(!is.finite(dat$elev)) || any(!is.finite(dat$ghg))) stop("Invalid covariates")
  if (length(dat$alpha0) != 1 || !is.finite(dat$alpha0) || dat$alpha0 <= 0)
    stop("alpha0 must be a positive fixed scalar")
  if (!is.null(dat$radii) && (length(dat$radii)!=K ||
      any(!is.finite(dat$radii)) || any(dat$radii!=6)))
    stop("Support radii must be fixed at 6")
  dat$L <- L; dat$W <- W; dat$n <- n; dat$nt <- nt; dat$alpha <- 0.5
  dat$A <- power_weights(W)
  dat
}

# L-BFGS-B may evaluate a coordinate a few ulps below its lower bound.
# Treat only roundoff-sized negative eta as zero, leaving positive values exact.
# The optimizer objective, gradient, cache, and accepted state use this same rule.
.maxid_eta <- function(eta,K) {
  if (!is.numeric(eta) || is.complex(eta) || length(eta)!=K)
    stop("Invalid sqrt(theta): expected ",K," numeric components; received ",length(eta))
  bad <- which(!is.finite(eta))
  if (length(bad)) stop("Invalid sqrt(theta): non-finite components at indices ",
                        paste(bad,collapse=", "))
  tolerance <- 64*.Machine$double.eps*pmax(1,abs(eta))
  bad <- which(eta < -tolerance)
  if (length(bad)) stop("Invalid sqrt(theta): negative components beyond roundoff at indices ",
    paste(bad,collapse=", ")," (values ",paste(sprintf("%.17g",eta[bad]),collapse=", "),")")
  pmax(as.numeric(eta),0)
}

# On eta=sqrt(theta), use u*A=u*W^(1/alpha), not u*W.
# Zero basis weights have ZERO derivative, including at eta=0.
# Derivatives below are for the same finite-range quadrature as the objective.
make_center_cache <- function(W, ng=601L, vlim=c(-36, 30), capacity=16L,
                              alpha=0.5, alpha0=0.1) {
  W <- power_weights(W,alpha); K <- ncol(W); n <- nrow(W)
  if (length(alpha0)!=1 || !is.finite(alpha0) || alpha0<=0) stop("Invalid alpha0 scale")
  if (ng < 3 || capacity < 1 || length(vlim) != 2 || vlim[1] >= vlim[2])
    stop("Invalid centering quadrature/cache controls")
  v <- seq(vlim[1], vlim[2], length.out=ng); u <- exp(v)
  quadrature <- rep(v[2]-v[1], ng); quadrature[c(1,ng)] <- quadrature[c(1,ng)]/2
  Wu <- lapply(seq_len(K), function(k) outer(W[,k],u))
  state <- new.env(parent=emptyenv()); state$keys <- character(); state$values <- list()
  state$hits <- 0L; state$misses <- 0L
  evaluate <- function(eta) {
    eta <- .maxid_eta(eta,K)
    key <- paste(sprintf("%.17g",eta),collapse=",") # no rounded/approximate memoization
    at <- match(key,state$keys)
    if (!is.na(at)) { state$hits <- state$hits+1L; return(state$values[[at]]) }
    state$misses <- state$misses+1L
    roots <- lapply(seq_len(K),function(k) sqrt(eta[k]^2+Wu[[k]]))
    a <- matrix(0,n,ng)
    for (k in seq_len(K)) a <- a+eta[k]-roots[[k]]
    B <- sweep(exp(a),2,quadrature,"*")
    # Direct difference avoids cancellation between two separately integrated constants.
    m <- rowSums(sweep(matrix(exp(-u),n,ng,byrow=TRUE)-exp(a),2,quadrature,"*"))
    dm <- matrix(0,n,K)
    for (k in seq_len(K)) {
      active <- which(W[,k]>0)
      if (length(active)) dm[active,k] <- -rowSums(
        B[active,,drop=FALSE]*(1-eta[k]/roots[[k]][active,,drop=FALSE]))
    }
    # S=B, so both centering views coincide and are independent of alpha0.
    value <- list(mean=m,gradient=dm,base_mean=m,base_gradient=dm)
    if (length(state$keys) >= capacity) {
      state$keys <- state$keys[-1L]; state$values <- state$values[-1L]
    }
    state$keys <- c(state$keys,key); state$values[[length(state$values)+1L]] <- value
    value
  }
  list(evaluate=evaluate,
       stats=function() c(hits=state$hits,misses=state$misses),
       ng=ng,vlim=vlim,alpha=alpha,alpha0=alpha0)
}

# One row per site/time, one column per posterior draw. Sort only ONCE per E-step.
.maxid_offset_cache <- function(logs,scale=1,row_shift=NULL,block_size=512L) {
  dm <- dim(logs)
  if (length(dm)!=3 || any(dm<1)) stop("Invalid latent log-field draws")
  N <- dm[1]*dm[2]; M <- dm[3]
  if (length(scale)!=1L || !is.finite(scale) || scale<=0 ||
      length(block_size)!=1L || !is.finite(block_size) || block_size<1 ||
      block_size!=as.integer(block_size)) stop("Invalid offset cache controls")
  if (is.null(row_shift)) row_shift <- numeric(N)
  if (length(row_shift)!=N || any(!is.finite(row_shift))) stop("Invalid row shifts")
  column_offsets <- (seq_len(M)-1)*N
  sorted <- matrix(0,N,M)
  # Avoid apply() collecting every sorted row and a full-size transpose.
  # Shift/standardize small blocks instead of copying the whole log-field array.
  for (first in seq.int(1L,N,by=block_size)) {
    rows <- seq.int(first,min(N,first+block_size-1L))
    block <- matrix(logs[as.vector(outer(rows,column_offsets,"+"))],length(rows),M)
    block <- sweep(block,1,row_shift[rows],"+")/scale
    if (any(!is.finite(block))) stop("Invalid latent log-field draws")
    sorted[rows,] <- t(vapply(seq_len(nrow(block)),function(j)sort(block[j,]),numeric(M)))
  }
  cumulative <- sorted
  if (M>1) for (j in 2:M) cumulative[,j] <- cumulative[,j]+cumulative[,j-1L]
  list(sorted=sorted,cumulative=cumulative,total=cumulative[,M],N=N,M=M,
       row=seq_len(N))
}

# E_draw |target-logB| and a subgradient, using row-wise binary searches.
.maxid_abs_stats <- function(offsets,target,value=TRUE) {
  N <- offsets$N; M <- offsets$M
  if (length(target)!=N || any(!is.finite(target))) stop("Invalid loss target")
  lo <- integer(N); hi <- rep.int(M,N)
  repeat {
    active <- which(lo<hi)
    if (!length(active)) break
    mid <- (lo[active]+hi[active]+1L)%/%2L
    yes <- offsets$sorted[active+(mid-1L)*N] <= target[active]
    lo[active[yes]] <- mid[yes]; hi[active[!yes]] <- mid[!yes]-1L
  }
  gradient <- 2*lo/M-1
  if (!value) return(list(gradient=gradient))
  prefix <- numeric(N); positive <- which(lo>0)
  prefix[positive] <- offsets$cumulative[positive+(lo[positive]-1L)*N]
  value <- target*gradient+(offsets$total-2*prefix)/M
  # Tiny negative values may arise from floating-point cancellation.
  list(value=pmax(value,0),gradient=gradient)
}

# At fixed beta/theta, the scale objective is convex in lambda=1/(c*alpha0):
# -log(lambda)+mean_draw,observation |lambda*(L-X beta)-log B/alpha0|.
# center is zero for the uncentered model; the helper accepts an explicit offset.
# The -log(lambda) term is the conditional Laplace normalizer, not optional.
# Bracket its monotone subgradient, then bisect; reuse the sorted offsets.
.maxid_scale_update <- function(offsets,residual,center,cscale,alpha0,
                                 tol=1e-8,maxit=80L) {
  if (length(residual)!=offsets$N || length(center)!=offsets$N ||
      any(!is.finite(c(residual,center))) || !any(residual!=0))
    stop("Scale update needs finite, nonzero regression residuals")
  lambda0 <- 1/(cscale*alpha0)
  evaluate <- function(lambda) {
    stats <- .maxid_abs_stats(offsets,lambda*residual+center)
    list(loss=-log(lambda)+mean(stats$value),
         gradient=-1/lambda+mean(residual*stats$gradient))
  }
  initial <- evaluate(lambda0); lo <- hi <- lambda0
  if (initial$gradient>0) {
    for (j in seq_len(60L)) {
      lo <- lo/2
      if (evaluate(lo)$gradient<=0) break
    }
  } else if (initial$gradient<0) {
    for (j in seq_len(60L)) {
      hi <- hi*2
      if (evaluate(hi)$gradient>=0) break
    }
  }
  if (evaluate(lo)$gradient>0 || evaluate(hi)$gradient<0)
    stop("Could not bracket a finite positive conditional scale estimate")
  iterations <- 0L
  for (j in seq_len(maxit)) {
    if (hi-lo<=tol*(lo+hi)/2) break
    mid <- (lo+hi)/2
    if (evaluate(mid)$gradient<0) lo <- mid else hi <- mid
    iterations <- j
  }
  candidates <- c(lambda0,lo,(lo+hi)/2,hi)
  losses <- vapply(candidates,function(x)evaluate(x)$loss,numeric(1))
  at <- which.min(losses)
  list(cscale=1/(alpha0*candidates[at]),iterations=iterations,
       loss_decrease=initial$loss-losses[at])
}

# Proximal map of mean_draw |z-(base-d*logB_draw/alpha0)|. Offsets are decreasing
# on the pseudo-response scale. Solve the piecewise-linear derivative exactly.
.maxid_lad_prox <- function(offsets,base,cscale,point,lambda) {
  N <- offsets$N; M <- offsets$M
  lo <- rep.int(1L,N); hi <- rep.int(M+1L,N)
  repeat {
    active <- which(lo<hi)
    if (!length(active)) break
    mid <- (lo[active]+hi[active])%/%2L
    knot <- base[active]-cscale*offsets$sorted[active+(M-mid)*N]
    yes <- knot+lambda*(2*mid/M-1) >= point[active]
    hi[active[yes]] <- mid[yes]; lo[active[!yes]] <- mid[!yes]+1L
  }
  k <- lo; z <- point+lambda*(1-2*(k-1L)/M)
  right <- which(k<=M)
  z[right] <- pmin(z[right],base[right]-cscale*offsets$sorted[right+(M-k[right])*N])
  left <- which(k>1L)
  z[left] <- pmax(z[left],base[left]-cscale*offsets$sorted[left+(M-k[left]+1L)*N])
  z
}

# Convex LAD regression, without constructing N*M repeated design rows.
# ADMM updates the fitted values with exact proximal maps and projects onto X.
# Returning convergence diagnostics is essential: a finite maxit need not solve LAD exactly.
.maxid_lad_update <- function(X, XtX_chol, offsets, L, center, cscale, beta,
                             maxit=200L, abstol=1e-5, reltol=1e-5,
                             rho=5, relaxation=1.6, state=NULL) {
  N <- nrow(X); p <- ncol(X)
  # The fixed-c location augmentation has only one LAD coefficient.
  # Solve its convex objective directly instead of iterating ADMM.
  if (p==1L) return(.maxid_lad_scalar(as.numeric(X),offsets,L,center,
                                     cscale,beta,maxit=maxit))
  base <- L+cscale*center
  fitted <- as.numeric(X%*%beta); z <- fitted; dual <- numeric(N)
  if (!is.null(state)) {
    if (length(state$z)!=N || length(state$dual)!=N || any(!is.finite(state$z)) ||
        any(!is.finite(state$dual)) || !is.finite(state$rho) || state$rho<=0)
      stop("Invalid LAD warm-start state")
    z <- state$z; dual <- state$dual; rho <- state$rho
  }
  project <- function(v) as.numeric(backsolve(XtX_chol,
    forwardsolve(t(XtX_chol),crossprod(X,v))))
  loss <- function(b) sum(.maxid_abs_stats(offsets,(L-as.numeric(X%*%b))/cscale+center)$value)
  initial <- loss(beta); best <- beta; best_loss <- initial
  converged <- FALSE; primal <- dual_res <- Inf
  for (it in seq_len(maxit)) {
    beta <- project(z-dual); fitted <- as.numeric(X%*%beta); previous <- z
    relaxed <- relaxation*fitted+(1-relaxation)*previous
    z <- .maxid_lad_prox(offsets,base,cscale,relaxed+dual,1/rho)
    dual <- dual+relaxed-z
    # RMS feasibility and stationarity for the PER-OBSERVATION LAD objective.
    # Using an unnormalized summed gradient makes the tolerance depend on sample size.
    primal <- sqrt(mean((fitted-z)^2))
    dual_res <- rho*sqrt(sum(crossprod(X,z-previous)^2))/N
    eps_primal <- abstol+reltol*max(sqrt(mean(fitted^2)),sqrt(mean(z^2)))
    eps_dual <- sqrt(p)*abstol+reltol*sqrt(sum(crossprod(X,rho*dual)^2))/N
    if (it%%10L==0L || it==1L || it==maxit ||
        (primal<=eps_primal && dual_res<=eps_dual)) {
      current <- loss(beta)
      if (is.finite(current) && current<best_loss) {best<-beta;best_loss<-current}
    }
    if (primal<=eps_primal && dual_res<=eps_dual) {converged<-TRUE;break}
    # Residual balancing; rescale the scaled dual to preserve the unscaled multiplier.
    if (it%%25L==0L) {
      s <- rho*sqrt(mean((z-previous)^2)); next_rho <- rho
      if (primal>10*s) next_rho <- min(1e4,2*rho)
      if (s>10*primal) next_rho <- max(1e-4,rho/2)
      dual <- dual*(rho/next_rho); rho <- next_rho
    }
  }
  list(beta=best,converged=converged,iterations=it,
       loss_decrease=initial-best_loss,primal=primal,dual=dual_res,
       state=list(z=z,dual=dual,rho=rho))
}

.maxid_lad_scalar <- function(x,offsets,L,center,cscale,beta,maxit=200L,tol=1e-8) {
  use <- which(x!=0)
  if (!length(use)) stop("The scalar LAD covariate must vary")
  base <- L+cscale*center
  a <- (base[use]-cscale*offsets$sorted[use,1])/x[use]
  b <- (base[use]-cscale*offsets$sorted[use,offsets$M])/x[use]
  lo <- min(a,b); hi <- max(a,b)
  evaluate <- function(coef,value=TRUE) {
    ans <- .maxid_abs_stats(offsets,(L-x*coef)/cscale+center,value=value)
    list(loss=if(value)sum(ans$value) else NA_real_,gradient=-sum(x*ans$gradient)/cscale)
  }
  initial <- evaluate(beta)$loss
  for (it in seq_len(maxit)) {
    if (hi-lo<=tol*(1+abs(lo)+abs(hi))) break
    mid <- (lo+hi)/2
    if (evaluate(mid,FALSE)$gradient<0) lo <- mid else hi <- mid
  }
  candidates <- c(beta,lo,(lo+hi)/2,hi)
  losses <- vapply(candidates,function(v)evaluate(v)$loss,numeric(1))
  take <- which.min(losses)
  list(beta=candidates[take],converged=hi-lo<=tol*(1+abs(lo)+abs(hi)),
       iterations=it,loss_decrease=initial-losses[take],
       primal=hi-lo,dual=NA_real_,state=NULL)
}

.maxid_spd <- function(A) {
  E <- eigen((A+t(A))/2,symmetric=TRUE)
  values <- pmax(E$values,max(1e-8,max(E$values)*1e-6))
  tcrossprod(sweep(E$vectors,2,sqrt(values),"*"))
}

# Laplace location Fisher information + expPS prior curvature provides a
# covariance reflecting overlap between basis functions, even before adaptation.
.maxid_joint_metric <- function(W,q,theta,alpha0=1) {
  K <- nrow(q); nt <- ncol(q); metric <- array(0,c(K,K,nt))
  for (t in seq_len(nt)) {
    z <- exp(q[,t]); B <- as.numeric(W%*%z)
    D <- sweep(sweep(W,2,z,"*"),1,B,"/")
    # Conditional location information for Laplace scale alpha0 is 1/alpha0^2.
    H <- crossprod(D)/alpha0^2+diag(0.25*exp(-q[,t])+theta*z+1e-8,K)
    metric[,,t] <- .maxid_spd(chol2inv(chol(H)))
  }
  metric
}

# Includes the Jacobian for q=log Z. Theta-dependent prior constants cancel in MH.
.maxid_joint_logposterior <- function(q,W,target,theta,alpha0=1) {
  nt <- ncol(q); result <- rep(-Inf,nt); logB <- matrix(NA_real_,nrow(W),nt)
  valid <- colSums(!is.finite(q)|q>700|q< -700)==0L
  if (any(valid)) {
    z <- exp(q[,valid,drop=FALSE]); B <- W%*%z
    ok <- colSums(!is.finite(B)|B<=0)==0L
    times <- which(valid)[ok]
    if (length(times)) {
      logB[,times] <- log(B[,ok,drop=FALSE])
      result[times] <- -colSums(abs(target[,times,drop=FALSE]-logB[,times,drop=FALSE]))/alpha0-
        0.5*colSums(q[,times,drop=FALSE])-0.25*colSums(exp(-q[,times,drop=FALSE]))-
        as.numeric(crossprod(theta,z[,ok,drop=FALSE]))
    }
  }
  list(value=result,logB=logB)
}

# A posterior-mode search is ONLY an initialization; E-steps still sample/integrate Z.
.maxid_initial_q <- function(W,target,theta,maxit=60L,alpha0=1) {
  K <- ncol(W); nt <- ncol(target); q <- matrix(0,K,nt)
  for (t in seq_len(nt)) {
    evaluate <- function(x,gradient=FALSE) {
      z <- exp(x); B <- as.numeric(W%*%z); residual <- log(B)-target[,t]
      if (!gradient) return(sum(abs(residual))/alpha0+0.5*sum(x)+
        0.25*sum(exp(-x))+sum(theta*z))
      D <- sweep(sweep(W,2,z,"*"),1,B,"/")
      as.numeric(crossprod(D,sign(residual)))/alpha0+0.5-0.25*exp(-x)+theta*z
    }
    opt <- optim(rep(0,K),function(x)evaluate(x),function(x)evaluate(x,TRUE),
      method="L-BFGS-B",lower=rep(-50,K),upper=rep(50,K),control=list(maxit=maxit))
    if (is.finite(opt$value)) q[,t] <- opt$par
  }
  q
}

.maxid_ess <- function(x) {
  M <- length(x)
  if (M<8L || !is.finite(var(x)) || var(x)==0) return(NA_real_)
  a <- as.numeric(acf(x,lag.max=min(100L,M-1L),plot=FALSE)$acf)[-1L]
  pairs <- a[seq.int(1L,length(a)-1L,2L)]+a[seq.int(2L,length(a),2L)]
  bad <- which(pairs<=0)
  if (length(bad)) pairs <- if(bad[1]==1L) numeric() else pairs[seq_len(bad[1]-1L)]
  min(M,M/max(1,1+2*sum(pairs)))
}

# Returned logs contain unscaled log B. center is retained only for call compatibility.
# A nonzero center is rejected, preventing accidental use of the previous likelihood.
e_step_joint <- function(dat,beta,theta,q,joint_state=NULL,center=NULL,
                         ndraw=100L,burn=200L,thin=1L,adapt_every=50L,
                         coordinate_every=5L) {
  W <- power_weights(dat$W); n <- dat$n; nt <- dat$nt; K <- ncol(W)
  if (!is.null(center) && (length(center)!=n || any(!is.finite(center)) || any(center!=0)))
    stop("The uncentered model requires a zero center")
  target <- (dat$L-outer(beta[1]+beta[2]*dat$elev,rep(1,nt))-
    outer(rep(1,n),beta[3]*dat$ghg))/dat$cscale
  reference <- .maxid_joint_metric(W,q,theta,dat$alpha0)
  if (is.null(joint_state)) {
    covariance <- reference; scales <- rep(2.38/sqrt(K),nt)
  } else {
    covariance <- joint_state$covariance; scales <- joint_state$scales
    if (!identical(dim(covariance),c(K,K,nt)) || length(scales)!=nt ||
        any(!is.finite(covariance)) || any(!is.finite(scales)) || any(scales<=0))
      stop("Invalid joint proposal state")
    for (t in seq_len(nt)) covariance[,,t] <- .maxid_spd(
      0.8*covariance[,,t]+0.2*reference[,,t])
  }
  cholesky <- lapply(seq_len(nt),function(t)t(chol(covariance[,,t])))
  current <- .maxid_joint_logposterior(q,W,target,theta,dat$alpha0)
  if (any(!is.finite(current$value))) stop("Nonfinite initial latent posterior")
  qdraws <- array(0,c(K,nt,ndraw)); logs <- array(0,c(n,nt,ndraw))
  warm_window <- array(0,c(K,nt,adapt_every)); window_count <- 0L
  warm_accept <- numeric(nt); retained_accept <- numeric(nt); retained_moves <- 0L
  zsum <- numeric(K); draw <- 0L
  for (it in seq_len(burn+ndraw*thin)) {
    proposed <- q
    for (t in seq_len(nt)) proposed[,t] <- q[,t]+scales[t]*as.numeric(cholesky[[t]]%*%rnorm(K))
    candidate <- .maxid_joint_logposterior(proposed,W,target,theta,dat$alpha0)
    accept <- log(runif(nt)) < candidate$value-current$value
    accept[is.na(accept)] <- FALSE
    q[,accept] <- proposed[,accept,drop=FALSE]
    current$value[accept] <- candidate$value[accept]
    current$logB[,accept] <- candidate$logB[,accept,drop=FALSE]
    if (it<=burn) warm_accept <- warm_accept+accept
    else {retained_accept<-retained_accept+accept;retained_moves<-retained_moves+1L}
    # Coordinate moves supplement the joint kernel; both kernels leave the same target invariant.
    if (coordinate_every>0L && it%%coordinate_every==0L) {
      for (k in seq_len(K)) {
        proposed <- q
        sd <- vapply(seq_len(nt),function(t)sqrt(covariance[k,k,t]),numeric(1))
        proposed[k,] <- proposed[k,]+rnorm(nt,sd=sd)
        candidate <- .maxid_joint_logposterior(proposed,W,target,theta,dat$alpha0)
        accept <- log(runif(nt)) < candidate$value-current$value
        accept[is.na(accept)] <- FALSE
        q[,accept] <- proposed[,accept,drop=FALSE]
        current$value[accept] <- candidate$value[accept]
        current$logB[,accept] <- candidate$logB[,accept,drop=FALSE]
      }
    }
    if (it<=burn) {
      window_count <- window_count+1L; warm_window[,,window_count] <- q
      if (window_count==adapt_every) {
        scales <- pmax(0.05,pmin(5,scales*exp(warm_accept/window_count-0.234)))
        if (window_count>max(K+2L,10L)) {
          for (t in seq_len(nt)) {
            empirical <- cov(t(matrix(warm_window[,t,],K,adapt_every)))
            covariance[,,t] <- .maxid_spd(0.5*empirical+0.5*reference[,,t])
            cholesky[[t]] <- t(chol(covariance[,,t]))
          }
        }
        warm_accept[] <- 0; window_count <- 0L
      }
    } else if ((it-burn)%%thin==0L) {
      draw <- draw+1L; qdraws[,,draw] <- q; logs[,,draw] <- current$logB
      zsum <- zsum+rowSums(exp(q))
    }
  }
  ess <- matrix(apply(matrix(qdraws,K*nt,ndraw),1,.maxid_ess),K,nt)
  list(q=q,logs=logs,qdraws=qdraws,meanZ=zsum/(ndraw*nt),ess=ess,
       acceptance=retained_accept/retained_moves,
       joint_state=list(covariance=covariance,scales=scales))
}

# At fixed c, use V_kt=exp(h_t) Z_kt, h_t=(beta0+beta_GHG*g_t)/c.
# The observation location becomes beta_elev*e_i+c*log(sum A_ik V_kt).
# On r=log V, prior log density is
# -.5*(r-h)-.25*exp(h-r)-theta*exp(r-h)+sqrt(theta)+constant.
# These moments are computed at the OLD parameters, so optimize a shift delta h.
# Profile eta and the intercept shift analytically, leaving a one-dimensional
# GHG-shift optimization. Log-sum-exp avoids extreme exponentials.
.maxid_location_update <- function(meanZ_time,meanInvZ_time,ghg,beta,cscale,
                                    eta,maxit=60L) {
  K <- nrow(meanZ_time); T <- ncol(meanZ_time)
  if (!identical(dim(meanInvZ_time),c(K,T)) || length(ghg)!=T ||
      any(!is.finite(c(meanZ_time,meanInvZ_time))) ||
      any(meanZ_time<=0) || any(meanInvZ_time<=0)) stop("Invalid latent moments")
  lse <- function(x) { m<-max(x); m+log(sum(exp(x-m))) }
  logZ <- log(meanZ_time); logU <- log(colSums(meanInvZ_time))
  profile <- function(gamma) {
    H <- sweep(logZ,2,gamma*ghg,"-")
    logH <- apply(H,1,lse)
    Hweights <- exp(H-logH)
    dlogH <- -as.numeric(Hweights%*%ghg)
    v <- logU+gamma*ghg; logI <- lse(v)
    dlogI <- sum(exp(v-logI)*ghg)
    reciprocal <- 2*log(T)-logH; logD <- lse(reciprocal)
    dlogD <- sum(exp(reciprocal-logD)*(-dlogH))
    difference <- logD-logI
    if (!is.finite(difference) || difference>=0)
      stop("Degenerate latent scale profile: nonpositive C")
    ratio <- exp(difference)
    logC <- logI-log(4)+log(-expm1(difference))
    derivative <- (dlogI-ratio*dlogD)/(1-ratio)
    a <- log(K*T/2)-logC
    new_eta <- exp(a+log(T/2)-logH)
    list(loss=K*T*logC/2-K*gamma*sum(ghg)/2,
         gradient=K*T*derivative/2-K*sum(ghg)/2,
         a=a,eta=new_eta)
  }
  opt <- optim(0,function(g)profile(g)$loss,function(g)profile(g)$gradient,
               method="BFGS",control=list(maxit=maxit,reltol=1e-11))
  candidate <- profile(opt$par)
  shift <- candidate$a+opt$par*ghg
  prior_loss <- function(et,delta) {
    .25*sum(sweep(meanInvZ_time,2,exp(delta),"*"))+
      sum(et^2*rowSums(sweep(meanZ_time,2,exp(-delta),"*")))-
      K*sum(delta)/2-T*sum(et)
  }
  before <- prior_loss(eta,rep(0,T))
  after <- prior_loss(candidate$eta,shift)
  if (!is.finite(after) || after>before+1e-8) {
    # An optimizer failure must never be accepted as an EM update.
    return(list(beta=beta,eta=eta,shift=rep(0,T),code=max(1L,opt$convergence),
                loss_decrease=0))
  }
  out <- beta
  out[c(1,3)] <- beta[c(1,3)]+cscale*c(candidate$a,opt$par)
  list(beta=out,eta=candidate$eta,shift=shift,code=opt$convergence,
       loss_decrease=before-after)
}

# Exact leave-one-replicate-out jackknife bias correction for theta=eta^2.
# eta enters both M-steps through a profiled reciprocal-type transform of
# time-averaged posterior Z moments; since expPS(1/2,theta) is heavy-tailed
# (Var(Z) ~ theta^-1.5), this transform is convex enough that Jensen's
# inequality gives E[theta_hat]>theta_true, worse for smaller true theta.
# Verified against true (noiseless) simulated Z: naive relative bias of
# 0.5-25% at T=600, dropping to under ~3.4% after this correction, with no
# variance inflation. beta is held fixed throughout (it already shows no
# meaningful bias at convergence); only eta/theta is refit per leave-one-out
# fold. Costs T+1 calls to .maxid_location_update (location augmentation) or
# a closed-form recomputation (classical augmentation) -- a few seconds to
# under a minute at T~600, negligible next to a 35-iteration MCEM fit.
.maxid_jackknife_theta_locaug <- function(meanZ_time,meanInvZ_time,ghg,beta,cscale,
                                          eta,maxit=60L) {
  K <- nrow(meanZ_time); T <- ncol(meanZ_time)
  full <- .maxid_location_update(meanZ_time,meanInvZ_time,ghg,beta,cscale,eta,maxit)
  theta_full <- full$eta^2
  theta_loo_sum <- numeric(K)
  for (t in seq_len(T)) {
    loo <- .maxid_location_update(meanZ_time[,-t,drop=FALSE],meanInvZ_time[,-t,drop=FALSE],
                                  ghg[-t],beta,cscale,eta,maxit)
    theta_loo_sum <- theta_loo_sum+loo$eta^2
  }
  theta_loo_mean <- theta_loo_sum/T
  pmax(T*theta_full-(T-1)*theta_loo_mean,0)
}

# Classical Z-augmentation counterpart: eta=1/(2*mean(Z)) has a closed-form
# leave-one-out mean, so no refit loop is needed.
.maxid_jackknife_theta_classical <- function(meanZ_time) {
  K <- nrow(meanZ_time); T <- ncol(meanZ_time)
  zbar <- rowMeans(meanZ_time)
  theta_full <- 1/(4*zbar^2)
  zbar_loo <- sweep(-meanZ_time,1,T*zbar,"+")/(T-1)
  theta_loo_mean <- rowMeans(1/(4*zbar_loo^2))
  pmax(T*theta_full-(T-1)*theta_loo_mean,0)
}

# Estimate c by default; hold it fixed with estimate_c=FALSE and c_fixed.
# Fixed-c fits default to the equivalent location augmentation above.
# Set location_augmentation=FALSE for the classical Z augmentation.
# center_ng is a compatibility control; fitting does not evaluate centering.
fit_mcem_fast <- function(dat,iterations=35L,draws_start=80L,draws_end=600L,
                          burn_first=800L,burn=150L,thin=1L,m_cycles=1L,
                          lad_maxit=200L,lad_abstol=1e-5,lad_reltol=1e-5,
                          theta_maxit=60L,center_ng=601L,
                          start=NULL,seed=91L,verbose=TRUE,c_init=NULL,
                          estimate_c=TRUE,c_fixed=NULL,location_augmentation=NULL,
                          uncertainty=c("none","louis"),uncertainty_level=.95,
                          uncertainty_bandwidth=.1,uncertainty_final=TRUE,
                          theta_jackknife=TRUE) {
  dat <- .maxid_check(dat); n <- dat$n; nt <- dat$nt; K <- ncol(dat$W); N <- n*nt
  uncertainty <- match.arg(uncertainty)
  if (uncertainty=="louis" && (length(uncertainty_level)!=1L ||
      !is.finite(uncertainty_level) || uncertainty_level<=0 || uncertainty_level>=1 ||
      length(uncertainty_bandwidth)!=1L || !is.finite(uncertainty_bandwidth) ||
      uncertainty_bandwidth<=0 || !is.logical(uncertainty_final) ||
      length(uncertainty_final)!=1L || is.na(uncertainty_final)))
    stop("Invalid uncertainty controls")
  if (!is.logical(theta_jackknife) || length(theta_jackknife)!=1L || is.na(theta_jackknife))
    stop("theta_jackknife must be TRUE or FALSE")
  if (!is.null(start) && !identical(start$model_spec,maxid_model_spec()))
    stop("Cannot resume a fit from the previous model design")
  if (!is.null(start) && (is.null(start$cscale) || is.null(start$alpha0) ||
      start$alpha0!=dat$alpha0)) stop("Resume needs cscale and the same fixed alpha0")
  if (missing(estimate_c) && !is.null(start$estimate_c)) estimate_c <- start$estimate_c
  if (!is.logical(estimate_c) || length(estimate_c)!=1L || is.na(estimate_c))
    stop("estimate_c must be TRUE or FALSE")
  if (estimate_c && !is.null(c_fixed))
    stop("Use estimate_c=FALSE with c_fixed; use c_init to initialize an estimated c")
  if (uncertainty=="louis" && estimate_c)
    stop("Iteration-level Louis diagnostics currently require fixed c")
  if (uncertainty=="louis" && min(draws_start,draws_end)<25L)
    stop("Louis diagnostics need at least 25 retained draws per E-step")
  if (!estimate_c) {
    if (!is.null(c_init)) stop("Use c_fixed, not c_init, when estimate_c=FALSE")
    if (is.null(c_fixed) && !is.null(start)) c_fixed <- start$cscale
    if (!is.numeric(c_fixed) || length(c_fixed)!=1L || !is.finite(c_fixed) || c_fixed<=0)
      stop("Fixed-c fitting needs a positive finite scalar c_fixed (or start$cscale)")
    c_fixed <- as.numeric(c_fixed)
  }
  if (is.null(location_augmentation)) {
    location_augmentation <- if(!is.null(start$location_augmentation))
      start$location_augmentation else !estimate_c
  }
  if (!is.logical(location_augmentation) || length(location_augmentation)!=1L ||
      is.na(location_augmentation)) stop("location_augmentation must be TRUE or FALSE")
  if (location_augmentation && estimate_c)
    stop("Location augmentation requires fixed c; set location_augmentation=FALSE to estimate c")
  integer_controls <- c(iterations,draws_start,draws_end,thin,m_cycles,lad_maxit,theta_maxit)
  if (any(!is.finite(integer_controls)) || any(integer_controls<1) ||
      any(integer_controls!=as.integer(integer_controls)) ||
      any(!is.finite(c(burn_first,burn))) || burn_first<0 || burn<0 ||
      burn_first!=as.integer(burn_first) || burn!=as.integer(burn) ||
      !is.finite(lad_abstol) || lad_abstol<=0 || !is.finite(lad_reltol) || lad_reltol<=0)
    stop("Invalid iteration/draw controls")
  set.seed(seed)
  X <- cbind(intercept=1,elevation=rep(dat$elev,nt),GHG=rep(dat$ghg,each=n))
  if (qr(X)$rank<ncol(X)) stop("Regression design is rank deficient")
  L <- as.vector(dat$L)
  X_lad <- if(location_augmentation) X[,2,drop=FALSE] else X
  XtX_chol <- chol(crossprod(X_lad)); zero <- numeric(N)
  beta <- if(is.null(start$beta)) as.numeric(qr.solve(X,L)) else as.numeric(start$beta)
  theta <- if(is.null(start$theta)) rep(0.002,K) else as.numeric(start$theta)
  if (length(beta)!=3L || any(!is.finite(beta)) || length(theta)!=K ||
      any(!is.finite(theta)) || any(theta<0)) stop("Invalid starting parameters")
  cscale <- if(estimate_c) .maxid_initial_scale(L,X,beta,dat$alpha0,
                            if(is.null(start)) c_init else start$cscale) else c_fixed
  dscale <- cscale*dat$alpha0
  if (!is.finite(dscale) || dscale<=0) stop("c * alpha0 must be positive and finite")
  eta <- sqrt(theta)
  q <- start$q; joint_state <- start$joint_state; lad_state <- start$lad_state
  if (!is.null(start) && !identical(start$location_augmentation,location_augmentation))
    lad_state <- NULL
  if (is.null(q)) {
    target <- matrix((L-as.numeric(X%*%beta))/cscale,n,nt)
    q <- .maxid_initial_q(dat$A,target,theta,alpha0=dat$alpha0)
  }
  if (!identical(dim(q),c(K,nt)) || any(!is.finite(q))) stop("Invalid starting log-latents")
  rows <- vector("list",iterations); last_posterior <- NULL
  uncertainty_rows <- list(); uncertainty_details <- list()
  for (it in seq_len(iterations)) {
    # Release the previous sampled fields/cache before allocating the next
    # E-step. This matters for 500 sites with several thousand retained draws.
    e <- NULL; offsets <- NULL
    gc(verbose=FALSE)
    nd <- round(draws_start*(draws_end/draws_start)^((it-1)/max(1,iterations-1)))
    begin <- proc.time()[["elapsed"]]
    dat$cscale <- cscale
    e <- e_step_joint(dat,beta,eta^2,q,joint_state,ndraw=nd,
      burn=if(it==1L)burn_first else burn,thin=thin)
    q <- e$q; joint_state <- e$joint_state
    e_seconds <- proc.time()[["elapsed"]]-begin
    if (uncertainty=="louis") {
      u <- .maxid_louis_uncertainty(dat,beta,eta^2,e,it-1L,
                  level=uncertainty_level,bandwidth=uncertainty_bandwidth)
      uncertainty_rows[[it]] <- u$table
      uncertainty_details[[it]] <- u[setdiff(names(u),"table")]
    }
    begin_m <- proc.time()[["elapsed"]]
    previous_beta <- beta; previous_eta <- eta; previous_c <- cscale
    if (location_augmentation) {
      meanZ_time <- apply(exp(e$qdraws),c(1,2),mean)
      meanInvZ_time <- apply(exp(-e$qdraws),c(1,2),mean)
      h_old <- (beta[1]+beta[3]*dat$ghg)/cscale
      # The cache shifts/standardizes small blocks, without materializing V
      # or another full-size array of transformed log fields.
    }
    offsets <- .maxid_offset_cache(e$logs,scale=dat$alpha0,
                  row_shift=if(location_augmentation)rep(h_old,each=n) else NULL)
    e$logs <- NULL
    sampled_loss <- function(b,et,cs=cscale) {
      et <- .maxid_eta(et,K); d <- cs*dat$alpha0
      if (location_augmentation) {
        a <- (L-as.numeric(X_lad%*%b[2]))/d
        delta <- ((b[1]-previous_beta[1])+(b[3]-previous_beta[3])*dat$ghg)/cs
        prior <- .25*sum(sweep(meanInvZ_time,2,exp(delta),"*"))+
          sum(et^2*rowSums(sweep(meanZ_time,2,exp(-delta),"*")))-
          K*sum(delta)/2-nt*sum(et)
      } else {
        a <- (L-as.numeric(X%*%b))/d
        prior <- nt*sum(et^2*e$meanZ-et)
      }
      log(d)+(sum(.maxid_abs_stats(offsets,a)$value)+prior)/N
    }
    old_loss <- sampled_loss(beta,eta)
    lad_iterations <- 0L; lad_all_converged <- TRUE; scale_iterations <- 0L
    theta_code <- 0L
    update_beta <- function(b) {
      initial <- if(location_augmentation) b[2] else b
      out <- .maxid_lad_update(X_lad,XtX_chol,offsets,L,zero,dscale,initial,
                       lad_maxit,lad_abstol,lad_reltol,state=lad_state)
      lad_state <<- out$state
      if (location_augmentation) { full<-b;full[2]<-out$beta;out$beta<-full }
      out
    }
    update_c <- function() {
      if (!estimate_c) return(invisible(NULL))
      scale <- .maxid_scale_update(offsets,L-as.numeric(X%*%beta),zero,cscale,dat$alpha0)
      scale_iterations <<- scale_iterations+scale$iterations
      if (sampled_loss(beta,eta,scale$cscale)<=sampled_loss(beta,eta)) {
        cscale <<- scale$cscale; dscale <<- cscale*dat$alpha0
      }
    }
    if (location_augmentation) {
      loc <- .maxid_location_update(meanZ_time,meanInvZ_time,dat$ghg,
                                   beta,cscale,eta,maxit=theta_maxit)
      theta_code <- loc$code
      if (sampled_loss(loc$beta,loc$eta)<=sampled_loss(beta,eta)+1e-12) {
        beta <- loc$beta; eta <- loc$eta
      }
    } else {
      # The observation likelihood has no theta dependence conditional on Z.
      eta <- 1/(2*e$meanZ)
    }
    for (cycle in seq_len(m_cycles+1L)) {
      lad <- update_beta(beta); lad_iterations <- lad_iterations+lad$iterations
      lad_all_converged <- lad_all_converged && lad$converged
      if (sampled_loss(lad$beta,eta)<=sampled_loss(beta,eta)) beta <- lad$beta
      update_c()
      # With fixed-c location augmentation the scalar LAD solution is already
      # complete and no other conditional update changes its objective.
      if (location_augmentation && lad$converged) break
    }
    gain <- old_loss-sampled_loss(beta,eta)
    if (!is.finite(gain) || gain< -1e-10) stop("Sampled Q ascent check failed")
    if (location_augmentation) {
      # Warm-start on the SAME retained V scale after the M-step changes h.
      delta <- ((beta[1]-previous_beta[1])+(beta[3]-previous_beta[3])*dat$ghg)/cscale
      q <- sweep(q,2,delta,"-")
    }
    m_seconds <- proc.time()[["elapsed"]]-begin_m
    ess_values <- as.numeric(e$ess); ess_values <- ess_values[is.finite(ess_values)]
    rows[[it]] <- c(beta,eta^2,cscale=cscale,effective_scale=dscale,Qgain=gain,ndraw=nd,
      delta_beta=max(abs(beta-previous_beta)),delta_eta=max(abs(eta-previous_eta)),
      delta_c=abs(cscale-previous_c),scale_iterations=scale_iterations,
      joint_acceptance=mean(e$acceptance),
      ess_min=if(length(ess_values))min(ess_values) else NA_real_,
      ess_median=if(length(ess_values))median(ess_values) else NA_real_,
      lad_iterations=lad_iterations,lad_converged=as.numeric(lad_all_converged),
      theta_code=theta_code,E_seconds=e_seconds,M_seconds=m_seconds,
      center_hits=0,center_misses=0)
    names(rows[[it]])[seq_len(3+K)] <- c(colnames(X),paste0("theta",seq_len(K)))
    if (verbose) cat(sprintf(
      "MCECM %02d/%02d: beta=(%.3f, %.3f, %.3f), c=%.5f (%s), Q gain=%.3g, accept=%.2f, median ESS=%.1f, E/M=%.1f/%.1fs\n",
      it,iterations,beta[1],beta[2],beta[3],cscale,if(estimate_c)"estimated" else "fixed",
      gain,mean(e$acceptance),rows[[it]][["ess_median"]],e_seconds,m_seconds))
    last_posterior <- list(meanZ=e$meanZ,ess=e$ess,acceptance=e$acceptance)
  }
  # Standard post-fit step: an exact leave-one-replicate-out jackknife bias
  # correction for theta, reusing the last iteration's own E-step (no extra
  # MCMC draws) since it already reflects the posterior at (near-)converged
  # parameters. This only affects theta_jackknife below; theta itself stays
  # the raw EM fixed-point estimate used for resuming and for the Louis bands.
  theta_jackknife_estimate <- NULL
  if (theta_jackknife) {
    meanZ_time_final <- apply(exp(e$qdraws),c(1,2),mean)
    theta_jackknife_estimate <- if (location_augmentation) {
      meanInvZ_time_final <- apply(exp(-e$qdraws),c(1,2),mean)
      .maxid_jackknife_theta_locaug(meanZ_time_final,meanInvZ_time_final,dat$ghg,
                                    beta,cscale,eta,maxit=theta_maxit)
    } else .maxid_jackknife_theta_classical(meanZ_time_final)
    names(theta_jackknife_estimate) <- paste0("theta",seq_len(K))
  }
  uncertainty_final_seconds <- 0
  if (uncertainty=="louis" && uncertainty_final) {
    e <- NULL; offsets <- NULL; gc(verbose=FALSE)
    begin <- proc.time()[["elapsed"]]; dat$cscale <- cscale
    e <- e_step_joint(dat,beta,eta^2,q,joint_state,ndraw=nd,burn=burn,thin=thin)
    u <- .maxid_louis_uncertainty(dat,beta,eta^2,e,iterations,
                level=uncertainty_level,bandwidth=uncertainty_bandwidth)
    uncertainty_rows[[iterations+1L]] <- u$table
    uncertainty_details[[iterations+1L]] <- u[setdiff(names(u),"table")]
    uncertainty_final_seconds <- proc.time()[["elapsed"]]-begin
  }
  list(beta=setNames(beta,colnames(X)),theta=setNames(eta^2,paste0("theta",seq_len(K))),
       theta_jackknife=theta_jackknife_estimate,
       cscale=cscale,estimate_c=estimate_c,c_fixed=if(estimate_c)NULL else c_fixed,
       location_augmentation=location_augmentation,
       history=as.data.frame(do.call(rbind,rows)),q=q,joint_state=joint_state,lad_state=lad_state,
       uncertainty=if(length(uncertainty_rows))do.call(rbind,uncertainty_rows) else NULL,
       uncertainty_details=uncertainty_details,
       uncertainty_final_seconds=uncertainty_final_seconds,
       posterior=last_posterior,center_cache_stats=c(hits=0L,misses=0L),
       controls=list(iterations=iterations,draws_start=draws_start,draws_end=draws_end,
         burn_first=burn_first,burn=burn,thin=thin,m_cycles=m_cycles,center_ng=center_ng,
         c_init=c_init,estimate_c=estimate_c,c_fixed=if(estimate_c)NULL else c_fixed,
         location_augmentation=location_augmentation,uncertainty=uncertainty,
         uncertainty_level=uncertainty_level,uncertainty_bandwidth=uncertainty_bandwidth,
         uncertainty_final=uncertainty_final,theta_jackknife=theta_jackknife),
       model_spec=maxid_model_spec(),effective_scale=dscale,alpha=dat$alpha,alpha0=dat$alpha0,
       method=if(location_augmentation) "MCEM with uncentered field and location augmentation" else
              "Monte Carlo ECM with uncentered field and classical Z augmentation",
       note=paste("Qgain is an improvement of the same sampled complete-data objective.",
                  "ESS estimates are approximate; no observed-likelihood convergence is certified.",
                  "posterior diagnostics describe the E-step preceding the final M-step.",
                  "Finite latent means give positive theta updates even for true zero tilts.",
                  "theta is the raw EM fixed-point estimate (used for resuming fits and",
                  "for the Louis uncertainty bands); theta_jackknife applies an exact",
                  "leave-one-replicate-out bias correction on top of it and is the",
                  "recommended point estimate to report. It corrects the average bias of",
                  "the estimator across repeated sampling, not necessarily the gap to",
                  "truth in any one fit -- assess that via repeated-simulation recovery,",
                  "not proximity to truth on a single dataset."))
}

# Independent chains at the RETURNED parameters (after the final M-step).
posterior_check_fast <- function(dat,fit,ndraw=1000L,burn=500L,seeds=c(101L,202L)) {
  dat <- .maxid_check(dat)
  if (!identical(fit$model_spec,maxid_model_spec())) stop("Incompatible fitted model")
  if (is.null(fit$cscale) || is.null(fit$alpha0) || fit$alpha0!=dat$alpha0)
    stop("Posterior check needs the fitted cscale and the same fixed alpha0")
  dat$cscale <- fit$cscale
  lapply(seeds,function(seed) {
    set.seed(seed)
    e <- e_step_joint(dat,fit$beta,fit$theta,fit$q,fit$joint_state,
                      ndraw=ndraw,burn=burn,thin=1L)
    list(seed=seed,meanZ=e$meanZ,ess=e$ess,acceptance=e$acceptance,
         q=e$q,joint_state=e$joint_state)
  })
}

# Iteration-level observed-information diagnostic for FIXED c and alpha0.
# Use the smooth r=log V augmentation for beta0, beta_GHG and eta=sqrt(theta).
# The Laplace elevation Hessian is a delta function: approximate its posterior
# expectation by a Gaussian kernel of residuals, bandwidth* c * alpha0.
# This smoothing affects diagnostics only, never the fitted likelihood.
.maxid_louis_uncertainty <- function(dat,beta,theta,e,iteration=0L,
                                      level=.95,bandwidth=.1) {
  K <- ncol(dat$W); T <- dat$nt; n <- dat$n; M <- dim(e$qdraws)[3]
  P <- K+3L; eta <- sqrt(theta); cscale <- dat$cscale; d <- cscale*dat$alpha0
  params <- c("intercept","elevation","GHG",paste0("theta",seq_len(K)))
  root_names <- c(params[1:3],paste0("eta",seq_len(K)))
  complete <- missing <- matrix(0,P,P)
  score_draws <- matrix(0,P,M)
  kernel_width <- bandwidth*d
  for (t in seq_len(T)) {
    z <- exp(matrix(e$qdraws[,t,],K,M)); invz <- 1/z
    residual <- matrix(dat$L[,t]-beta[1]-beta[2]*dat$elev-
                       beta[3]*dat$ghg[t],n,M)-cscale*matrix(e$logs[,t,],n,M)
    hscore <- K/2-.25*colSums(invz)+as.numeric(crossprod(theta,z))
    s <- matrix(0,P,M)
    s[1,] <- hscore/cscale
    s[2,] <- as.numeric(crossprod(dat$elev,sign(residual)))/d
    s[3,] <- dat$ghg[t]*hscore/cscale
    s[3L+seq_len(K),] <- 1-2*eta*z
    score_draws <- score_draws+s
    # Times are conditionally independent. Sum each time's posterior covariance
    # to avoid estimating unnecessary cross-time score covariances.
    centered <- s-rowMeans(s)
    missing <- missing+tcrossprod(centered)/(M-1)
    mz <- rowMeans(z); mi <- rowMeans(invz)
    curvature <- (.25*sum(mi)+sum(theta*mz))/cscale^2
    v <- c(1,dat$ghg[t]); ix <- c(1L,3L); iz <- 3L+seq_len(K)
    complete[ix,ix] <- complete[ix,ix]+curvature*tcrossprod(v)
    cross <- -2*eta*mz/cscale
    complete[ix,iz] <- complete[ix,iz]+outer(v,cross)
    complete[iz,ix] <- t(complete[ix,iz])
    complete[cbind(iz,iz)] <- complete[cbind(iz,iz)]+2*mz
    density_zero <- rowMeans(matrix(dnorm(residual/kernel_width)/kernel_width,n,M))
    complete[2,2] <- complete[2,2]+2*sum(dat$elev^2*density_zero)/d
  }
  information <- (complete-missing+t(complete-missing))/2
  scaling <- sqrt(diag(complete))
  scaled <- information/outer(scaling,scaling)
  eig <- eigen(scaled,symmetric=TRUE,only.values=TRUE)$values
  ok <- all(is.finite(eig)) && min(eig)>1e-8
  covariance <- mc_covariance <- matrix(NA_real_,P,P)
  batch_size <- max(2L,floor(sqrt(M))); batches <- M%/%batch_size
  score_mc_cov <- matrix(NA_real_,P,P)
  if (batches>=5L) {
    averages <- vapply(seq_len(batches),function(b)
      rowMeans(score_draws[,(b-1L)*batch_size+seq_len(batch_size),drop=FALSE]),
      numeric(P))
    score_mc_cov <- cov(t(averages))/batches
  }
  if (ok) {
    covariance <- chol2inv(chol(scaled))/outer(scaling,scaling)
    if (all(is.finite(score_mc_cov)))
      mc_covariance <- covariance%*%score_mc_cov%*%covariance
  }
  dimnames(information) <- dimnames(covariance) <- dimnames(mc_covariance) <-
    list(root_names,root_names)
  root_estimate <- c(beta,eta)
  root_se <- sqrt(pmax(diag(covariance),0))
  root_mcse <- sqrt(pmax(diag(mc_covariance),0))
  quantile <- qnorm((1+level)/2)
  lower <- root_estimate-quantile*root_se; upper <- root_estimate+quantile*root_se
  mc_lower <- root_estimate-quantile*root_mcse
  mc_upper <- root_estimate+quantile*root_mcse
  iz <- 3L+seq_len(K)
  lower[iz] <- pmax(lower[iz],0)^2; upper[iz] <- pmax(upper[iz],0)^2
  mc_lower[iz] <- pmax(mc_lower[iz],0)^2; mc_upper[iz] <- pmax(mc_upper[iz],0)^2
  jacobian <- c(rep(1,3),2*eta)
  table <- data.frame(iteration=iteration,parameter=params,estimate=c(beta,theta),
    se=jacobian*root_se,lower=lower,upper=upper,
    mcse=jacobian*root_mcse,mc_lower=mc_lower,mc_upper=mc_upper,
    score=rowMeans(score_draws),score_mcse=sqrt(pmax(diag(score_mc_cov),0)),
    root_se=root_se,information_ok=ok,information_min_eigen=min(eig),
    bandwidth=bandwidth,ndraw=M,level=level,row.names=NULL)
  list(table=table,information=information,covariance_root=covariance,
       mc_covariance_root=mc_covariance,score_mc_covariance=score_mc_cov,
       note=paste("Approximate Louis observed information at the displayed iterate.",
         "Laplace elevation curvature uses a residual-density kernel; it is a diagnostic approximation.",
         "Before convergence these are local-curvature bands, not confidence intervals for an MLE.",
         "Root-theta Wald bands are nonregular near zero; no simultaneous or boundary coverage is claimed.",
         "mcse linearizes posterior score sampling error with batch means; it excludes curvature error",
         "and the accumulated Monte Carlo error of previous EM updates.",
         "Non-positive information returns NA, without repairing eigenvalues."))
}

# Both panels retain all coefficients together and all local tilts together.
# Statistical bands use likelihood curvature; Monte Carlo bands describe only
# the local score sampling error. Fixed c and alpha0 have no fitted uncertainty.
plot_mcem_traces <- function(fit,truth_beta=NULL,truth_theta=NULL,
                             bands=c("statistical","monte_carlo","none"),
                             main="MCEM parameter traces") {
  bands <- match.arg(bands); h <- fit$history
  band_level <- if(is.null(fit$uncertainty))95 else 100*fit$uncertainty$level[1]
  beta_names <- c("intercept","elevation","GHG")
  theta_names <- grep("^theta[0-9]+$",names(h),value=TRUE)
  groups <- list(beta_names,theta_names)
  colors <- list(c("#193c50","#107d84","#bd5737"),hcl.colors(length(theta_names),"Dark 3"))
  u <- fit$uncertainty
  if (bands!="none" && is.null(u)) stop("Fit with uncertainty='louis' to plot bands")
  oldpar <- par(mfrow=c(1,2),mar=c(4,4,3,1),oma=c(3,0,2,0))
  on.exit(par(oldpar))
  show_jackknife <- !is.null(fit$theta_jackknife)
  for (panel in 1:2) {
    keys <- groups[[panel]]; col <- colors[[panel]]
    truth <- if(panel==1)truth_beta else truth_theta
    relevant <- if(is.null(u))NULL else u[u$parameter%in%keys,]
    lo_col <- if(bands=="monte_carlo")"mc_lower" else "lower"
    hi_col <- if(bands=="monte_carlo")"mc_upper" else "upper"
    panel_jackknife <- panel==2 && show_jackknife
    values <- c(as.matrix(h[,keys,drop=FALSE]),truth)
    if (bands!="none") values <- c(values,relevant[[lo_col]],relevant[[hi_col]])
    if (panel_jackknife) values <- c(values,fit$theta_jackknife[keys])
    values <- values[is.finite(values)]
    if (!length(values)) stop("No finite trace values")
    xmax <- nrow(h); xright <- xmax+if(panel_jackknife)1 else 0
    plot(NA,xlim=c(if(is.null(u))1 else 0,max(2,xright)),ylim=range(values),
      xlab="MCEM iteration",ylab=if(panel==1)"Beta" else "Theta",
      main=if(panel==1)"Regression coefficients" else "Local expPS tilts")
    for (j in seq_along(keys)) {
      if (bands!="none") {
        z <- relevant[relevant$parameter==keys[j],]; z <- z[order(z$iteration),]
        good <- is.finite(z[[lo_col]]) & is.finite(z[[hi_col]])
        runs <- split(seq_len(nrow(z)),cumsum(c(TRUE,diff(good)!=0)))
        for (rr in runs) if(all(good[rr]) && length(rr)>1L)
          polygon(c(z$iteration[rr],rev(z$iteration[rr])),
            c(z[[lo_col]][rr],rev(z[[hi_col]][rr])),
            border=NA,col=adjustcolor(col[j],alpha.f=.12))
      }
      lines(seq_len(nrow(h)),h[[keys[j]]],col=col[j],lwd=1.7)
      if(!is.null(truth))abline(h=truth[j],col=col[j],lty=2,lwd=.8)
    }
    if (panel_jackknife)
      points(rep(xmax+.6,length(keys)),fit$theta_jackknife[keys],col=col,pch=16,cex=1.1,xpd=NA)
    legend("topright",keys,col=col,lty=1,bty="n",ncol=if(panel==1)1 else 3,
           cex=if(panel==1).8 else .65)
  }
  mtext(main,outer=TRUE,side=3,line=.5)
  jackknife_caption <- if(show_jackknife) "; filled points: jackknife-corrected theta (recommended estimate)" else ""
  mtext(switch(bands,statistical=sprintf("Shading: approximate %g%% local likelihood bands; dashed lines: simulation truths%s",band_level,jackknife_caption),
      monte_carlo=sprintf("Shading: approximate %g%% local Monte Carlo error bands; dashed lines: simulation truths%s",band_level,jackknife_caption),
      none=sprintf("Dashed lines: simulation truths%s",jackknife_caption)),outer=TRUE,side=1,line=1,cex=.75)
  invisible(u)
}
