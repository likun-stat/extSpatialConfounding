# Source the three model modules before calling these functions.
# Batch-means MCSE describes posterior Monte Carlo error, not parameter SE.
.maxid_mcse_batches <- function(draws) {
  draws <- as.matrix(draws) # posterior functionals by retained sweep
  M <- ncol(draws)
  batch_size <- max(2L,floor(sqrt(M)))
  batches <- M%/%batch_size
  if (batches<5L) return(rep(NA_real_,nrow(draws)))
  averages <- vapply(seq_len(batches),function(b)
    rowMeans(draws[,(b-1L)*batch_size+seq_len(batch_size),drop=FALSE]),
    numeric(nrow(draws)))
  averages <- matrix(averages,nrow(draws),batches)
  apply(averages,1,sd)/sqrt(batches)
}

check_final_mcem <- function(sim,fit,ndraw=2000L,burn=800L,seeds=c(101L,202L)) {
  if (length(seeds)<2L) stop("Supply at least two independent sampling seeds")
  dat <- .maxid_check(sim)
  if (!identical(fit$model_spec,maxid_model_spec()) || fit$alpha0!=dat$alpha0)
    stop("Use a fitted object with the same model and fixed alpha0")
  dat$cscale <- fit$cscale
  n <- dat$n; T <- dat$nt; K <- ncol(dat$W); N <- n*T
  X <- cbind(intercept=1,elevation=rep(dat$elev,T),GHG=rep(dat$ghg,each=n))
  dscale <- fit$cscale*dat$alpha0
  eta <- sqrt(fit$theta)
  target <- matrix((as.vector(dat$L)-as.numeric(X%*%fit$beta))/fit$cscale,n,T)
  out <- lapply(seeds,function(seed) {
    # Release caches from a preceding high-draw fit/chain before allocating
    # another full log-field draw array. This does not change the RNG or target.
    gc(verbose=FALSE)
    set.seed(seed)
    e <- e_step_joint(dat,fit$beta,fit$theta,fit$q,fit$joint_state,
                      ndraw=ndraw,burn=burn,thin=1L)
    meanZ_draws <- vapply(seq_len(ndraw),function(r)
      rowMeans(exp(matrix(e$qdraws[,,r],K,T))),numeric(K))
    meanZ_draws <- matrix(meanZ_draws,K,ndraw)
    # Fisher's identity: posterior-average complete-data scores estimate
    # the observed marginal likelihood score at the fixed returned parameters.
    score <- vapply(seq_len(ndraw),function(r) {
      g <- sign(target-e$logs[,,r])
      beta_score <- as.numeric(crossprod(X,as.vector(g)))/dscale
      # Equivalent location augmentation gives much less noisy scores for
      # the intercept and temporal trend when the white-noise scale is small.
      # Include the r=log V prior Jacobian: its h score has the +K/2 term.
      q_r <- matrix(e$qdraws[,,r],K,T)
      h_score <- K/2-.25*colSums(exp(-q_r))+
                 as.numeric(crossprod(fit$theta,exp(q_r)))
      beta_score[c(1,3)] <- c(sum(h_score),sum(h_score*dat$ghg))/fit$cscale
      eta_score <- T*(1-2*eta*meanZ_draws[,r])
      c(beta_score,eta_score)/N
    },numeric(3+K))
    rownames(score) <- c(names(fit$beta),paste0("eta",seq_len(K)))
    result <- list(seed=seed,meanZ=rowMeans(meanZ_draws),meanZ_mcse=.maxid_mcse_batches(meanZ_draws),
         score=rowMeans(score),score_mcse=.maxid_mcse_batches(score),
         ess=e$ess,acceptance=e$acceptance)
    e <- NULL
    gc(verbose=FALSE)
    result
  })
  means <- data.frame(k=seq_len(K),meanZ_1=out[[1]]$meanZ,meanZ_2=out[[2]]$meanZ,
              mcse_1=out[[1]]$meanZ_mcse,mcse_2=out[[2]]$meanZ_mcse)
  means$standardized_difference <- abs(means$meanZ_1-means$meanZ_2)/
                                  sqrt(means$mcse_1^2+means$mcse_2^2)
  scores <- data.frame(parameter=names(out[[1]]$score),score_1=out[[1]]$score,
                       score_2=out[[2]]$score,mcse_1=out[[1]]$score_mcse,
                       mcse_2=out[[2]]$score_mcse,row.names=NULL)
  list(chains=out,meanZ_comparison=means,score_comparison=scores,
       note=paste("Scores are divided by n*T; eta=sqrt(theta).",
                  "Intercept/GHG scores use the equivalent location prior to reduce Monte Carlo variance;",
                  "exact-zero tilts are nonregular and this is not regular boundary inference.",
                  "Interior eta scores should approach zero at a stationary fit;",
                  "at eta=0 a nonpositive one-sided score is compatible with a boundary maximum.",
                  "MCSE uses batch means and is approximate, not a parameter standard error."))
}

run_mcem_performance <- function(iterations=35L,draws_start=200L,draws_end=600L,
                                 outdir="mcem_results",seed=91L,alpha0=0.1,nt=400L,
                                 uncertainty="louis") {
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  terrain <- read.csv("elevation_grid.csv")
  d <- make_spatial_design(nx=25,ny=20,nt=nt,radii=6,elevation=terrain$elevation_m)
  sim <- simulate_spatial(d,beta=c(intercept=1,elevation=0.5,GHG=0.3),
                          cscale=0.12,alpha0=alpha0,seed=43)
  timing <- system.time(fit <- fit_mcem_fast(sim,iterations=iterations,
                    draws_start=draws_start,draws_end=draws_end,burn_first=800,burn=150,
                    m_cycles=2,seed=seed,estimate_c=FALSE,c_fixed=0.12,
                    uncertainty=uncertainty))
  saveRDS(sim,file.path(outdir,"simulation.rds"))
  saveRDS(fit,file.path(outdir,"fit.rds"))
  write.csv(fit$history,file.path(outdir,"history.csv"),row.names=FALSE)
  if(!is.null(fit$uncertainty))
    write.csv(fit$uncertainty,file.path(outdir,"uncertainty.csv"),row.names=FALSE)
  checks <- check_final_mcem(sim,fit)
  saveRDS(checks,file.path(outdir,"posterior_checks.rds"))
  write.csv(checks$meanZ_comparison,file.path(outdir,"posterior_meanZ.csv"),row.names=FALSE)
  write.csv(checks$score_comparison,file.path(outdir,"posterior_scores.csv"),row.names=FALSE)
  list(sim=sim,fit=fit,checks=checks,timing=timing)
}
