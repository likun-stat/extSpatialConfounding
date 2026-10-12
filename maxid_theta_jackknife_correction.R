# Jackknife theta correction is now automatic: fit_mcem_fast(..., theta_jackknife=TRUE,
# the default) computes fit$theta_jackknife as part of every fit, at no extra MCMC
# cost (it reuses the last iteration's own E-step). See maxid_regression_fast_fit.R
# for the diagnosis, the .maxid_jackknife_theta_locaug/.maxid_jackknife_theta_classical
# helpers, and the integration -- this file no longer needs to duplicate that math.
#
# What's left here is a standalone re-check for two cases the automatic path doesn't
# cover:
#   (a) a fit saved before this change, or fit with theta_jackknife=FALSE, that has
#       no fit$theta_jackknife of its own;
#   (b) wanting an independent-seed re-verification: the automatic correction reuses
#       the last EM iteration's own E-step draws, which is cheap but ties the
#       correction to that one posterior draw. correct_theta_jackknife() instead
#       draws a brand-new E-step at the fit's converged parameters (same mechanism
#       as check_final_mcem() in check_mcem_performance.R), at the cost of one extra
#       full E-step (tens of seconds at T=600).
#
# Source maxid_regression_mcem.R and maxid_regression_fast_fit.R first.

correct_theta_jackknife <- function(sim,fit,ndraw=800L,burn=800L,seed=505L,maxit=60L) {
  dat <- .maxid_check(sim)
  if (!identical(fit$model_spec,maxid_model_spec()) || fit$alpha0!=dat$alpha0)
    stop("Use a fitted object with the same model and fixed alpha0")
  dat$cscale <- fit$cscale
  K <- ncol(dat$W)
  set.seed(seed)
  e <- e_step_joint(dat,fit$beta,fit$theta,fit$q,fit$joint_state,
                    ndraw=ndraw,burn=burn,thin=1L)
  meanZ_time <- apply(exp(e$qdraws),c(1,2),mean)
  theta_naive <- if (fit$location_augmentation) {
    meanInvZ_time <- apply(exp(-e$qdraws),c(1,2),mean)
    e <- NULL; gc(verbose=FALSE)
    .maxid_location_update(meanZ_time,meanInvZ_time,dat$ghg,fit$beta,fit$cscale,
                           sqrt(fit$theta),maxit)$eta^2
  } else {
    e <- NULL; gc(verbose=FALSE)
    1/(4*rowMeans(meanZ_time)^2)
  }
  theta_jack <- if (fit$location_augmentation) {
    .maxid_jackknife_theta_locaug(meanZ_time,meanInvZ_time,dat$ghg,fit$beta,fit$cscale,
                                  sqrt(fit$theta),maxit)
  } else .maxid_jackknife_theta_classical(meanZ_time)
  data.frame(parameter=names(fit$theta),fit_theta=as.numeric(fit$theta),
             naive_resampled=theta_naive,jackknife_corrected=theta_jack)
}
