if (!require("pacman")) install.packages("pacman")
pacman::p_load(optparse, pracma, data.table, matrixStats, dplyr, irlba, ggplot2)

option_list = list(
  make_option(c("--n_train"), type="integer", default="8", 
              help="training sample size", metavar="number"),
  make_option(c("--n_test"), type="integer", default="0", 
              help="test sample size", metavar="number"),
  make_option(c("--d"), type="character", default="0", 
              help="rank", metavar="number"),
  make_option(c("--beta0"), type="character", default='1-3-5', 
              help="betax", metavar="number"),
  make_option(c("--betau"), type="numeric", default="2", 
              help="betau", metavar="number"),
  make_option(c("--betav"), type="numeric", default="0", 
              help="betav", metavar="number"),
  make_option(c("--epsa"), type="numeric", default="0", 
              help="sigma_A", metavar="number"),
  make_option(c("--epsy"), type="numeric", default="-6", 
              help="sigma_y", metavar="number"),
  make_option(c("--nsim"), type="integer", default="100", 
              help="number of sim", metavar="number"),
  make_option(c("--seed"), type="integer", default="1", 
              help="seed", metavar="number"),
  make_option(c("--max_iter"), type="integer", default="1000", 
              help="max iteration number", metavar="number"),
  make_option(c("--tol"), type="numeric", default="1e-4", 
              help="tolerate", metavar="number"),
  make_option(c("--xmat"), type="character", default="a", 
              help="X matrix type", metavar="number"),
  make_option(c("--jobid"), type="character", default="1", 
              help="jobid", metavar="number"),
  make_option(c("--power_beta"), type="numeric", default="1", 
              help="power opt with betau betav", metavar="number"),
  make_option(c("--confint"), type="numeric", default="0", 
              help="compute confint", metavar="number"),
  make_option(c("--confint_A"), type="numeric", default="0", 
              help="compute confint for A", metavar="number"),
  make_option(c("--hulc"), type="numeric", default="0", 
              help="compute hulc confint", metavar="number"),
  make_option(c("--observed_ratio"), type="numeric", default="1", 
              help="Ratio of observed X, y", metavar="number"),
  make_option(c("--ibatch"), type="numeric", default="1", 
              help="batch index", metavar="number"),
  make_option(c("--batch_size"), type="numeric", default="100", 
              help="batch size", metavar="number"),
  make_option(c("--sym"), type="numeric", default="0", 
              help="Symmetric A or not", metavar="number")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt

# Rcpp::sourceCpp("utilities.cpp")
# source("utilities.R")
if(!require("SuperCENT")) devtools::install_github("cccfran/SuperCENT")
library("SuperCENT")
source("utils_sim.R")

new_dir <- paste0("../hpcc_output/confint_", opt$jobid, "/")
dir.create(new_dir)

# Setup -------------------------------------------------------------------

# 2^c("n_train", "n_test", "d", "epsa", "epsy", "betau", "betav")
if(opt$power_beta) {
  opt <- power_opt(opt, which = c("n_train", "n_test", "d", "epsa", "epsy", "betau", "betav"))  
} else {
  opt <- power_opt(opt, which = c("n_train", "n_test", "d", "epsa", "epsy"))
}
opt <- beta0_opt(opt)
if(opt$sym) opt$betav <- 0
beta0vec <- c(opt$beta0, opt$betau, opt$betav)
max_iter <- opt$max_iter

alpha <- .05


# Sim -------------------------------------------------------------------

job_id <- opt$jobid
# possible methods: c("ols", "oracle", "two_stage", "lr_oracle", "lr_plugin", "lr_cv")
# possible oracles CIs: c("two_stage_oracle", "lr_oracle_oracle", "lr_plugin_oracle", "lr_cv_oracle")
methods <- c("ols", "two_stage", "lr_oracle", "lr_cv")
ci_oracles <- c("two_stage_oracle", "lr_oracle_oracle", "lr_cv_oracle")
if(opt$confint) methods <- c(methods, ci_oracles)
if(opt$observed_ratio < 1) methods <- c(methods, "lr_oracle_sub", "two_stage_sub")

# ss ----------------------------------------------------------------------

for(ss in 1:opt$batch_size) {
  # set seed for this batch
  opt$seed <- (opt$ibatch-1)*opt$batch_size+ss
  if(ss %% 10 == 0) {print(ss); print(ret_df)}
  
  # generate fixed sample
  sample <- generate_sample(opt = opt, beta0 = opt$beta0, seed = 0, fixed_only = T)
  
  D = sample$D
  d = sample$d
  p = sample$p
  U_train = sample$U_train
  U_test = sample$U_test
  V_train = sample$V_train
  V_test = sample$V_test
  X_train = sample$X_train
  X_test = sample$X_test
  
  U <- rbind(U_train, U_test)
  V <- rbind(V_train, V_test)
  X <- rbind(X_train, X_test)
  
  observed <- sample$observed
  obs <- which(sample$observed == 1)
  unobs <- which(sample$observed != 1)
  n_obs <- sum(sample$observed == 1)
  
  sample_A_y <- generate_A_y(opt = opt, beta0 = opt$beta0, U = U, V = V, X = X)
  A_train <- sample_A_y$A_train
  A <- sample_A_y$A
  y_train <- sample_A_y$y_train
  y_test <- sample_A_y$y_test
  
  if(opt$observed_ratio < 1) {
    X_unobs <- X_train[unobs,,drop=F]
    X_train <- X_train[obs,,drop=F]
    y_unobs <- y_train[unobs]
    y_train <- y_train[obs]
  }
  
  betahat_two_stage <- betahat_lr_oracle <- betahat_lr_cv <- betahat_lr_plugin <- NULL
  A_hat_two_stage <- A_hat_lr_oracle <- A_hat_lr_cv <- A_hat_lr_plugin <- NULL
  
  # main --------------------------------------------------------------------
  
  ret_df <- NULL
  ret_ <- NULL
  ret_list <- list()
  iter <- 1
  system.time({
    for(m in methods) {
      
      print(m)
      
      df <- opt$n_train - length(beta0vec)
      lopt <- l_optimal_lr(d[1], beta0vec, opt$epsy^2, opt$epsa^2, opt$n_train, weights = observed)
      U_train_ <- U_train; V_train_ <- V_train; A_train_ <- A_train
      
      if(m == "ols") {
        ret <- two_stage(A_train, X_train, y_train, weights = observed)
        ret <- ret_sign(ret, beta0vec, U_train_, V_train_)
        fit <- lm(y_train ~ cbind(X_train, ret$u[obs], ret$v[obs]) - 1)
        interval <- confint.default(fit)
        sduv <- summary(fit)$coef[(p+1):(p+2),2]
        ci_a <- NULL
      }
      
      if(m == "oracle") {
        lopt <- l_optimal_lr(d[1], beta0vec, opt$epsy^2, opt$epsa^2, opt$n_train, weights = observed)
        ret <- oracle(A_train, X_train, y_train, d, U_train, V_train, 
                      beta0vec = beta0vec, beta_hat = beta0vec, 
                      epsa = opt$epsa, epsy = opt$epsy, l = lopt,
                      A_hat = U_train %*% D %*% t(V_train),
                      weights = observed, method = m)
      }
      
      if(m == "two_stage") {
        ret <- two_stage(A_train, X_train, y_train, weights = observed)
        betahat_two_stage <- ret$beta
        A_hat_two_stage <- ret$d * ret$u %*% t(ret$v)
        
        if(opt$hulc) interval_hulc <- hulc(A_train, X_train, y_train, FUN = two_stage, weights = observed)
      }
      
      if(m == "two_stage_oracle") {
        ret <- oracle(A_train, X_train, y_train, d, U_train, V_train, 
                      beta0vec = beta0vec, beta_hat = betahat_two_stage, 
                      epsa = opt$epsa, epsy = opt$epsy, 
                      A_hat = A_hat_two_stage,
                      weights = observed, method = m)
      }
      
      
      if(m == "lr_cv") {
        lopt = l_optimal_lr(d[1], beta0vec, opt$epsy^2, opt$epsa^2, opt$n_train, weights = observed)
        ret <- cv.supercent(A_train, X_train, y_train, l = lopt, lrange = 2^4, 
                            gap = 2, folds = 10, max_iter = max_iter, weights = observed)
        betahat_lr_cv <- ret$beta
        A_hat_lr_cv <- ret$d * ret$u %*% t(ret$v)
        
        if(opt$hulc) interval_hulc <- hulc(A_train, X_train, y_train, FUN = cv.supercent, l = lopt,  lrange = 2^4, 
                                           gap = 2^2, folds = 10, max_iter = max_iter)
      }
      
      if(m == "lr_cv_oracle") {
        lopt = l_optimal_lr(d[1], beta0vec, opt$epsy^2, opt$epsa^2, opt$n_train, weights = observed)
        ret <- oracle(A_train, X_train, y_train, d, U_train, V_train, 
                      beta0vec = beta0vec, beta_hat = betahat_lr_cv, 
                      epsa = opt$epsa, epsy = opt$epsy, l = lopt,
                      A_hat = A_hat_lr_cv, weights = observed, method = m)
      }
      
      if(m == "lr_plugin") {
        ret <- supercent(A_train, X_train, y_train, max_iter = max_iter, weights = observed)
        betahat_lr_plugin <- ret$beta
        A_hat_lr_plugin <- ret$d * ret$u %*% t(ret$v)
        
        if(opt$hulc) interval_hulc <- hulc(A_train, X_train, y_train, FUN = supercent, max_iter = max_iter, weights = observed)
      }
      
      if(m == "lr_plugin_oracle") {
        lopt = l_optimal_lr(d[1], beta0vec, opt$epsy^2, opt$epsa^2, opt$n_train, weights = observed)
        ret <- oracle(A_train, X_train, y_train, d, U_train, V_train, 
                      beta0vec = beta0vec, beta_hat = betahat_lr_plugin, 
                      epsa = opt$epsa, epsy = opt$epsy, l = lopt,
                      A_hat = A_hat_lr_plugin, weights = observed, method = m)
      }
      
      if(m == "lr_oracle") {
        lopt = l_optimal_lr(d[1], beta0vec, opt$epsy^2, opt$epsa^2, opt$n_train, weights = observed)
        ret <- supercent(A_train, X_train, y_train, l = lopt, max_iter = max_iter, weights = observed, verbose = 0)
        betahat_lr_oracle <- ret$beta
        A_hat_lr_oracle <- ret$d * ret$u %*% t(ret$v)
        
        if(opt$hulc) interval_hulc <- hulc(A_train, X_train, y_train, FUN = supercent, l = lopt, max_iter = max_iter, weights = observed)
      }
      
      if(m == "lr_oracle_oracle") {
        lopt = l_optimal_lr(d[1], beta0vec, opt$epsy^2, opt$epsa^2, opt$n_train, weights = observed)
        ret <- oracle(A_train, X_train, y_train, d, U_train, V_train, 
                      beta0vec = beta0vec, beta_hat = betahat_lr_oracle, 
                      epsa = opt$epsa, epsy = opt$epsy, l = lopt,
                      A_hat = A_hat_lr_oracle, weights = observed, method = m)
      }
      
      if(m == "lr_oracle_sub") {
        A_train_ <- A_train[obs, obs]
        U_train_ <- U_train[obs,,drop=F]
        U_train_ <- U_train_/sqrt(sum(U_train_^2)/n_obs)
        V_train_ <- V_train[obs,,drop=F]
        V_train_ <- V_train_/sqrt(sum(V_train_^2)/n_obs)
        lopt = l_optimal_lr(d[1], beta0vec, opt$epsy^2, opt$epsa^2, n_obs)
        ret <- supercent(A_train_, X_train, y_train, l = lopt, max_iter = max_iter)
        
        if(opt$hulc) interval_hulc <- hulc(A_train_, X_train, y_train, FUN = supercent, l = lopt, max_iter = max_iter)
      }
      
      if(m == "two_stage_sub") {
        A_train_ <- A_train[obs, obs]
        U_train_ <- U_train[obs,,drop=F]
        U_train_ <- U_train_/sqrt(sum(U_train_^2)/n_obs)
        V_train_ <- V_train[obs,,drop=F]
        V_train_ <- V_train_/sqrt(sum(V_train_^2)/n_obs)
        ret <- two_stage(A_train_, X_train, y_train)
        
        if(opt$hulc) interval_hulc <- hulc(A_train_, X_train, y_train, FUN = two_stage)
      }
      
      ret <- ret_sign(ret, beta0vec, U_train_, V_train_)
      
      sub_flag <- grepl("sub", m)
      
      ## testing
      if(opt$observed_ratio < 1 & !sub_flag) {
        test_err <- y_unobs - cbind(X_unobs, ret$u[unobs], ret$v[unobs]) %*% ret$beta
      } else {
        y_test_hat <- predict_supervised(ret, A, X_test, weights = observed)
        test_err <- y_test - y_test_hat
      }
      
      ## Confint
      if(opt$confint) {
        if(m != "ols") {
          interval <- confint(ret, ci = T)
          sduv <- sqrt(rate_betauv_two_stage(X = ret$X, 
                                             u = ret$u, 
                                             v = ret$v, 
                                             beta0 = ret$beta, 
                                             d = ret$d, 
                                             sigmay2 = ret$epsy^2,
                                             sigmaa2 = ret$epsa^2, 
                                             n = n, 
                                             output = "uv"))
        }
        
        
        # print(interval)
        # interval
        coverage <- matrix( ((interval[,1]-beta0vec) * (interval[,2]-beta0vec)) < 0, nrow = 1)
        width <- matrix(interval[,2] - interval[,1], nrow = 1)
      }
      
      if(opt$confint_A & m != "ols") {
        ci_a <- confint_A(ret, ci = T, A0 = U_train %*% D %*% t(V_train))
      }
      
      ## hulc
      if(opt$hulc & !(m %in% c("ols", "oracle"))) {
        coverage_hulc <- matrix( ((interval_hulc[,1]-beta0vec) * (interval_hulc[,2]-beta0vec)) < 0, nrow = 1)
        width_hulc <- matrix(interval_hulc[,2] - interval_hulc[,1], nrow = 1) 
      }
      
      # testing 
      test_ret <- predict_svd(A, ret$beta, X_test)
      test_ret <- ret_sign(test_ret, beta0vec, U_test, V_test)
      test_ret$y <- cbind(X_test, test_ret$u, test_ret$v) %*% ret$beta
      test_err_svd <- sum((y_test-test_ret$y)^2)
      
      # output
      tmp <- data.table(method = m,
                        alpha = alpha,
                        epsy = opt$epsy,
                        epsa = opt$epsa, 
                        n = ifelse(grepl("sub", m), n_obs, opt$n_train),
                        n_test = opt$n_test,
                        observed_ratio = opt$observed_ratio,
                        beta0 = paste(opt$beta0, collapse = "-"),
                        betau = opt$betau,
                        betav = opt$betav,
                        d = opt$d,
                        s = opt$seed,
                        d_hat = ret$d,
                        betau_hat = ret$beta[p+1],
                        betav_hat = ret$beta[p+2],
                        lopt = lopt,
                        iter = ifelse(is.null(ret$iter), NA, ret$iter),
                        residual = sum(ret$residuals^2),
                        test_err = sum(test_err^2),
                        u = ifelse(is.null(ret$u)|sub_flag, NA, spec_norm_diff(ret$u, U_train)),
                        v = ifelse(is.null(ret$u)|sub_flag, NA, spec_norm_diff(ret$v, V_train)),
                        u_f = ifelse(is.null(ret$u)|sub_flag, NA, norm(ret$u - U_train, "f")^2 ),
                        v_f = ifelse(is.null(ret$v)|sub_flag, NA, norm(ret$v - V_train, "f")^2 ),
                        u_s = ifelse(sub_flag, spec_norm_diff(ret$u, U_train_), spec_norm_diff(ret$u[obs], U_train[obs,])),
                        v_s = ifelse(sub_flag, spec_norm_diff(ret$v, V_train_), spec_norm_diff(ret$v[obs], V_train[obs,])),
                        u_s_f = ifelse(sub_flag, sum((ret$u - U_train_)^2), sum((ret$u[obs] - U_train[obs,])^2) ),
                        v_s_f = ifelse(sub_flag, sum((ret$v - V_train_)^2), sum((ret$v[obs] - V_train[obs,])^2) ),
                        u_sc = ifelse(sub_flag, NA, spec_norm_diff(ret$u[unobs], U_train[unobs,])),
                        v_sc = ifelse(sub_flag, NA, spec_norm_diff(ret$v[unobs], V_train[unobs,])),
                        u_sc_f = ifelse(sub_flag, NA, sum((ret$u[unobs] - U_train[unobs,])^2) ),
                        v_sc_f = ifelse(sub_flag, NA, sum((ret$v[unobs] - V_train[unobs,])^2) ),
                        A_err = ifelse(is.null(ret$u)|sub_flag, NA, norm(ret$d*ret$u%*%t(ret$v) - U_train %*% D %*% t(V_train), "f")^2),
                        A_s_err = ifelse(sub_flag, norm(ret$d*ret$u%*%t(ret$v) - U_train_ %*% D %*% t(V_train_), "f")^2, 
                                         norm(ret$d*ret$u[obs]%*%t(ret$v[obs]) - U_train[obs,,drop=F] %*% D %*% t(V_train[obs,,drop=F]), "f")^2),
                        A0_2 = norm(U_train %*% D %*% t(V_train), "f")^2,
                        A0_s_2 = norm(U_train[obs,,drop=F] %*% D %*% t(V_train[obs,,drop=F]), "f")^2,
                        rss_y = sum((ret$residuals)^2),
                        rss_a = ifelse(sub_flag, NA, sum((ret$d * ret$u %*% t(ret$v) - A_train)^2)),
                        hbetau = ret$beta[p+1], 
                        hbetav = ret$beta[p+2], 
                        betau_rate = (ret$beta[p+1] - opt$betau),
                        betav_rate = (ret$beta[p+2] - opt$betav),
                        hat_sigmaa2 = ret$epsa^2,
                        hat_sigmay2 = ret$epsy^2,
                        hsd_betau = sduv[1],
                        hsd_betav = sduv[2],
                        l = ret$l,
                        # u_angle = ifelse(is.null(ret$u)|sub_flag, NA, acos(c(t(ret$u) %*% U_train/opt$n_train))*180/pi),
                        # v_angle = ifelse(is.null(ret$v)|sub_flag, NA, acos(c(t(ret$v) %*% V_train/opt$n_train))*180/pi),
                        test_err_svd = test_err_svd,
                        cov_A_FWER = ifelse(opt$confint_A, all(ci_a$covered_FWER), NA),
                        width_A_FWER = ifelse(opt$confint_A, mean(ci_a$upper_FWER - ci_a$lower_FWER), NA),
                        cov_A = ifelse(opt$confint_A, mean(ci_a$covered), NA),
                        width_A = ifelse(opt$confint_A, mean(ci_a$upper - ci_a$lower), NA)
                        
      )
      
      # output beta hat and CIs
      tmp[, paste0("hbetax_", 1:p) := as.list(ret$beta[1:p])]
      new_ret_df <- tmp
      
      if(opt$confint) {
        coverage <- as.data.table(coverage)
        names(coverage) <- paste0("cov_", c(1:p, "u", "v"))
        width <- as.data.table(width)
        names(width) <- paste0("width_", c(1:p, "u", "v"))
        
        lower <- as.data.table(t(interval[,1]))
        names(lower) <- paste0("lower_", c(1:p, "u", "v"))
        upper <- as.data.table(t(interval[,2]))
        names(upper) <- paste0("upper_", c(1:p, "u", "v"))
        
        new_ret_df <- cbind(tmp, coverage, width, lower, upper) 
        
        
      }
      
      if(opt$hulc & !(m %in% c("ols", "oracle"))) {
        coverage_hulc <- as.data.table(coverage_hulc)
        names(coverage_hulc) <- paste0("cov_hulc_", c(1:p, "u", "v"))
        width_hulc <- as.data.table(width_hulc)
        names(width_hulc) <- paste0("width_hulc_", c(1:p, "u", "v"))
        
        new_ret_df <- c(new_ret_df, coverage_hulc, width_hulc)
      }
      
      ret_df <- rbind(ret_df, 
                      new_ret_df,
                      fill = T
      )
      
      print(ret$iter)
      
      # write this round
      try({
        fwrite(ret_df, paste0(new_dir,
                              "confint",
                              "_n", opt$n_train,
                              "_ntest", opt$n_test,
                              "_xmat", opt$xmat,
                              "_d", paste(opt$d, collapse = "-"),
                              "_beta", paste(opt$beta0, collapse = "-"),
                              "_betau", opt$betau,
                              "_betav", opt$betav,
                              "_epsa", opt$epsa,
                              "_epsy", opt$epsy,
                              "_s", opt$seed,
                              "_obs", opt$observed_ratio,
                              ".csv"))
      })
      
    }
    # end loop for methods 
  } )
  # end system.time()
} 
# end loop for s

print(ret_df)

# Output ------------------------------------------------------------------

fwrite(ret_df, paste0(new_dir,
                      "confint",
                      "_n", opt$n_train,
                      "_ntest", opt$n_test,
                      "_xmat", opt$xmat,
                      "_d", paste(opt$d, collapse = "-"),
                      "_beta", paste(opt$beta0, collapse = "-"),
                      "_betau", opt$betau,
                      "_betav", opt$betav,
                      "_epsa", opt$epsa,
                      "_epsy", opt$epsy,
                      "_s", opt$seed,
                      "_obs", opt$observed_ratio,
                      ".csv"))
