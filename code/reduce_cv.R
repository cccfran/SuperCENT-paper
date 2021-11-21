#!/usr/bin/env Rscript
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table)

args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Job id", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "confint"
}

job_id <- args[1]
job_name <- "confint"

sim_setup <- fread(paste0("../output_simulation/", job_name, "_", job_id,".log"), nrows=2)

nsim = sim_setup$nsim

beta0 = sim_setup$beta0
beta_u = sim_setup$betau
beta_v = sim_setup$betav
beta0vec <- c(beta0, beta_u, beta_v)

xmat = sim_setup$x

ns = seq(sim_setup$n_min, sim_setup$n_max, sim_setup$n_gap)
n_tests = seq(sim_setup$n_test_min, sim_setup$n_test_max, sim_setup$n_test_gap)
ds = seq(sim_setup$d_min, sim_setup$d_max, sim_setup$d_gap)
betaus = seq(sim_setup$betau_min, sim_setup$betau_max, sim_setup$betau_gap)
betavs = seq(sim_setup$betav_min, sim_setup$betav_max, sim_setup$betav_gap)
epsas = seq(sim_setup$epsa_min, sim_setup$epsa_max, sim_setup$epsa_gap)
epsys = seq(sim_setup$epsy_min, sim_setup$epsy_max, sim_setup$epsy_gap)
obs_ratios = seq(sim_setup$obs_ratio_min, sim_setup$obs_ratio_max, sim_setup$obs_ratio_gap)

ss = 1:nsim

ret_list <- lapply(Sys.glob(paste0("../hpcc_output/", job_name, "_", job_id, "/*")), fread)
ret <- rbindlist(ret_list, fill = T)

outfile <- paste0("../output_simulation/", job_name, "_", job_id,  
                                   "_n", paste(ns, collapse = "-"),
                                   "_ntest", paste(n_tests, collapse = "-"),
                                   "_xmat", xmat,
                                   "_d", paste(ds, collapse = "-"),
                                   "_beta", paste(beta0, collapse = "-"),
                                   "_betau", paste(betaus, collapse = "-"),
                                   "_betav", paste(betavs, collapse = "-"),
                                   "_epsa", paste(epsas, collapse = "-"),
                                   "_epsy", paste(epsys, collapse = "-"), 
                                   "_obs", paste(obs_ratios, collapse = "-"), 
                                   "_nsim", nsim,
                                   ".csv")
fwrite(ret, outfile, verbose = T)



saveRDS(list(job_id = job_id,
             nsim = nsim,
             ss = ss,
             epsas = epsas,
             epsys = epsys,
             beta0 = beta0,
             betaus = betaus,
             betavs = betavs,
             ds = ds,
             ns = ns,
             n_tests = n_tests,
             obs_ratios = obs_ratios,
             xmat = sim_setup$xmat,
             scaleuvtest = sim_setup$scaleuvtest), 
        paste0("../output_simulation/", job_name, "_", job_id, ".RDS"))


# render_report(c(1601513), "cv")