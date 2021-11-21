if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, matrixStats, dplyr, ggplot2, igraph,
               latex2exp, tidyverse,
               irlba, xtable, stargazer, circlize)

if(!require("SuperCENT")) devtools::install_github("cccfran/SuperCENT")
library("SuperCENT")

# Paramter setting
args <- commandArgs(trailingOnly = TRUE)

# ## time windows
gap <- as.numeric(ifelse(is.null(args[1]), 1000, args[1]))
## maximum iteration
max_iter <- as.numeric(ifelse(is.null(args[2]), 1000, args[2]))

print(paste0("Running trade premium with max iteration ",
             max_iter, " with rolling windows ", gap))

# read data
# fx
fx_raw <- fread("../data/FX_Sub.csv")
setnames(fx_raw,
         old = c("iso3"),
         new = c("code"))

fx_raw[, year := year(Date)]
fx_raw[, rx := exp(logrx)]
fx_raw_year <- fx_raw[, lapply(.SD, mean, na.rm=T), 
                      by = .(code, year), 
                      .SDcols = c("rx")]
fx_raw_year[, logrx := log(rx)]

# gdp 
if(!file.exists("../data/real_gdp_long.csv")) source("construct_gdp_data.R")
gdp_raw <- fread("../data/real_gdp_long.csv")

# trade
if(!file.exists("../data/trade_data_sub.csv")) source("construct_trade_data.R")
trade_raw <- fread("../data/trade_data_sub.csv")

for(yy in 1999:(2012 - gap + 1)) {
  
  print(yy)
  ymax <- min(yy+gap-1, 2013)
  years <- yy:ymax
  
  # prepare trade, gdp and rx as a rolling average
  ## trade
  trade <- trade_raw[year %in% years][, .(value = mean(value,na.rm = T)),
                                      by=.(from, to)]
  ## gdp
  gdp <- gdp_raw[year %in% years][, .(gdp = mean(gdp)), by = .(code)]
  # gdp_old <- merge(gdp_old, gdp, by = "code")
  # gdp <- gdp_raw[year %in% years][, .(gdp = mean(gdp)), by = .(code)]
  ## risk premium
  vars <- c("rx")
  risk_premium <- fx_raw_year
  premium <- fx_raw_year[year %in% years][, lapply(.SD, mean, na.rm=T), 
                                          by=.(code),.SDcols = vars]
  premium[, logrx := log(rx)]
  
  # merge
  premium <- merge(premium, gdp, by = "code", all.x = T)
  premium_exist <- premium[!is.na(logrx) & code %in% unique(trade$from) & !is.infinite(logrx),]
  premium_exist[, gdp_share := gdp / sum(premium_exist$gdp)]
  setorderv(premium_exist, "code")
  
  # trade
  trade_sub <- merge(premium_exist[, .(to = code, gdp_to = gdp)], 
                     trade, by = "to")
  trade_sub <- merge(premium_exist[, .(from = code, gdp_from = gdp)], 
                     trade_sub, by = "from")
  
  # standardize the unit: gdp in billion and trade in million
  trade_sub[, weight := value / (gdp_to + gdp_from) / 1e3]
  graph <- graph_from_data_frame(trade_sub[,.(from, to, weight)], directed = T)
  mat <- as_adjacency_matrix(graph, attr = "weight", sparse = F)
  
  # main
  y <- premium_exist$logrx
  X <- model.matrix(logrx ~ 1 + gdp_share, premium_exist)
  
  # two-stage
  ret <- two_stage(mat, X, y)
  # SuperCENT-lambdahat
  ret_lr <- supercent(mat, X, y, max_iter = max_iter)
  # SuperCENT-CV
  ret_cv <- cv.supercent(mat, X, y, lrange = 2^30, 
                         l = length(y)*(ret$epsy)^2*(ret$epsa)^2,
                         gap = 2^1, folds = 10, max_iter = max_iter)
  
  # centrality ranking
  ret_rank <- data.table(code = rownames(ret$A),
                         u = ret$u, v = ret$v)
  ret_lr_rank <- data.frame(code = rownames(ret_lr$A),
                            u = ret_lr$u, v = ret_lr$v)
  setDT(ret_lr_rank)
  
  ret_cv_rank <- data.frame(code = rownames(ret_cv$A),
                            u = ret_cv$u, v = ret_cv$v)
  setDT(ret_cv_rank)
  
  # save
  ret$cent <- ret_rank
  ret_lr$cent <- ret_lr_rank
  ret_cv$cent <- ret_cv_rank
  
  filepath <- paste0("../output_trade_premium/",
                     yy, "_gap", gap, "_miter", max_iter, ".RDS")
  saveRDS(list(ret = ret, ret_lr = ret_lr, ret_cv = ret_cv),
          filepath)
}
