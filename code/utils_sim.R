library(ggplot2)


generate_sample <- function(opt, beta0 = NULL, seed = 0, 
                            fixed_only = FALSE)
{
  d = as.numeric(strsplit(as.character(opt$d), "-")[[1]])
  r = length(d)
  # if(is.null(beta0)) beta0 <- as.numeric(strsplit(as.character(opt$beta0), "-")[[1]])
  beta0
  beta_u = opt$betau
  beta_v = opt$betav
  epsy = opt$epsy
  epsa = opt$epsa
  N = opt$n_train + opt$n_test
  
  # Fixed
  ## U, V
  set.seed(seed)
  D <- diag(d, nrow = r)
  
  U_train <- matrix(rnorm(r * opt$n_train), nrow = opt$n_train)
  U_norm <- sqrt(colSums(U_train^2))
  U_train <- sweep(U_train, 2, U_norm, "/") * sqrt(opt$n_train)
  
  V_train <- matrix(rnorm(r * opt$n_train), nrow = opt$n_train)
  V_norm <- sqrt(colSums(V_train^2))
  V_train <- sweep(V_train, 2, V_norm, "/") * sqrt(opt$n_train)
  
  U_test <- matrix(rnorm(opt$n_test), nrow = opt$n_test)
  V_test <- matrix(rnorm(opt$n_test), nrow = opt$n_test)
  U_test <- sweep(U_test, 2, U_norm, "/") * sqrt(opt$n_train) 
  V_test <- sweep(V_test, 2, V_norm, "/") * sqrt(opt$n_train) 
  
  if(opt$xmat == "c") {
    V_train <- 0.5*U_train + matrix(rnorm(r * opt$n_train), nrow = opt$n_train)
    V_norm <- sqrt(colSums(V_train^2))
    V_train <- sweep(V_train, 2, V_norm, "/") * sqrt(opt$n_train)
    
    V_test <- 0.5*U_test + matrix(rnorm(opt$n_test), nrow = opt$n_test)
    V_test <- sweep(V_test, 2, V_norm, "/") * sqrt(opt$n_train) 
  }
  
  
  U = rbind(U_train, U_test)
  V = rbind(V_train, V_test)
  
  beta0vec = c(beta0, beta_u, beta_v)
  
  ret <- list(d = d, D = d,
              U_train = U_train, U_test = U_test, 
              V_train = V_train, V_test = V_test,
              beta0vec = beta0vec, 
              n_test = opt$n_test)
  
  if(!is.null(beta0)) {
    ## X
    p = length(beta0) 
    X <- matrix(rnorm(N*(p-1)), nrow=N, ncol = p - 1, byrow = F)
    X <- pracma::gramSchmidt(X)$Q
    ## orthogonal design
    if(opt$xmat %in% c("a", "c") ) {
      X <- X
    } else if(opt$xmat == "b") { ## X depends on U, V
      # tmp <- opt$betau/(p-1) * U[, 1:(p-1)] + opt$betav/(p-1) * V[, 1:(p-1)]
      tmp <- opt$betau/(p-1) * U[, 1] + opt$betav/(p-1) * V[, 1]
      # if(r < (p-1)) {tmp <- cbind(tmp, rep(0, p-1-r))}
      X <- sweep(X, 1, tmp, "*")
      X <- sweep(X, 2, sqrt(colSums(X^2)), "/") * sqrt(N)
    }
    X <- cbind(1, X)
    X_train <- X[1:opt$n_train, , drop=F]
    X_test <- X[(opt$n_train+1):N, , drop=F]
    
    ret <- list(d = d, p = p, D = d,
                U_train = U_train, U_test = U_test, 
                V_train = V_train, V_test = V_test,
                X_train = X_train, X_test = X_test,
                beta0vec = beta0vec, 
                n_test = opt$n_test)
  }
  
  
  # Random
  ## A, y
  if(!fixed_only) {
    set.seed(opt$seed)
    y <- X %*% beta0 + beta_u * U[,1] + beta_v * V[,1] + rnorm(N, sd = epsy)
    A <- matrix(rnorm(N^2, sd = epsa), nrow = N) 
    if(opt$sym == 0) {
      A <- A + U %*% D %*% t(V)
    } else {
      A[lower.tri(A)] <- t(A)[lower.tri(A)]
      A <- A + U %*% D %*% t(U) 
    }
    
    A_train <- A[1:opt$n_train, 1:opt$n_train]
    y_train <- y[1:opt$n_train]
    y_test <- y[(opt$n_train+1):N]
    
    ret <- list(d = d, p = p, D = d,
                A = A, A_train = A_train, 
                U_train = U_train, U_test = U_test, 
                V_train = V_train, V_test = V_test,
                X_train = X_train, X_test = X_test,
                y_train = y_train, y_test = y_test,
                beta0vec = beta0vec, 
                epsy = epsy, epsa = epsa,
                n_test = opt$n_test)
  }
  
  # observed <- rep(1, opt$n_train)
  # if(opt$observed_ratio < 1) {
  #   set.seed(opt$seed)
  #   observed[sample(1:opt$n_train, opt$n_train*(1-opt$observed_ratio), replace = F)] <- 0
  # }
  
  n_obs <- opt$n_train*opt$observed_ratio
  observed <- c(rep(1, n_obs), 
                rep(0, opt$n_train - n_obs))
  
  ret$observed <- observed
  
  return(ret) 
}


generate_A_y <- function(opt, beta0 = NULL, U, V, X)
{
  set.seed(opt$seed)
  # if(is.null(beta0)) beta0 <- as.numeric(strsplit(as.character(opt$beta0), "-")[[1]])
  d = as.numeric(strsplit(as.character(opt$d), "-")[[1]])
  r = length(d)
  D <- diag(d, nrow = r)
  
  N <- opt$n_train + opt$n_test
  epsy_p2 = opt$epsy
  epsa_p2 = opt$epsa
  
  beta0vec = c(beta0, opt$betau, opt$betav)

  if(is.null(X)) {
    y <- opt$betau * U[,1] + opt$betav * V[,1] + rnorm(N, sd = epsy_p2)
  } else {
    y <- X %*% beta0 + opt$betau * U[,1] + opt$betav * V[,1] + rnorm(N, sd = epsy_p2)
  }
  
  A <- matrix(rnorm(N^2, sd = epsa_p2), nrow = N) 
  if(opt$sym == 0) {
    A <- A + U %*% D %*% t(V)
  } else {
    A[lower.tri(A)] <- t(A)[lower.tri(A)]
    A <- A + U %*% D %*% t(U) 
  }
  # A <- U %*% D %*% t(V) + matrix(rnorm(N^2, sd = epsa_p2), nrow = N)

  A_train <- A[1:opt$n_train, 1:opt$n_train]
  y_train <- y[1:opt$n_train]
  y_test <- y[(opt$n_train+1):N]

  ret <- list(A = A, A_train = A_train, 
            y_train = y_train, y_test = y_test,
            beta0vec = beta0vec)
}

generate_UV <- function(n, r, seed = 0)
{
  set.seed(seed)
  
  U_train <- matrix(rnorm(r * n), nrow = n)
  U_train <- sweep(U_train, 2, sqrt(colSums(U_train^2)), "/") * sqrt(n)
  V_train <- matrix(rnorm(r * n), nrow = n)
  V_train <- sweep(V_train, 2, sqrt(colSums(V_train^2)), "/") * sqrt(n)

  return(list(U_train = U_train, V_train = V_train))
}

power_opt <- function(opt, which = c("n_train","n_test", "d", "epsa", "epsy", "betau", "betav")) {
  if("n_train" %in% which) opt$n_train <- 2^opt$n_train
  if("n_test" %in% which) opt$n_test <- 2^opt$n_test
  if("d" %in% which) opt$d <- 2^as.numeric(opt$d)
  if("epsa" %in% which) opt$epsa <- 2^opt$epsa
  if("epsy" %in% which) opt$epsy <- 2^opt$epsy
  if("betau" %in% which) opt$betau <- 2^opt$betau
  if("betav" %in% which) opt$betav <- 2^opt$betav

  return(opt)
}

beta0_opt <- function(opt) {
  beta0 <- NULL
  if(opt$beta0 != '0') {
    beta0 <- as.numeric(strsplit(as.character(opt$beta0), "-")[[1]])
  }
  opt$beta0 <- beta0

  return(opt)
}

add_general_label <- function(p, labelR, labelT) {
  
  # Get the ggplot grob
  z <- ggplotGrob(p)
  
  # Get the positions of the strips in the gtable: t = top, l = left, ...
  posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
  posT <- subset(z$layout, grepl("strip-t", name), select = t:r)
  
  # Add a new column to the right of current right strips, 
  # and a new row on top of current top strips
  width <- z$widths[max(posR$r)]    # width of current right strips
  height <- z$heights[min(posT$t)]  # height of current top strips
  
  z <- gtable_add_cols(z, width, max(posR$r))  
  z <- gtable_add_rows(z, height, min(posT$t)-1)
  
  # Construct the new strip grobs
  stripR <- gTree(name = "Strip_right", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelR, rot = -90, gp = gpar(fontsize = 15, col = "grey10"))))
  
  stripT <- gTree(name = "Strip_top", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelT, gp = gpar(fontsize = 15, col = "grey10"))))
  
  # Position the grobs in the gtable
  z <- gtable_add_grob(z, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
  z <- gtable_add_grob(z, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
  
  # Add small gaps between strips
  z <- gtable_add_cols(z, unit(1/5, "line"), max(posR$r))
  z <- gtable_add_rows(z, unit(1/5, "line"), min(posT$t))
  
  # Draw it
  # grid.newpage()
  # grid.draw(z)
  
  z
}

file_name <- function(sim_setup) {
  ds <- as.numeric(strsplit(as.character(sim_setup$d), "-")[[1]])
  beta0vec <- c(beta0, sim_setup$betau, sim_setup$betav)
  lopt = sim_setup$epsy_max^2/ sim_setup$epsa_max^2
  lmin = floor(lopt-10); lmax = floor(lopt+10)
  tol = 1e-4
  
  paste0(
    "_n", sim_setup$n_train,
    "_ntest", sim_setup$n_test,
    "_xmat", sim_setup$xmat,
    "_d", min(ds), "-", max(ds),
    "_beta", paste(beta0vec, collapse = "-"),
    "_epsa", sim_setup$epsa_min, "-", sim_setup$epsa_max, 
    "_epsy", sim_setup$epsy_min, "-", sim_setup$epsy_max, 
    "_l", paste(c(lmin, lmax, sim_setup$lgap), collapse = "-"),
    "_tol", tol,
    "_maxiter", sim_setup$max_iter, 
    "_seed", 1, "-", sim_setup$nsim
  )
}

tab_fig_path <- function(sim_setup, desc, params, ext = "pdf")
{
  paste0("../table_figure/", 
         params$job_name, "_", params$job_id, 
         "_", desc,
         "_n", paste(sim_setup$ns, collapse = "-"),
         "_ntest", paste(sim_setup$n_tests, collapse = "-"),
         "_xmat", sim_setup$xmat,
         "_d", paste(sim_setup$ds, collapse = "-"),
         "_beta", paste(sim_setup$beta0, collapse = "-"),
         "_betau", paste(sim_setup$betaus, collapse = "-"),
         "_betav", paste(sim_setup$betavs, collapse = "-"),
         "_epsa", paste(sim_setup$epsas, collapse = "-"),
         "_epsy", paste(sim_setup$epsys, collapse = "-"), 
         "_nsim", max(sim_setup$ss),
         ".", ext)
}


my_theme <- 
  theme_bw() +
  theme(plot.subtitle=element_text(size = 15),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        # axis.text.x = element_text(size = 15, angle = -90, vjust = 0.5, hjust = 1),
        panel.border = element_rect(colour = "gray75", fill=NA, size=2),
        # panel.grid.major = element_blank(),
        strip.text = element_text(size = 15)
  )


TeX_label <- function(var, pre = NULL, post = NULL)
{
  match_tbl <- data.table(vars = c("epsy", "epsa", "betau", "betav", "d",
                                   "observed_ratio"),
                          TeX = paste0(pre, c("$\\log_2(\\sigma_y)$", 
                                  "$\\log_2(\\sigma_a)$",
                                  "$\\log_2(\\beta_u)$",
                                  "$\\log_2(\\beta_v)$",
                                  "$\\log_2(d)$",
                                  "Observed ratio "), post))
  TeX(match_tbl[vars == var, TeX])
}


plot_ratio_mean <- function(ret, methods, value, 
                            base = "two_stage", 
                            FUN = mean,
                            title = NULL, subtitle = NULL,
                            xaxis = "epsy", hpanel = "epsa", vpanel = "d",
                            scale_y_percent = F, error_bar = F,
                            save = F, plotname = value,
                            width = 9, height = 12) {
  
  ret <- ret[method %in% methods]
  
  labs_cols <- dplyr::inner_join(method_labels_tbl(), 
                                 data.table(method = unique(ret$method)),
                                 by = "method")
  
  ret$method <- factor(ret$method, levels = labs_cols$method)
  
  # get ratio
  ret_tmp <- dcast(ret, 
                   epsy + epsa + betau + betav + d + s ~ method, value.var = value)
  ret_tmp[, (methods) := lapply(.SD, function(x) (x-get(base))/get(base)), 
          .SDcols = methods]
  ret_tmp <- melt(ret_tmp, 
                  id.vars = c("epsy", "epsa", "betau", "betav", "d", "s"),
                  variable.name = "method",
                  value.name = "ratio")
  
  # get title
  if(is.null(title)) title <- value
  
  # plot
  p <- ret_tmp[,.(mean_ratio = (FUN(ratio, na.rm = T)),
                  sd_ratio = (sd(ratio, na.rm = T))), 
               by=c("method", "epsy", "epsa", "betau", "betav", "d")] %>%
    ggplot(aes(x = log2(get(xaxis)), y = mean_ratio, 
               group = method, col = method, shape = method)) +
    # geom_line(position = position_dodge(width=0.9)) +
    # geom_point(position = position_dodge(width=0.9)) +
    geom_point(size = 3) + geom_line() +
    xlab(TeX_label(xaxis)) +
    scale_x_continuous(breaks = sim_setup$epsys) +
    # scale_y_continuous(labels = scales::percent) +
    ylab("") +
    ggtitle(title) +
    labs(subtitle = subtitle) +
    my_theme + 
    scale_color_manual(values = labs_cols$cols, labels = labs_cols$labels) +
    scale_shape_manual(values = labs_cols$shapes, labels = labs_cols$labels) +
    # scale_color_discrete(labels = c("lr", "oracle", "svd")) +
    # coord_cartesian(ylim = c(0,2)) +
    facet_grid(log2(get(hpanel)) ~ log2(get(vpanel)), scales = "free")
  
  if(error_bar) p <- p + geom_errorbar(aes(ymax = mean_ratio + sd_ratio/sqrt(sim_setup$nsim),
                                           ymin = mean_ratio - sd_ratio/sqrt(sim_setup$nsim)),
                                       position = position_dodge(width=0.9))
  if(scale_y_percent) p <- p + scale_y_continuous(labels = scales::percent)
  z <- add_general_label(p, TeX_label(hpanel), TeX_label(vpanel))
  
  if(save) {
    ggsave(z, 
             filename = paste0("../output_simulation/plot/",
                               paste(sim_setup$job_id, collapse = "-"),
                               "_ratio_", plotname, ".pdf"),
             width = width, height = height)
  }
  
  grid.newpage()
  grid.draw(z)
}

plot_boxplot <- function(ret, methods, value, 
                         FUN = mean,
                         scale_y_transform = log10,
                         title = NULL, subtitle = NULL,
                         xaxis = "epsy", hpanel = "epsa", vpanel = "d",
                         save = F, plotname = value,
                         width = 9, height = 12, semi = F, gray = F,  value0 = NULL,
                         ylim = NULL, x_breaks = NULL,
                         y_scales = NULL) {
  
  # get title
  if(is.null(title)) title <- value
  
  ret <- ret[method %in% methods]
  
  labs_cols <- dplyr::inner_join(method_labels_tbl(semi = semi, gray = gray), 
                                 data.table(method = unique(ret$method)),
                                 by = "method")
  
  dot_col <- ifelse(gray, "red", "red")
  
  ret$method <- factor(ret$method, levels = labs_cols$method, ordered = T)
  
  param_groups <- c("method", "epsy", "epsa", "betau", "betav", "d")
  if(semi) param_groups <- c(param_groups, "observed_ratio")
  
  if(is.null(x_breaks)) x_breaks <- sim_setup[[xaxis]]
  if(is.null(x_breaks)) x_breaks <- log2(unique(ret[[xaxis]]))
    
  # theoretical
  # value0 <- paste0(value, '0')
  
  # plot
  p <- ret %>%
    ggplot(aes(x = log2(get(xaxis)), 
               y = scale_y_transform(get(value)), 
               group = interaction(get(xaxis), method), 
               # col = method,
               fill = method 
               )) +
    geom_boxplot(coef = 5,
                 col = "grey30",
                 outlier.size = 0.1) +
    # geom_line(data = ret[,
    #                      .(y = (FUN(get(value), na.rm = T))),
    #                      by=param_groups],
    #           aes(x = log2(get(xaxis)),
    #               y = scale_y_transform(y),
    #               group = method,
    #               col = method),
    #           position = position_dodge(width=0.1)
    # ) +
    facet_grid(log2(get(hpanel)) ~ log2(get(vpanel)), scales = "free") +
    xlab(TeX_label(xaxis)) +
    scale_x_continuous(breaks = x_breaks) +
    scale_color_manual(values = labs_cols$cols, labels = labs_cols$labels) +
    scale_fill_manual(values = labs_cols$cols, labels = labs_cols$labels) +
    scale_shape_manual(values = labs_cols$shape, labels = labs_cols$labels) +
    # scale_y_continuous(labels = scales::percent) +
    ylab("") +
    ggtitle(title) +
    labs(subtitle = subtitle) +
    my_theme +
    theme(legend.key.width = unit(1, "cm"),
          legend.key.height = unit(1.5, "cm"),
          legend.margin =margin(t = -10))
  # scale_color_discrete(labels = c("lr", "oracle", "svd")) +
  # coord_cartesian(ylim = c(0,2)) +
  
  if(!is.null(y_scales)) {
    p <- p + facet_grid_sc(rows = vars(log2(get(hpanel))),
                           cols = vars(log2(get(vpanel))),
                           scales = list(y = y_scales)) 
  }
  
  if(!is.null(value0)) p <- p +
    geom_point(data = unique(ret, by = c(param_groups)),
               aes(y = scale_y_transform(get(value0)),
                   shape = method),
               col = dot_col, size = 4, stroke = 1.1,
               alpha = 1,
               position = position_dodge(1.5)) 
  
  if(!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  
  if(length(methods) >4 ) p <- p + guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  if(grepl("diff", value)) p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  if(grepl("bias", value)) p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  
  z <- add_general_label(p, TeX_label(hpanel), TeX_label(vpanel))
  
  if(save) {
    ggsave(z, 
             filename = paste0("../output_simulation/plot/",
                               paste(sim_setup$job_id, collapse = "-"),
                               "_boxplot_", plotname, ".pdf"),
             width = width, height = height)
  }
  
  grid.newpage()
  grid.draw(z)
}

plot_mean <- function(ret, methods, value, 
                      FUN = mean,
                      error_bar = T,
                      scale_y_transform = log2,
                      title = NULL, subtitle = NULL,
                      xaxis = "epsy", hpanel = "epsa", vpanel = "d",
                      save = F, plotname = value,
                      width = 9, height = 12, scales = "free",
                      ylim = NULL) {
  
  ret <- ret[method %in% methods]
  
  labs_cols <- dplyr::inner_join(method_labels_tbl(gray = F), 
                                 data.table(method = unique(ret$method)),
                                 by = "method")
  
  ret$method <- factor(ret$method, levels = labs_cols$method)
  
  # get title
  if(is.null(title)) title <- value
  
  # plot
  p <- ret[,
           .(mean_ratio = (FUN(get(value), na.rm = T)),
             sd_ratio = (sd(get(value), na.rm = T))), 
           by=c("method", "epsy", "epsa", "betau", "betav", "d")] %>%
    ggplot(aes(x = log2(get(xaxis)), 
               y = scale_y_transform(mean_ratio), 
               group = method, 
               fill = method,
               col = method,
               shape = method)) +
    # geom_line(position = position_dodge(width=0.9)) +
    # geom_point(position = position_dodge(width=0.9)) +
    geom_point(size = 3, stroke = 1.5,
               position = position_dodge(width = 1.2)) + 
    # geom_line(position = position_dodge(width = 1)) +
    facet_grid(log2(get(hpanel)) ~ log2(get(vpanel)), scales = scales) +
    xlab(TeX_label(xaxis)) +
    scale_x_continuous(breaks = sim_setup$epsys) +
    scale_color_manual(values = labs_cols$cols, labels = labs_cols$labels) +
    scale_fill_manual(values = labs_cols$cols, labels = labs_cols$labels) +
    scale_shape_manual(values = labs_cols$shapes, labels = labs_cols$labels) +
    ylab("") +
    ggtitle(title) +
    labs(subtitle = subtitle) +
    my_theme 
  # coord_cartesian(ylim = c(0,2)) +
  
  if(error_bar) p <- p +
    geom_errorbar(aes(ymax = log2(mean_ratio + sd_ratio/sqrt(sim_setup$nsim)),
                      ymin = log2(mean_ratio - sd_ratio/sqrt(sim_setup$nsim))),
                  position = position_dodge(width=0.9)) 
  
  if(grepl("cov", value)) {
    p <- p + geom_hline(yintercept = 0.95, linetype = "dashed")
    p <- p + scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))
  } 
  
  if(length(methods) > 5) p <- p + guides(col=guide_legend(nrow=2,byrow=TRUE))
  if(!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  
  z <- add_general_label(p, TeX_label(hpanel), TeX_label(vpanel))
  
  if(save) {
    ggsave(z, 
             filename = paste0("../output_simulation/plot/",
                               paste(sim_setup$job_id, collapse = "-"),
                               "_", plotname, ".pdf"),
             width = width, height = height)
  }
  
  grid.newpage()
  grid.draw(z)
}


plot_line_mean <- function(ret, methods, value, 
                      FUN = mean,
                      error_bar = T,
                      scale_y_transform = log2,
                      title = NULL, subtitle = NULL,
                      xaxis = "epsy", hpanel = "epsa", vpanel = "d",
                      save = F, plotname = value,
                      width = 9, height = 12, scales = "free",
                      ylim = NULL, x_breaks = NULL,
                      x_lab_pre = NULL) {
  
  ret <- ret[method %in% methods]
  
  labs_cols <- dplyr::inner_join(method_labels_tbl(gray = F), 
                                 data.table(method = unique(ret$method)),
                                 by = "method")
  
  labs_cols[grepl("lr", method), labels := "SuperCENT"]
  labs_cols[!grepl("lr", method), labels := "Two-stage"]
  
  ret$method <- factor(ret$method, levels = labs_cols$method)
  
  if(is.null(x_breaks)) x_breaks <- sim_setup[[xaxis]]
  if(is.null(x_breaks)) x_breaks <- log2(unique(ret[[xaxis]]))
  
  # get title
  if(is.null(title)) title <- value
  
  # plot
  p <- ret[,
           .(mean_ratio = (FUN(get(value), na.rm = T)),
             sd_ratio = (sd(get(value), na.rm = T))), 
           by=c("method", "epsy", "epsa", "betau", "betav", "d")] %>%
    ggplot(aes(x = log2(get(xaxis)), 
               y = scale_y_transform(mean_ratio), 
               group = method, 
               fill = method,
               col = method,
               linetype = method,
               shape = method)) +
    geom_line(size = 2) +
    # geom_point(position = position_dodge(width=0.9)) +
    # geom_point(size = 3, stroke = 1.5,
    #            position = position_dodge(width = 1.2)) + 
    # geom_line(position = position_dodge(width = 1)) +
    # facet_grid(log2(get(hpanel)) ~ log2(get(vpanel)), scales = scales) +
    xlab(TeX_label(xaxis, x_lab_pre)) +
    scale_x_continuous(breaks = x_breaks) +
    scale_color_manual(values = labs_cols$cols, labels = labs_cols$labels) +
    scale_fill_manual(values = labs_cols$cols, labels = labs_cols$labels) +
    scale_shape_manual(values = labs_cols$shapes, labels = labs_cols$labels) +
    scale_linetype_manual(values = labs_cols$linetypes, labels = labs_cols$labels) +
    ylab("") +
    ggtitle(title) +
    labs(subtitle = subtitle) +
    my_theme +
    theme(legend.key.width = unit(2,"cm"))
  # coord_cartesian(ylim = c(0,2)) +
  
  if(grepl("cov", value)) {
    p <- p + geom_hline(yintercept = 0.95, linetype = "dashed")
    p <- p + scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))
  } 
  
  if(length(methods) > 5) p <- p + guides(col=guide_legend(nrow=2,byrow=TRUE))
  if(!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  
  # z <- add_general_label(p, TeX_label(hpanel), TeX_label(vpanel))
  # 
  # if(save) {
  #   ggsave(z, 
  #          filename = paste0("../output_simulation/plot/",
  #                            paste(sim_setup$job_id, collapse = "-"),
  #                            "_", plotname, ".pdf"),
  #          width = width, height = height)
  # }
  # 
  # grid.newpage()
  # grid.draw(z)
  p
}


render_report <- function(job_id, job_name, type = "pdf") {
  rmarkdown::render(
    paste0(job_name, ".Rmd"), params = list(
      job_id = job_id,
      job_name = job_name
    ),
    output_file = paste0("../write-ups/simulation/Simulation-", job_name, "_", job_id, ".", type)
  )
}


render_case_study <- function(case, setting, gap, max_iter, excess_return, type = "html") {
  params_list = list(
    setting = setting, gap = gap,
    max_iter = max_iter
  )
  
  if(!is.null(excess_return)) params_list = append(params_list, 
                                                  excess_return = excess_return)
  
  rmarkdown::render(
    paste0(case, ".Rmd"), params = params_list,
    output_file = paste0("../output_simulation/trade_premium/trade_premium_gap", gap, "_",
                         setting, "_miter", max_iter, 
                         ifelse(!is.null(excess_return) && excess_return,
                                "_excessreturn", ""), ".", type)
  )
}


# read sim setup
read_setup <- function(filepath = "../output_simulation/", params) {
  rds_list <- lapply(params$job_id, 
                     function(job) readRDS(paste0(filepath, 
                                                  params$job_name, "_", 
                                                  job, 
                                                  ".RDS")))
  
  # get the union of multiple sim setup
  unionFun <- function(n, obj) {
    unique(unlist(lapply(obj, `[[`, n)))
  }
  
  unions <- lapply(seq_along(rds_list[[1]]), FUN = unionFun, obj = rds_list)
  names(unions) <- names(rds_list[[1]])
  
  unions
  
}

# read multiple sim output
read_output <- function(filepath = "../output_simulation/", params) 
{
  filepaths <- unlist(lapply(params$job_id,
                             function(job_id) Sys.glob(paste0(filepath, 
                                                              params$job_name, "_", 
                                                              job_id, 
                                                              "_[!aux|!theo]*.csv"))))
  ret_list <- lapply(filepaths, fread)
  
  rbindlist(ret_list, fill = T)
}

# read aux
read_aux <- function(filepath = "../output_simulation/", params) 
{
  filepaths <- unlist(lapply(params$job_id,
                             function(job_id) Sys.glob(paste0(filepath, 
                                                              params$job_name, "_", 
                                                              job_id, 
                                                              "_aux.csv"))))
  ret_list <- lapply(filepaths, fread)
  
  rbindlist(ret_list)
}

# read aux
read_theorectical <- function(filepath = "../output_simulation/", params) 
{
  filepaths <- unlist(lapply(params$job_id,
                             function(job_id) Sys.glob(paste0(filepath, 
                                                              params$job_name, "_", 
                                                              job_id, 
                                                              "_theorectical.csv"))))
  
  ret_list <- lapply(filepaths, fread)
  
  rbindlist(ret_list)
}

# Generate stuffs that were missed during simulation
generate_sample_from_seteup <- function(filepath = "../output_simulation/", params) {
  
  ret_all <- NULL
  
  for(job_id in params$job_id) {
    
    sim_setup <- readRDS(paste0(filepath, 
                                params$job_name, "_", 
                                job_id, 
                                ".RDS"))
    
    print(sim_setup)
    
    opt <- NULL
    opt$beta0 <- "1-3-5"
    opt$power_beta <- 1
    opt$n_test <- sim_setup$n_tests
    opt$xmat <- sim_setup$xmat
    opt$observed_ratio <- 1
    
    ret <- NULL
    
    for(n_train in sim_setup$ns) {
      opt$n_train <- n_train
      for(betauu in sim_setup$betaus) {
        opt$betau <- betauu
        for(betavv in sim_setup$betavs) {
          opt$betav <- betavv
          for(d in sim_setup$d) {
            opt$d <- d
            for(epsa in sim_setup$epsas) {
              opt$epsa <- epsa
              for(epsy in sim_setup$epsys) {
                opt$epsy <- epsy
                for(ss in sim_setup$ss) {
                  opt$seed <- ss
                  
                  if(opt$power_beta) {
                    optt <- power_opt(opt, which = c("n_train", "n_test", "d", "epsa", "epsy", "betau", "betav"))  
                  } else {
                    optt <- power_opt(opt, which = c("n_train", "n_test", "d", "epsa", "epsy"))
                  }
                  
                  optt <- beta0_opt(optt)
                  beta0vec <- c(optt$beta0, optt$betau, optt$betav)
                  
                  
                  # generate fixed sample
                  sample <- generate_sample(opt = optt, beta0 = optt$beta0, seed = 0, fixed_only = T)
                  
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
                  
                  alpha <- .05
                  
                  sample_A_y <- generate_A_y(opt = optt, beta0 = optt$beta0, U = U, V = V, X = X)
                  A_train <- sample_A_y$A_train
                  A <- sample_A_y$A
                  y_train <- sample_A_y$y_train
                  y_test <- sample_A_y$y_test  
                  
                  tss <- sum((y_train - mean(y_train))^2)
                  
                  tmp <- data.table(epsy = optt$epsy,
                                    epsa = optt$epsa, 
                                    n = optt$n_train,
                                    # n_test = optt$n_test,
                                    beta0 = paste(optt$beta0, collapse = "-"),
                                    betau = optt$betau,
                                    betav = optt$betav,
                                    d = optt$d,
                                    s = optt$seed, 
                                    tss = tss,
                                    A2 = sum(A_train^2)
                  )
                  
                  ret <- rbind(ret, tmp)
                  
                }
              }
            }
          }
        }
      }
    }
    
    ret <- rbind(ret_all, ret)
    
    fwrite(ret, paste0(filepath, params$job_name, "_", job_id, "_aux.csv"))
  }
  
  ret_all
}



# Generate stuffs that were missed during simulation
theorectical_property <- function(filepath = "../output_simulation/", params) {
  
  ret_all <- NULL
  
  for(job_id in params$job_id) {
    
    sim_setup <- readRDS(paste0(filepath, 
                                params$job_name, "_", 
                                job_id, 
                                ".RDS"))
    
    print(sim_setup)
    
    outpath <- paste0(filepath, params$job_name, "_", job_id, "_theorectical.csv")
    if(file.exists(outpath)) next
    
    opt <- NULL
    opt$beta0 <- "1-3-5"
    opt$power_beta <- 1
    opt$n_test <- sim_setup$n_tests
    opt$xmat <- sim_setup$xmat
    opt$observed_ratio <- 1
    
    ret <- NULL
    
    for(n_train in sim_setup$ns) {
      opt$n_train <- n_train
      for(betauu in sim_setup$betaus) {
        opt$betau <- betauu
        for(betavv in sim_setup$betavs) {
          opt$betav <- betavv
          for(d in sim_setup$d) {
            opt$d <- d
            for(epsa in sim_setup$epsas) {
              opt$epsa <- epsa
              for(epsy in sim_setup$epsys) {
                opt$epsy <- epsy
                
                if(opt$power_beta) {
                  optt <- power_opt(opt, which = c("n_train", "n_test", "d", "epsa", "epsy", "betau", "betav"))  
                } else {
                  optt <- power_opt(opt, which = c("n_train", "n_test", "d", "epsa", "epsy"))
                }
                
                optt <- beta0_opt(optt)
                beta0vec <- c(optt$beta0, optt$betau, optt$betav)
                p <- length(beta0vec)
                
                # generate fixed sample
                sample <- generate_sample(opt = optt, beta0 = optt$beta0, seed = 0, fixed_only = T)
                
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
                
                # rate for ts
                lu_ts <- lv_ts <- optt$epsa^2*(optt$n_train-1)/optt$d^2/optt$n_train^2
                pu_ts <- pv_ts <- lu_ts*(1-lu_ts/4)
                betauv <- rate_betauv_two_stage(X = X_train, u = U_train, v = V_train, beta0 = beta0vec, 
                                      d = optt$d, sigmay2 = optt$epsy^2, sigmaa2 = optt$epsa^2, n = optt$n_train)
                betax <- rate_betax_two_stage(X = X_train, u = U_train, v = V_train, beta0 = beta0vec, 
                                     d = optt$d, sigmay2 = optt$epsy^2, sigmaa2 = optt$epsa^2, n = optt$n_train)
                lA_ts <- optt$epsa^2*(2*optt$n_train-1)/optt$d^2/optt$n_train^2
                MSE_ts <- (optt$n_train-p)/optt$n_train*
                  (optt$epsy^2 + 
                     optt$epsa^2/optt$d^2/optt$n_train*(optt$betau^2 + optt$betav^2) )
                
                # delta 
                l <- l_optimal_lr(optt$d, beta0vec, optt$epsy^2, optt$epsa^2, optt$n_train)
                delta <- (2*l*optt$d^2 + optt$betau^2 + optt$betav^2)*optt$epsa^2/optt$d^2/optt$n_train - optt$epsy^2
                delta <- delta / (l*optt$d^2 + optt$betau^2 + optt$betav^2)^2
                  
                # rate for supercent
                n_scale <- (optt$n_train - p)/optt$n_train
                lu_supercent <- lu_ts - n_scale * optt$betau^2 * delta
                lv_supercent <- lv_ts - n_scale * optt$betav^2 * delta
                pu_supercent <- lu_supercent*(1-lu_supercent/4)
                pv_supercent <- lv_supercent*(1-lv_supercent/4)
                lA_supercent <- lA_ts - n_scale * (optt$betau^2+optt$betav^2) * delta
                MSE_supercent <- MSE_ts - 
                  n_scale * (optt$betau^2+optt$betav^2) * 
                  (2*l*optt$d^2 + optt$betau^2 + optt$betav^2) /
                  (l*optt$d^2 + optt$betau^2 + optt$betav^2)^2 *
                  (optt$epsy^2 + (optt$betau^2+optt$betav^2)*optt$epsa^2/optt$d^2/optt$n_train)
                
                methods <- c("two_stage", "lr_oracle", "cv_lr", "lr_approx")
                
                tmp_ts <- data.table(method = "two_stage",
                                  epsy = optt$epsy,
                                  epsa = optt$epsa,
                                  n = optt$n_train,
                                  # n_test = optt$n_test,
                                  beta0 = paste(optt$beta0, collapse = "-"),
                                  betau = optt$betau,
                                  betav = optt$betav,
                                  d = optt$d,
                                  lopt = l,
                                  u_0 = pu_ts,
                                  v_0 = pv_ts, 
                                  u_f_0 = lu_ts,
                                  v_f_0 = lv_ts,
                                  A_norm_0 = lA_ts,
                                  MSE_0 = MSE_ts,
                                  betau_0 = betauv[1]/optt$betau^2,
                                  betav_0 = betauv[2]/optt$betav^2
                )
                tmp_betax <- as.data.table(t(betax))
                names(tmp_betax) <- paste0("betax", 1:length(optt$beta0), "_0")
                tmp_ts <- cbind(tmp_ts, tmp_betax)
                
                tmp_supercent <- data.table(method = "supercent",
                                     epsy = optt$epsy,
                                     epsa = optt$epsa,
                                     n = optt$n_train,
                                     # n_test = optt$n_test,
                                     beta0 = paste(optt$beta0, collapse = "-"),
                                     betau = optt$betau,
                                     betav = optt$betav,
                                     d = optt$d,
                                     lopt = l,
                                     u_0 = pu_supercent,
                                     v_0 = pv_supercent, 
                                     u_f_0 = lu_supercent,
                                     v_f_0 = lv_supercent,
                                     A_norm_0 = lA_supercent,
                                     MSE_0 = MSE_supercent,
                                     betau_0 = betauv[1]/optt$betau^2,
                                     betav_0 = betauv[2]/optt$betav^2
                )
                tmp_supercent <- cbind(tmp_supercent, tmp_betax)
                
                ret <- rbind(ret, tmp_ts, tmp_supercent)
                
                
              }
            }
          }
        }
      }
    }
    
    ret <- rbind(ret_all, ret)
    
    fwrite(ret, outpath)
  }
  
  ret_all
}


method_labels_tbl <- function(semi = F, gray = T) {

  
  if(gray) {
    cols <- c("grey100", "grey80",
              "grey60", "grey60", 
              "grey40", "grey40", 
              "grey40", 
              "grey20", "grey0",
              "grey100", "grey100",
              "grey100", "grey100",
              "grey70", "grey10")
  } else {
    cols <- c(scales::hue_pal()(10)[c(5, 5, 
                                      7, 7, 
                                      3, 3, 3)],
              scales::hue_pal()(7)[c(6,1,7)], # ts, ols, ts oracle
              scales::hue_pal()(10)[c(5)], #lr oracle oracle
              scales::hue_pal()(7)[c(3,4)],
              scales::hue_pal()(10)[c(2, 8)])
  }
  shapes <- c(3, 1, 
              16, 16, 
              4, 4,
              4,
              17, 25,
              2, 1, 
              3, 13,
              4, 2)
  linetypes <- c("solid", "solid", 
                 "longdash", "longdash",
                 "dotted", "dotted", "dotted",
                 "dashed", "dashed",
                 "dotdash", "solid",
                 "dashed", "longdash",
                 "solid", "dotted")
  
  methods <- c("oracle", "lr_oracle", 
               "lr_cv", "cv_lr", 
               "lr", "lr_approx", 
               "lr_plugin",
               "two_stage", "ols",
               "two_stage_oracle", "lr_oracle_oracle", 
               "lr_plugin_oracle", "lr_cv_oracle",
               "lr_oracle_sub", "two_stage_sub")
  
  
  labels_preTeX <-  c("Oracle", "SuperCENT-$\\lambda_0$",
               "SuperCENT-$\\widehat{\\lambda}_{cv}$", "SuperCENT-$\\widehat{\\lambda}_{cv}$",
               "SuperCENT-$\\widehat{\\lambda}_{0}$", "SuperCENT-$\\widehat{\\lambda}_{0}$",
               "SuperCENT-$\\widehat{\\lambda}_{0}$",
               "Two-stage", "Two-stage-adhoc",
               "Two-stage-oracle", "SuperCENT-$\\lambda_{0}$-oracle",
               "SuperCENT-$\\widehat{\\lambda}_{0}$-oracle", "SuperCENT-$\\widehat{\\lambda}_{cv}$-oracle",
               "SuperCENT-oracle-Sub", "Two-stage-Sub"
  )
  
  if(semi) labels_preTeX[1:5] <- paste0("Semi-", labels_preTeX[1:5])
  labels <- TeX(labels_preTeX)
  
  ret_tbl <- data.table(method = methods,
                        labels = labels,
                        shapes = shapes,
                        cols = cols,
                        linetypes = linetypes)
  
  # order
  methods_order <- data.table(method = c("oracle", 
               "lr_oracle_oracle", 
               "lr_plugin_oracle", "lr_cv_oracle",
               "lr_oracle", 
               "lr_cv", "cv_lr", 
               "lr", "lr_approx", 
               "lr_plugin",
               "two_stage_oracle", 
               "ols",
               "two_stage", 
               "lr_oracle_sub", "two_stage_sub"))
  ret_tbl <- dplyr::left_join(methods_order, ret_tbl, by = "method")
  
  ret_tbl$method <- factor(ret_tbl$method, levels = methods, ordered = T)
  
  ret_tbl
}
