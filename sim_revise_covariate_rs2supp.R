# Adding a covariate
# TODO: Boundary avoiding EBMC

library(SimDesign)
library(bootmlm)
library(lme4)
library(blme)

DESIGNFACTOR <- createDesign(
    ave_clus_size = c(5, 25),
    num_clus = c(20, 50, 100),
    gamma0 = 1,
    gamma1 = c(0.4, 0.1),
    gamma2 = -0.2,
    gamma3 = 0.5,
    gamma4 = 0.3,
    gamma5 = 0.2,
    tau0sq = c(.10, .40),
    tauxsq = c(.05, .25, 1),
    balanced = c(TRUE, FALSE)
    # subset = !(gamma1 == 0.2 & tauxsq == 1 & tau0sq == .10)
)

DESIGNFACTOR$gamma2 <- ifelse(DESIGNFACTOR$gamma1 == 0.1,
    0.3, -0.2
)
DESIGNFACTOR$tau0sq <- with(
    DESIGNFACTOR,
    tau0sq - tauxsq * gamma2^2
)

#-------------------------------------------------------------------

# Finite X (Truncated normal)
# With lb = -1, ub = -1,
# variance = 1 + (-1 * dnorm(-1) - (1 * dnorm(1))) / (pnorm(1) - pnorm(-1))
# VAR_TRUNC_NORMAL <- 1 + (-1 * dnorm(-1) - (1 * dnorm(1))) /
#   (pnorm(1) - pnorm(-1))
sd_trunc_norm <- function(orig_sd, range = c(-2, 2)) {
    al <- range[1]
    be <- range[2]
    Z <- pnorm(be) - pnorm(al)
    orig_sd * sqrt(1 + (al * dnorm(al) - be * dnorm(be)) / Z -
        ((dnorm(al) - dnorm(be)) / Z)^2)
}

rtruncnorm <- function(n) {
    x <- NULL
    while (length(x) != n) {
        n_to_sim <- n - length(x)
        z <- runif(n_to_sim, min = -1, max = 1)
        Qz <- exp(-z^2 / 2)
        u <- runif(n_to_sim)
        x <- c(x, z[u < Qz])
    }
    x
}

generate_x <- function(j, jpop, cid,
                       range = c(-2, 2), target_icc = 0.1) {
    meanx_pop <- rtruncate(jpop, rnorm,
        range = range,
        sd = target_icc / (1 - target_icc)
    )
    # Simulate lv-1 finite sample
    x_sam <- rtruncate(length(cid), rnorm,
        range = range,
        mean = meanx_sam[cid]
    )
}

generate_mux <- function(j, jpop, target_icc = 0.2,
                         range = c(-2, 2)) {
    meanx_pop <- rtruncate(jpop, rnorm,
        range = range,
        sd = sqrt(target_icc / (1 - target_icc))
    )
    # Draw finite sample
    sample(meanx_pop, size = j)
}

Generate <- function(condition, fixed_objects = NULL) {
    ave_clus_size <- condition$ave_clus_size
    num_clus <- condition$num_clus
    num_obs <- ave_clus_size * num_clus
    gamma0 <- condition$gamma0
    gamma1 <- condition$gamma1
    gamma2 <- condition$gamma2
    gamma3 <- condition$gamma3
    gamma4 <- condition$gamma4
    gamma5 <- condition$gamma5
    tau0 <- sqrt(condition$tau0sq)
    taux <- sqrt(condition$tauxsq)
    if (condition$balanced) {
        clus_size <- rep(ave_clus_size, 5)
    } else {
        clus_size <- seq(ave_clus_size * 1 / 5,
            ave_clus_size * 9 / 5,
            length.out = 5
        )
    }
    clus_id <- rep.int(seq_len(num_clus), rep_len(clus_size, num_clus))
    # Design matrix
    # mux <- generate_mux(num_clus, 20000)
    z <- rnorm(num_clus, sd = taux)
    mux <- 0.5 + 0.3 * z + rnorm(num_clus, sd = sqrt(1 - 0.3^2) * taux)
    z <- z[clus_id]
    mux <- mux[clus_id] # expand to lv 1
    # x <- rtruncate(num_obs, rnorm, range = c(-2, 2),
    #                mean = mux)
    # xc <- x - mux
    wc <- rnorm(num_obs)
    wc <- wc - ave(wc, clus_id, FUN = mean)
    xc <- 0.5 * wc + rnorm(num_obs, sd = sqrt(1 - .5^2))
    x <- mux + xc
    # y <- gamma0 + gamma1 * xc + gamma2 * mux +
    #   rnorm(num_clus, sd = sqrt(.2 - (sd_trunc_norm(.1) * gamma2)^2))[clus_id] +
    #   rnorm(num_obs, sd = sqrt(1 - (sd_trunc_norm(1) * gamma1)^2))
    y <- gamma0 + gamma1 * xc + gamma2 * mux +
        gamma3 * z + gamma4 * wc + gamma5 * mux * wc +
        rnorm(num_clus, sd = tau0)[clus_id] +
        rnorm(num_clus, sd = sqrt(.05))[clus_id] * wc +
        rnorm(num_obs, sd = sqrt(1 - gamma1^2))
    dat <- data.frame(
        y = y, x = x, mux = mux, xc = xc, z = z,
        wc = wc, cid = clus_id
    )
    dat
}
# Test:
# test_dat <- Generate(DESIGNFACTOR[1, ])

# Helper functions for analyses
get_rel <- function(object) {
    nj <- ave(object@resp$y, object@flist[[1]], FUN = length)
    theta2x <- object@theta[1]^2
    theta2x / (theta2x + 1 / nj)
}

get_waldci_tau <- function(x) {
    vc_est <- data.frame(VarCorr(x))
    var_sig2 <- diag(vcov_vc(x, sd_cor = FALSE))[c(1, 3)]
    vc_ci <- cbind(
        "2.5 %" = vc_est$vcov[c(1, 2)] - qnorm(.975) * sqrt(var_sig2),
        "97.5 %" = vc_est$vcov[c(1, 2)] + qnorm(.975) * sqrt(var_sig2)
    )
    rownames(vc_ci) <- c(".sig01", ".sig03")
    vc_ci
}

run_mlm <- function(data) {
    fit <- lmer(y ~ x + z + x * wc + (wc | cid), data = data)
    est <- fixef(fit)[c("x", "z", "wc", "x:wc")]
    vc_est <- data.frame(VarCorr(fit))
    ci <- confint(fit, parm = c("x", "z", "wc", "x:wc", ".sig01", ".sig03"))
    c(
        est[1], ci["x", ], NA, 0, 0.1,
        est[2], ci["z", ], est[3], ci["wc", ],
        est["x:wc"], ci["x:wc", ],
        vc_est[vc_est$grp == "cid" & vc_est$var1 == "(Intercept)" &
            is.na(vc_est$var2), "vcov"], ci[".sig01", ]^2,
        vc_est[vc_est$grp == "cid" & vc_est$var1 == "wc" &
            is.na(vc_est$var2), "vcov"], ci[".sig03", ]^2
    )
}

run_cmc <- function(data) {
    fit <- lmer(y ~ x_cmc + x_cm + z + x_cm * wc + (wc | cid), data = data,
                REML = FALSE)
    est <- fixef(fit)[c("x_cmc", "x_cm", "z", "wc", "x_cm:wc")]
    vc_est <- data.frame(VarCorr(fit))
    ci <- confint(fit,
        parm = c("x_cmc", "x_cm", "z", "wc", "x_cm:wc"),
        method = "Wald"
    )
    ci <- rbind(ci, get_waldci_tau(fit))
    c(
        est["x_cmc"], ci["x_cmc", ], est["x_cm"], ci["x_cm", ],
        est["z"], ci["z", ], est["wc"], ci["wc", ],
        est["x_cm:wc"], ci["x_cm:wc", ],
        vc_est[vc_est$grp == "cid" & vc_est$var1 == "(Intercept)" &
            is.na(vc_est$var2), "vcov"], ci[".sig01", ],
        vc_est[vc_est$grp == "cid" & vc_est$var1 == "wc" &
            is.na(vc_est$var2), "vcov"], ci[".sig03", ]
    )
}

# run_acmc <- function(data) {
#   sample_mlmx <- lmer(x ~ (1 | cid), data = data)
#   data$x_cmw <- data$x_cm * sqrt(get_rel(sample_mlmx))
#   fit <- lmer(y ~ x_cmc + x_cmw + (1 | cid), data = data)
#   est <- fixef(fit)[c("x_cmc", "x_cmw")]
#   ci <- confint(fit, parm = c("x_cmc", "x_cmw", ".sig01"))
#   if (!("x_cmw" %in% rownames(ci))) ci <- rbind(ci, "x_cmw" = c(NA, NA))
#   c(est["x_cmc"], ci["x_cmc", ], est["x_cmw"], ci["x_cmw", ],
#     (fit@theta * sigma(fit))^2, ci[".sig01", ]^2)
# }

run_acmc <- function(data, objectx) {
    data$x_cmw <- (data$x_cm - mean(data$x_cm)) * get_rel(objectx) +
        mean(data$x_cm)
    data$x_cmwc <- data$x - data$x_cmw
    fit <- lmer(y ~ x_cmwc + x_cmw + z + x_cmw * wc + (wc | cid), data = data)
    est <- fixef(fit)[c("x_cmwc", "x_cmw", "z", "wc", "x_cmw:wc")]
    vc_est <- data.frame(VarCorr(fit))
    ci <- confint(fit, parm = c(
        "x_cmwc", "x_cmw", "z", "wc", "x_cmw:wc",
        ".sig01", ".sig03"
    ))
    if (!("x_cmw" %in% rownames(ci))) ci <- rbind(ci, "x_cmw" = c(NA, NA))
    c(
        est["x_cmwc"], ci["x_cmwc", ], est["x_cmw"], ci["x_cmw", ],
        est["z"], ci["z", ], est["wc"], ci["wc", ],
        est["x_cmw:wc"], ci["x_cmw:wc", ],
        vc_est[vc_est$grp == "cid" & vc_est$var1 == "(Intercept)" &
            is.na(vc_est$var2), "vcov"], ci[".sig01", ],
        vc_est[vc_est$grp == "cid" & vc_est$var1 == "wc" &
            is.na(vc_est$var2), "vcov"], ci[".sig03", ]
    )
}

# run_acmc3 <- function(data) {
#   sample_mlmx <- lmer(x ~ (1 | cid), data = data)
#   data$x_cmw <- data$x_cm * sqrt(get_rel(sample_mlmx))
#   fit <- lmer(y ~ x + x_cmw + (1 | cid), data = data)
#   est <- fixef(fit)[c("x", "x_cmw")]
#   ci <- confint(fit, parm = c("x", "x_cmw", ".sig01"))
#   if (!("x_cmw" %in% rownames(ci))) ci <- rbind(ci, "x_cmw" = c(NA, NA))
#   c(est["x"], ci["x", ], est["x_cmw"], ci["x_cmw", ],
#     (fit@theta * sigma(fit))^2, ci[".sig01", ]^2)
# }

run_ebmc <- function(data, objectx, ci_method = "profile", REML = TRUE) {
    relx_cm <- get_rel(objectx)
    data$x_ebm <- predict(objectx, re.form = NA) * (1 - relx_cm) +
        relx_cm * (data$x_cm)
    data$x_ebmc <- data$x - data$x_ebm
    fit <- lmer(y ~ x_ebmc + x_ebm + z + x_ebm * wc + (wc | cid), data = data,
                REML = REML)
    est <- fixef(fit)[c("x_ebmc", "x_ebm", "z", "wc", "x_ebm:wc")]
    vc_est <- data.frame(VarCorr(fit))
    est_tau0sq <- vc_est[vc_est$grp == "cid" & vc_est$var1 == "(Intercept)" &
        is.na(vc_est$var2), "vcov"]
    if (ci_method == "profile") {
        ci <- confint(fit,
            parm = c("x_ebmc", "x_ebm", "z", "wc", "x_ebm:wc", ".sig01", ".sig03"),
            method = ci_method
        )
        ci[c(".sig01", ".sig03"), ] <- ci[c(".sig01", ".sig03"), ]^2
    } else {
        ci <- confint(fit,
            parm = c("x_ebmc", "x_ebm", "z", "wc", "x_ebm:wc"),
            method = ci_method
        )
        ci <- rbind(ci, get_waldci_tau(fit))
    }
    tau0sq_bias <- (1 - mean(relx_cm)) * (est["x_ebm"] - est["x_ebmc"])^2 *
        (objectx@theta * sigma(objectx))^2
    tau0sq_bias <- min(tau0sq_bias, est_tau0sq)
    if (!("x_ebm" %in% rownames(ci))) ci <- rbind(ci, "x_ebm" = c(NA, NA))
    c(
        est["x_ebmc"], ci["x_ebmc", ], est["x_ebm"], ci["x_ebm", ],
        est["z"], ci["z", ], est["wc"], ci["wc", ],
        est["x_ebm:wc"], ci["x_ebm:wc", ],
        est_tau0sq - tau0sq_bias, ci[".sig01", ] - tau0sq_bias,
        vc_est[vc_est$grp == "cid" & vc_est$var1 == "wc" &
            is.na(vc_est$var2), "vcov"], ci[".sig03", ]
    )
}

run_lmc_mplus <- function(data, rep_num, npar = 10,
                          sel_par = c(1, 6, 7, 3, 5, 10, 8),
                          clean = TRUE) {
    dat_name <- sprintf("sim%i.dat", rep_num)
    write.table(data, file.path("sim3_mplus", dat_name),
        row.names = FALSE, col.names = FALSE
    )
    inp_name <- sprintf("lmc_cov%i.inp", rep_num)
    savedata_name <- sprintf("res_lmc%i.dat", rep_num)
    inp <- paste(c(
        "TITLE:  Covariate Model",
        sprintf("DATA:   FILE = %s;", dat_name),
        "VARIABLE:
        NAMES = y x mux xc z wc clus;
        USEVAR = y x z wc;
        WITHIN = wc;
        BETWEEN = z;
        CLUSTER = clus;
!DEFINE: CENTER x (GROUPMEAN);
ANALYSIS:
        TYPE = TWOLEVEL RANDOM;
        ESTIMATOR = MLR;
MODEL:  %WITHIN%
        y ON x;
        s | y ON wc;
        %BETWEEN%
        y ON x z;
        s ON x;
        y WITH s;
        !x WITH s;
SAVEDATA:",
        sprintf("        RESULTS = %s;", savedata_name)
    ), collapse = "\n")
    writeLines(inp, file.path("sim3_mplus", inp_name))
    system(
        sprintf(
            "cd sim3_mplus && mplus %s",
            inp_name
        ),
        ignore.stdout = TRUE
    )
    pars <- scan(file.path("sim3_mplus", savedata_name),
        quiet = TRUE
    )
    pars <- matrix(pars[c(sel_par, sel_par + npar)], ncol = 2)
    pars <- cbind(
        pars[, 1], pars[, 1] - qnorm(.975) * pars[, 2],
        pars[, 1] + qnorm(.975) * pars[, 2]
    )
    if (clean) {
        file.remove(
            file.path("sim3_mplus", dat_name),
            file.path("sim3_mplus", inp_name),
            file.path("sim3_mplus", savedata_name),
            file.path("sim3_mplus", sprintf("lmc_cov%i.out", rep_num))
        )
    }
    c(t(pars))
}

run_lmc_bayes <- function(data, rep_num,
                          sel_par = c(1, 8, 9, 4, 7, 12, 10),
                          clean = TRUE) {
    dat_name <- sprintf("sim%i.dat", rep_num)
    write.table(data, file.path("sim3_mplus", dat_name),
        row.names = FALSE, col.names = FALSE
    )
    inp_name <- sprintf("lmc_cov_bayes%i.inp", rep_num)
    savedata_name <- sprintf("draws%i.dat", rep_num)
    inp <- paste(c(
        "TITLE:  Covariate Model (MCMC)",
        sprintf("DATA:   FILE = %s;", dat_name),
        "VARIABLE:
        NAMES = y x mux xc z wc clus;
        USEVAR = y x z wc;
        WITHIN = wc;
        BETWEEN = z;
        CLUSTER = clus;
!DEFINE: CENTER x (GROUPMEAN);
ANALYSIS:
        TYPE = TWOLEVEL RANDOM;
        ESTIMATOR = BAYES;
MODEL:  %WITHIN%
        y ON x;
        s | y ON wc;
        %BETWEEN%
        y ON x z;
        s ON x;
        y WITH s;
        !x WITH s;
SAVEDATA:",
        sprintf("        BPARAMETERS = %s;", savedata_name)
    ), collapse = "\n")
    writeLines(inp, file.path("sim3_mplus", inp_name))
    system(
        sprintf(
            "cd sim3_mplus && mplus %s",
            inp_name
        ),
        ignore.stdout = TRUE
    )
    pars <- read.table(file.path("sim3_mplus", savedata_name))
    ndraws <- nrow(pars) / 4
    pars <- pars[c(
        seq_len(ndraws) + ndraws,
        seq_len(ndraws) + 3 * ndraws
    ), 2 + sel_par]
    if (clean) {
        file.remove(
            file.path("sim3_mplus", dat_name),
            file.path("sim3_mplus", inp_name),
            file.path("sim3_mplus", savedata_name),
            file.path("sim3_mplus", sprintf("lmc_cov_bayes%i.out", rep_num))
        )
    }
    c(
        apply(pars,
            MARGIN = 2,
            FUN = function(x) c(median(x), quantile(x, c(.025, .975)))
        )
    )
}

Analyse_cmc <- function(condition, dat, fixed_objects = NULL) {
    stats <- c("est", "ll", "ul")
    pars <- c(
        "gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "tau0sq",
        "tau1sq"
    )
    dat$x_cm <- ave(dat$x, dat$cid)
    dat$x_cmc <- dat$x - dat$x_cm
    out <- run_cmc(dat)
    names(out) <- outer(stats, pars, paste, sep = "_")
    out
}

Analyse_ebmc <- function(condition, dat, fixed_objects = NULL,
                         ci_method = "wald",
                         REML = FALSE) {
    stats <- c("est", "ll", "ul")
    pars <- c(
        "gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "tau0sq",
        "tau1sq"
    )
    dat$x_cm <- ave(dat$x, dat$cid)
    sample_mlmx <- lmer(x ~ z + (1 | cid),
        data = dat,
        REML = REML,
        control = lmerControl(calc.derivs = FALSE)
    )
    out <- run_ebmc(dat,
        objectx = sample_mlmx
    )
    names(out) <- outer(stats, pars, paste, sep = "_")
    out
}

Analyse_ebmc_rs <- function(condition, dat, fixed_objects = NULL,
                            ci_method = "profile",
                            REML = TRUE) {
    stats <- c("est", "ll", "ul")
    pars <- c(
        "gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "tau0sq",
        "tau1sq"
    )
    dat$x_cm <- ave(dat$x, dat$cid)
    sample_mlmx <- lmer(x ~ z + (wc | cid),
        data = dat,
        REML = REML,
        control = lmerControl(calc.derivs = FALSE)
    )
    out <- run_ebmc(dat,
        objectx = sample_mlmx,
        REML = REML,
        ci_method = ci_method
    )
    names(out) <- outer(stats, pars, paste, sep = "_")
    out
}

Analyse_ebmc_rs_wald <- function(condition, dat, fixed_objects = NULL) {
    Analyse_ebmc_rs(condition,
        dat = dat, fixed_objects = fixed_objects,
        ci_method = "Wald", REML = FALSE
    )
}

# Boundary avoiding
Analyse_ebmc_ba <- function(condition, dat, fixed_objects = NULL) {
    stats <- c("est", "ll", "ul")
    pars <- c(
        "gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "tau0sq",
        "tau1sq"
    )
    dat$x_cm <- ave(dat$x, dat$cid)
    sample_mlmx <- blmer(x ~ z + (wc | cid),
        data = dat,
        cov.prior = wishart,
        control = lmerControl(calc.derivs = FALSE)
    )
    out <- run_ebmc(dat,
        objectx = sample_mlmx
    )
    names(out) <- outer(stats, pars, paste, sep = "_")
    out
}

Analyse_mplus <- function(condition, dat, fixed_objects = NULL) {
    stats <- c("est", "ll", "ul")
    pars <- c(
        "gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "tau0sq",
        "tau1sq"
    )
    out <- run_lmc_mplus(dat, rep = condition$REPLICATION)
    names(out) <- outer(stats, pars, paste, sep = "_")
    out
}

Analyse_bayes <- function(condition, dat, fixed_objects = NULL) {
    stats <- c("est", "ll", "ul")
    pars <- c(
        "gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "tau0sq",
        "tau1sq"
    )
    out <- run_lmc_bayes(dat, rep = condition$REPLICATION)
    names(out) <- outer(stats, pars, paste, sep = "_")
    out
}

Analyse <- function(condition, dat, fixed_objects = NULL) {
    # Initialize output
    methods <- c("mlm", "cmc", "acmc", "ebmc")
    stats <- c("est", "ll", "ul")
    pars <- c(
        "gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "tau0sq",
        "tau1sq"
    )
    ret <- vector("list", length(methods))
    names(ret) <- methods
    error_msgs <- list()
    # No centering
    # ret$mlm <- run_mlm(dat)
    ret$mlm <- tryCatch(
        run_mlm(dat),
        error = function(cnd) {
            error_msgs <<- append(error_msgs, "run_mlm() failed")
            NULL
        }
    )
    # group-mean centering
    dat$x_cm <- ave(dat$x, dat$cid)
    dat$x_cmc <- dat$x - dat$x_cm
    # ret$cmc <- run_cmc(dat)
    ret$cmc <- tryCatch(
        run_cmc(dat),
        error = function(cnd) {
            error_msgs <<- append(error_msgs, "run_cmc() failed")
            NULL
        }
    )
    # Adjusted group-mean centering
    # sample_mlmx <- lmer(x ~ (1 | cid), data = dat,
    #                     control = lmerControl(calc.derivs = TRUE))
    sample_mlmxcov <- lmer(x ~ z + (1 | cid),
        data = dat,
        control = lmerControl(calc.derivs = TRUE)
    )
    # ret$acmc <- run_acmc(dat, objectx = sample_mlmxcov)
    ret$acmc <- tryCatch(
        run_acmc(dat, objectx = sample_mlmxcov),
        error = function(cnd) {
            error_msgs <<- append(error_msgs, "run_acmc() failed")
            NULL
        }
    )
    # Adjusted EB-mean centering
    # ret$ebmc <- run_ebmc(dat, objectx = sample_mlmxcov)
    ret$ebmc <- tryCatch(
        run_ebmc(dat, objectx = sample_mlmxcov),
        error = function(cnd) {
            error_msgs <<- append(error_msgs, "run_ebmc() failed")
            NULL
        }
    )
    # Latent-mean centering
    # ret$lmc <- run_lmc(dat)
    # Reformat the output
    if (length(error_msgs) > 0) {
        rlang::abort(paste(error_msgs, collapse = ", "))
    }
    ret <- unlist(ret)
    names(ret) <- outer(as.vector(outer(stats, pars, paste, sep = "_")),
        methods, paste,
        sep = "_"
    )
    ret
}
# Test:
# Analyse(DESIGNFACTOR[1, ], dat = test_dat)

# Helper function for computing trimmed bias
tr_bias <- function(est, par, trim = .1, na.rm = TRUE) {
    est <- as.matrix(est)
    apply(sweep(est, 2, par, "-"), 2, mean,
        trim = trim,
        na.rm = na.rm
    )
}

# Helper function for computing trimmed RMSE
tr_rmse <- function(est, par, trim = .1) {
    est <- as.matrix(est)
    if (trim == 0) {
        emp_sd <- apply(est, 2L, sd, na.rm = TRUE)
    } else {
        emp_sd <- apply(est, 2L, mad, na.rm = TRUE)
    }
    sqrt(apply(sweep(est, 2, par, "-"), 2, mean, trim = trim, na.rm = TRUE)^2 +
        emp_sd^2)
}

tr_rmse2 <- function(est, par, trim = .1) {
    est <- as.matrix(est)
    apply(sweep(est, MARGIN = 2, STAT = par, "-")^2,
        MARGIN = 2, FUN = mean, trim = trim, na.rm = TRUE
    )
}

# Helper function for computing trimmed SE bias
tr_rsebias <- function(est_se, est, trim = .1) {
    est_se <- as.matrix(est_se)
    est <- as.matrix(est)
    est_se <- apply(est_se, 2, mean, trim = trim, na.rm = TRUE)
    if (trim == 0) {
        emp_sd <- apply(est, 2L, sd, na.rm = TRUE)
    } else {
        emp_sd <- apply(est, 2L, mad, na.rm = TRUE)
    }
    est_se / emp_sd - 1
}

# Coverage accounting for missing data
coverage <- function(cis, par, ci_width = FALSE, names = NULL,
                     trim = 0) {
    nc <- ncol(cis)
    stopifnot(ncol(cis) %% 2 == 0L)
    oddcs <- seq.int(1, nc, by = 2)
    evencs <- seq.int(2, nc, by = 2)
    if (ci_width & trim == 0) {
        ret <- colMeans(cis[, evencs, drop = FALSE] -
            cis[, oddcs, drop = FALSE], na.rm = TRUE)
    } else if (ci_width & trim == .1) {
        ret <- apply(
            cis[, evencs, drop = FALSE] -
                cis[, oddcs, drop = FALSE], 2,
            mean,
            trim = trim, na.rm = TRUE
        )
    } else {
        ret <- colMeans(cis[, oddcs, drop = FALSE] < 0 &
            cis[, evencs, drop = FALSE] > 0, na.rm = TRUE)
    }
    if (!is.null(names)) {
        names(ret) <- names
    }
    ret
}

# Proportion of outliers
prop_outliers <- function(est) {
    apply(est, 2, function(x) {
        length(boxplot.stats(x)$out) / length(x)
    })
}

Evaluate <- function(condition, results, fixed_objects = NULL) {
    methods <- c("cmc", "ebmc_rs", "lmc_mplus")
    # stats <- c("est", "ll", "ul")
    pars <- c(
        "gamma1", "gamma2", "gamma3", "gamma4", "gamma5",
        "tau0sq", "tau1sq"
    )
    pop <- c(
        condition$gamma1, condition$gamma2,
        condition$gamma3, condition$gamma4,
        condition$gamma5,
        condition$tau0sq, .05
    )
    results_est <- results[, grep("est\\_", colnames(results))]
    results_ci <- results[, grep("(ll|ul)_", colnames(results))]
    results_ci_centered <- sweep(results_ci,
        MARGIN = 2,
        STATS = rep(pop, each = 2)
    )
    ret <- c(
        bias = tr_bias(results_est, pop, trim = 0),
        rmse = tr_rmse2(results_est, pop, trim = 0),
        trimmed_bias = tr_bias(results_est, pop),
        trimmed_rmse = tr_rmse(results_est, pop),
        # coverage = SimDesign::ECR(
        #   results_ci_centered, parameter = 0,
        #   names = paste0("ci_",
        #                  outer(pars, methods, paste, sep = "_"))
        #   ),
        coverage = coverage(
            results_ci_centered,
            par = 0,
            names = t(outer(paste0(methods, ".ci"), pars, paste, sep = "_"))
        ),
        # ci_width = SimDesign::ECR(
        #   results_ci_centered, parameter = 0, CI_width = TRUE,
        #   names = paste0("ci_",
        #                  outer(pars, methods, paste, sep = "_"))
        #   )
        ciwidth = coverage(
            results_ci_centered,
            par = 0, ci_width = TRUE,
            names = t(outer(paste0(methods, ".ci"), pars, paste, sep = "_"))
        ),
        trimmed_ciwidth = coverage(
            results_ci_centered,
            par = 0, ci_width = TRUE, trim = 0.1,
            names = t(outer(paste0(methods, ".ci"), pars, paste, sep = "_"))
        )
    )
    ret
}

#-------------------------------------------------------------------

res <- runSimulation(
    design = DESIGNFACTOR,
    replications = 2000,
    generate = Generate,
    analyse = list(
        # cmc = Analyse_cmc,
        # ebmc = Analyse_ebmc,
        ebmc_rs = Analyse_ebmc_rs,
        ebmc_reml_rs = Analyse_ebmc_rs,
        ebmc_ba = Analyse_ebmc_ba,
        lmc_mplus_mlr = Analyse_mplus,
        lmc_bayes = Analyse_bayes
    ),
    summarise = Evaluate,
    packages = c("lme4", "bootmlm", "blme"),
    filename = "sim_revise_cov-rs2supp-results",
    seed = rep(23506100, nrow(DESIGNFACTOR)),
    extra_options = list(
        allow_na = TRUE,
        include_replication_index = TRUE
    ),
    parallel = TRUE,
    ncores = 20L,
    save = TRUE,
    save_results = TRUE,
    save_details = list(
        save_results_dirname = "sim_revise_cov-rs2supp-results"
    )
)
# saveRDS(res, "res-sim_trial_covariate3.RDS")
