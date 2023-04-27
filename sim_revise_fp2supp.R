# Last update: Use raw tau0sq estimate (upward bias)
# TODO: Replicate https://www.statmodel.com/download/CenteredMediation.pdf
#       Write to handle multiple processors

library(SimDesign)
library(lme4)
library(bootmlm)
library(blme)
# library(OpenMx)
# mxOption(key = "Calculate Hessian", value = "No")
# mxOption(key = "Standard Errors", value = "No")

DESIGNFACTOR <- createDesign(
    finite_pop_num = 1:100,
    ave_clus_size = c(5, 25),
    num_clus = c(20, 50, 100),
    gamma0 = 0,
    gamma1 = c(0.7),
    gamma2 = -0.3,
    tau0sq = c(.10, .40),
    tauxsq = c(.05, .25, 1),
    samp_frac = c(0, .2, .5),
    balanced = c(FALSE)
)

DESIGNFACTOR$gamma2 <- ifelse(DESIGNFACTOR$gamma1 == 0.2,
    0.5, -0.3
)

#-------------------------------------------------------------------

# Finite sample X (npop = 100)
XCPOP <- replicate(max(DESIGNFACTOR$finite_pop_num),
    expr = lapply(
        seq_len(max(DESIGNFACTOR$num_clus)),
        function(i) {
            xc <- rnorm(250)
        }
    ), simplify = FALSE
)

generate_x <- function(j, jpop, cid, target_icc,
                       range = c(-2, 2)) {
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

generate_mux <- function(j, jpop, target_sd,
                         range = c(-2, 2)) {
    meanx_pop <- rtruncate(jpop, rnorm,
        range = range,
        sd = target_sd
    )
    # Draw finite sample
    sample(meanx_pop, size = j)
}

Generate <- function(condition, fixed_objects = NULL) {
    ave_clus_size <- condition$ave_clus_size
    num_clus <- condition$num_clus
    num_obs <- ave_clus_size * num_clus
    samp_frac <- condition$samp_frac
    gamma0 <- condition$gamma0
    gamma1 <- condition$gamma1
    gamma2 <- condition$gamma2
    tau0 <- sqrt(condition$tau0sq)
    taux <- sqrt(condition$tauxsq)
    tau1 <- tau0 / 2
    if (condition$balanced) {
        clus_size <- rep(ave_clus_size, 5)
    } else {
        clus_size <- seq(ave_clus_size * 1 / 5,
            ave_clus_size * 9 / 5,
            length.out = 5
        )
    }
    clus_size <- rep_len(clus_size, num_clus)
    clus_id <- rep.int(seq_len(num_clus), clus_size)
    # Design matrix
    mux <- rnorm(num_clus, sd = taux)
    mux <- mux[clus_id] # expand to lv 1
    xcpop <- fixed_objects$xcpop[[condition$finite_pop_num]]
    if (samp_frac == 0) {
        xc <- rnorm(num_obs)
    } else {
        pop_clus_size <- ave_clus_size / samp_frac
        xc <- lapply(1:num_clus, function(j) {
            nj <- clus_size[j]
            xcfinite <- xcpop[[j]][seq_len(pop_clus_size)]
            xcfinite <- (xcfinite - mean(xcfinite)) /
                sd(xcfinite) * (pop_clus_size - 1) / pop_clus_size
            xcfinite[sample.int(pop_clus_size, size = nj)]
        })
        xc <- unlist(xc)
    }
    x <- mux + xc
    y <- gamma0 + gamma1 * xc + gamma2 * mux +
        rnorm(num_clus, sd = tau0)[clus_id] +
        rnorm(num_clus, sd = tau1)[clus_id] * xc +
        rnorm(num_obs, sd = sqrt(1 - gamma1^2))
    dat <- data.frame(y = y, x = x, mux = mux, xc = xc, cid = clus_id)
    dat
}
# Test:
# test_dat <- Generate(DESIGNFACTOR[3601, ], fixed_objects = list(xcpop = XCPOP))

# Helper functions for analyses
get_rel <- function(object, popsize = NULL) {
    nj <- ave(object@frame[[2]], object@frame[[2]], FUN = length)
    theta2x <- object@theta^2
    if (is.null(popsize)) {
        fpc <- 1
    } else {
        fpc <- (popsize - nj) / (popsize - 1)
    }
    theta2x / (theta2x + 1 / nj * fpc)
}

run_mlm <- function(data) {
    fit <- lmer(y ~ x + (x | cid), data = data)
    est <- fixef(fit)
    vc_est <- data.frame(VarCorr(fit))
    ci <- confint(fit, parm = c("x", ".sig01", ".sig03"))
    c(
        est["x"], ci["x", ], NA, 0, 0.1,
        vc_est[vc_est$grp == "cid" & vc_est$var1 == "(Intercept)" &
            is.na(vc_est$var2), "vcov"], ci[".sig01", ]^2,
        vc_est[vc_est$grp == "cid" & vc_est$var1 == "x" &
            is.na(vc_est$var2), "vcov"], ci[".sig03", ]^2
    )
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

run_cmc <- function(data) {
    fit <- lmer(y ~ x_cmc + x_cm + (x_cmc | cid), data = data,
                REML = FALSE)
    est <- fixef(fit)[c("x_cmc", "x_cm")]
    vc_est <- data.frame(VarCorr(fit))
    ci <- confint(fit, parm = c("x_cmc", "x_cm"), method = "Wald")
    ci <- rbind(ci, get_waldci_tau(fit))
    c(
        est["x_cmc"], ci["x_cmc", ],
        est["x_cm"], ci["x_cm", ],
        vc_est[vc_est$grp == "cid" & vc_est$var1 == "(Intercept)" &
            is.na(vc_est$var2), "vcov"], ci[".sig01", ],
        vc_est[vc_est$grp == "cid" & vc_est$var1 == "x_cmc" &
            is.na(vc_est$var2), "vcov"], ci[".sig03", ]
    )
}

run_acmc <- function(data, objectx, popsize) {
    relx_cm <- get_rel(objectx, popsize = popsize)
    data$x_cmw <- (data$x_cm - mean(data$x_cm)) * relx_cm +
        mean(data$x_cm)
    data$x_cmcw <- data$x - data$x_cmw
    fit <- lmer(y ~ x_cmcw + x_cmw + (x_cmcw | cid), data = data)
    est <- fixef(fit)[c("x_cmcw", "x_cmw")]
    vc_est <- data.frame(VarCorr(fit))
    est_tau0sq <- vc_est[vc_est$grp == "cid" & vc_est$var1 == "(Intercept)" &
        is.na(vc_est$var2), "vcov"]
    ci <- confint(fit, parm = c("x_cmcw", "x_cmw", ".sig01", ".sig03"))
    est_sigmax <- sigma(objectx)
    est_taux <- objectx@theta * est_sigmax
    nj <- ave(objectx@frame[[2]], objectx@frame[[2]], FUN = length)
    tau0sq_bias <- mean(est["x_cmw"]^2 * (1 - relx_cm^2) * est_taux^2 +
        est_sigmax^2 / nj * est["x_cmcw"]^2)
    if (!("x_cmw" %in% rownames(ci))) ci <- rbind(ci, "x_cmw" = c(NA, NA))
    c(
        est["x_cmcw"], ci["x_cmcw", ],
        est["x_cmw"], ci["x_cmw", ],
        est_tau0sq - tau0sq_bias, ci[".sig01", ]^2 - tau0sq_bias,
        vc_est[vc_est$grp == "cid" & vc_est$var1 == "x_cmcw" &
            is.na(vc_est$var2), "vcov"], ci[".sig03", ]^2
    )
}

run_ebmc <- function(data, objectx, popsize, ci_method = "profile", REML = TRUE) {
    relx_cm <- get_rel(objectx, popsize = popsize)
    data$x_ebm <- predict(objectx, re.form = NA) * (1 - relx_cm) +
        relx_cm * (data$x_cm)
    data$x_ebmc <- data$x - data$x_ebm
    fit <- lmer(y ~ x_ebmc + x_ebm + (x_ebmc | cid), data = data, REML = REML)
    est <- fixef(fit)[c("x_ebmc", "x_ebm")]
    vc_est <- data.frame(VarCorr(fit))
    est_tau0sq <- vc_est[vc_est$grp == "cid" & vc_est$var1 == "(Intercept)" &
        is.na(vc_est$var2), "vcov"]
    if (ci_method == "profile") {
        ci <- confint(fit,
            parm = c("x_ebmc", "x_ebm", ".sig01", ".sig03"),
            method = ci_method
        )
        ci[c(".sig01", ".sig03"), ] <- ci[c(".sig01", ".sig03"), ]^2
    } else {
        ci <- confint(fit,
            parm = c("x_ebmc", "x_ebm"),
            method = ci_method
        )
        ci <- rbind(ci, get_waldci_tau(fit))
    }
    if (!("x_ebm" %in% names(est))) {
        tau0sq_bias <- 0
    } else {
        tau0sq_bias <- (1 - mean(relx_cm)) * (est["x_ebm"] - est["x_ebmc"])^2 *
            (objectx@theta * sigma(objectx))^2
        tau0sq_bias <- min(tau0sq_bias, est_tau0sq)
    }
    if (!("x_ebm" %in% rownames(ci))) ci <- rbind(ci, "x_ebm" = c(NA, NA))
    c(
        est["x_ebmc"], ci["x_ebmc", ],
        est["x_ebm"], ci["x_ebm", ],
        est_tau0sq - tau0sq_bias, ci[".sig01", ] - tau0sq_bias,
        vc_est[vc_est$grp == "cid" & vc_est$var1 == "x_ebmc" &
            is.na(vc_est$var2), "vcov"], ci[".sig03", ]
    )
}

# # Define OpenMx Model
# BYCLUS <- mxModel(
#     model = "byClus", type = "RAM",
#     latentVars = c("slope", "intercept", "intercept_x"),
#     # mxData(data.frame(cid = unique(test_dat$cid)),
#     #        type = 'raw', primaryKey = 'cid'),
#     mxPath(
#         from = c("intercept", "slope", "intercept_x"),
#         arrows = 2, values = .2,
#         labels = c("tau0sq", "tau1sq", NA)
#     ),
#     mxPath(from = "intercept", to = "slope", arrows = 2, values = 0),
#     mxPath(
#         from = "intercept_x", to = "intercept", values = 1,
#         labels = "contextual"
#     )
# )

# HYBRID_MX <- mxModel(
#     model = "within", type = "RAM", BYCLUS,
#     manifestVars = c("y", "x"), # latentVars = c("xd"),
#     # mxData(test_dat[c("cid", "x", "y")], type = 'raw', sort = FALSE),
#     mxPath(from = "one", to = c("y", "x"), arrows = 1, free = TRUE),
#     # mxPath(
#     #     from = "one", to = "xd", arrows = 1, free = FALSE,
#     #     labels = "data.x"
#     # ),
#     # mxPath(
#     #     from = "xd", to = "y", arrows = 1, free = TRUE,
#     #     labels = "gamma1"
#     # ),
#     mxPath(
#         from = "x", to = "y", arrows = 1, free = TRUE,
#         labels = "gamma1"
#     ),
#     mxPath(from = c("y", "x"), arrows = 2, values = c(1, 1)),
#     mxPath(paste0("byClus.", c("intercept", "slope", "intercept_x")),
#         c("y", "y", "x"),
#         arrows = 1, free = FALSE,
#         values = c(1, NA, 1),
#         labels = c(NA, "data.x", NA), joinKey = "cid"
#     ),
#     mxAlgebra(contextual + gamma1, name = "gamma2"),
#     mxCI(c("gamma1", "gamma2", "tau0sq", "tau1sq"))
# )

# run_lmc <- function(data, hybrid_mx, startvals) {
#     hybrid_mx$byClus$data <- mxData(data.frame(cid = unique(data$cid)),
#         type = "raw", primaryKey = "cid"
#     )
#     hybrid_mx$data <- mxData(data[c("cid", "x", "y")],
#         type = "raw", sort = FALSE
#     )
#     hybrid_mx$A$values["y", "x"] <- startvals$gamma1
#     hybrid_mx$S$values["y", "y"] <- 1 - startvals$gamma1^2
#     hybrid_mx$byClus$A$values["intercept", "intercept_x"] <-
#         startvals$gamma2 - startvals$gamma1
#     # hybrid_mx$byClus$S$values["slope", "slope"] <- startvals$tau0sq / 4
#     # hybrid_mx$byClus$S$values["intercept", "intercept"] <- startvals$tau0sq
#     hybrid_mx$byClus$S$values["intercept_x", "intercept_x"] <- startvals$tauxsq
#     fit <- mxRun(hybrid_mx, intervals = TRUE, silent = TRUE)
#     c(t(fit@output$confidenceIntervals[, c("estimate", "lbound", "ubound")]))
# }

# Add Mplus code
run_lmc_mplus <- function(data, rep_num, npar = 8, sel_par = c(2, 8, 7, 5),
                          clean = TRUE) {
    dat_name <- sprintf("sim%i.dat", rep_num)
    write.table(data, file.path("sim2_mplus", dat_name),
        row.names = FALSE, col.names = FALSE
    )
    inp_name <- sprintf("lmc_ran_slp%i.inp", rep_num)
    savedata_name <- sprintf("res_lmc%i.dat", rep_num)
    inp <- paste(c(
        "TITLE:  Random Slope Model",
        sprintf("DATA:   FILE = %s;", dat_name), "
VARIABLE:
        NAMES = y x mux xc clus;
        USEVAR = y x;
        !WITHIN = xc;
        !BETWEEN = x;
        CLUSTER = clus;
!DEFINE: CENTER x (GROUPMEAN);
ANALYSIS:
        TYPE = TWOLEVEL RANDOM;
MODEL:  %WITHIN%
        s | y ON x;
        %BETWEEN%
        y ON x (gamma1);
        [s] (cont);
        y WITH s;
        !x WITH s;
MODEL CONSTRAINT:
        NEW(gamma2);
        gamma2 = gamma1 + cont;
SAVEDATA:",
        sprintf("        RESULTS = %s;", savedata_name)
    ), collapse = "\n")
    writeLines(inp, file.path("sim2_mplus", inp_name))
    system(sprintf(
        "cd sim2_mplus && mplus %s",
        inp_name
    ),
    ignore.stdout = TRUE
    )
    pars <- scan(file.path("sim2_mplus", savedata_name),
        quiet = TRUE
    )
    pars <- matrix(pars[c(sel_par, sel_par + npar)], ncol = 2)
    pars <- cbind(
        pars[, 1], pars[, 1] - qnorm(.975) * pars[, 2],
        pars[, 1] + qnorm(.975) * pars[, 2]
    )
    if (clean) {
        file.remove(
            file.path("sim2_mplus", dat_name),
            file.path("sim2_mplus", inp_name),
            file.path("sim2_mplus", savedata_name),
            file.path("sim2_mplus", sprintf("lmc_ran_slp%i.out", rep_num))
        )
    }
    c(t(pars))
}

# Add Mplus Bayes
run_lmc_bayes <- function(data, rep_num, sel_par = c(3, 6, 9, 7),
                          clean = TRUE) {
    dat_name <- sprintf("sim%i.dat", rep_num)
    write.table(data, file.path("sim2_mplus", dat_name),
        row.names = FALSE, col.names = FALSE
    )
    inp_name <- sprintf("lmc_ran_slp_bayes%i.inp", rep_num)
    savedata_name <- sprintf("draws%i.dat", rep_num)
    inp <- paste(c(
        "
TITLE:  Random Slope Model (MCMC)",
        sprintf("DATA:   FILE = %s;", dat_name),
        "VARIABLE:
        NAMES = y x mux xc clus;
        USEVAR = y x;
        !WITHIN = xc;
        !BETWEEN = x;
        CLUSTER = clus;
!DEFINE: CENTER x (GROUPMEAN);
ANALYSIS:
        TYPE = TWOLEVEL RANDOM;
        ESTIMATOR = BAYES;
MODEL:  %WITHIN%
        s | y ON x;
        %BETWEEN%
        y ON x (gamma1);
        [s] (cont);
        y WITH s;
        !x WITH s;
SAVEDATA:",
        sprintf("        BPARAMETERS = %s;", savedata_name)
    ), collapse = "\n")
    writeLines(inp, file.path("sim2_mplus", inp_name))
    system(sprintf(
        "cd sim2_mplus && mplus %s",
        inp_name
    ),
    ignore.stdout = TRUE
    )
    pars <- read.table(file.path("sim2_mplus", savedata_name))
    ndraws <- nrow(pars) / 4
    pars <- pars[c(
        seq_len(ndraws) + ndraws,
        seq_len(ndraws) + 3 * ndraws
    ), 2 + sel_par]
    if (clean) {
        file.remove(
            file.path("sim2_mplus", dat_name),
            file.path("sim2_mplus", inp_name),
            file.path("sim2_mplus", savedata_name),
            file.path("sim2_mplus", sprintf("lmc_ran_slp_bayes%i.out", rep_num))
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
    pars <- c("gamma1", "gamma2", "tau0sq", "tau1sq")
    dat$x_cm <- ave(dat$x, dat$cid)
    dat$x_cmc <- dat$x - dat$x_cm
    out <- run_cmc(dat)
    names(out) <- outer(stats, pars, paste, sep = "_")
    out
}

Analyse_ebmc <- function(condition, dat, fixed_objects = NULL,
                         ci_method = "profile",
                         REML = TRUE, fp = FALSE) {
    stats <- c("est", "ll", "ul")
    pars <- c("gamma1", "gamma2", "tau0sq", "tau1sq")
    dat$x_cm <- ave(dat$x, dat$cid)
    sample_mlmx <- lmer(x ~ (1 | cid),
        data = dat,
        REML = REML,
        control = lmerControl(calc.derivs = FALSE)
    )
    samp_frac <- condition$samp_frac
    if (!fp || samp_frac == 0) {
        pop_clus_size <- NULL
    } else {
        pop_clus_size <- condition$ave_clus_size / samp_frac
    }
    out <- run_ebmc(dat,
        objectx = sample_mlmx, popsize = pop_clus_size,
        ci_method = ci_method, REML = REML
    )
    names(out) <- outer(stats, pars, paste, sep = "_")
    out
}

Analyse_ebmc_fp <- function(condition, dat, fixed_objects = NULL) {
    Analyse_ebmc(condition, dat, fp = TRUE)
}

Analyse_ebmc_wald <- function(condition, dat, fixed_objects = NULL) {
    Analyse_ebmc(condition, dat,
        fixed_objects = NULL,
        ci_method = "Wald",
        REML = FALSE
    )
}

Analyse_ebmc_wald_fp <- function(condition, dat, fixed_objects = NULL) {
    Analyse_ebmc(condition, dat,
        fixed_objects = NULL,
        ci_method = "Wald",
        REML = FALSE,
        fp = TRUE
    )
}

# Boundary avoiding
Analyse_ebmc_ba <- function(condition, dat, fixed_objects = NULL,
                            ci_method = "profile") {
    stats <- c("est", "ll", "ul")
    pars <- c("gamma1", "gamma2", "tau0sq", "tau1sq")
    dat$x_cm <- ave(dat$x, dat$cid)
    sample_mlmx <- blmer(x ~ (1 | cid),
        data = dat,
        cov.prior = gamma(shape = 2, rate = 0.1, posterior.scale = "sd"),
        control = lmerControl(calc.derivs = FALSE)
    )
    samp_frac <- condition$samp_frac
    if (samp_frac == 0) {
        pop_clus_size <- NULL
    } else {
        pop_clus_size <- condition$ave_clus_size / samp_frac
    }
    out <- run_ebmc(dat,
        objectx = sample_mlmx, popsize = pop_clus_size,
        ci_method = ci_method
    )
    names(out) <- outer(stats, pars, paste, sep = "_")
    out
}

Analyse_lmc <- function(condition, dat, fixed_objects = NULL) {
    stats <- c("est", "ll", "ul")
    pars <- c("gamma1", "gamma2", "tau0sq", "tau1sq")
    hybrid_mx <- fixed_objects$hybrid_mx
    startvals <- condition[c("gamma1", "gamma2", "tau0sq", "tauxsq")]
    out <- run_lmc(dat, hybrid_mx = hybrid_mx, startvals = startvals)
    names(out) <- outer(stats, pars, paste, sep = "_")
    out
}

Analyse_mplus <- function(condition, dat, fixed_objects = NULL) {
    stats <- c("est", "ll", "ul")
    pars <- c("gamma1", "gamma2", "tau0sq", "tau1sq")
    out <- run_lmc_mplus(dat, rep = condition$REPLICATION)
    names(out) <- outer(stats, pars, paste, sep = "_")
    out
}

Analyse_bayes <- function(condition, dat, fixed_objects = NULL) {
    stats <- c("est", "ll", "ul")
    pars <- c("gamma1", "gamma2", "tau0sq", "tau1sq")
    out <- run_lmc_bayes(dat, rep = condition$REPLICATION)
    names(out) <- outer(stats, pars, paste, sep = "_")
    out
}

#-------------------------------------------------------------------

res <- runSimulation(
    design = DESIGNFACTOR,
    replications = 20,
    generate = Generate,
    analyse = list(
        # cmc = Analyse_cmc,
        ebmc_fp = Analyse_ebmc_wald_fp,
        ebmc_reml_fp = Analyse_ebmc_fp,
        # ebmc_wald = Analyse_ebmc_wald,
        ebmc_ba = Analyse_ebmc_ba,
        lmc_mplus_mlr = Analyse_mplus,
        lmc_bayes = Analyse_bayes
    ),
    summarise = NA,
    packages = c("lme4", "bootmlm", "blme"),
    fixed_objects = list(
        xcpop = XCPOP
    ),
    filename = "sim_revise_100fp2supp-results",
    seed = rep(
        23506100:23506199,
        nrow(DESIGNFACTOR) / 100
    ),
    extra_options = list(
        allow_na = TRUE,
        include_replication_index = TRUE
    ),
    parallel = TRUE,
    ncores = 20L,
    save = TRUE,
    save_results = TRUE,
    save_details = list(
        save_results_dirname = "sim_revise_100fp2supp-results_"
    )
)
