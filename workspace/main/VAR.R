# ================================================================
#  TER M1 Actuariat — ISFA
#  Impact of Excess-of-Loss Reinsurance on Catastrophe Risk
#  Complete R Code
# ================================================================

library(ggplot2)
library(scales)

# ================================================================
#  1. PARAMETERS
# ================================================================

set.seed(123)
n_sim  <- 10000
lambda <- 5        # Poisson frequency
mu     <- 10       # Lognormal location parameter
sigma  <- 1        # Lognormal scale parameter
d      <- 300000   # XL retention (priority)
L      <- 800000   # XL treaty limit

# ================================================================
#  2. SEVERITY DISTRIBUTIONS
# ================================================================

# -- Lognormal (baseline) ----------------------------------------
sim_lnorm <- function(n) rlnorm(n, meanlog = mu, sdlog = sigma)

# -- Pareto (heavy tail) -----------------------------------------
# Scale theta calibrated so E[Pareto] = E[Lognormal]
# E[X] = theta / (alpha - 1)  =>  theta = E[X] * (alpha - 1)
alpha_p <- 2
theta_p <- exp(mu + sigma^2 / 2) * (alpha_p - 1)
sim_pareto <- function(n) theta_p / runif(n)^(1 / alpha_p)

# -- Burr XII (intermediate tail) --------------------------------
# F(x) = 1 - (1 + (x/theta)^alpha)^(-gamma)
# Simulated via quantile inversion: X = theta*(U^(-1/gamma)-1)^(1/alpha)
alpha_b <- 2
gamma_b <- 1
theta_b <- exp(mu + sigma^2 / 2) * 0.9
sim_burr <- function(n) {
  u <- runif(n)
  theta_b * (u^(-1 / gamma_b) - 1)^(1 / alpha_b)
}

# ================================================================
#  3. GENERIC SIMULATION FUNCTION
# ================================================================

simulate_losses <- function(sim_sev, n_sim, lambda, d, L) {
  S     <- numeric(n_sim)
  S_net <- numeric(n_sim)
  for (i in seq_len(n_sim)) {
    N <- rpois(1, lambda)
    if (N > 0) {
      X        <- sim_sev(N)
      X_net    <- pmin(X, d) + pmax(0, X - L)
      S[i]     <- sum(X)
      S_net[i] <- sum(X_net)
    }
  }
  VaR      <- quantile(S,     0.99)
  TVaR     <- mean(S[S > VaR])
  VaR_net  <- quantile(S_net, 0.99)
  TVaR_net <- mean(S_net[S_net > VaR_net])
  list(
    S        = S,
    S_net    = S_net,
    VaR      = VaR,
    TVaR     = TVaR,
    VaR_net  = VaR_net,
    TVaR_net = TVaR_net,
    red_VaR  = (1 - VaR_net  / VaR)  * 100,
    red_TVaR = (1 - TVaR_net / TVaR) * 100
  )
}

# -- Run simulations for all three distributions -----------------
res_ln  <- simulate_losses(sim_lnorm,  n_sim, lambda, d, L)
res_bur <- simulate_losses(sim_burr,   n_sim, lambda, d, L)
res_par <- simulate_losses(sim_pareto, n_sim, lambda, d, L)

# -- Print results -----------------------------------------------
cat("===== LOGNORMAL =====\n")
cat("VaR  99% gross:", round(res_ln$VaR,      0), "\n")
cat("TVaR 99% gross:", round(res_ln$TVaR,     0), "\n")
cat("VaR  99% net  :", round(res_ln$VaR_net,  0), "\n")
cat("TVaR 99% net  :", round(res_ln$TVaR_net, 0), "\n")
cat("Red. VaR      :", round(res_ln$red_VaR,  1), "%\n")
cat("Red. TVaR     :", round(res_ln$red_TVaR, 1), "%\n\n")

cat("===== BURR XII =====\n")
cat("VaR  99% gross:", round(res_bur$VaR,      0), "\n")
cat("TVaR 99% gross:", round(res_bur$TVaR,     0), "\n")
cat("VaR  99% net  :", round(res_bur$VaR_net,  0), "\n")
cat("TVaR 99% net  :", round(res_bur$TVaR_net, 0), "\n")
cat("Red. VaR      :", round(res_bur$red_VaR,  1), "%\n")
cat("Red. TVaR     :", round(res_bur$red_TVaR, 1), "%\n\n")

cat("===== PARETO =====\n")
cat("VaR  99% gross:", round(res_par$VaR,      0), "\n")
cat("TVaR 99% gross:", round(res_par$TVaR,     0), "\n")
cat("VaR  99% net  :", round(res_par$VaR_net,  0), "\n")
cat("TVaR 99% net  :", round(res_par$TVaR_net, 0), "\n")
cat("Red. VaR      :", round(res_par$red_VaR,  1), "%\n")
cat("Red. TVaR     :", round(res_par$red_TVaR, 1), "%\n\n")

# -- Summary statistics ------------------------------------------
cat("===== SUMMARY S (gross) =====\n")
cat("Lognormal:\n") ; print(summary(res_ln$S))
cat("Burr XII:\n")  ; print(summary(res_bur$S))
cat("Pareto:\n")    ; print(summary(res_par$S))

# ================================================================
#  4. BOOTSTRAP CONFIDENCE INTERVALS (Lognormal)
# ================================================================

set.seed(42)
B <- 1000
boot_stats <- replicate(B, {
  idx    <- sample(n_sim, n_sim, replace = TRUE)
  s_b    <- res_ln$S[idx]
  var_b  <- quantile(s_b, 0.99)
  tvar_b <- mean(s_b[s_b > var_b])
  c(VaR = var_b, TVaR = tvar_b)
})

ci_var  <- quantile(boot_stats["VaR",],  c(0.025, 0.975))
ci_tvar <- quantile(boot_stats["TVaR",], c(0.025, 0.975))

cat("===== BOOTSTRAP CI (Lognormal, B=1000) =====\n")
cat("VaR  95% CI: [", round(ci_var[1],0),  ";", round(ci_var[2],0),  "]\n")
cat("TVaR 95% CI: [", round(ci_tvar[1],0), ";", round(ci_tvar[2],0), "]\n\n")

# ================================================================
#  5. CONVERGENCE ANALYSIS
# ================================================================

steps <- seq(500, n_sim, by = 500)

conv_all <- do.call(rbind, lapply(
  list(list(res_ln,  "Lognormal"),
       list(res_bur, "Burr XII"),
       list(res_par, "Pareto")),
  function(x) {
    r <- x[[1]] ; nom <- x[[2]]
    data.frame(
      n    = steps,
      VaR  = sapply(steps, function(k) quantile(r$S[1:k], 0.99)),
      TVaR = sapply(steps, function(k) {
        v <- quantile(r$S[1:k], 0.99)
        mean(r$S[1:k][r$S[1:k] > v])
      }),
      loi  = nom
    )
  }
))
conv_all$loi <- factor(conv_all$loi,
                       levels = c("Lognormal", "Burr XII", "Pareto"))

# ================================================================
#  6. SENSITIVITY ANALYSIS ON RETENTION LEVEL d
# ================================================================

d_grid <- seq(100000, 700000, by = 50000)

sens_all <- do.call(rbind, lapply(
  list(list(sim_lnorm,  "Lognormal"),
       list(sim_burr,   "Burr XII"),
       list(sim_pareto, "Pareto")),
  function(x) {
    sim_sev <- x[[1]] ; nom <- x[[2]]
    set.seed(123)
    S_gross <- numeric(n_sim)
    for (i in seq_len(n_sim)) {
      N <- rpois(1, lambda)
      if (N > 0) S_gross[i] <- sum(sim_sev(N))
    }
    TVaR_gross <- mean(S_gross[S_gross > quantile(S_gross, 0.99)])
    do.call(rbind, lapply(d_grid, function(dv) {
      set.seed(123)
      S_net <- numeric(n_sim)
      for (i in seq_len(n_sim)) {
        N <- rpois(1, lambda)
        if (N > 0) {
          X        <- sim_sev(N)
          S_net[i] <- sum(pmin(X, dv) + pmax(0, X - L))
        }
      }
      TVaR_net <- mean(S_net[S_net > quantile(S_net, 0.99)])
      premium  <- 1.10 * mean(S_gross - S_net)
      data.frame(
        d        = dv,
        loi      = nom,
        TVaR_net = TVaR_net,
        premium  = premium,
        red_TVaR = (1 - TVaR_net / TVaR_gross) * 100
      )
    }))
  }
))
sens_all$loi <- factor(sens_all$loi,
                       levels = c("Lognormal", "Burr XII", "Pareto"))

# ================================================================
#  7. FIGURES (exported as PDF for LaTeX)
# ================================================================

# -- Figure 1 : aggregate loss densities (gross) -----------------
df <- data.frame(
  losses = c(res_ln$S, res_bur$S, res_par$S),
  dist   = factor(
    rep(c("Lognormal", "Burr XII", "Pareto"), each = n_sim),
    levels = c("Lognormal", "Burr XII", "Pareto"))
)

p_densite <- ggplot(df, aes(x = losses, fill = dist, colour = dist)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  coord_cartesian(xlim = c(0, 2e6)) +
  scale_x_continuous(labels = comma) +
  scale_fill_manual(values = c("#185FA5", "#3B6D11", "#993C1D")) +
  scale_colour_manual(values = c("#185FA5", "#3B6D11", "#993C1D")) +
  labs(x = "Annual aggregate losses (euros)", y = "Density",
       fill = "Distribution", colour = "Distribution") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

pdf("fig_densite.pdf", width = 7, height = 4)
print(p_densite)
dev.off()

# -- Figure 2 : gross vs net by distribution ---------------------
df2 <- data.frame(
  losses = c(res_ln$S,  res_ln$S_net,
             res_bur$S, res_bur$S_net,
             res_par$S, res_par$S_net),
  type   = rep(c("Gross", "Net"), times = 3, each = n_sim),
  dist   = factor(
    rep(c("Lognormal", "Burr XII", "Pareto"), each = 2 * n_sim),
    levels = c("Lognormal", "Burr XII", "Pareto"))
)

p_brutnet <- ggplot(df2, aes(x = losses, colour = type, linetype = type)) +
  geom_density(linewidth = 0.7) +
  facet_wrap(~dist, scales = "free", ncol = 1) +
  coord_cartesian(xlim = c(0, 2e6)) +
  scale_x_continuous(labels = comma) +
  scale_colour_manual(values = c("Gross" = "#993C1D", "Net" = "#185FA5")) +
  scale_linetype_manual(values = c("Gross" = "solid", "Net" = "dashed")) +
  labs(x = "Annual aggregate losses (euros)", y = "Density",
       colour = "", linetype = "") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

pdf("fig_brutnet.pdf", width = 7, height = 7)
print(p_brutnet)
dev.off()

# -- Figure 3 : VaR convergence ----------------------------------
p_conv_var <- ggplot(conv_all, aes(x = n, y = VaR, colour = loi)) +
  geom_line(linewidth = 0.7) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_colour_manual(values = c("#185FA5", "#3B6D11", "#993C1D")) +
  labs(x = "Number of simulations", y = "VaR 99% (euros)",
       colour = "Distribution") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

pdf("fig_conv_var.pdf", width = 7, height = 4)
print(p_conv_var)
dev.off()

# -- Figure 4 : TVaR convergence ---------------------------------
p_conv_tvar <- ggplot(conv_all, aes(x = n, y = TVaR, colour = loi)) +
  geom_line(linewidth = 0.7) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_colour_manual(values = c("#185FA5", "#3B6D11", "#993C1D")) +
  labs(x = "Number of simulations", y = "TVaR 99% (euros)",
       colour = "Distribution") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

pdf("fig_conv_tvar.pdf", width = 7, height = 4)
print(p_conv_tvar)
dev.off()

# -- Figure 5 : sensitivity TVaR ---------------------------------
p_sens_tvar <- ggplot(sens_all, aes(x = d, y = TVaR_net, colour = loi)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 300000, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_colour_manual(values = c("#185FA5", "#3B6D11", "#993C1D")) +
  labs(x = "Retention level d (euros)", y = "Net TVaR 99% (euros)",
       colour = "Distribution") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

pdf("fig_sens_tvar.pdf", width = 7, height = 4)
print(p_sens_tvar)
dev.off()

# -- Figure 6 : sensitivity premium ------------------------------
p_sens_premium <- ggplot(sens_all, aes(x = d, y = premium, colour = loi)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 300000, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_colour_manual(values = c("#185FA5", "#3B6D11", "#993C1D")) +
  labs(x = "Retention level d (euros)",
       y = "Reinsurance premium (euros)",
       colour = "Distribution") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

pdf("fig_sens_premium.pdf", width = 7, height = 4)
print(p_sens_premium)
dev.off()

cat("All figures exported successfully.\n")
