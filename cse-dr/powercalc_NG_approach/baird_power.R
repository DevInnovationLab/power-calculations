# Baird, Bohren, McIntosh & Ozler (2016) randomized-saturation power functions.
# Source: PDEL "R Implementation.R" (base R, no packages).
#   https://pdel.ucsd.edu/_files/R%20Implementation.R
# One bug fixed below: in power_ind, `zeros(1,length(pi))` -> `rep(0,length(pi))`.
#
# Functions give analytic MDEs (Theorems 1-3):
#   power_pooled : pooled treatment (MDE_T), pooled spillover (MDE_S),
#                  treatment-only-with-pure-control (MDE_Tonly)
#   power_slope  : slope MDEs (how spillover changes with saturation): MDSE_T / MDSE_S
#   power_ind    : individual-saturation treatment / spillover MDEs
# Args: n=cluster size, C=#clusters, alpha, gamma=power, tau=between-cluster var,
#       sigma=within-cluster var, pi=saturations, f=fraction of clusters at each.
#
# NOTE: written for your R/RStudio (this was not executed in-session; the
# numbers in the comments were verified via the identical Python implementation).

power_pooled = function(n,C,alpha,gamma,tau,sigma,pi,f){
  t_alpha = qt(1 - alpha/2,n*C-3); t_gamma = qt(gamma,n*C-3)
  mu_ind = rep(0,length(pi)); eta_ind = rep(0,length(pi))
  for(i in 1:length(pi)){ mu_ind[i]=pi[i]*f[i]; eta_ind[i]=pi[i]^2*f[i] }
  mu  = sum(mu_ind); eta = sum(eta_ind)
  psi = 0; if(pi[1]==0) psi = f[1]
  muS = 1 - mu - psi
  etaT = (eta-mu^2)/(1-psi)-(psi/(1-psi)^2)*mu^2
  varN = tau+sigma; varCo = (n-1)*tau
  Var  = 1/(n*C)*(varCo*(1/(psi*(1-psi))+(1-psi)/(mu ^2)*etaT)+varN*(psi+mu )/(mu *psi))
  VarS = 1/(n*C)*(varCo*(1/(psi*(1-psi))+(1-psi)/(muS^2)*etaT)+varN*(psi+muS)/(muS*psi))
  MDE_T = (t_alpha+t_gamma)*Var^0.5
  MDE_S = (t_alpha+t_gamma)*VarS^0.5
  Var_T = 1/(n*C)*(varCo*(eta-mu^2)/(mu^2*(1-mu)^2)+varN/(mu*(1-mu)))
  MDE_Tonly = (t_alpha+t_gamma)*Var_T^0.5
  return(list(MDE_T=MDE_T, MDE_S=MDE_S, MDE_Tonly=MDE_Tonly))
}

power_slope = function(n,C,alpha,gamma,tau,sigma,pi,f,j,k){
  t_alpha = qt(1 - alpha/2,n*C-3); t_gamma = qt(gamma,n*C-3)
  mu_ind = rep(0,length(pi)); p_ind = rep(0,length(pi))
  for(i in 1:length(pi)){ mu_ind[i]=pi[i]*f[i]; p_ind[i]=(1-pi[i])*f[i] }
  varN = tau+sigma; varCo = (n-1)*tau
  Var_T = (varCo*(1/f[j]+1/f[k])+varN*(1/mu_ind[j]+1/mu_ind[k]))/(n*C)
  Var_S = (varCo*(1/f[j]+1/f[k])+varN*(1/p_ind[j] +1/p_ind[k]))/(n*C)
  MDSE_T = ((t_alpha+t_gamma)/(pi[k]-pi[j]))*Var_T^0.5
  MDSE_S = ((t_alpha+t_gamma)/(pi[k]-pi[j]))*Var_S^0.5
  return(list(MDSE_T=MDSE_T, MDSE_S=MDSE_S))
}

power_ind = function(n,C,alpha,gamma,tau,sigma,pi,f,p){
  t_alpha = qt(1 - alpha/2,n*C-3); t_gamma = qt(gamma,n*C-3)
  mu_ind = rep(0,length(pi)); p_ind = rep(0,length(pi))   # fixed: was zeros(1,length(pi))
  for(i in 1:length(pi)){ mu_ind[i]=pi[i]*f[i]; p_ind[i]=(1-pi[i])*f[i] }
  psi = 0; if(pi[1]==0) psi = f[1]
  varN = tau+sigma; varCo = (n-1)*tau
  MDE_ind_T = (t_alpha+t_gamma)*(1/(n*C)*(varCo*(1/f[p]+1/psi)+varN*(1/mu_ind[p]+1/psi)))^0.5
  MDE_ind_S = (t_alpha+t_gamma)*(1/(n*C)*(varCo*(1/f[p]+1/psi)+varN*(1/p_ind[p] +1/psi)))^0.5
  return(list(MDE_ind_T=MDE_ind_T, MDE_ind_S=MDE_ind_S))
}

# ============================================================================
# Wiring to our study (ESI en Valores chatbot): spillover-SLOPE MDE per outcome,
# comparing equal-thirds {25/50/75} vs extreme-weighted {20/50/80} @ 40/20/40.
# (The slope is the relevant spillover estimand when there is no pure-control
#  saturation; it is identified from variation across saturations.)
# ============================================================================
if (sys.nframe() == 0) {              # run only when sourced/executed directly

  # (a) paper-replication check (Baird et al. Table 2, col 5):
  r <- power_pooled(10,100,.05,.8,.1,.9, c(0,.25,.5,.75,1), rep(.2,5))
  s <- power_slope (10,100,.05,.8,.1,.9, c(0,.25,.5,.75,1), rep(.2,5), 2, 4)
  cat(sprintf("paper check  MDE_T=%.4f  MDE_S=%.4f  slope MDSE_T=%.4f\n",
              r$MDE_T, r$MDE_S, s$MDSE_T))
  # expected (from validated Python run): MDE_T=0.3179  MDE_S=0.3387  slope=1.0592

  # (b) our outcomes
  n <- 90; C <- 450
  outcomes <- list(c("Pregnancy",          0.050),
                   c("Overconfidence",      0.065),
                   c("Depression/Knowledge",0.100),
                   c("Test scores",         0.200))
  cat(sprintf("\nSpillover-slope MDE (SD), 80%% power, a=0.05, %d schools x %d students:\n", C, n))
  cat(sprintf("%-22s %8s %14s %16s\n","Outcome","ICC","equal thirds","extreme 40/20/40"))
  for (o in outcomes) {
    name <- o[1]; icc <- as.numeric(o[2]); tau <- icc; sig <- 1 - icc
    eq <- power_slope(n,C,.05,.8,tau,sig, c(.25,.50,.75), c(1/3,1/3,1/3), 1, 3)
    ex <- power_slope(n,C,.05,.8,tau,sig, c(.20,.50,.80), c(.40,.20,.40), 1, 3)
    cat(sprintf("%-22s %6.3f %12.4f %15.4f\n", name, icc, eq$MDSE_S, ex$MDSE_S))
  }
  # expected (validated via Python):
  #   Pregnancy            0.050   0.1819   0.1428
  #   Overconfidence       0.065   0.1983   0.1549
  #   Depression/Knowledge 0.100   0.2320   0.1799
  #   Test scores          0.200   0.3085   0.2373
}
