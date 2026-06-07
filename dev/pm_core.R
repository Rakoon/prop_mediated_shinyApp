# dev/pm_core.R
# Core math for Proportion of Total Effect Mediated (PM) and Absolute PM.
# Vectorized: scalars in -> point estimate; vectors (length nsims) in -> MC draws.
#
# Conventions (MediationModels.pdf, Figure 1):
#   total effect c  = (sum of all indirect terms) + c'(direct)
#   PM    = mediated effect / total effect
#   AbsPM = |mediated terms| / (|mediated terms| + |c'|)   -> stays in [0, 1]
#
# Input mode for the direct/total effect (the "third" argument):
#   mode = "direct" : `third` is c' (direct).  total = sum(indirect) + c'
#   mode = "total"  : `third` is the total effect. c' = total - sum(indirect)
# Either way the missing quantity (and its CI) is derived; you only need the SE
# of whichever one you supply.

# ---- helpers -----------------------------------------------------------

draw_pair <- function(n, m1, s1, m2, s2, r = 0) {
  z1 <- rnorm(n); z2 <- rnorm(n)
  list(x = m1 + s1 * z1,
       y = m2 + s2 * (r * z1 + sqrt(max(0, 1 - r^2)) * z2))
}

# Equicorrelated normal draws (shared-factor model): every pair of the returned
# vectors has correlation r. Used for the b-paths of a parallel model, which all
# come from the same Y equation. Valid for r in [0, 1).
draw_equicorr <- function(n, means, sds, r = 0) {
  if (r <= 0) return(Map(function(m, s) rnorm(n, m, s), means, sds))
  z0 <- rnorm(n)
  Map(function(m, s) m + s * (sqrt(r) * z0 + sqrt(1 - r) * rnorm(n)), means, sds)
}

mc_ci <- function(x, level = 0.95) {
  x <- x[is.finite(x)]
  a <- (1 - level) / 2
  unname(quantile(x, probs = c(a, 1 - a), na.rm = TRUE))
}

# Combine goal-specific numerator + the full set of indirect `terms` + the
# direct/total `third` input into PM, AbsPM, total effect, and direct effect.
finalize <- function(num, absnum, terms, third, mode = c("direct", "total"),
                     extra = list()) {
  mode     <- match.arg(mode)
  full_med <- Reduce(`+`, terms)
  abs_full <- Reduce(`+`, lapply(terms, abs))
  if (mode == "direct") { cp <- third;        total <- full_med + cp }
  else                  { total <- third;     cp    <- total - full_med }
  c(list(pm       = num / total,
         abspm    = absnum / (abs_full + abs(cp)),
         num      = num,
         denom    = total,
         total    = total,
         direct   = cp,
         full_med = full_med), extra)
}

# ---- Model 1: simple mediation ----------------------------------------
pm_simple <- function(a, b, third, mode = "direct") {
  ie <- a * b
  finalize(ie, abs(ie), list(ie), third, mode, list(ie = ie))
}

# ---- Model 2: first- or second-stage moderated mediation --------------
# first  stage: a(z) = a1 + a3*z ; IE = a(z)*b ; index = a3*b
# second stage: b(z) = b1 + b3*z ; IE = a *b(z); index = a *b3
pm_modmed <- function(third, z, stage = c("first", "second"),
                      a1 = NULL, a3 = NULL, b = NULL, b1 = NULL, b3 = NULL,
                      mode = "direct") {
  stage <- match.arg(stage)
  if (stage == "first") { cond <- a1 + a3 * z; ie <- cond * b; index <- a3 * b }
  else                  { cond <- b1 + b3 * z; ie <- a1 * cond; index <- a1 * b3 }
  finalize(ie, abs(ie), list(ie), third, mode,
           list(ie = ie, index = index, cond = cond))
}

# ---- Model 3: parallel multiple mediation (N mediators) ----------------
# a, b: equal-length lists (or vectors) of the a_i and b_i paths.
# which = NULL -> total mediated (all mediators); which = k -> isolate M_k.
pm_parallel <- function(a, b, third, which = NULL, mode = "direct") {
  if (!is.list(a)) a <- as.list(a)
  if (!is.list(b)) b <- as.list(b)
  terms <- Map(function(ai, bi) ai * bi, a, b)
  if (is.null(which)) {
    num    <- Reduce(`+`, terms)
    absnum <- Reduce(`+`, lapply(terms, abs))
  } else {
    num <- terms[[which]]; absnum <- abs(terms[[which]])
  }
  finalize(num, absnum, terms, third, mode, list(terms = terms))
}

# ---- Model 4: serial mediation (M1 -> M2) ------------------------------
pm_serial <- function(a1, a2, a3, b1, b2, third,
                      goal = c("isolate", "total"), mode = "direct") {
  goal <- match.arg(goal)
  ie_m1 <- a1 * b1; ie_m2 <- a2 * b2; ie_serial <- a1 * a3 * b2
  terms <- list(ie_serial, ie_m1, ie_m2)
  if (goal == "isolate") { num <- ie_serial + ie_m2; absnum <- abs(ie_serial) + abs(ie_m2) }
  else                   { num <- ie_serial + ie_m1 + ie_m2
                           absnum <- abs(ie_serial) + abs(ie_m1) + abs(ie_m2) }
  finalize(num, absnum, terms, third, mode,
           list(ie_m1 = ie_m1, ie_m2 = ie_m2, ie_serial = ie_serial))
}
