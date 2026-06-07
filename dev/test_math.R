# dev/test_math.R
# Standalone tests of the PM math (no shiny dependency).

source("dev/pm_core.R")

approx <- function(x, y, tol = 1e-9) abs(x - y) < tol
ok <- function(label, cond) cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))

## --- Model 1: simple ----------------------------------------------------
r1 <- pm_simple(a = .3, b = .4, third = .2)            # direct mode
ok("M1 PM = .375",            approx(r1$pm, .375))
ok("M1 absPM = .375",         approx(r1$abspm, .375))
ok("M1 total effect = .32",   approx(r1$total, .32))
ok("M1 direct effect = .20",  approx(r1$direct, .20))

# total-input mode: supply total = .32, derive c' = .20, same PM
r1t <- pm_simple(a = .3, b = .4, third = .32, mode = "total")
ok("M1 total-mode PM = .375",        approx(r1t$pm, .375))
ok("M1 total-mode derives c' = .20", approx(r1t$direct, .20))
ok("M1 total-mode total = .32",      approx(r1t$total, .32))

# inconsistent: a=.3 b=-.4 c'=.5
r1b <- pm_simple(a = .3, b = -.4, third = .5)
ok("M1 inconsistent PM = -0.3158",   approx(r1b$pm, -.12/.38, 1e-6))
ok("M1 inconsistent absPM = 0.1935", approx(r1b$abspm, .12/.62, 1e-6))

## --- Model 2: first-stage moderated -------------------------------------
r2 <- pm_modmed(a1 = .3, a3 = .1, b = .4, third = .2, z = 1, stage = "first")
ok("M2 first-stage PM = .4444",      approx(r2$pm, .16/.36, 1e-9))
ok("M2 index = a3*b = .04",          approx(r2$index, .04))
ok("M2 conditional total = .36",     approx(r2$total, .36))
r2z0 <- pm_modmed(a1 = .3, a3 = .1, b = .4, third = .2, z = 0, stage = "first")
ok("M2 first-stage z=0 == simple",   approx(r2z0$pm, pm_simple(.3,.4,.2)$pm))

## --- Model 2: second-stage moderated ------------------------------------
r2s <- pm_modmed(a1 = .3, b1 = .4, b3 = .1, third = .2, z = 1, stage = "second")
ok("M2 second-stage PM = .4286",     approx(r2s$pm, .15/.35, 1e-9))
ok("M2 second-stage index = .03",    approx(r2s$index, .03))

## --- Model 3: parallel, 2 mediators (back-compat) -----------------------
r3i <- pm_parallel(a = list(.3,.2), b = list(.4,.3), third = .2, which = 1)
ok("M3 isolate M1 PM = .3158",       approx(r3i$pm, .12/.38, 1e-9))
r3t <- pm_parallel(a = list(.3,.2), b = list(.4,.3), third = .2)        # total
ok("M3 total PM = .4737",            approx(r3t$pm, .18/.38, 1e-9))
ok("M3 total effect = .38",          approx(r3t$total, .38))

## --- Model 3: parallel, 5 mediators -------------------------------------
av <- list(.3,.2,.1,.15,.05); bv <- list(.4,.3,.5,.2,.6)
# terms: .12,.06,.05,.03,.03 = sum .29 ; total = .29+.2 = .49
r3n <- pm_parallel(a = av, b = bv, third = .2)
ok("M3 5-med total mediated = .29",  approx(r3n$full_med, .29, 1e-9))
ok("M3 5-med total effect = .49",    approx(r3n$total, .49, 1e-9))
ok("M3 5-med total PM = .5918",      approx(r3n$pm, .29/.49, 1e-9))
r3n3 <- pm_parallel(a = av, b = bv, third = .2, which = 3)   # isolate M3 = .05
ok("M3 5-med isolate M3 PM = .1020", approx(r3n3$pm, .05/.49, 1e-9))

## --- Model 4: serial ----------------------------------------------------
# a1a3b2=.0225 a1b1=.12 a2b2=.06 denom=.4025
r4m2 <- pm_serial(a1=.3,a2=.2,a3=.25,b1=.4,b2=.3,third=.2, goal="m2")
ok("M4 isolate M2 PM = .2050",       approx(r4m2$pm, .0825/.4025, 1e-9))
r4m1 <- pm_serial(a1=.3,a2=.2,a3=.25,b1=.4,b2=.3,third=.2, goal="m1")
ok("M4 isolate M1 PM = .3540",       approx(r4m1$pm, .1425/.4025, 1e-9))
r4t <- pm_serial(a1=.3,a2=.2,a3=.25,b1=.4,b2=.3,third=.2, goal="total")
ok("M4 total PM = .5031",            approx(r4t$pm, .2025/.4025, 1e-9))
ok("M4 total effect = .4025",        approx(r4t$total, .4025, 1e-9))
ok("M4 M1+M2 isolations overlap (> total num)",
   (r4m1$num + r4m2$num) > r4t$num)   # serial path double-counted

## --- Monte Carlo CI sanity ---------------------------------------------
set.seed(1)
n  <- 20000
A  <- rnorm(n, .3, .05); B <- rnorm(n, .4, .05); CP <- rnorm(n, .2, .05)
mc <- mc_ci(pm_simple(A, B, CP)$pm, level = .95)
ok("MC CI brackets point .375", mc[1] < .375 && mc[2] > .375)
cat(sprintf("    MC 95%% CI = [%.3f, %.3f]\n", mc[1], mc[2]))

set.seed(2)
pr <- draw_pair(50000, 0, 1, 0, 1, r = 0.6)
ok("draw_pair correlation ~ 0.6", abs(cor(pr$x, pr$y) - 0.6) < 0.02)

cat("\nDone.\n")
