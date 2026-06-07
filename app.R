# app.R
# Proportion of Total Effect Mediated (PM) and Absolute PM, with Monte Carlo
# confidence intervals, for simple, moderated, parallel, and serial mediation.
#
# Models (MediationModels.pdf, Figure 1):
#   1. Simple        X -> M -> Y
#   2. Moderated     first- OR second-stage moderation by z
#   3. Parallel      2 to 5 mediators
#   4. Serial        M1 -> M2
#
# PM    = mediated effect / total effect (= mediated + c')
# AbsPM = |mediated terms| / (|mediated terms| + |c'|)   -> stays in [0, 1].
#
# Direct/total input modes:
#   "direct" : user gives c' (direct);  total effect c is derived (= med + c').
#   "total"  : user gives total effect c; direct effect c' is derived (= c - med).
#   Only the supplied quantity needs an SE; the other (and its CI) is derived.
#
# Deployed as a static shinylive (webR) site; `shiny` is supplied by the runtime.

library(shiny)

# =======================================================================
# Core math (vectorized: scalars -> point estimate, vectors -> MC draws)
# =======================================================================

draw_pair <- function(n, m1, s1, m2, s2, r = 0) {
  z1 <- rnorm(n); z2 <- rnorm(n)
  list(x = m1 + s1 * z1,
       y = m2 + s2 * (r * z1 + sqrt(max(0, 1 - r^2)) * z2))
}

# Equicorrelated draws (shared factor): every pair has correlation r, r in [0,1).
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

finalize <- function(num, absnum, terms, third, mode = c("direct", "total"),
                     extra = list()) {
  mode     <- match.arg(mode)
  full_med <- Reduce(`+`, terms)
  abs_full <- Reduce(`+`, lapply(terms, abs))
  if (mode == "direct") { cp <- third;    total <- full_med + cp }
  else                  { total <- third; cp    <- total - full_med }
  c(list(pm = num / total, abspm = absnum / (abs_full + abs(cp)),
         num = num, denom = total, total = total, direct = cp,
         full_med = full_med), extra)
}

pm_simple <- function(a, b, third, mode = "direct") {
  ie <- a * b
  finalize(ie, abs(ie), list(ie), third, mode, list(ie = ie))
}

pm_modmed <- function(third, z, stage = c("first", "second"),
                      a1 = NULL, a3 = NULL, b = NULL, b1 = NULL, b3 = NULL,
                      mode = "direct") {
  stage <- match.arg(stage)
  if (stage == "first") { cond <- a1 + a3 * z; ie <- cond * b; index <- a3 * b }
  else                  { cond <- b1 + b3 * z; ie <- a1 * cond; index <- a1 * b3 }
  finalize(ie, abs(ie), list(ie), third, mode,
           list(ie = ie, index = index, cond = cond))
}

pm_parallel <- function(a, b, third, which = NULL, mode = "direct") {
  if (!is.list(a)) a <- as.list(a)
  if (!is.list(b)) b <- as.list(b)
  terms <- Map(function(ai, bi) ai * bi, a, b)
  if (is.null(which)) { num <- Reduce(`+`, terms); absnum <- Reduce(`+`, lapply(terms, abs)) }
  else                { num <- terms[[which]];     absnum <- abs(terms[[which]]) }
  finalize(num, absnum, terms, third, mode, list(terms = terms))
}

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

# =======================================================================
# Path diagrams (inline SVG)
# =======================================================================

SUB <- c("₁", "₂", "₃", "₄", "₅")  # subscripts 1..5

.svg_wrap <- function(inner, h = 210)
  sprintf('<svg viewBox="0 0 560 %g" width="100%%" style="max-width:560px;font-family:sans-serif">%s</svg>', h, inner)

.box <- function(x, y, txt, w = 70, h = 40)
  sprintf('<rect x="%g" y="%g" width="%g" height="%g" rx="4" fill="#f4f6fb" stroke="#33415c"/><text x="%g" y="%g" text-anchor="middle" font-size="16">%s</text>',
          x, y, w, h, x + w/2, y + h/2 + 5, txt)

.arrow <- function(x1, y1, x2, y2, lab, lx = NULL, ly = NULL, dashed = FALSE) {
  if (is.null(lx)) lx <- (x1 + x2)/2
  if (is.null(ly)) ly <- (y1 + y2)/2 - 6
  dash <- if (dashed) ' stroke-dasharray="5,4"' else ''
  sprintf('<line x1="%g" y1="%g" x2="%g" y2="%g" stroke="#33415c" marker-end="url(#ah)"%s/><text x="%g" y="%g" text-anchor="middle" font-size="15" font-style="italic">%s</text>',
          x1, y1, x2, y2, dash, lx, ly, lab)
}

.defs <- '<defs><marker id="ah" markerWidth="9" markerHeight="9" refX="7" refY="3" orient="auto"><path d="M0,0 L7,3 L0,6 Z" fill="#33415c"/></marker></defs>'

diagram_parallel <- function(n = 2) {
  top <- 15; gap <- 58
  ys  <- top + (seq_len(n) - 1) * gap
  span <- (n - 1) * gap + 40
  mid <- top + span / 2 - 20
  h   <- top + (n - 1) * gap + 50
  parts <- c(.box(20, mid, "X"), .box(470, mid, "Y"),
             .arrow(92, mid + 18, 470, mid + 18, "c'"))
  for (i in seq_len(n)) {
    yi <- ys[i]
    parts <- c(parts,
      .box(235, yi, paste0("M", SUB[i])),
      .arrow(92, mid + 12, 235, yi + 20, paste0("a", SUB[i])),
      .arrow(305, yi + 20, 470, mid + 12, paste0("b", SUB[i])))
  }
  HTML(.svg_wrap(paste0(.defs, paste(parts, collapse = "")), h))
}

diagram_svg <- function(model, stage = "first", n = 2) {
  if (model == "parallel") return(diagram_parallel(n))
  s <- switch(model,
    "simple" = paste0(
      .box(20, 90, "X"), .box(245, 20, "M"), .box(470, 90, "Y"),
      .arrow(92, 95, 245, 50, "a"),
      .arrow(315, 50, 470, 95, "b"),
      .arrow(92, 120, 470, 120, "c'")),
    "moderated" = paste0(
      .box(20, 110, "X"), .box(245, 30, "M"), .box(470, 110, "Y"),
      .box(20, 20, "z", 50, 30),
      .arrow(70, 50, if (stage == "first") 250 else 360,
                     if (stage == "first") 60 else 130, "", dashed = TRUE),
      .arrow(92, 115, 245, 65, if (stage == "first") "a×z" else "a"),
      .arrow(315, 65, 470, 115, if (stage == "first") "b" else "b×z"),
      .arrow(92, 140, 470, 140, "c'")),
    "serial" = paste0(
      .box(20, 95, "X"), .box(190, 15, "M₁"), .box(360, 15, "M₂"), .box(470, 95, "Y"),
      .arrow(80, 100, 190, 45, "a₁", lx = 120, ly = 60),
      .arrow(260, 35, 360, 35, "a₃"),
      .arrow(90, 110, 360, 50, "a₂", lx = 250, ly = 120),
      .arrow(225, 55, 470, 105, "b₁", lx = 300, ly = 95),
      .arrow(430, 50, 478, 95, "b₂"),
      .arrow(92, 125, 470, 125, "c'"))
  )
  HTML(.svg_wrap(paste0(.defs, s)))
}

# =======================================================================
# UI
# =======================================================================

num2 <- function(id, lab, val, step = 0.01)
  numericInput(id, lab, value = val, step = step)

coefSE <- function(cid, clab, cval, sid, sval = 0.05)
  fluidRow(column(6, num2(cid, clab, cval)),
           column(6, num2(sid, paste0("SE(", clab, ")"), sval, step = 0.005)))

# the direct/total effect input pair, depending on input mode
thirdSE <- function(mode) {
  if (identical(mode, "total"))
    coefSE("tot", "Total effect", 0.32, "se_tot", 0.06)
  else
    coefSE("cp", "c'", 0.20, "se_cp")
}

ui <- fluidPage(
  tags$head(tags$style(HTML(
    ".small-note{color:#5b6470;font-size:13px} .res-table td,.res-table th{padding:4px 10px}"))),
  titlePanel("Proportion Mediated (PM) — Monte Carlo Confidence Intervals"),
  sidebarLayout(
    sidebarPanel(
      width = 5,
      selectInput("model", "Mediation model",
        c("Simple mediation"              = "simple",
          "Moderated mediation"           = "moderated",
          "Parallel (multiple) mediation" = "parallel",
          "Serial mediation"              = "serial")),

      conditionalPanel("input.model == 'moderated'",
        radioButtons("stage", "Moderated stage",
          c("First stage (z moderates a-path, X→M)"  = "first",
            "Second stage (z moderates b-path, M→Y)" = "second"))),

      conditionalPanel("input.model == 'parallel'",
        sliderInput("n_med", "Number of mediators", min = 2, max = 5, value = 2, step = 1),
        radioButtons("goal_par", "Research goal",
          c("Isolate each mediator (Mₖ individually)" = "isolate",
            "Total mediated effect (all mediators)"        = "total"))),

      conditionalPanel("input.model == 'serial'",
        radioButtons("goal_ser", "Research goal",
          c("Isolate mediator M₂ (effects through M₂)" = "isolate",
            "Total mediated effect (M₁+M₂)"            = "total"))),

      tags$hr(),
      conditionalPanel("input.model != 'moderated'",
        radioButtons("effect_mode", "Specify the direct vs. total effect",
          c("Direct effect c' (+ SE)"    = "direct",
            "Total effect c (+ SE)"      = "total"), inline = TRUE),
        helpText(class = "small-note",
          "Provide whichever you have: c' (direct) or c (total). The other is derived—its SE is not needed.")),

      tags$b("Path coefficients & standard errors"),
      uiOutput("coef_inputs"),

      conditionalPanel("input.model == 'moderated'",
        tags$hr(),
        tags$b("Moderator z"),
        fluidRow(column(6, num2("z_mean", "Mean of z", 0, step = 0.1)),
                 column(6, num2("z_sd",   "SD of z",   1, step = 0.1))),
        helpText(class = "small-note",
          "PM is evaluated at z = Mean and ±1 SD. Add a custom value below if needed."),
        checkboxInput("use_customz", "Add a custom z value", FALSE),
        conditionalPanel("input.use_customz",
          num2("customz", "Custom z", 0, step = 0.1))),

      tags$hr(),
      selectInput("level", "Confidence level",
        c("90%" = "0.90", "95%" = "0.95", "99%" = "0.99"), selected = "0.95"),

      checkboxInput("adv", "Advanced options (correlations, draws, seed)", FALSE),
      conditionalPanel("input.adv",
        uiOutput("cor_inputs"),
        num2("nsims", "Monte Carlo draws", 100000, step = 10000),
        num2("seed",  "Random seed", 123, step = 1))
    ),
    mainPanel(
      width = 7,
      uiOutput("diagram"),
      tags$hr(),
      h4("Results"),
      tableOutput("results"),
      uiOutput("stability"),
      tags$hr(),
      h4("Copyable summary"),
      verbatimTextOutput("summary_txt")
    )
  )
)

# =======================================================================
# Server
# =======================================================================

server <- function(input, output, session) {

  # effect input mode is forced to "direct" for moderated (total effect is
  # conditional on z, so a single total is not well defined there)
  emode <- reactive(if (identical(input$model, "moderated")) "direct"
                    else if (is.null(input$effect_mode)) "direct" else input$effect_mode)

  # ---- dynamic coefficient inputs --------------------------------------
  output$coef_inputs <- renderUI({
    third <- thirdSE(emode())
    switch(input$model,
      "simple" = tagList(
        coefSE("a", "a", 0.30, "se_a"),
        coefSE("b", "b", 0.40, "se_b"), third),
      "moderated" = if (identical(input$stage, "second")) tagList(
        coefSE("a",  "a",  0.30, "se_a"),
        coefSE("b1", "b1", 0.40, "se_b1"),
        coefSE("b3", "b3 (b×z)", 0.10, "se_b3"), third
      ) else tagList(
        coefSE("a1", "a1", 0.30, "se_a1"),
        coefSE("a3", "a3 (a×z)", 0.10, "se_a3"),
        coefSE("b",  "b",  0.40, "se_b"), third),
      "parallel" = {
        n <- input$n_med %||% 2
        ad <- c(.30, .20, .15, .10, .05); bd <- c(.40, .30, .25, .20, .15)
        rows <- lapply(seq_len(n), function(i) tagList(
          coefSE(paste0("a", i), paste0("a", SUB[i]), ad[i], paste0("se_a", i)),
          coefSE(paste0("b", i), paste0("b", SUB[i]), bd[i], paste0("se_b", i))))
        tagList(rows, third)
      },
      "serial" = tagList(
        coefSE("a1", "a1", 0.30, "se_a1"),
        coefSE("a2", "a2", 0.20, "se_a2"),
        coefSE("a3", "a3", 0.25, "se_a3"),
        coefSE("b1", "b1", 0.40, "se_b1"),
        coefSE("b2", "b2", 0.30, "se_b2"), third)
    )
  })

  # ---- optional same-equation correlation inputs -----------------------
  output$cor_inputs <- renderUI({
    cor1 <- function(id, lab, lo = -0.95) sliderInput(id, lab, min = lo, max = 0.95, value = 0, step = 0.05)
    switch(input$model,
      "moderated" = if (identical(input$stage, "second"))
          cor1("r_b1b3", "cor(b1, b3)  [same Y equation]")
        else
          cor1("r_a1a3", "cor(a1, a3)  [same M equation]"),
      "parallel" = cor1("r_bpar", "cor among b-paths  [same Y equation]", lo = 0),
      "serial"   = tagList(
          cor1("r_a2a3", "cor(a2, a3)  [same M2 equation]"),
          cor1("r_b1b2", "cor(b1, b2)  [same Y equation]")),
      helpText(class = "small-note", "No within-equation pairs to correlate here.")
    )
  })

  output$diagram <- renderUI(
    diagram_svg(input$model,
                if (is.null(input$stage)) "first" else input$stage,
                input$n_med %||% 2))

  g <- function(id, default = 0) { v <- input[[id]]; if (is.null(v) || is.na(v)) default else v }

  # third (direct or total) mean + SE, per current mode
  third_in <- function() if (identical(emode(), "total"))
      list(m = g("tot", .32), s = g("se_tot", .06)) else list(m = g("cp", .2), s = g("se_cp", .05))

  # ---- main computation -------------------------------------------------
  compute <- reactive({
    level <- as.numeric(input$level)
    n     <- max(2000, g("nsims", 100000))
    mode  <- emode()
    set.seed(g("seed", 123))
    th    <- third_in()
    TH    <- rnorm(n, th$m, th$s)            # direct or total draws
    model <- input$model

    if (model == "simple") {
      A <- rnorm(n, g("a",.3), g("se_a",.05)); B <- rnorm(n, g("b",.4), g("se_b",.05))
      d <- pm_simple(A, B, TH, mode)
      p <- pm_simple(g("a",.3), g("b",.4), th$m, mode)
      build_single("Simple mediation", p, d, level, mode)

    } else if (model == "moderated") {
      stage <- input$stage
      if (identical(stage, "second")) {
        A  <- rnorm(n, g("a",.3), g("se_a",.05))
        pr <- draw_pair(n, g("b1",.4), g("se_b1",.05), g("b3",.1), g("se_b3",.05), g("r_b1b3",0))
        f_d <- function(z) pm_modmed(third = TH, z = z, stage = "second", a1 = A, b1 = pr$x, b3 = pr$y, mode = mode)
        f_p <- function(z) pm_modmed(third = th$m, z = z, stage = "second", a1 = g("a",.3), b1 = g("b1",.4), b3 = g("b3",.1), mode = mode)
        ttl <- "Second-stage moderated mediation"
      } else {
        pr <- draw_pair(n, g("a1",.3), g("se_a1",.05), g("a3",.1), g("se_a3",.05), g("r_a1a3",0))
        B  <- rnorm(n, g("b",.4), g("se_b",.05))
        f_d <- function(z) pm_modmed(third = TH, z = z, stage = "first", a1 = pr$x, a3 = pr$y, b = B, mode = mode)
        f_p <- function(z) pm_modmed(third = th$m, z = z, stage = "first", a1 = g("a1",.3), a3 = g("a3",.1), b = g("b",.4), mode = mode)
        ttl <- "First-stage moderated mediation"
      }
      build_moderated(ttl, f_p, f_d, level, mode, g("z_mean",0), g("z_sd",1),
                      if (isTRUE(input$use_customz)) g("customz", NA) else NA)

    } else if (model == "parallel") {
      nm  <- input$n_med %||% 2
      ad  <- c(.30,.20,.15,.10,.05); bd <- c(.40,.30,.25,.20,.15)
      A   <- lapply(seq_len(nm), function(i) rnorm(n, g(paste0("a",i), ad[i]), g(paste0("se_a",i), .05)))
      Bm  <- sapply(seq_len(nm), function(i) g(paste0("b",i), bd[i]))
      Bs  <- sapply(seq_len(nm), function(i) g(paste0("se_b",i), .05))
      B   <- draw_equicorr(n, as.list(Bm), as.list(Bs), g("r_bpar", 0))
      ap  <- lapply(seq_len(nm), function(i) g(paste0("a",i), ad[i]))
      bp  <- lapply(seq_len(nm), function(i) g(paste0("b",i), bd[i]))
      build_parallel(p_n = nm, A = A, B = B, ap = ap, bp = bp, TH = TH, thm = th$m,
                     goal = input$goal_par, level = level, mode = mode)

    } else { # serial
      goal <- input$goal_ser
      A1  <- rnorm(n, g("a1",.3), g("se_a1",.05))
      pra <- draw_pair(n, g("a2",.2), g("se_a2",.05), g("a3",.25), g("se_a3",.05), g("r_a2a3",0))
      prb <- draw_pair(n, g("b1",.4), g("se_b1",.05), g("b2",.3),  g("se_b2",.05), g("r_b1b2",0))
      d <- pm_serial(A1, pra$x, pra$y, prb$x, prb$y, TH, goal, mode)
      p <- pm_serial(g("a1",.3), g("a2",.2), g("a3",.25), g("b1",.4), g("b2",.3), th$m, goal, mode)
      gl <- if (goal == "isolate") "isolating M₂" else "total M₁+M₂"
      build_single(paste0("Serial mediation (", gl, ")"), p, d, level, mode)
    }
  })

  # ---- row/builder helpers ---------------------------------------------
  mkrow <- function(name, est, draws, level) {
    ci <- mc_ci(draws, level)
    data.frame(Quantity = name, Estimate = est,
               `CI lower` = ci[1], `CI upper` = ci[2], check.names = FALSE)
  }
  eff_rows <- function(p, d, level, mode) {
    r <- mkrow("Total effect", p$total, d$total, level)
    if (identical(mode, "total")) r <- rbind(r, mkrow("Direct effect (derived)", p$direct, d$direct, level))
    r
  }
  pmline <- function(lbl, p, d, level) {
    pct <- paste0(round(level*100), "%")
    ci1 <- mc_ci(d$pm, level); ci2 <- mc_ci(d$abspm, level)
    sprintf("%sPM = %.3f, %s CI [%.3f, %.3f]; Absolute PM = %.3f [%.3f, %.3f]",
            lbl, p$pm, pct, ci1[1], ci1[2], p$abspm, ci2[1], ci2[2])
  }
  effline <- function(p, d, level, mode) {
    pct <- paste0(round(level*100), "%"); ci <- mc_ci(d$total, level)
    s <- sprintf("Total effect = %.3f, %s CI [%.3f, %.3f]", p$total, pct, ci[1], ci[2])
    if (identical(mode, "total")) {
      cid <- mc_ci(d$direct, level)
      s <- c(s, sprintf("Direct effect (derived) = %.3f, %s CI [%.3f, %.3f]", p$direct, pct, cid[1], cid[2]))
    }
    s
  }

  build_single <- function(title, p, d, level, mode) {
    rows <- rbind(
      mkrow("PM", p$pm, d$pm, level),
      mkrow("Absolute PM", p$abspm, d$abspm, level),
      eff_rows(p, d, level, mode))
    list(rows = rows, point = p, level = level,
         txt = c(title, pmline("", p, d, level), effline(p, d, level, mode)))
  }

  build_parallel <- function(p_n, A, B, ap, bp, TH, thm, goal, level, mode) {
    title <- sprintf("Parallel mediation, %d mediators (%s)", p_n,
                     if (goal == "isolate") "each mediator" else "total")
    rows <- NULL; txt <- title
    if (goal == "isolate") {
      for (k in seq_len(p_n)) {
        d <- pm_parallel(A, B, TH, which = k, mode = mode)
        p <- pm_parallel(ap, bp, thm, which = k, mode = mode)
        lbl <- paste0("M", SUB[k], ": ")
        rows <- rbind(rows,
          mkrow(paste0(lbl, "PM"), p$pm, d$pm, level),
          mkrow(paste0(lbl, "Absolute PM"), p$abspm, d$abspm, level))
        txt <- c(txt, pmline(lbl, p, d, level))
      }
      d0 <- pm_parallel(A, B, TH, mode = mode); p0 <- pm_parallel(ap, bp, thm, mode = mode)
    } else {
      d0 <- pm_parallel(A, B, TH, mode = mode); p0 <- pm_parallel(ap, bp, thm, mode = mode)
      rows <- rbind(mkrow("PM", p0$pm, d0$pm, level),
                    mkrow("Absolute PM", p0$abspm, d0$abspm, level))
      txt <- c(txt, pmline("", p0, d0, level))
    }
    rows <- rbind(rows, eff_rows(p0, d0, level, mode))
    list(rows = rows, point = p0, level = level, txt = c(txt, effline(p0, d0, level, mode)))
  }

  build_moderated <- function(title, f_p, f_d, level, mode, zmean, zsd, customz) {
    levs <- list(c("z = −1 SD", zmean - zsd), c("z = Mean", zmean), c("z = +1 SD", zmean + zsd))
    if (!is.na(customz)) levs <- c(levs, list(c(sprintf("z = %.3g (custom)", customz), customz)))
    rows <- NULL; txt <- title
    for (L in levs) {
      lab <- L[1]; z <- as.numeric(L[2])
      p <- f_p(z); d <- f_d(z)
      rows <- rbind(rows,
        mkrow(paste0(lab, ": PM"), p$pm, d$pm, level),
        mkrow(paste0(lab, ": Absolute PM"), p$abspm, d$abspm, level),
        mkrow(paste0(lab, ": Total effect"), p$total, d$total, level))
      txt <- c(txt, sprintf("%s (z=%.3g): %s; %s",
                            lab, z, pmline("", p, d, level), effline(p, d, level, mode)[1]))
    }
    p0 <- f_p(0); d0 <- f_d(0); pct <- paste0(round(level*100), "%")
    ci <- mc_ci(d0$index, level)
    txt <- c(txt, sprintf("Index of moderated mediation = %.3f, %s CI [%.3f, %.3f]",
                          p0$index, pct, ci[1], ci[2]))
    list(rows = rows, point = p0, level = level, txt = txt)
  }

  # ---- outputs ----------------------------------------------------------
  output$results <- renderTable(compute()$rows, digits = 3, na = "—",
                                striped = TRUE, width = "100%")
  output$summary_txt <- renderText(paste(compute()$txt, collapse = "\n"))

  output$stability <- renderUI({
    p <- compute()$point; msg <- NULL
    if (is.finite(p$pm) && (p$pm < 0 || p$pm > 1))
      msg <- "PM falls outside [0, 1]: the indirect and direct effects have opposing signs (inconsistent mediation / suppression). Ordinary PM is unstable here—rely on Absolute PM."
    else if (is.finite(p$denom) && abs(p$denom) < 0.02)
      msg <- "The total effect (denominator) is near zero, so PM and its CI are unstable. Interpret with caution and prefer Absolute PM."
    if (is.null(msg)) return(NULL)
    div(class = "small-note",
        style = "background:#fff6e5;border:1px solid #e6c47a;border-radius:6px;padding:8px 10px;margin-top:8px",
        tags$b("⚠ Note: "), msg)
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a

shinyApp(ui, server)
