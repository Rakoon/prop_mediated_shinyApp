# Proportion Mediated (PM) Calculator

A web app for calculating the **Proportion of Total Effect Mediated (PM)** and
**Absolute PM**, with **Monte Carlo confidence intervals**, across common
mediation designs.

**▶ Live app: https://rakoon.github.io/prop_mediated_shinyApp/**

It runs entirely in your browser — an R Shiny app compiled to WebAssembly with
[shinylive](https://posit-dev.github.io/r-shinylive/). Nothing is uploaded; no
server is involved.

## What it does

You enter the path coefficients and their standard errors from your fitted
model. The app reports, for the quantity of interest:

- **PM** = mediated effect ÷ total effect
- **Absolute PM** = |mediated terms| ÷ (|mediated terms| + |c′|), which stays in
  [0, 1] even under inconsistent mediation / suppression (opposite-sign paths)
- The **total effect** and (when relevant) the **direct effect**
- A **Monte Carlo confidence interval** at your chosen level (90 / 95 / 99%),
  obtained by simulating from the coefficients' sampling distributions

Results appear as a table plus a copyable text summary, alongside a labeled path
diagram of the selected model.

## Supported models

| Model | Description | Reported quantities |
|---|---|---|
| **Simple** | X → M → Y | PM |
| **Moderated** | First- **or** second-stage moderation by z | PM at z = Mean and ±1 SD (plus an optional custom z) and the **index of moderated mediation** |
| **Parallel** | 2–5 mediators | PM for each mediator individually, or the total mediated effect |
| **Serial** | X → M₁ → M₂ → Y | PM isolating M₂, or the total mediated effect |

### Specifying the direct vs. total effect

For the simple, parallel, and serial models you can supply **either**:

- the **direct effect** *c′* (with its SE) — the app derives the total effect, or
- the **total effect** *c* (with its SE) — the app derives the direct effect
  *c′ = c − Σ(indirect)*.

You only need the standard error of whichever one you provide; the other
quantity and its confidence interval are derived. (Moderated mediation uses the
direct-effect input only, because the total effect there is conditional on *z*.)

### Optional coefficient correlations

By default, Monte Carlo draws are independent. Under **Advanced options** you can
supply correlations for coefficients estimated in the *same* regression (e.g.
`a1` & `a3` in moderation, the b-paths in parallel/serial mediation), which
yields more accurate intervals when you have them. Cross-equation covariances
are assumed ≈ 0.

## Notes on interpretation

- PM is most informative when the indirect and direct effects share the same
  sign and the total effect is not near zero. When they don't, ordinary PM can
  fall outside [0, 1] or explode; the app flags this and you should rely on
  **Absolute PM**.
- Confidence intervals are **Monte Carlo** intervals (a.k.a. the "distribution of
  the product" / parametric-bootstrap approach), based on the point estimates and
  standard errors you enter.

## Development

The app source of truth is [`app.R`](app.R). The deployed payload is
[`app.json`](app.json) (read by the shinylive runtime in `shinylive/`), which is
generated from `app.R`.

```sh
# regenerate app.json after editing app.R
python dev/build_appjson.py

# run the (pure-R) math and diagram checks — no shiny install required
Rscript dev/test_math.R
Rscript dev/test_diagram.R
```

To preview locally, serve the repository root over HTTP and open it in a browser
(the included service worker supplies the headers webR needs):

```sh
python -m http.server 8753
# then visit http://127.0.0.1:8753/
```

The model equations follow `MediationModels.pdf` in this repository.
