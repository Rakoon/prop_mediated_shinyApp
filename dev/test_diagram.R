# dev/test_diagram.R -- smoke-test the SVG diagram generators outside shiny
HTML <- function(x) x
SUB <- c("1","2","3","4","5")
src <- readLines("app.R")
i1 <- grep("^.svg_wrap", src, fixed = FALSE)[1]
ends <- grep("^# =====", src)
i2 <- min(ends[ends > i1])
eval(parse(text = paste(src[i1:(i2 - 1)], collapse = "\n")))

for (m in c("simple", "serial"))
  cat(sprintf("%-10s OK  %d chars\n", m, nchar(diagram_svg(m))))
for (s in c("first", "second"))
  cat(sprintf("moderated/%-6s OK  %d chars\n", s, nchar(diagram_svg("moderated", s))))
for (k in 2:5)
  cat(sprintf("parallel n=%d OK  %d chars\n", k, nchar(diagram_svg("parallel", n = k))))

# custom names + XML escaping
nmlist <- list(X = "Stress", Y = "Well-being", Z = "Gender",
               med = c("Coping", "Support", "M3", "M4", "M5"))
svg <- diagram_svg("serial", nm = nmlist)
cat(sprintf("serial+names: has 'Coping'=%s, has 'Stress'=%s\n",
            grepl("Coping", svg), grepl("Stress", svg)))
esc <- diagram_svg("simple", nm = list(X = "A&B", Y = "x<y"))
cat(sprintf("xml-escape: '&amp;'=%s, '&lt;'=%s\n", grepl("&amp;", esc), grepl("&lt;", esc)))
cat("all diagrams generated without error\n")
