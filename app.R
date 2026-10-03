# app.R
# Expected Counts in Contingency Tables
# Why expected counts are R x C / N, and not simply N / 4
#
# Builds on the Multiplication Rule from the Probability Rules app:
#   if two traits are independent, P(A and B) = P(A) x P(B)
#
# Example data are HYPOTHETICAL (frogs, chytrid infection, and pond habitat)
# and are simulated, not drawn from a real study.
#
# Base R only - no tidyverse / ggplot / MASS / boot dependencies.
# Avoids |>, \(x), and other R >= 4.1 syntax for older Shiny servers.

library(shiny)

# ---- shared palette ---------------------------------------------------------
blue   <- "#195190"
teal   <- "#009499"
orange <- "#E07B39"
grey   <- "gray60"

col_inf <- orange    # infected (top row)
col_un  <- teal      # not infected (bottom row)

row_lab <- c("Infected", "Not infected")
col_lab <- c("Forest", "Meadow")

# ---- small helpers ----------------------------------------------------------
info_box <- function(..., border_col = teal) {
  div(style = paste0(
    "border-left: 4px solid ", border_col, ";",
    "background-color: #f0f6ff;",
    "border-radius: 4px;",
    "padding: 10px 14px;",
    "margin-bottom: 10px;",
    "font-size: 16px; line-height: 1.5;"
  ), ...)
}

formula_box <- function(...) {
  div(style = paste0(
    "background-color: #ffffff;",
    "border: 1px solid #d0d7de;",
    "border-radius: 4px;",
    "padding: 8px 12px;",
    "margin: 8px 0;",
    "font-size: 16px;"
  ), ...)
}

# Build a labelled 2 x 2 table from the four cell counts
make_tab <- function(a, b, c, d) {
  # a = infected & forest,     b = infected & meadow
  # c = not infected & forest, d = not infected & meadow
  matrix(c(a, c, b, d), nrow = 2,
         dimnames = list(row_lab, col_lab))
}

# Expected counts under independence: row total x column total / grand total
expected_tab <- function(tab) {
  N <- sum(tab)
  outer(rowSums(tab), colSums(tab)) / N
}

# Pearson chi-square statistic (no continuity correction)
chisq_stat <- function(O, E) {
  ok <- E > 0
  sum((O[ok] - E[ok])^2 / E[ok])
}

# Count dots laid out in a tidy grid inside a rectangle
dot_grid <- function(n, xl, xr, yb, yt, col) {
  if (n <= 0 || (xr - xl) <= 0 || (yt - yb) <= 0) return(invisible())
  w  <- xr - xl
  h  <- yt - yb
  nc <- max(1, ceiling(sqrt(n * w / h)))
  nr <- ceiling(n / nc)
  xs <- xl + (seq_len(nc) - 0.5) * w / nc
  ys <- yt - (seq_len(nr) - 0.5) * h / nr
  g  <- expand.grid(x = xs, y = ys)[seq_len(n), ]
  spacing <- min(w / nc, h / nr)
  cex <- max(0.35, min(1.8, 0.7 * spacing / 0.025))
  points(g$x, g$y, pch = 19, col = adjustcolor(col, alpha.f = 0.75), cex = cex)
}

# Number label on a white backing box so it stays readable over dots
label_box <- function(x, y, lab, col, cex = 1.3) {
  w <- strwidth(lab,  cex = cex, font = 2) * 1.35
  h <- strheight(lab, cex = cex, font = 2) * 1.9
  rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2,
       col = "white", border = col, lwd = 1.5)
  text(x, y, lab, col = col, cex = cex, font = 2)
}

# ---- mosaic plot ------------------------------------------------------------
# Column widths are proportional to column totals; within each column the
# height of each coloured block is the proportion infected / not infected.
# Area of each block is therefore proportional to its count.
draw_mosaic <- function(tab, main, dots = FALSE, digits = 0,
                        indep_line = FALSE, quarters = FALSE) {
  N    <- sum(tab)
  ctot <- colSums(tab)
  wcol <- ctot / N
  gap  <- 0.03
  xl   <- c(0, wcol[1] + gap)
  xr   <- c(wcol[1], 1 + gap)
  p_inf <- ifelse(ctot > 0, tab[1, ] / ctot, 0)

  op <- par(mar = c(4.5, 8.5, 3, 1))
  on.exit(par(op))
  plot(NA, xlim = c(0, 1 + gap), ylim = c(0, 1), axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  title(main = main, cex.main = 1.3, font.main = 2, line = 1.6)

  fmt <- function(v) formatC(v, format = "f", digits = digits)
  fmt_tot <- function(v) {
    if (abs(v - round(v)) < 1e-9) formatC(v, format = "d")
    else formatC(v, format = "f", digits = 1)
  }

  for (j in 1:2) {
    if (wcol[j] <= 0) next
    split <- 1 - p_inf[j]
    # top block: infected
    rect(xl[j], split, xr[j], 1,
         col = adjustcolor(col_inf, alpha.f = 0.25), border = "white", lwd = 2)
    # bottom block: not infected
    rect(xl[j], 0, xr[j], split,
         col = adjustcolor(col_un, alpha.f = 0.25), border = "white", lwd = 2)
    if (dots) {
      pad <- 0.012
      dot_grid(tab[1, j], xl[j] + pad, xr[j] - pad, split + pad, 1 - pad, col_inf)
      dot_grid(tab[2, j], xl[j] + pad, xr[j] - pad, pad, split - pad, col_un)
    }
    xm <- (xl[j] + xr[j]) / 2
    label_box(xm, (split + 1) / 2, fmt(tab[1, j]), col_inf)
    label_box(xm, split / 2,       fmt(tab[2, j]), col_un)
    # column label underneath
    mtext(paste0(col_lab[j], "\n(", fmt_tot(ctot[j]), " frogs)"),
          side = 1, at = xm, line = 2.2, cex = 1.1, font = 2)
  }

  # row labels on the left, lined up with the first non-empty column
  jref  <- if (wcol[1] > 0) 1 else 2
  split <- 1 - p_inf[jref]
  mtext("Infected",     side = 2, at = (split + 1) / 2, las = 1,
        line = 0.5, col = col_inf, font = 2, cex = 1.1)
  mtext("Not\ninfected", side = 2, at = split / 2, las = 1,
        line = 0.5, col = col_un, font = 2, cex = 1.1)

  if (indep_line) {
    y_ind <- 1 - sum(tab[1, ]) / N
    abline(h = y_ind, lty = 2, lwd = 2.5, col = blue)
    mtext("Dashed line = where 'no link' would split the frogs",
          side = 3, line = 0.2, cex = 0.95, col = blue)
  }

  if (quarters) {
    abline(v = (1 + gap) / 2, lty = 3, lwd = 3, col = "gray30")
    abline(h = 0.5,           lty = 3, lwd = 3, col = "gray30")
    mtext("Dotted lines = the equal-split guess (N / 4)",
          side = 3, line = 0.2, cex = 0.95, col = "gray30")
  }
  box(col = "gray80")
}

# ---- natural-frequency tree -------------------------------------------------
draw_tree <- function(N, pF, pI) {
  op <- par(mar = c(0.5, 0.5, 3, 0.5))
  on.exit(par(op))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  title(main = "Splitting the frogs, one trait at a time",
        cex.main = 1.3, font.main = 2)

  nF  <- N * pF
  nM  <- N * (1 - pF)
  fmt <- function(v) {
    if (abs(v - round(v)) < 1e-9) formatC(round(v), format = "d")
    else formatC(v, format = "f", digits = 1)
  }
  fmt_p <- function(v) formatC(v, format = "f", digits = 2)

  x0 <- 0.10; x1 <- 0.42; x2 <- 0.80
  yF <- 0.73; yM <- 0.27
  leaves_y <- c(0.89, 0.60, 0.40, 0.11)

  # branch label on a white backing so it never sits on top of a line
  blab <- function(x, y, lab, col) {
    w <- strwidth(lab,  cex = 1.0, font = 2) * 1.15
    h <- strheight(lab, cex = 1.0, font = 2) * 1.8
    rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = "white", border = NA)
    text(x, y, lab, col = col, cex = 1.0, font = 2)
  }

  # branches: root to habitat (labels sit on the middle of each branch)
  segments(x0 + 0.08, 0.5, x1 - 0.08, c(yF, yM), lwd = 2, col = "gray40")
  blab((x0 + x1) / 2, (0.5 + yF) / 2, paste0("\u00d7 ", fmt_p(pF)),     "gray20")
  blab((x0 + x1) / 2, (0.5 + yM) / 2, paste0("\u00d7 ", fmt_p(1 - pF)), "gray20")

  # branches: habitat to infection status
  leaf_left <- x2 + 0.06 - 0.12
  segments(x1 + 0.08, yF, leaf_left, leaves_y[1:2], lwd = 2, col = "gray40")
  segments(x1 + 0.08, yM, leaf_left, leaves_y[3:4], lwd = 2, col = "gray40")
  bx <- (x1 + 0.08 + leaf_left) / 2
  blab(bx, (yF + leaves_y[1]) / 2, paste0("\u00d7 ", fmt_p(pI)),     col_inf)
  blab(bx, (yF + leaves_y[2]) / 2, paste0("\u00d7 ", fmt_p(1 - pI)), col_un)
  blab(bx, (yM + leaves_y[3]) / 2, paste0("\u00d7 ", fmt_p(pI)),     col_inf)
  blab(bx, (yM + leaves_y[4]) / 2, paste0("\u00d7 ", fmt_p(1 - pI)), col_un)

  # nodes
  node <- function(x, y, lab, col, bg = "white") {
    w <- 0.15; h <- 0.10
    rect(x - w / 2 - 0.005, y - h / 2, x + w / 2 + 0.005, y + h / 2,
         col = bg, border = col, lwd = 2)
    text(x, y, lab, col = col, cex = 1.0, font = 2)
  }
  node(x0, 0.5, paste0("All frogs\n", N), "gray20")
  node(x1, yF, paste0("Forest\n", fmt(nF)), "gray20")
  node(x1, yM, paste0("Meadow\n", fmt(nM)), "gray20")

  leaf_lab <- c(paste0("Infected, forest\n", fmt(nF * pI)),
                paste0("Not inf., forest\n", fmt(nF * (1 - pI))),
                paste0("Infected, meadow\n", fmt(nM * pI)),
                paste0("Not inf., meadow\n", fmt(nM * (1 - pI))))
  leaf_col <- c(col_inf, col_un, col_inf, col_un)
  for (k in 1:4) {
    w <- 0.24; h <- 0.12
    rect(x2 - w / 2 + 0.06, leaves_y[k] - h / 2, x2 + w / 2 + 0.06, leaves_y[k] + h / 2,
         col = adjustcolor(leaf_col[k], alpha.f = 0.15),
         border = leaf_col[k], lwd = 2)
    text(x2 + 0.06, leaves_y[k], leaf_lab[k], col = leaf_col[k], cex = 0.95, font = 2)
  }
}

# ---- HTML contingency table -------------------------------------------------
# working: optional 2 x 2 matrix of strings shown in small text under each cell
html_ct <- function(tab, caption, digits = 0, working = NULL,
                    margins = NULL, accent = "gray20") {
  fmt <- function(v) formatC(v, format = "f", digits = digits)
  if (is.null(margins)) {
    margins <- list(R = rowSums(tab), C = colSums(tab), N = sum(tab))
  }
  fmt_m <- function(v) {
    if (all(abs(v - round(v)) < 1e-9)) formatC(v, format = "d", big.mark = "")
    else formatC(v, format = "f", digits = 1)
  }
  td  <- "padding:6px 10px; border:1px solid #d0d7de; text-align:center;"
  th  <- paste0(td, " background-color:#f6f8fa;")
  rc  <- c(col_inf, col_un)

  cell <- function(i, j) {
    val <- paste0("<span style='font-size:18px; font-weight:700; color:", rc[i], ";'>",
                  fmt(tab[i, j]), "</span>")
    if (!is.null(working)) {
      val <- paste0(val, "<br><span style='font-size:12px; color:#555;'>",
                    working[i, j], "</span>")
    }
    paste0("<td style='", td, "'>", val, "</td>")
  }

  rows <- character(2)
  for (i in 1:2) {
    rows[i] <- paste0(
      "<tr><th style='", th, " color:", rc[i], "; text-align:left;'>", row_lab[i], "</th>",
      cell(i, 1), cell(i, 2),
      "<td style='", td, " color:#555;'>", fmt_m(margins$R[i]), "</td></tr>")
  }

  paste0(
    "<div style='margin-bottom:6px; font-weight:700; font-size:16px; color:", accent, ";'>",
    caption, "</div>",
    "<table style='border-collapse:collapse; width:100%; font-size:15px;'>",
    "<tr><th style='", th, "'></th>",
    "<th style='", th, "'>", col_lab[1], "</th>",
    "<th style='", th, "'>", col_lab[2], "</th>",
    "<th style='", th, " color:#555;'>Row total</th></tr>",
    rows[1], rows[2],
    "<tr><th style='", th, " color:#555; text-align:left;'>Column total</th>",
    "<td style='", td, " color:#555;'>", fmt_m(margins$C[1]), "</td>",
    "<td style='", td, " color:#555;'>", fmt_m(margins$C[2]), "</td>",
    "<td style='", td, " color:#555; font-weight:700;'>", fmt_m(margins$N), "</td></tr>",
    "</table>")
}

# ---- shared process controls (tabs 2 and 3 use the same layout) ------------
process_controls <- function(id) {
  tagList(
    radioButtons(paste0(id, "_proc"), "How is infection related to habitat?",
                 choices = c("No link (independent)" = "indep",
                             "Linked (infection rate differs)" = "assoc"),
                 selected = "indep"),
    sliderInput(paste0(id, "_n"), "Number of frogs sampled (N):",
                value = 100, min = 20, max = 400, step = 10),
    sliderInput(paste0(id, "_pf"), "% of frogs living in forest ponds:",
                value = 70, min = 5, max = 95, step = 5),
    conditionalPanel(
      condition = paste0("input.", id, "_proc == 'indep'"),
      sliderInput(paste0(id, "_pi"), "% infected (same in both habitats):",
                  value = 30, min = 5, max = 95, step = 5)
    ),
    conditionalPanel(
      condition = paste0("input.", id, "_proc == 'assoc'"),
      sliderInput(paste0(id, "_pif"), "% infected in forest ponds:",
                  value = 40, min = 5, max = 95, step = 5),
      sliderInput(paste0(id, "_pim"), "% infected in meadow ponds:",
                  value = 15, min = 5, max = 95, step = 5)
    )
  )
}

# ================================ UI =========================================
ui <- fluidPage(
  titlePanel("Expected Counts in Contingency Tables",
             windowTitle = "Expected Counts"),

  tags$head(tags$style(HTML(
    ".action-button { color:#fff; background-color:#569BBD; border:none; }
     .action-button:hover { color:#fff; background-color:#3E7C99; }
     .action-button:active { transform:scale(0.97); }
     .verdict p { font-size:16px; line-height:1.6; }"
  ))),

  tabsetPanel(
    type = "tabs",

    # ======= 1. WHERE EXPECTED COUNTS COME FROM ==============================
    tabPanel(
      "1. Where expected counts come from",
      br(),
      fluidRow(
        column(
          width = 4,
          wellPanel(
            p(strong("The question")),
            p("A researcher swabs frogs from forest ponds and meadow ponds and records ",
              "whether each frog is infected with a skin fungus. ",
              em("(Hypothetical example - numbers are invented for teaching.)")),
            p("If habitat has ", strong("nothing to do"), " with infection, how many frogs ",
              "should we expect in each of the four boxes of the table?"),
            p("It is tempting to say 'N / 4 in every box'. But that would only be true if ",
              "half the frogs came from each habitat ", em("and"), " half of all frogs ",
              "were infected. Move the sliders to see why."),
            hr(),
            sliderInput("t1_n", "Number of frogs sampled (N):",
                        value = 100, min = 20, max = 400, step = 10),
            sliderInput("t1_pf", "% of frogs living in forest ponds:",
                        value = 70, min = 5, max = 95, step = 5),
            sliderInput("t1_pi", "% of all frogs that are infected:",
                        value = 30, min = 5, max = 95, step = 5),
            checkboxInput("t1_quarters",
                          "Show the 'equal split' guess (N / 4) as dotted lines",
                          value = FALSE),
            hr(),
            helpText("Glenn Tattersall, PhD"),
            helpText("For use in BIOL 3P96 - Biostatistics")
          )
        ),
        column(
          width = 8,
          fluidRow(
            column(6, plotOutput("t1_tree",   height = "380px")),
            column(6, plotOutput("t1_mosaic", height = "380px"))
          ),
          fluidRow(
            column(6, wellPanel(htmlOutput("t1_table"))),
            column(6, wellPanel(htmlOutput("t1_naive")))
          ),
          info_box(
            strong("The rule behind every expected count"),
            p("If two traits are independent, the Multiplication Rule says ",
              "P(A and B) = P(A) \u00d7 P(B). Multiply by the number of frogs to turn ",
              "a probability into a count:"),
            formula_box(
              p("Expected count = N \u00d7 P(row) \u00d7 P(column)"),
              p("= N \u00d7 (Row total / N) \u00d7 (Column total / N)"),
              p("= (Row total \u00d7 Column total \u00d7 N) / (N \u00d7 N)"),
              p(strong("= Row total \u00d7 Column total / N"))
            ),
            p("One N cancels, which is why the familiar shortcut R \u00d7 C / N is the ",
              "same thing as multiplying the two probabilities.")
          ),
          wellPanel(class = "verdict", div(uiOutput("t1_verdict"), align = "justify"))
        )
      )
    ),

    # ======= 2. ONE SAMPLE ===================================================
    tabPanel(
      "2. One sample",
      br(),
      fluidRow(
        column(
          width = 4,
          wellPanel(
            p(strong("Real data are messy")),
            p("Now we actually go out and sample frogs. You decide the ",
              strong("true process"), " in nature: either infection has no link to ",
              "habitat, or infection is more common in one habitat."),
            p("The computer then catches N frogs at random. Each dot on the left ",
              "plot is one frog. The expected table is built only from the ",
              strong("totals"), " of what we caught - exactly as you would with ",
              "real data."),
            hr(),
            process_controls("t2"),
            checkboxInput("t2_quarters",
                          "Show the 'equal split' guess (N / 4) as dotted lines",
                          value = FALSE),
            actionButton("t2_new", "Resample")
          )
        ),
        column(
          width = 8,
          fluidRow(
            column(6, plotOutput("t2_obs", height = "380px")),
            column(6, plotOutput("t2_exp", height = "380px"))
          ),
          fluidRow(
            column(4, wellPanel(htmlOutput("t2_tab_obs"))),
            column(4, wellPanel(htmlOutput("t2_tab_exp"))),
            column(4, wellPanel(htmlOutput("t2_tab_naive")))
          ),
          wellPanel(class = "verdict", div(uiOutput("t2_verdict"), align = "justify"))
        )
      )
    ),

    # ======= 3. MANY SAMPLES =================================================
    tabPanel(
      "3. Many samples",
      br(),
      fluidRow(
        column(
          width = 4,
          wellPanel(
            p(strong("Which guess is right in the long run?")),
            p("A single sample bounces around by chance. So here we repeat the whole ",
              "study many times and keep track of one box: ",
              span(style = paste0("color:", col_inf, "; font-weight:700;"),
                   "infected frogs from forest ponds"), "."),
            p("We also run a chi-square test on every sample twice - once with the ",
              "correct expected counts (R \u00d7 C / N) and once with the equal-split ",
              "guess (N / 4) - and count how often each one claims that infection ",
              "is linked to habitat."),
            hr(),
            process_controls("t3"),
            sliderInput("t3_sims", "Number of repeated studies:",
                        value = 1000, min = 200, max = 3000, step = 100),
            actionButton("t3_new", "Resample")
          )
        ),
        column(
          width = 8,
          fluidRow(
            column(6, plotOutput("t3_hist", height = "380px")),
            column(6, plotOutput("t3_bars", height = "380px"))
          ),
          wellPanel(class = "verdict", div(uiOutput("t3_verdict"), align = "justify"))
        )
      )
    )
  )
)

# ================================ SERVER =====================================
server <- function(input, output, session) {

  # ---------------- shared: read the true process for a tab -----------------
  get_process <- function(id) {
    type <- input[[paste0(id, "_proc")]]
    pF   <- input[[paste0(id, "_pf")]] / 100
    if (type == "indep") {
      piF <- input[[paste0(id, "_pi")]] / 100
      piM <- piF
    } else {
      piF <- input[[paste0(id, "_pif")]] / 100
      piM <- input[[paste0(id, "_pim")]] / 100
    }
    list(type = type, N = input[[paste0(id, "_n")]], pF = pF, piF = piF, piM = piM)
  }

  # ======================= TAB 1 ============================================
  t1 <- reactive({
    N  <- input$t1_n
    pF <- input$t1_pf / 100
    pI <- input$t1_pi / 100
    E  <- outer(c(pI, 1 - pI), c(pF, 1 - pF)) * N
    dimnames(E) <- list(row_lab, col_lab)
    list(N = N, pF = pF, pI = pI, E = E)
  })

  output$t1_tree <- renderPlot({
    d <- t1()
    draw_tree(d$N, d$pF, d$pI)
  })

  output$t1_mosaic <- renderPlot({
    d <- t1()
    draw_mosaic(d$E, "Expected counts (no link)", digits = 1,
                quarters = input$t1_quarters)
  })

  output$t1_table <- renderUI({
    d <- t1()
    pr <- c(d$pI, 1 - d$pI)
    pc <- c(d$pF, 1 - d$pF)
    w  <- matrix("", 2, 2)
    for (i in 1:2) for (j in 1:2) {
      w[i, j] <- paste0(d$N, " \u00d7 ", formatC(pr[i], format = "f", digits = 2),
                        " \u00d7 ", formatC(pc[j], format = "f", digits = 2))
    }
    HTML(html_ct(d$E, "Expected counts = N \u00d7 P(row) \u00d7 P(column)",
                 digits = 1, working = w, accent = blue))
  })

  output$t1_naive <- renderUI({
    d <- t1()
    naive <- matrix(d$N / 4, 2, 2, dimnames = list(row_lab, col_lab))
    w <- matrix(paste0(d$N, " / 4"), 2, 2)
    HTML(html_ct(naive, "The 'equal split' guess = N / 4",
                 digits = 1, working = w, accent = "gray40"))
  })

  output$t1_verdict <- renderUI({
    d   <- t1()
    N   <- d$N
    nF  <- N * d$pF
    nM  <- N * (1 - d$pF)
    eIF <- d$E[1, 1]
    eIM <- d$E[1, 2]
    f1  <- function(v) {
      if (abs(v - round(v)) < 1e-9) formatC(round(v), format = "d")
      else formatC(v, format = "f", digits = 1)
    }
    pct <- function(v) paste0(round(100 * v), "%")

    s1 <- paste0("Out of ", N, " frogs, ", pct(d$pF), " come from forest ponds (",
                 f1(nF), " frogs) and ", pct(1 - d$pF), " from meadow ponds (",
                 f1(nM), " frogs).")
    s2 <- paste0("If habitat has no effect on infection, then ", pct(d$pI),
                 " of the forest frogs should be infected (", f1(nF), " \u00d7 ",
                 formatC(d$pI, format = "f", digits = 2), " = ", f1(eIF),
                 ") and ", pct(d$pI), " of the meadow frogs should be infected (",
                 f1(nM), " \u00d7 ", formatC(d$pI, format = "f", digits = 2),
                 " = ", f1(eIM), ").")

    if (abs(d$pF - 0.5) < 1e-9 && abs(d$pI - 0.5) < 1e-9) {
      s3 <- paste0("Here both splits are exactly 50:50, so every box expects ",
                   f1(N / 4), " frogs and the equal-split guess happens to be right. ",
                   "This is the ", strong("only"), " situation where N / 4 works - ",
                   "move either slider away from 50% and it breaks.")
    } else {
      s3 <- paste0("The equal-split guess would put ", f1(N / 4),
                   " frogs in every box. That ignores the fact that ",
                   if (d$pF != 0.5) paste0("there are ",
                     if (d$pF > 0.5) "more forest frogs than meadow frogs"
                     else "more meadow frogs than forest frogs") else "",
                   if (d$pF != 0.5 && d$pI != 0.5) " and that " else "",
                   if (d$pI != 0.5) paste0(
                     if (d$pI < 0.5) "most frogs are not infected"
                     else "most frogs are infected") else "",
                   ". A box can only fill up with frogs that exist in that row ",
                   strong("and"), " that column.")
    }

    s4 <- paste0("Notice in the right-hand plot that the line between infected and ",
                 "not infected sits at the ", strong("same height"), " in both habitats. ",
                 "That flat line is what 'no link' looks like: the infection rate does not ",
                 "change when you switch habitat.")
    tagList(p(HTML(s1)), p(HTML(s2)), p(HTML(s3)), p(HTML(s4)))
  })

  # ======================= TAB 2 ============================================
  t2 <- reactive({
    input$t2_new
    pr <- get_process("t2")
    set.seed(2000 + input$t2_new)
    N       <- pr$N
    habitat <- rbinom(N, 1, pr$pF)                    # 1 = forest
    p_inf   <- ifelse(habitat == 1, pr$piF, pr$piM)
    inf     <- rbinom(N, 1, p_inf)                    # 1 = infected
    O <- make_tab(sum(inf == 1 & habitat == 1), sum(inf == 1 & habitat == 0),
                  sum(inf == 0 & habitat == 1), sum(inf == 0 & habitat == 0))
    E <- expected_tab(O)
    E[!is.finite(E)] <- 0
    naive <- matrix(N / 4, 2, 2, dimnames = list(row_lab, col_lab))
    list(type = pr$type, pr = pr, O = O, E = E, naive = naive, N = N,
         X2_E = chisq_stat(O, E), X2_naive = chisq_stat(O, naive))
  })

  output$t2_obs <- renderPlot({
    d <- t2()
    draw_mosaic(d$O, "What we caught (observed)", dots = TRUE,
                indep_line = TRUE)
  })

  output$t2_exp <- renderPlot({
    d <- t2()
    draw_mosaic(d$E, "Expected if no link (R \u00d7 C / N)", digits = 1,
                quarters = input$t2_quarters)
  })

  output$t2_tab_obs <- renderUI({
    d <- t2()
    HTML(html_ct(d$O, "Observed (O)", digits = 0, accent = "gray20"))
  })

  output$t2_tab_exp <- renderUI({
    d <- t2()
    R <- rowSums(d$O); C <- colSums(d$O)
    w <- matrix("", 2, 2)
    for (i in 1:2) for (j in 1:2) {
      w[i, j] <- paste0(R[i], " \u00d7 ", C[j], " / ", d$N)
    }
    HTML(html_ct(d$E, "Expected (E) = R \u00d7 C / N", digits = 1,
                 working = w, margins = list(R = R, C = C, N = d$N), accent = blue))
  })

  output$t2_tab_naive <- renderUI({
    d <- t2()
    HTML(html_ct(d$naive, "Equal-split guess = N / 4", digits = 1,
                 accent = "gray40"))
  })

  output$t2_verdict <- renderUI({
    d  <- t2()
    O  <- d$O; E <- d$E; N <- d$N
    R  <- rowSums(O); C <- colSums(O)
    f1 <- function(v) formatC(v, format = "f", digits = 1)
    f2 <- function(v) formatC(v, format = "f", digits = 2)
    pct <- function(v) paste0(round(100 * v), "%")

    if (any(R == 0) || any(C == 0)) {
      return(p("One row or column of this sample is empty, so expected counts cannot ",
               "be split between the two habitats. Press Resample or increase N."))
    }

    rate_F <- O[1, 1] / C[1]
    rate_M <- O[1, 2] / C[2]

    s1 <- paste0("We caught ", C[1], " forest frogs and ", C[2], " meadow frogs. In total ",
                 R[1], " of ", N, " frogs (", pct(R[1] / N), ") were infected.")
    s2 <- paste0("If infection had nothing to do with habitat, ", pct(R[1] / N),
                 " of the forest frogs should be infected: ", R[1], " \u00d7 ", C[1],
                 " / ", N, " = ", f1(E[1, 1]), " frogs. We actually saw ", O[1, 1],
                 " (", pct(rate_F), " of forest frogs), compared with ", O[1, 2],
                 " of the meadow frogs (", pct(rate_M), ").")
    s3 <- paste0("The expected table keeps exactly the same row and column totals as ",
                 "the observed table. It does not invent new frogs - it only asks how ",
                 "those same frogs would be shared out among the four boxes if there ",
                 "were no link.")

    if (d$type == "indep") {
      s4 <- paste0("In this simulation there truly is ", strong("no link"), " (both habitats ",
                   "have a ", pct(d$pr$piF), " infection rate), so the observed counts sit ",
                   "close to the expected counts and any gaps are just chance. ",
                   "Press Resample a few times to see how much they wobble.")
    } else {
      s4 <- paste0("In this simulation there ", strong("is"), " a real link (", pct(d$pr$piF),
                   " infected in forest ponds vs ", pct(d$pr$piM), " in meadow ponds), ",
                   "so the observed split lines are staggered and sit away from the dashed ",
                   "'no link' line.")
    }

    s5 <- paste0("The chi-square statistic adds up (O \u2212 E)\u00b2 / E over the four boxes. ",
                 "Using the correct expected counts gives \u03c7\u00b2 = ", f2(d$X2_E),
                 ". Using the equal-split guess instead gives \u03c7\u00b2 = ", f2(d$X2_naive),
                 if (d$X2_naive > 2 * d$X2_E + 1)
                   paste0(" - much larger, because N / 4 mixes up 'unequal numbers of frogs ",
                          "in each habitat' with 'infection depends on habitat'.")
                 else ".")
    tagList(p(HTML(s1)), p(HTML(s2)), p(HTML(s3)), p(HTML(s4)), p(HTML(s5)))
  })

  # ======================= TAB 3 ============================================
  t3 <- reactive({
    input$t3_new
    pr <- get_process("t3")
    S  <- input$t3_sims
    N  <- pr$N
    set.seed(3000 + input$t3_new)
    nF <- rbinom(S, N, pr$pF)
    nM <- N - nF
    a  <- rbinom(S, nF, pr$piF)          # infected, forest
    b  <- rbinom(S, nM, pr$piM)          # infected, meadow
    c_ <- nF - a                          # not infected, forest
    d_ <- nM - b                          # not infected, meadow
    R1 <- a + b; R2 <- c_ + d_

    E11 <- R1 * nF / N; E12 <- R1 * nM / N
    E21 <- R2 * nF / N; E22 <- R2 * nM / N
    X2  <- (a - E11)^2 / E11 + (b - E12)^2 / E12 +
           (c_ - E21)^2 / E21 + (d_ - E22)^2 / E22
    ok  <- is.finite(X2)
    p_correct <- pchisq(X2[ok], df = 1, lower.tail = FALSE)

    q  <- N / 4
    X2n <- ((a - q)^2 + (b - q)^2 + (c_ - q)^2 + (d_ - q)^2) / q
    p_naive <- pchisq(X2n, df = 3, lower.tail = FALSE)

    list(type = pr$type, pr = pr, S = S, N = N, a = a, E11 = E11,
         sig_correct = mean(p_correct < 0.05), sig_naive = mean(p_naive < 0.05),
         n_ok = sum(ok))
  })

  output$t3_hist <- renderPlot({
    d <- t3()
    par(mar = c(5, 5, 3, 1))
    br_ <- seq(min(d$a, floor(d$N / 4)) - 0.5,
               max(d$a, ceiling(d$N / 4)) + 0.5, by = 1)
    if (length(br_) > 60) br_ <- pretty(c(br_[1], br_[length(br_)]), 40)
    h <- hist(d$a, breaks = br_, plot = FALSE)
    plot(h, col = adjustcolor(col_inf, alpha.f = 0.35), border = "white",
         main = "Infected forest frogs, across all studies",
         xlab = "Number of infected forest frogs in a study",
         ylab = "Number of studies", cex.main = 1.2, cex.lab = 1.1,
         xlim = range(br_), ylim = c(0, max(h$counts) * 1.3),
         axes = FALSE)
    axis(1); axis(2, las = 1)
    box(bty = "l")
    m_obs <- mean(d$a)
    m_exp <- mean(d$E11, na.rm = TRUE)
    abline(v = m_obs,     col = col_inf, lwd = 3)
    abline(v = m_exp,     col = blue,    lwd = 3, lty = 2)
    abline(v = d$N / 4,   col = "gray30", lwd = 3, lty = 3)
    legend("topright", bg = "white", cex = 0.95,
           legend = c(paste0("Average observed: ", formatC(m_obs, format = "f", digits = 1)),
                      paste0("Average R \u00d7 C / N: ", formatC(m_exp, format = "f", digits = 1)),
                      paste0("N / 4 guess: ", formatC(d$N / 4, format = "f", digits = 1))),
           col = c(col_inf, blue, "gray30"), lty = c(1, 2, 3), lwd = 3)
  })

  output$t3_bars <- renderPlot({
    d <- t3()
    par(mar = c(5, 5, 3, 1))
    vals <- 100 * c(d$sig_correct, d$sig_naive)
    bp <- barplot(vals, names.arg = c("Correct E\n(R \u00d7 C / N)", "Equal split\n(N / 4)"),
                  col = c(adjustcolor(blue, alpha.f = 0.6), "gray70"), border = NA,
                  ylim = c(0, 110), las = 1, cex.names = 1.05,
                  ylab = "% of studies claiming a link (p < 0.05)",
                  main = "How often does the test say 'linked'?", cex.main = 1.2)
    text(bp, vals + 5, paste0(round(vals, 1), "%"), font = 2, cex = 1.2)
    if (d$type == "indep") {
      abline(h = 5, lty = 2, lwd = 2, col = col_inf)
      text(par("usr")[2], 9, "5% expected by chance", adj = c(1, 0),
           col = col_inf, cex = 0.95)
    }
    box(bty = "l")
  })

  output$t3_verdict <- renderUI({
    d   <- t3()
    f1  <- function(v) formatC(v, format = "f", digits = 1)
    pct <- function(v) paste0(round(100 * v, 1), "%")
    m_obs <- mean(d$a)
    m_exp <- mean(d$E11, na.rm = TRUE)

    s1 <- paste0("Across ", d$S, " repeated studies of ", d$N, " frogs, the number of ",
                 "infected forest frogs averaged ", f1(m_obs), ".")

    if (d$type == "indep") {
      s2 <- paste0("Because there is truly no link here, the observed average matches the ",
                   "R \u00d7 C / N expected count (", f1(m_exp), "). The equal-split guess (",
                   f1(d$N / 4), ") ",
                   if (abs(d$N / 4 - m_exp) < 0.5) "happens to land close by only because both splits are near 50:50."
                   else "misses, even though nothing interesting is going on.")
      s3 <- paste0("Using the correct expected counts, ", pct(d$sig_correct),
                   " of studies wrongly claimed a link - close to the 5% we accept as a ",
                   "false-alarm rate. Using N / 4, ", pct(d$sig_naive),
                   " of studies claimed a link.")
      s4 <- paste0("The N / 4 version is really answering a different question: 'are all ",
                   "four combinations equally common?' It finds 'something going on' ",
                   "whenever the habitats hold different numbers of frogs or most frogs share ",
                   "one infection status, even when infection and habitat are completely ",
                   "unrelated.")
    } else {
      s2 <- paste0("Here there is a real link, so the observed average (", f1(m_obs),
                   ") differs from the R \u00d7 C / N expected count (", f1(m_exp),
                   "). That gap is exactly what a chi-square test looks for.")
      s3 <- paste0("Using the correct expected counts, ", pct(d$sig_correct),
                   " of studies detected the link (this is the test's ", em("power"), "). ",
                   "Using N / 4, ", pct(d$sig_naive), " flagged 'something' - but that ",
                   "number cannot be trusted, because N / 4 also flags differences in ",
                   "the totals that have nothing to do with a link.")
      s4 <- paste0("Switch back to 'No link' to see how often the N / 4 version raises a ",
                   "false alarm.")
    }
    note <- p(style = "font-size:13px; color:#666;",
              "Chi-square tests here use no continuity correction, so results can differ ",
              "slightly from chisq.test(), which applies Yates' correction to 2 \u00d7 2 tables ",
              "by default.")
    tagList(p(HTML(s1)), p(HTML(s2)), p(HTML(s3)), p(HTML(s4)), note)
  })
}

shinyApp(ui, server)
