# evalinger Shiny App
# Design Calculator, Monitoring Dashboard, Method Comparison
#
# Launch with: shiny::runApp(system.file("shiny", package = "evalinger"))

library(shiny)
library(bslib)

# Source package functions if not installed
pkg_loaded <- requireNamespace("evalinger", quietly = TRUE)
if (!pkg_loaded) {
  pkg_dir <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."),
                                     "..", ".."), mustWork = FALSE)
  r_files <- list.files(file.path(pkg_dir, "R"), full.names = TRUE)
  if (length(r_files) > 0) {
    for (f in sort(r_files)) source(f, local = TRUE)
  }
}

# ---------- UI ----------
ui <- page_navbar(
  title = "evalinger: E-Values for Clinical Trials",
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2c3e50",
    "navbar-bg" = "#2c3e50"
  ),

  # Tab 1: Design Calculator
  nav_panel(
    "Design",
    layout_sidebar(
      sidebar = sidebar(
        title = "Trial Parameters",
        numericInput("p_C", "Control rate (p_C)", 0.30, min = 0.01, max = 0.99,
                     step = 0.01),
        numericInput("delta", "Treatment effect (delta)", 0.15, min = 0.01,
                     max = 0.50, step = 0.01),
        numericInput("alpha_d", "Alpha", 0.025, min = 0.001, max = 0.10,
                     step = 0.005),
        numericInput("power_d", "Target power", 0.80, min = 0.50, max = 0.99,
                     step = 0.05),
        numericInput("nrep_d", "MC replications", 2000, min = 500, max = 50000,
                     step = 500),
        actionButton("run_design", "Compute Design", class = "btn-primary")
      ),
      card(
        card_header("Design Results"),
        verbatimTextOutput("design_output")
      ),
      card(
        card_header("Lambda Sensitivity"),
        plotOutput("lambda_plot", height = "350px")
      )
    )
  ),

  # Tab 2: Monitoring Dashboard
  nav_panel(
    "Monitor",
    layout_sidebar(
      sidebar = sidebar(
        title = "Monitoring Setup",
        numericInput("mon_lambda", "Betting fraction (lambda)", 0.31,
                     min = 0.01, max = 0.99, step = 0.01),
        numericInput("mon_alpha", "Alpha", 0.025, min = 0.001, max = 0.10,
                     step = 0.005),
        hr(),
        h5("Simulate interim data"),
        numericInput("mon_pT", "True p_T (for simulation)", 0.45,
                     min = 0.01, max = 0.99, step = 0.01),
        numericInput("mon_pC", "True p_C (for simulation)", 0.30,
                     min = 0.01, max = 0.99, step = 0.01),
        numericInput("batch_size", "Patients per batch (per arm)", 25,
                     min = 5, max = 200, step = 5),
        actionButton("add_batch", "Add Batch", class = "btn-primary"),
        actionButton("reset_mon", "Reset Monitor", class = "btn-outline-danger")
      ),
      card(
        card_header("Monitor Status"),
        verbatimTextOutput("monitor_status")
      ),
      card(
        card_header("E-Process Path"),
        plotOutput("monitor_plot", height = "350px")
      ),
      card(
        card_header("Confidence Sequence"),
        plotOutput("cs_plot", height = "300px")
      )
    )
  ),

  # Tab 3: Method Comparison
  nav_panel(
    "Compare",
    layout_sidebar(
      sidebar = sidebar(
        title = "Comparison Parameters",
        numericInput("cmp_pC", "Control rate (p_C)", 0.30, min = 0.01,
                     max = 0.99, step = 0.01),
        numericInput("cmp_pT", "Alternative (p_T)", 0.45, min = 0.02,
                     max = 0.99, step = 0.01),
        numericInput("cmp_Nmax", "Max N per arm", 200, min = 50, max = 2000,
                     step = 50),
        numericInput("cmp_looks", "Number of looks", 4, min = 1, max = 20,
                     step = 1),
        numericInput("cmp_nrep", "MC replications", 2000, min = 200,
                     max = 50000, step = 500),
        numericInput("cmp_alpha", "Alpha", 0.025, min = 0.001, max = 0.10,
                     step = 0.005),
        actionButton("run_compare", "Run Comparison", class = "btn-primary")
      ),
      card(
        card_header("Comparison Results"),
        verbatimTextOutput("compare_output")
      ),
      card(
        card_header("Visualization"),
        plotOutput("compare_plot", height = "400px")
      )
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {

  # --- Tab 1: Design ---
  design_result <- reactiveVal(NULL)

  observeEvent(input$run_design, {
    withProgress(message = "Computing design...", {
      des <- edesign_binary(
        p_C = input$p_C,
        delta = input$delta,
        alpha = input$alpha_d,
        power = input$power_d,
        nrep = input$nrep_d,
        seed = 42
      )
      design_result(des)
    })
  })

  output$design_output <- renderPrint({
    req(design_result())
    print(design_result())
  })

  output$lambda_plot <- renderPlot({
    req(design_result())
    des <- design_result()
    grid <- grow_lambda_grid(p_C = des$p_C,
                             delta_grid = des$delta,
                             lambda_grid = seq(0.05, 0.95, by = 0.02),
                             alpha = des$alpha)
    plot(grid$lambda, grid$expected_n, type = "l", lwd = 2, col = "steelblue",
         xlab = expression(lambda), ylab = "Expected N per arm",
         main = "Expected Stopping Time vs Betting Fraction")
    abline(v = des$lambda, col = "red", lty = 2)
    legend("topright", legend = sprintf("GROW optimal: %.3f", des$lambda),
           col = "red", lty = 2, bty = "n")
  })

  # --- Tab 2: Monitor ---
  mon_state <- reactiveVal(NULL)
  all_x_T <- reactiveVal(integer(0))
  all_x_C <- reactiveVal(integer(0))

  observeEvent(input$reset_mon, {
    mon_state(emonitor(alpha = input$mon_alpha, lambda = input$mon_lambda))
    all_x_T(integer(0))
    all_x_C(integer(0))
  })

  observeEvent(input$add_batch, {
    mon <- mon_state()
    if (is.null(mon)) {
      mon <- emonitor(alpha = input$mon_alpha, lambda = input$mon_lambda)
    }
    set.seed(as.integer(Sys.time()))
    new_T <- rbinom(input$batch_size, 1, input$mon_pT)
    new_C <- rbinom(input$batch_size, 1, input$mon_pC)
    mon <- update(mon, x_T = new_T, x_C = new_C)
    mon_state(mon)
    all_x_T(c(all_x_T(), new_T))
    all_x_C(c(all_x_C(), new_C))
  })

  output$monitor_status <- renderPrint({
    req(mon_state())
    print(mon_state())
  })

  output$monitor_plot <- renderPlot({
    req(mon_state())
    mon <- mon_state()
    if (mon$n == 0) return(NULL)
    nn <- seq_len(mon$n)
    plot(nn, mon$history_log_e, type = "l", col = "steelblue", lwd = 1.5,
         xlab = "Observations per arm", ylab = "log(E)",
         main = "E-Process Path")
    abline(h = log(1 / mon$alpha), col = "red", lty = 2, lwd = 1.5)
    abline(h = 0, col = "gray", lty = 3)
    if (mon$rejected) {
      abline(v = mon$rejection_time, col = "darkgreen", lty = 3)
    }
  })

  output$cs_plot <- renderPlot({
    req(mon_state())
    mon <- mon_state()
    if (mon$n < 10) return(NULL)
    xT <- all_x_T()
    xC <- all_x_C()
    cs <- confseq_binary(xT, xC, alpha = 1 - 0.95)
    nn <- seq_len(cs$n)
    plot(nn, cs$delta_hat, type = "l", ylim = range(c(cs$lower, cs$upper)),
         xlab = "Observations per arm", ylab = "Treatment effect",
         main = "Confidence Sequence (95%)")
    polygon(c(nn, rev(nn)), c(cs$lower, rev(cs$upper)),
            col = rgb(0.27, 0.51, 0.71, 0.2), border = NA)
    lines(nn, cs$lower, col = "steelblue", lty = 2)
    lines(nn, cs$upper, col = "steelblue", lty = 2)
    abline(h = 0, col = "gray50", lty = 3)
  })

  # --- Tab 3: Compare ---
  compare_result <- reactiveVal(NULL)

  observeEvent(input$run_compare, {
    withProgress(message = "Running simulation...", {
      cmp <- simulate_comparison(
        p_C = input$cmp_pC,
        p_T_alt = input$cmp_pT,
        Nmax = input$cmp_Nmax,
        n_looks = input$cmp_looks,
        alpha = input$cmp_alpha,
        nrep = input$cmp_nrep,
        seed = 42
      )
      compare_result(cmp)
    })
  })

  output$compare_output <- renderPrint({
    req(compare_result())
    print(compare_result())
  })

  output$compare_plot <- renderPlot({
    req(compare_result())
    cmp <- compare_result()
    res <- cmp$results

    par(mfrow = c(1, 3), mar = c(7, 4, 3, 1))
    cols <- c(evalue = "steelblue", gs_obf = "darkorange",
              naive_p = "salmon", cal_bayes = "seagreen")
    bar_cols <- cols[res$method]

    barplot(res$null_rej, names.arg = res$method, main = "Type I Error",
            col = bar_cols, las = 2, ylab = "Rejection rate", ylim = c(0, max(0.1, max(res$null_rej) * 1.2)))
    abline(h = cmp$design$alpha, col = "red", lty = 2)

    barplot(res$alt_rej, names.arg = res$method, main = "Power",
            col = bar_cols, las = 2, ylab = "Rejection rate", ylim = c(0, 1))

    barplot(res$avg_n_alt, names.arg = res$method,
            main = "Avg N (alternative)", col = bar_cols, las = 2,
            ylab = "Average N per arm")
  })
}

shinyApp(ui, server)
