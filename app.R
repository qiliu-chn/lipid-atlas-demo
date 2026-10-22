# Lipid Atlas Demo — an R Shiny app
# Explore a SIMULATED quantitative lipidomics dataset:
#   1) Explorer   – search a lipid, view concentrations by group / sex
#   2) Volcano    – univariate Welch t-tests + fold change, multiple-testing correction
#   3) PCA        – log transform + Pareto scaling, scores and loadings
#   4) Classifier – logistic regression, k-fold cross-validated ROC/AUC,
#                   and a label-permutation test of the AUC
#
# Run locally:  shiny::runApp()

library(shiny)
library(ggplot2)
library(DT)

# ---- Data ----------------------------------------------------------------
dat <- read.csv("data/demo_lipids.csv", check.names = FALSE, stringsAsFactors = FALSE)
dat$group <- factor(dat$group, levels = c("Control", "Cancer"))
lipid_cols <- setdiff(names(dat), c("sample_id", "group", "sex"))
lipid_class <- sub(" .*", "", lipid_cols)
log_mat <- log2(as.matrix(dat[, lipid_cols]))

grp_colors <- c(Control = "#4C72B0", Cancer = "#DD8452")

# ---- Helper functions ----------------------------------------------------
pareto_scale <- function(m) {
  m <- scale(m, center = TRUE, scale = FALSE)
  sweep(m, 2, sqrt(apply(m, 2, sd)), "/")
}

auc_rank <- function(score, y) {          # y: 1 = Cancer, 0 = Control
  r <- rank(score)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

roc_points <- function(score, y) {
  th <- sort(unique(c(-Inf, score, Inf)), decreasing = TRUE)
  data.frame(
    fpr = sapply(th, function(t) mean(score[y == 0] >= t)),
    tpr = sapply(th, function(t) mean(score[y == 1] >= t))
  )
}

# k-fold CV. Lipid selection happens INSIDE each training fold,
# so the held-out fold never influences which lipids are used.
cv_predict <- function(X, y, k_folds = 5, n_lipids = 5, seed = 1) {
  set.seed(seed)
  folds <- sample(rep(seq_len(k_folds), length.out = length(y)))
  pred <- numeric(length(y))
  for (f in seq_len(k_folds)) {
    tr <- folds != f
    p_tr <- apply(X[tr, , drop = FALSE], 2,
                  function(v) t.test(v[y[tr] == 1], v[y[tr] == 0])$p.value)
    keep <- names(sort(p_tr))[seq_len(n_lipids)]
    d_tr <- data.frame(y = y[tr], X[tr, keep, drop = FALSE], check.names = TRUE)
    d_te <- data.frame(X[!tr, keep, drop = FALSE], check.names = TRUE)
    fit <- suppressWarnings(glm(y ~ ., data = d_tr, family = binomial))
    pred[!tr] <- predict(fit, newdata = d_te, type = "response")
  }
  pred
}

# ---- UI --------------------------------------------------------------------
ui <- navbarPage(
  title = "Lipid Atlas Demo",
  header = div(style = "padding: 0 15px; color: #a33;",
               strong("Demo only: all data are simulated, not real patient data.")),

  tabPanel("Explorer",
    sidebarLayout(
      sidebarPanel(
        selectInput("cls", "Lipid class", c("All", sort(unique(lipid_class)))),
        selectizeInput("lip", "Lipid species (type to search)", choices = lipid_cols),
        checkboxInput("by_sex", "Split by sex", FALSE),
        checkboxInput("log_y", "Log scale", TRUE),
        downloadButton("dl", "Download data for this lipid")
      ),
      mainPanel(plotOutput("box", height = 420), DTOutput("summary_tbl"))
    )
  ),

  tabPanel("Volcano",
    sidebarLayout(
      sidebarPanel(
        selectInput("padj", "Multiple-testing correction",
                    c("Bonferroni" = "bonferroni", "Benjamini-Hochberg (FDR)" = "BH", "None" = "none")),
        sliderInput("fc_thr", "Fold-change threshold (±%)", 0, 50, 20, step = 5),
        sliderInput("p_thr", "Adjusted p-value threshold", 0.001, 0.1, 0.05, step = 0.001)
      ),
      mainPanel(plotOutput("volcano", height = 450), DTOutput("uni_tbl"))
    )
  ),

  tabPanel("PCA",
    sidebarLayout(
      sidebarPanel(
        radioButtons("pca_col", "Colour samples by", c("group", "sex")),
        helpText("Preprocessing: log2 transform, mean-centring, Pareto scaling.")
      ),
      mainPanel(plotOutput("pca_scores", height = 420), plotOutput("pca_load", height = 350))
    )
  ),

  tabPanel("Classifier",
    sidebarLayout(
      sidebarPanel(
        sliderInput("k", "Number of CV folds", 3, 10, 5),
        sliderInput("nlip", "Lipids used in model (chosen inside each fold)", 1, 15, 5),
        sliderInput("nperm", "Number of label permutations", 20, 500, 100, step = 20),
        actionButton("run", "Run cross-validation + permutation test", class = "btn-primary"),
        helpText("Cross-validation estimates how well the model predicts unseen samples.",
                 "The permutation test asks whether that AUC is better than chance,",
                 "by repeating the whole CV on shuffled group labels.")
      ),
      mainPanel(
        fluidRow(column(6, plotOutput("roc", height = 380)),
                 column(6, plotOutput("perm", height = 380))),
        verbatimTextOutput("clf_txt")
      )
    )
  ),

  tabPanel("About",
    fluidPage(
      h3("About this demo"),
      p("A small prototype of an interactive lipid concentration database: search,",
        "explore, test and model quantitative lipidomics data."),
      p("Data: simulated with R/simulate_data.R (300 samples, 37 lipid species, 6 classes)."),
      p("Built with R, shiny, ggplot2 and DT.")
    )
  )
)

# ---- Server ----------------------------------------------------------------
server <- function(input, output, session) {

  # Explorer ---------------------------------------------------------------
  observeEvent(input$cls, {
    ch <- if (input$cls == "All") lipid_cols else lipid_cols[lipid_class == input$cls]
    updateSelectizeInput(session, "lip", choices = ch, selected = ch[1])
  })

  lipid_df <- reactive({
    req(input$lip)
    data.frame(sample_id = dat$sample_id, group = dat$group, sex = dat$sex,
               conc = dat[[input$lip]])
  })

  output$box <- renderPlot({
    d <- lipid_df()
    p <- ggplot(d, aes(group, conc, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
      scale_fill_manual(values = grp_colors) +
      labs(title = input$lip, x = NULL, y = "Concentration (nmol/mL)") +
      theme_minimal(base_size = 14) + theme(legend.position = "none")
    if (input$log_y) p <- p + scale_y_log10()
    if (input$by_sex) p <- p + facet_wrap(~ sex)
    p
  })

  output$summary_tbl <- renderDT({
    d <- lipid_df()
    keys <- if (input$by_sex) list(d$group, d$sex) else list(d$group)
    s <- aggregate(d$conc, keys, function(x)
      c(n = length(x), median = median(x), Q1 = quantile(x, .25), Q3 = quantile(x, .75)))
    s <- do.call(data.frame, s)
    names(s) <- c(if (input$by_sex) c("Group", "Sex") else "Group",
                  "n", "Median", "Q1", "Q3")
    datatable(s, rownames = FALSE, options = list(dom = "t")) |>
      formatRound(c("Median", "Q1", "Q3"), 3)
  })

  output$dl <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9]", "_", input$lip), ".csv"),
    content = function(file) write.csv(lipid_df(), file, row.names = FALSE)
  )

  # Volcano ----------------------------------------------------------------
  uni <- reactive({
    y <- dat$group == "Cancer"
    res <- data.frame(
      lipid = lipid_cols, class = lipid_class,
      fold_change = colMeans(dat[y, lipid_cols]) / colMeans(dat[!y, lipid_cols]),
      p = apply(log_mat, 2, function(v) t.test(v[y], v[!y])$p.value)  # Welch by default
    )
    res$p_adj <- p.adjust(res$p, method = input$padj)
    thr <- input$fc_thr / 100
    res$significant <- res$p_adj < input$p_thr &
      (res$fold_change >= 1 + thr | res$fold_change <= 1 - thr)
    res[order(res$p), ]
  })

  output$volcano <- renderPlot({
    r <- uni(); thr <- input$fc_thr / 100
    ggplot(r, aes(log2(fold_change), -log10(p_adj), colour = significant)) +
      geom_point(size = 2.5) +
      geom_text(data = r[r$significant, ], aes(label = lipid), vjust = -0.8, size = 3.5,
                show.legend = FALSE) +
      geom_hline(yintercept = -log10(input$p_thr), linetype = 2) +
      geom_vline(xintercept = log2(c(1 - thr, 1 + thr)), linetype = 2) +
      scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#C44E52")) +
      labs(x = "log2 fold change (Cancer / Control)", y = "-log10 adjusted p",
           colour = "Significant") +
      theme_minimal(base_size = 14)
  })

  output$uni_tbl <- renderDT({
    datatable(uni(), rownames = FALSE, options = list(pageLength = 10)) |>
      formatRound("fold_change", 3) |> formatSignif(c("p", "p_adj"), 3)
  })

  # PCA --------------------------------------------------------------------
  pca <- prcomp(pareto_scale(log_mat), center = FALSE, scale. = FALSE)
  ve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  output$pca_scores <- renderPlot({
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], colour = dat[[input$pca_col]])
    ggplot(d, aes(PC1, PC2, colour = colour)) +
      geom_point(alpha = 0.7, size = 2) + stat_ellipse() +
      labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
           colour = input$pca_col, title = "PCA scores") +
      theme_minimal(base_size = 14)
  })

  output$pca_load <- renderPlot({
    d <- data.frame(lipid = rownames(pca$rotation), PC1 = pca$rotation[, 1], class = lipid_class)
    d <- head(d[order(-abs(d$PC1)), ], 15)
    ggplot(d, aes(reorder(lipid, PC1), PC1, fill = class)) + geom_col() + coord_flip() +
      labs(x = NULL, y = "PC1 loading", title = "Top 15 PC1 loadings") +
      theme_minimal(base_size = 13)
  })

  # Classifier -------------------------------------------------------------
  clf <- eventReactive(input$run, {
    y <- as.integer(dat$group == "Cancer")
    withProgress(message = "Running", value = 0, {
      pred <- cv_predict(log_mat, y, input$k, input$nlip)
      obs_auc <- auc_rank(pred, y)
      perm_auc <- numeric(input$nperm)
      for (i in seq_len(input$nperm)) {
        y_perm <- sample(y)
        perm_auc[i] <- auc_rank(cv_predict(log_mat, y_perm, input$k, input$nlip, seed = i), y_perm)
        incProgress(1 / input$nperm)
      }
    })
    list(pred = pred, y = y, auc = obs_auc, perm = perm_auc,
         p_perm = (sum(perm_auc >= obs_auc) + 1) / (input$nperm + 1))
  })

  output$roc <- renderPlot({
    r <- clf(); pts <- roc_points(r$pred, r$y)
    ggplot(pts, aes(fpr, tpr)) + geom_path(linewidth = 1.2, colour = "#DD8452") +
      geom_abline(linetype = 2, colour = "grey50") + coord_equal() +
      labs(x = "1 - specificity", y = "Sensitivity",
           title = sprintf("Cross-validated ROC (AUC = %.3f)", r$auc)) +
      theme_minimal(base_size = 13)
  })

  output$perm <- renderPlot({
    r <- clf()
    ggplot(data.frame(auc = r$perm), aes(auc)) +
      geom_histogram(bins = 25, fill = "grey70", colour = "white") +
      geom_vline(xintercept = r$auc, colour = "#C44E52", linewidth = 1.2) +
      labs(x = "AUC with shuffled labels", y = "Count",
           title = sprintf("Permutation test (p = %.3f)", r$p_perm)) +
      theme_minimal(base_size = 13)
  })

  output$clf_txt <- renderText({
    r <- clf(); cls <- ifelse(r$pred > 0.5, 1, 0)
    sprintf(paste0("Cut-off 0.5 on cross-validated predictions\n",
                   "Sensitivity: %.1f%%   Specificity: %.1f%%   Accuracy: %.1f%%\n",
                   "Observed AUC: %.3f   Mean AUC under permutation: %.3f   Permutation p: %.3f"),
            100 * mean(cls[r$y == 1] == 1), 100 * mean(cls[r$y == 0] == 0),
            100 * mean(cls == r$y), r$auc, mean(r$perm), r$p_perm)
  })
}

shinyApp(ui, server)
