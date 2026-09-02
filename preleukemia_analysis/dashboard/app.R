# Shiny dashboard for the revised preleukemia analysis (stages 02-07).
# Reads only the small artifacts in dashboard/data/ produced by bin/prepare_dashboard_data.R.

library(shiny)
library(bslib)
library(plotly)
library(DT)
library(dplyr)
library(tidyr)

# --- Aequorea abyssal-marine palette (matches naraen.net) --------------------
aeq <- list(
  glow    = "#4BDDE6",
  deep    = c("#0B1D2A", "#123245", "#1B4A60", "#26637B", "#337D95", "#4198AE",
              "#50B3C5", "#4BDDE6", "#8FE8EE", "#C8F4F6"),
  accent2 = "#337D95",
  danger  = "#E67E4B"
)
pal_d <- function(n) grDevices::colorRampPalette(aeq$deep[3:10])(n)
cond_cols <- c(WT = "#337D95", Mutant = "#E67E4B")

dat <- new.env()
for (f in list.files("data", pattern = "\\.rds$", full.names = TRUE)) {
  assign(sub("\\.rds$", "", basename(f)), readRDS(f), envir = dat)
}

theme_aeq <- function(p) {
  p |> layout(
    paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
    font = list(family = "system-ui", color = "#123245")
  )
}

ui <- page_navbar(
  title = "Preleukemia HSPC atlas",
  theme = bs_theme(preset = "cosmo", primary = aeq$accent2),

  nav_panel("Overview",
    layout_columns(
      value_box("Samples", "38", "8 preleukemic mouse models", theme = "primary"),
      value_box("Cells after QC", textOutput("n_cells", inline = TRUE),
                "emptyDrops + adaptive QC + scDblFinder", theme = "primary"),
      value_box("TCGA patients", "151", "survival cohort (age + OS complete)",
                theme = "primary")
    ),
    card(card_header("About"),
      markdown(paste(
        "Revised analysis of Isobe et al. (*Cell Genomics* 2023) preleukemic",
        "mouse HSPC scRNA-seq (GSE227026). Key statistical revisions vs. the",
        "original pipeline: **emptyDrops** cell calling, **per-sample adaptive",
        "(MAD) QC**, **scDblFinder** doublet removal, **Harmony** integration,",
        "**SingleR** annotation against the Dahlin 2018 atlas, **propeller**",
        "for composition (sample-level, replaces pooled chi-square),",
        "**pseudobulk edgeR** DE (replaces cell-level tests), and age-adjusted",
        "**Cox PH** for TCGA survival. Rendered stage reports live in",
        "`analysis/*.html`."))
    )
  ),

  nav_panel("QC",
    card(card_header("Cells per sample: called vs. final"),
         plotlyOutput("qc_bars", height = 420)),
    card(card_header("Per-sample QC summary"), DTOutput("qc_table"))
  ),

  nav_panel("Atlas",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("color_by", "Color cells by",
                    c("celltype", "condition", "model", "seurat_clusters",
                      "dpt_pseudotime")),
        helpText("Random 60k-cell subsample of the integrated atlas.")
      ),
      card(plotlyOutput("umap", height = 620), full_screen = TRUE)
    )
  ),

  nav_panel("Composition",
    card(card_header("Cell-type proportions per sample (points = mice)"),
         plotlyOutput("prop_plot", height = 420)),
    card(card_header("propeller results (logit, model-adjusted, BH-FDR)"),
         DTOutput("prop_table"))
  ),

  nav_panel("Differential expression",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("de_celltype", "Cell type",
                    sort(unique(dat$de_results$celltype))),
        sliderInput("fdr", "FDR threshold", 0.01, 0.25, 0.05, step = 0.01)
      ),
      card(card_header("Pseudobulk volcano (Mutant vs WT, model-adjusted)"),
           plotlyOutput("volcano", height = 480)),
      card(DTOutput("de_table"))
    )
  ),

  nav_panel("Survival",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("signature", "Signature", unique(dat$km_data$signature)),
        helpText("KM curves at median split (display); inference is the",
                 "age-adjusted Cox model on the continuous score.")
      ),
      card(card_header("TCGA-LAML overall survival"),
           plotlyOutput("km_plot", height = 460)),
      card(card_header("Cox PH results (HR per 1 SD of score)"),
           DTOutput("cox_table"))
    )
  )
)

server <- function(input, output, session) {

  output$n_cells <- renderText(format(sum(dat$qc_summary$n_final), big.mark = ","))

  output$qc_bars <- renderPlotly({
    d <- dat$qc_summary |>
      select(title, model, called = n_called, final = n_final) |>
      pivot_longer(c(called, final))
    plot_ly(d, x = ~title, y = ~value, color = ~name,
            colors = c(aeq$deep[5], aeq$glow), type = "bar") |>
      layout(barmode = "group", xaxis = list(tickangle = -60, title = ""),
             yaxis = list(title = "cells")) |> theme_aeq()
  })

  output$qc_table <- renderDT(datatable(
    dat$qc_summary, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE))

  output$umap <- renderPlotly({
    d <- dat$umap_sample
    if (input$color_by == "dpt_pseudotime") {
      plot_ly(d, x = ~UMAP1, y = ~UMAP2, color = ~dpt_pseudotime,
              colors = aeq$deep[c(2, 8, 10)],
              type = "scattergl", mode = "markers",
              marker = list(size = 2), text = ~celltype, hoverinfo = "text") |>
        theme_aeq()
    } else {
      col <- d[[input$color_by]]
      plot_ly(d, x = ~UMAP1, y = ~UMAP2, color = ~col,
              colors = if (input$color_by == "condition") unname(cond_cols)
                       else pal_d(n_distinct(col)),
              type = "scattergl", mode = "markers",
              marker = list(size = 2), text = ~celltype, hoverinfo = "text") |>
        theme_aeq()
    }
  })

  output$prop_plot <- renderPlotly({
    plot_ly(dat$celltype_props, x = ~celltype, y = ~prop, color = ~condition,
            colors = cond_cols, type = "box", boxpoints = "all",
            jitter = 0.4, pointpos = 0, marker = list(size = 4)) |>
      layout(boxmode = "group", xaxis = list(tickangle = -45, title = ""),
             yaxis = list(title = "proportion of sample")) |> theme_aeq()
  })

  output$prop_table <- renderDT(datatable(
    dat$propeller |> mutate(across(where(is.numeric), \(x) signif(x, 3))),
    rownames = FALSE))

  output$volcano <- renderPlotly({
    d <- dat$de_results |> filter(celltype == input$de_celltype) |>
      mutate(sig = FDR < input$fdr)
    plot_ly(d, x = ~logFC, y = ~-log10(FDR), color = ~sig,
            colors = c("#B9C7CE", aeq$danger),
            type = "scattergl", mode = "markers", marker = list(size = 5),
            text = ~paste0(gene, "<br>logFC ", round(logFC, 2), " | FDR ",
                           signif(FDR, 2)), hoverinfo = "text") |>
      layout(showlegend = FALSE, xaxis = list(title = "log2 FC (Mutant vs WT)"),
             yaxis = list(title = "-log10 FDR")) |> theme_aeq()
  })

  output$de_table <- renderDT({
    datatable(
      dat$de_results |> filter(celltype == input$de_celltype, FDR < input$fdr) |>
        arrange(FDR) |> mutate(across(where(is.numeric), \(x) signif(x, 3))),
      options = list(pageLength = 10), rownames = FALSE)
  })

  output$km_plot <- renderPlotly({
    d <- dat$km_data |> filter(signature == input$signature)
    km_curve <- function(sub) {
      sub <- arrange(sub, time)
      n <- nrow(sub); at_risk <- n; surv <- 1
      steps <- tibble(time = 0, surv = 1)
      for (i in seq_len(n)) {
        if (sub$event[i] == 1) surv <- surv * (1 - 1 / at_risk)
        at_risk <- at_risk - 1
        steps <- add_row(steps, time = sub$time[i], surv = surv)
      }
      steps
    }
    p <- plot_ly()
    for (grp in c("high", "low")) {
      s <- km_curve(filter(d, group == grp))
      p <- add_lines(p, data = s, x = ~time, y = ~surv, name = paste(input$signature, grp),
                     line = list(shape = "hv",
                                 color = if (grp == "high") aeq$danger else aeq$accent2))
    }
    p |> layout(xaxis = list(title = "months"),
                yaxis = list(title = "overall survival", range = c(0, 1))) |> theme_aeq()
  })

  output$cox_table <- renderDT(datatable(
    dat$cox_results |> mutate(across(where(is.numeric), \(x) signif(x, 3))),
    rownames = FALSE))
}

shinyApp(ui, server)
