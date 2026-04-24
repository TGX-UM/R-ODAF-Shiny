
library(shiny)
library(DT)
library(data.table)
library(dplyr)
library(tidyr)
library(DESeq2)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(cowplot)
library(corrplot)
library(PCAtools)
library(plotly)
library(stringr)
library(grid)
library(RSQLite)
library(DBI)
library(enrichR)
library(AnnotationDbi)

options(shiny.maxRequestSize = 100 * 1024^2)

# Explicit static resource mapping for logos
if (dir.exists("www/logos")) {
  addResourcePath(
    prefix = "logos",
    directoryPath = normalizePath("www/logos", winslash = "/", mustWork = TRUE)
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# FIXED FILTER PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
MIN_COVERAGE <- 5e6
MIN_CPM <- 1
FDR_CUTOFF <- 0.01
VAR_THRESHOLD <- 20

# Pathway analysis parameters
PATHWAY_INPUT_FDR <- 0.01
PATHWAY_RESULT_FDR <- 0.05
MIN_PATHWAY_GENES <- 5

# Bundled example data
EXAMPLE_COUNTS_PATH <- file.path("www", "examples", "Organoids_NORM_human.txt")
EXAMPLE_META_PATH <- file.path("www", "examples", "Organoids_meta_human.txt")

############################ helper functions ##############################

.fmt_time <- function(t) {
  if (is.null(t) || is.na(t)) return(NA_character_)
  format(as.POSIXct(t), "%Y-%m-%d %H:%M:%S %Z", tz = Sys.timezone())
}

.make_run_stamp <- function(t = Sys.time()) {
  format(as.POSIXct(t), "%Y%m%d_%H%M%S", tz = Sys.timezone())
}

.normalize_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  toupper(x)
}

.has_real_filter <- function(x) {
  !is.null(x) && nzchar(trimws(x)) && trimws(x) != "None"
}

check_outliers_vst <- function(dds, meta_df, design_col, conds, var_thresh) {
  vsd_mat <- vst(dds, blind = TRUE) %>% assay()
  pca_obj <- pca(vsd_mat, metadata = meta_df, removeVar = 0.1)
  out <- character()
  pcs <- which(pca_obj$variance >= var_thresh)
  
  for (pc in pcs) {
    max_pc <- abs(max(pca_obj$rotated[, pc]) - min(pca_obj$rotated[, pc]))
    for (c in conds) {
      samp_ids <- rownames(meta_df)[meta_df[[design_col]] == c]
      if (length(samp_ids) < 2) next
      vals <- pca_obj$rotated[samp_ids, pc, drop = FALSE] %>%
        as.data.frame() %>%
        arrange(desc(.[, 1]))
      if (nrow(vals) < 2) next
      dmat <- dist(vals[, 1]) %>% as.matrix()
      diag_vals <- dmat[row(dmat) == col(dmat) + 1]
      flags <- ((pca_obj$variance[pc] * diag_vals) / max_pc) > var_thresh
      if (any(flags)) out <- c(out, rownames(vals)[which(flags)])
    }
  }
  unique(out)
}

.read_counts <- function(fileinfo) {
  req(fileinfo)
  ext <- tolower(tools::file_ext(fileinfo$name))
  df <- if (ext == "csv") {
    read.csv(fileinfo$datapath, row.names = 1, check.names = FALSE)
  } else {
    read.delim(fileinfo$datapath, row.names = 1, check.names = FALSE)
  }
  as.matrix(df)
}

.read_meta <- function(fileinfo) {
  req(fileinfo)
  ext <- tolower(tools::file_ext(fileinfo$name))
  if (ext == "csv") {
    read.csv(fileinfo$datapath, row.names = 1, check.names = FALSE)
  } else {
    read.delim(fileinfo$datapath, row.names = 1, check.names = FALSE)
  }
}

.get_orgdb <- function(organism = c("human", "mouse", "rat")) {
  organism <- match.arg(organism)
  pkg <- switch(
    organism,
    human = "org.Hs.eg.db",
    mouse = "org.Mm.eg.db",
    rat = "org.Rn.eg.db"
  )
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Missing required annotation package: ", pkg, ".\n",
        "Install it (Bioconductor) or switch Organism to one that is available."
      ),
      call. = FALSE
    )
  }
  get(pkg, envir = asNamespace(pkg))
}

map_to_symbol <- function(gene_ids, organism = c("human", "mouse", "rat")) {
  organism <- match.arg(organism)
  if (length(gene_ids) == 0) return(character(0))
  ids_clean <- sub("\\.\\d+$", "", gene_ids)
  orgdb <- .get_orgdb(organism)
  res <- suppressMessages(
    AnnotationDbi::select(
      orgdb,
      keys = unique(ids_clean),
      keytype = "ENSEMBL",
      columns = c("SYMBOL")
    )
  )
  res <- res[!is.na(res$SYMBOL) & res$SYMBOL != "", , drop = FALSE]
  if (nrow(res) == 0) {
    out <- rep(NA_character_, length(ids_clean))
    names(out) <- gene_ids
    return(out)
  }
  res <- res[!duplicated(res$ENSEMBL), , drop = FALSE]
  sym_map <- setNames(res$SYMBOL, res$ENSEMBL)
  out <- sym_map[ids_clean]
  names(out) <- gene_ids
  out
}

.sanitize_label <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (is.na(x) || x == "") "dataset" else x
}

.make_internal_fileinfo <- function(path) {
  list(
    name = basename(path),
    datapath = path
  )
}

.example_data_ready <- function() {
  file.exists(EXAMPLE_COUNTS_PATH) && file.exists(EXAMPLE_META_PATH)
}

.get_enrichr_pathway_libs <- function(organism = c("human", "mouse", "rat")) {
  organism <- match.arg(organism)
  switch(
    organism,
    human = c("KEGG_2021_Human", "WikiPathways_2024_Human", "Reactome_2022"),
    mouse = c("KEGG_2021_Mouse", "WikiPathways_2024_Mouse", "Reactome_2022"),
    rat   = c("KEGG_2021_Human", "WikiPathways_2024_Human", "Reactome_2022")
  )
}

.check_enrichr_libs <- function(libs) {
  tryCatch({
    dbs <- enrichR::listEnrichrDbs()
    available <- unique(dbs$libraryName)
    libs[libs %in% available]
  }, error = function(e) {
    character(0)
  })
}

.run_pathway_enrichr <- function(pathway_degs_hgnc,
                                 organism = c("human", "mouse", "rat"),
                                 input_fdr = PATHWAY_INPUT_FDR,
                                 result_fdr = PATHWAY_RESULT_FDR,
                                 min_genes = MIN_PATHWAY_GENES) {
  organism <- match.arg(organism)
  
  if (is.null(pathway_degs_hgnc) || nrow(pathway_degs_hgnc) == 0) {
    return(list(
      input_genes = character(0),
      results = data.frame(),
      libraries_used = character(0)
    ))
  }
  
  if (!all(c("GeneSymbol", "padj") %in% colnames(pathway_degs_hgnc))) {
    return(list(
      input_genes = character(0),
      results = data.frame(),
      libraries_used = character(0)
    ))
  }
  
  deg_genes <- pathway_degs_hgnc %>%
    dplyr::filter(!is.na(padj), padj < input_fdr, !is.na(GeneSymbol), GeneSymbol != "") %>%
    dplyr::pull(GeneSymbol) %>%
    unique()
  
  libs <- .get_enrichr_pathway_libs(organism)
  libs <- .check_enrichr_libs(libs)
  
  if (length(libs) == 0) {
    return(list(
      input_genes = deg_genes,
      results = data.frame(),
      libraries_used = character(0)
    ))
  }
  
  if (length(deg_genes) < min_genes) {
    return(list(
      input_genes = deg_genes,
      results = data.frame(),
      libraries_used = libs
    ))
  }
  
  enr_raw <- tryCatch(
    enrichR::enrichr(
      genes = deg_genes,
      databases = libs
    ),
    error = function(e) {
      warning("enrichR pathway analysis failed: ", conditionMessage(e))
      return(list())
    }
  )
  
  if (length(enr_raw) == 0) {
    return(list(
      input_genes = deg_genes,
      results = data.frame(),
      libraries_used = libs
    ))
  }
  
  enr_tbl <- dplyr::bind_rows(
    lapply(names(enr_raw), function(db) {
      x <- enr_raw[[db]]
      if (is.null(x) || nrow(x) == 0) return(NULL)
      x$Library <- db
      x$InputGeneCount <- length(deg_genes)
      x
    })
  )
  
  if (is.null(enr_tbl) || nrow(enr_tbl) == 0) {
    enr_tbl <- data.frame()
  } else if ("Adjusted.P.value" %in% colnames(enr_tbl)) {
    enr_tbl <- enr_tbl %>%
      dplyr::filter(!is.na(Adjusted.P.value), Adjusted.P.value < result_fdr) %>%
      dplyr::arrange(Adjusted.P.value, dplyr::desc(Combined.Score))
  }
  
  list(
    input_genes = deg_genes,
    results = enr_tbl,
    libraries_used = libs
  )
}

.combine_all_pathway_results <- function(res_list) {
  if (is.null(res_list) || length(res_list) == 0) return(data.frame())
  out <- lapply(res_list, function(x) {
    tbl <- x$pathway_results
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    tbl$dataset_name <- x$dataset_name
    tbl$dataset_tag <- x$dataset_tag
    tbl$organism <- x$organism_used
    tbl$input_fdr <- x$pathway_input_fdr
    tbl
  })
  dplyr::bind_rows(out)
}

.combine_all_pathway_input_genes <- function(res_list) {
  if (is.null(res_list) || length(res_list) == 0) return(data.frame())
  out <- lapply(res_list, function(x) {
    if (is.null(x$pathway_input_genes) || length(x$pathway_input_genes) == 0) return(NULL)
    data.frame(
      dataset_name = x$dataset_name,
      dataset_tag = x$dataset_tag,
      organism = x$organism_used,
      input_fdr = x$pathway_input_fdr,
      GeneSymbol = x$pathway_input_genes,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

.build_summary_text <- function(res) {
  contrast_line <- if (!is.null(res$contrast_desc) && nzchar(res$contrast_desc)) {
    paste0("Contrast: ", res$contrast_desc)
  } else {
    paste0("Conditions: ", res$condition1, " vs ", res$condition2)
  }
  
  run_lines <- c(
    "",
    "Run timing:",
    paste0("Run stamp: ", if (!is.null(res$run_stamp)) res$run_stamp else NA_character_),
    paste0("Run started: ", .fmt_time(res$run_start)),
    paste0("Run ended: ", .fmt_time(res$run_end)),
    paste0("Run duration: ", if (!is.null(res$run_duration_sec)) round(res$run_duration_sec, 1) else NA_real_, " sec"),
    paste0("Dataset duration: ", if (!is.null(res$dataset_duration_sec)) round(res$dataset_duration_sec, 1) else NA_real_, " sec")
  )
  
  txt <- c(
    "R-ODAF DESeq2 Pipeline summary",
    "",
    paste0("Dataset: ", res$dataset_name),
    paste0("Dataset tag: ", res$dataset_tag),
    paste0("Organism: ", res$organism_used),
    paste0("Design column: ", res$design_col),
    contrast_line
  )
  
  if (.has_real_filter(res$filter_desc)) {
    txt <- c(txt, paste0("Filter: ", res$filter_desc))
  }
  
  txt <- c(
    txt,
    "",
    paste0("Min coverage: ", MIN_COVERAGE),
    paste0("Min CPM thresh: ", MIN_CPM),
    paste0("FDR cutoff: ", FDR_CUTOFF),
    paste0("Outlier var %: ", VAR_THRESHOLD),
    "",
    paste0("Samples used: ", res$n_samples_used),
    paste0("Samples removed: ", length(res$removed)),
    paste0("Total genes tested: ", nrow(res$all_res)),
    paste0("Final DEGs: ", nrow(res$degs))
  )
  
  if (isTRUE(res$pathway_enabled)) {
    pathway_lines <- c(
      "",
      "Pathway analysis:",
      paste0("Input DEG FDR: ", if (!is.null(res$pathway_input_fdr)) res$pathway_input_fdr else NA_real_),
      paste0("Result FDR: ", PATHWAY_RESULT_FDR),
      paste0("Input genes: ", if (!is.null(res$pathway_input_genes)) length(res$pathway_input_genes) else 0),
      paste0("Pathway terms: ", if (!is.null(res$pathway_results)) nrow(res$pathway_results) else 0)
    )
    txt <- c(txt, pathway_lines)
  }
  
  txt <- c(
    txt,
    run_lines,
    "",
    "Please cite:",
    "Citation information will be available"
  )
  
  paste(txt, collapse = "\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# SQL helper
# ─────────────────────────────────────────────────────────────────────────────

.ensure_sql_column <- function(con, table, col_name, col_type) {
  info <- dbGetQuery(con, sprintf("PRAGMA table_info(%s);", table))
  existing <- info$name
  if (!(col_name %in% existing)) {
    dbExecute(con, sprintf("ALTER TABLE %s ADD COLUMN %s %s;", table, col_name, col_type))
  }
}

.init_sql_warehouse <- function(con) {
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS master_degs (
      GeneID TEXT,
      GeneSymbol TEXT,
      GeneSymbol_Key TEXT,
      baseMean REAL,
      log2FoldChange REAL,
      lfcSE REAL,
      stat REAL,
      pvalue REAL,
      padj REAL,
      dataset_tag TEXT,
      organism TEXT,
      run_stamp TEXT,
      timestamp TEXT
    );
  ")
  
  .ensure_sql_column(con, "master_degs", "GeneSymbol_Key", "TEXT")
  
  dbExecute(con, "
    UPDATE master_degs
    SET GeneSymbol_Key = UPPER(TRIM(GeneSymbol))
    WHERE (GeneSymbol_Key IS NULL OR GeneSymbol_Key = '')
      AND GeneSymbol IS NOT NULL
      AND TRIM(GeneSymbol) <> '';
  ")
  
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_master_degs_dataset ON master_degs(dataset_tag);")
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_master_degs_symbol ON master_degs(GeneSymbol);")
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_master_degs_symbol_key ON master_degs(GeneSymbol_Key);")
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_master_degs_runstamp ON master_degs(run_stamp);")
}

############################## UI ##########################################

ui <- navbarPage(
  title = "R-ODAF",
  id = "mainNav",
  
  tabPanel(
    "Home",
    fluidPage(
      tags$head(
        tags$style(HTML("
          body { overflow-y: scroll; }

          .home-hero {
            padding: 28px 10px 18px 10px;
          }

          .home-title {
            font-size: 44px;
            font-weight: 600;
            margin-bottom: 10px;
            color: #1f2d3d;
          }

          .home-subtitle {
            font-size: 20px;
            color: #4d4d4d;
            margin-bottom: 10px;
          }

          .home-text {
            font-size: 16px;
            color: #5a5a5a;
            max-width: 980px;
            margin-bottom: 22px;
            line-height: 1.7;
          }

          .workflow-wrap {
            display: flex;
            flex-wrap: wrap;
            gap: 18px;
            margin-top: 8px;
            margin-bottom: 24px;
          }

          .workflow-card {
            flex: 1 1 220px;
            min-width: 220px;
            max-width: 260px;
            background: #f8fbff;
            border: 1px solid #cfe2ff;
            border-radius: 12px;
            overflow: hidden;
            box-shadow: 0 1px 3px rgba(0,0,0,0.05);
          }

          .workflow-icon {
            height: 110px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 42px;
            color: #1f77b4;
            border-bottom: 1px solid #cfe2ff;
            background: #eef6ff;
          }

          .workflow-body {
            padding: 18px 18px 16px 18px;
          }

          .workflow-body h4 {
            margin-top: 0;
            margin-bottom: 10px;
            color: #1f2d3d;
            font-size: 18px;
            font-weight: 600;
          }

          .workflow-body p {
            margin-bottom: 0;
            color: #4d4d4d;
            font-size: 15px;
            line-height: 1.6;
          }

          .home-about-card {
            background: #cfe2ff;
            border: 1px solid #9ec5fe;
            border-radius: 12px;
            padding: 24px;
            margin-bottom: 20px;
          }

          .home-about-card h3 {
            margin-top: 0;
            margin-bottom: 12px;
            color: #1f2d3d;
          }

          .home-about-card p {
            margin-bottom: 0;
            font-size: 16px;
            color: #4d4d4d;
            line-height: 1.7;
            max-width: 1100px;
          }

          .home-card {
            background: #eaf4ff;
            border: 1px solid #cfe2ff;
            border-radius: 12px;
            padding: 20px;
            margin-bottom: 20px;
            min-height: 210px;
          }

          .home-card h4 {
            margin-top: 0;
            margin-bottom: 12px;
            color: #1f2d3d;
          }

          .home-card p {
            margin-bottom: 10px;
            line-height: 1.6;
          }

          .home-card ul {
            padding-left: 18px;
            margin-bottom: 0;
            line-height: 1.6;
          }

          .home-note {
            margin-top: 12px;
            font-size: 14px;
            color: #5a5a5a;
          }

          .home-button-row {
            margin-top: 6px;
            margin-bottom: 8px;
          }

          .supported-by-wrap {
            margin-top: 28px;
            padding-top: 20px;
            padding-bottom: 8px;
            border-top: 1px solid #d9e6f2;
          }

          .supported-by-title {
            text-align: center;
            font-size: 16px;
            font-weight: 600;
            color: #4d4d4d;
            margin-bottom: 16px;
          }

          .supported-by-logos {
            display: flex;
            justify-content: center;
            align-items: center;
            flex-wrap: wrap;
            gap: 24px;
          }

          .supported-logo-box {
            background: #ffffff;
            border: 1px solid #e3edf7;
            border-radius: 12px;
            padding: 18px 24px;
            min-width: 230px;
            min-height: 110px;
            display: flex;
            align-items: center;
            justify-content: center;
            box-shadow: 0 1px 3px rgba(0,0,0,0.04);
          }

          .supported-logo {
            max-height: 72px;
            max-width: 210px;
            width: auto;
            height: auto;
            object-fit: contain;
            display: block;
          }
        "))
      ),
      
      div(
        class = "home-hero",
        div(class = "home-title", "Welcome to R-ODAF"),
        div(class = "home-subtitle", "A standardized Shiny workflow for RNA-seq differential expression analysis."),
        
        div(
          class = "workflow-wrap",
          
          div(
            class = "workflow-card",
            div(class = "workflow-icon", icon("upload")),
            div(
              class = "workflow-body",
              h4("Data upload"),
              p("Upload counts and metadata files for one or more datasets.")
            )
          ),
          
          div(
            class = "workflow-card",
            div(class = "workflow-icon", icon("search")),
            div(
              class = "workflow-body",
              h4("Quality control"),
              p("Check sample coverage, inspect PCA, and review outlier detection.")
            )
          ),
          
          div(
            class = "workflow-card",
            div(class = "workflow-icon", icon("sliders-h")),
            div(
              class = "workflow-body",
              h4("Differential expression"),
              p("Run the fixed DESeq2 workflow for standard or timepoint-specific contrasts.")
            )
          ),
          
          div(
            class = "workflow-card",
            div(class = "workflow-icon", icon("table")),
            div(
              class = "workflow-body",
              h4("Results exploration"),
              p("Review DEGs, heatmaps, gene expression plots, volcano plots, comparisons, and optional pathway analysis.")
            )
          )
        ),
        
        div(
          class = "home-button-row",
          actionButton("goToAnalysis", "Start analysis", class = "btn btn-primary btn-lg")
        )
      ),
      
      hr(),
      
      fluidRow(
        column(
          width = 12,
          div(
            class = "home-about-card",
            h3("About R-ODAF"),
            p("R-ODAF-Shiny is designed to make standardised RNA-seq differential expression analysis accessible to a wide range of users. Built on the R-ODAF framework for regulatory application, it provides a fixed and transparent workflow for quality control, DEG analysis, visualisation, and downloadable results. By combining accessibility with reproducibility, the app supports more consistent transcriptomic analysis across datasets.")
          )
        )
      ),
      
      fluidRow(
        column(
          width = 12,
          div(
            class = "supported-by-wrap",
            div(class = "supported-by-title", "Supported by"),
            div(
              class = "supported-by-logos",
              div(
                class = "supported-logo-box",
                tags$img(src = "logos/logo_1.png", class = "supported-logo", alt = "Logo 1")
              ),
              div(
                class = "supported-logo-box",
                tags$img(src = "logos/logo_2.png", class = "supported-logo", alt = "Logo 2")
              ),
              div(
                class = "supported-logo-box",
                tags$img(src = "logos/logo_3.png", class = "supported-logo", alt = "Logo 3")
              )
            )
          )
        )
      ),
      
      br()
    )
  ),
  
  tabPanel(
    "Analysis",
    fluidPage(
      titlePanel("R-ODAF DESeq2 Pipeline"),
      sidebarLayout(
        sidebarPanel(
          numericInput(
            "n_datasets",
            "Number of datasets",
            value = 1,
            min = 1,
            max = 5,
            step = 1
          ),
          uiOutput("dataset_selector_ui"),
          tags$hr(),
          uiOutput("dataset_file_inputs"),
          tags$hr(),
          tags$h5("Fixed filter settings:"),
          tags$p(strong("Min total reads / sample:"), format(MIN_COVERAGE, scientific = FALSE)),
          tags$p(strong("Min CPM threshold:"), MIN_CPM),
          tags$p(strong("FDR cutoff (alpha):"), FDR_CUTOFF),
          tags$p(strong("Outlier variance cutoff (%):"), VAR_THRESHOLD),
          tags$hr(),
          checkboxInput(
            "enable_pathway_analysis",
            "Add pathway analysis",
            value = FALSE
          ),
          conditionalPanel(
            condition = "input.enable_pathway_analysis",
            tagList(
              tags$h5("Pathway analysis:"),
              tags$p(strong("Input DEG threshold:"), "padj < 0.01"),
              tags$p(strong("Pathway result filter:"), "padj < 0.05"),
              tags$hr()
            )
          ),
          actionButton("runAnalysis", "Run analysis", class = "btn btn-primary"),
          width = 3
        ),
        
        mainPanel(
          uiOutput("analysis_tabs_ui")
        )
      )
    )
  ),
  
  tabPanel(
    "Tutorial",
    fluidPage(
      br(),
      if (file.exists("www/tutorials.md")) {
        includeMarkdown("www/tutorials.md")
      } else {
        tagList(
          h3("Tutorials"),
          p("The tutorial file 'www/tutorials.md' was not found."),
          p("Place your markdown tutorial in a folder called 'www' next to this app script,"),
          p("and ensure the images are in 'www/tutorials/'.")
        )
      }
    )
  ),
  
  tabPanel(
    "Documentation",
    fluidPage(
      br(),
      h3("Documentation"),
      p("R-ODAF-Shiny is a standardised Shiny tool for RNA-seq differential expression analysis. It applies a fixed DESeq2-based workflow so users can move from counts and metadata to quality-controlled DEGs, visualisations, and downloadable results in a transparent and reproducible way. The app also supports analysis of multiple datasets within the same session."),
      
      h4("When to use this app"),
      tags$ul(
        tags$li("Use it for RNA-seq experiments with a counts matrix and a metadata table."),
        tags$li("It is suitable for standard two-group comparisons such as control versus treated."),
        tags$li("It can also be used for timepoint-specific dose contrasts by subsetting selected timepoints before analysis."),
        tags$li("It is designed for reproducible DEG analysis, visualisation, and downstream interpretation.")
      ),
      
      h4("What you need"),
      tags$ul(
        tags$li("A counts file in CSV or tab-delimited format with genes as rows and samples as columns."),
        tags$li("A metadata file in CSV or tab-delimited format with samples as rows and variables as columns."),
        tags$li("Matching sample IDs between the counts matrix and metadata table."),
        tags$li("A design column that defines the comparison, for example control and case, or time and dose columns for timepoint-specific contrasts.")
      ),
      
      h4("What the app does"),
      tags$ul(
        tags$li("Checks sample coverage."),
        tags$li("Performs PCA-based quality control and outlier detection."),
        tags$li("Filters low-expression genes using fixed criteria."),
        tags$li("Runs a locked DESeq2 differential expression workflow."),
        tags$li("Returns DEGs, normalised counts, visualisations, and downloadable output tables."),
        tags$li("Stores DEG results in a session-only in-memory Results table for within-session comparison."),
        tags$li("Optionally performs pathway analysis with enrichR on filtered DEG input."),
        tags$li("Includes bundled example data for quick testing of the workflow.")
      ),
      
      h4("Fixed analysis settings"),
      tags$ul(
        tags$li(paste("Minimum coverage per sample:", format(MIN_COVERAGE, scientific = FALSE))),
        tags$li(paste("Minimum CPM:", MIN_CPM)),
        tags$li(paste("Main DEG FDR cutoff:", FDR_CUTOFF)),
        tags$li(paste("Outlier variance cutoff (%):", VAR_THRESHOLD)),
        tags$li(paste("Pathway input DEG FDR cutoff:", PATHWAY_INPUT_FDR)),
        tags$li(paste("Pathway result FDR cutoff:", PATHWAY_RESULT_FDR))
      ),
      
      h4("Full documentation"),
      p(
        "For detailed step-by-step guidance, examples, and additional notes, see the Tutorial tab and ",
        tags$a(
          href = "https://docs.vhp4safety.nl/en/latest/tutorials/r_odaf/r_odaf_intro.html",
          "the full R-ODAF documentation page",
          target = "_blank"
        ),
        "."
      )
    )
  ),
  
  tabPanel(
    "Citation",
    fluidPage(
      br(),
      h3("How to cite"),
      p("Formal citation information for this app is currently in preparation and will be added in a future update.")
    )
  ),
  
  tabPanel(
    "License",
    fluidPage(
      br(),
      h3("License"),
      tags$pre(
        "MIT License\n\n",
        "Copyright (c) 2025 Saad Lodhi\n\n",
        "Permission is hereby granted, free of charge, to any person obtaining a copy\n",
        "of this software and associated documentation files (the \"Software\"), to deal\n",
        "in the Software without restriction, including without limitation the rights\n",
        "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n",
        "copies of the Software, and to permit persons to whom the Software is\n",
        "furnished to do so, subject to the following conditions:\n\n",
        "The above copyright notice and this permission notice shall be included in all\n",
        "copies or substantial portions of the Software.\n\n",
        "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n",
        "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n",
        "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n",
        "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n",
        "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,\n",
        "ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n",
        "SOFTWARE.",
        style = "white-space: pre-wrap;"
      )
    )
  ),
  
  tabPanel(
    "Contact / Help",
    fluidPage(
      br(),
      h3("Contact / Help"),
      p("For questions, suggestions, or bug reports:"),
      tags$ul(
        tags$li(
          "GitHub issues: ",
          tags$a(
            href = "https://github.com/TGX-UM/R-ODAF-Shiny/issues/new",
            "start a GitHub issue",
            target = "_blank"
          )
        ),
        tags$li(
          "View existing issues: ",
          tags$a(
            href = "https://github.com/TGX-UM/R-ODAF-Shiny/issues",
            "open issue list",
            target = "_blank"
          )
        ),
        tags$li(
          "Email: ",
          tags$a(href = "mailto:saad.lodhi@maastrichtuniversity.nl", "saad.lodhi@maastrichtuniversity.nl")
        )
      )
    )
  )
)

############################## server ######################################

server <- function(input, output, session) {
  
  WINDOW_SIZE <- 5
  
  .get_counts_source <- function(i) {
    if (isTRUE(input[[paste0("use_example_", i)]])) {
      validate(
        need(
          .example_data_ready(),
          paste0(
            "Bundled example files not found. Expected:\n",
            EXAMPLE_COUNTS_PATH, "\n",
            EXAMPLE_META_PATH
          )
        )
      )
      .make_internal_fileinfo(EXAMPLE_COUNTS_PATH)
    } else {
      input[[paste0("counts_file_", i)]]
    }
  }
  
  .get_meta_source <- function(i) {
    if (isTRUE(input[[paste0("use_example_", i)]])) {
      validate(
        need(
          .example_data_ready(),
          paste0(
            "Bundled example files not found. Expected:\n",
            EXAMPLE_COUNTS_PATH, "\n",
            EXAMPLE_META_PATH
          )
        )
      )
      .make_internal_fileinfo(EXAMPLE_META_PATH)
    } else {
      input[[paste0("meta_file_", i)]]
    }
  }
  
  .render_heatmap_png <- function(mat, title_text, file, res = 150) {
    width_px  <- min(2400, max(1400, 90 * ncol(mat)))
    height_px <- min(5000, max(900, 12 * nrow(mat)))
    
    png(file, width = width_px, height = height_px, res = res)
    on.exit(dev.off(), add = TRUE)
    
    pheatmap(
      mat,
      scale = "row",
      show_rownames = FALSE,
      main = title_text,
      border_color = NA,
      fontsize_col = 10
    )
    
    invisible(list(width_px = width_px, height_px = height_px))
  }
  
  observeEvent(input$goToAnalysis, {
    updateNavbarPage(session, "mainNav", selected = "Analysis")
  })
  
  con_sql <- dbConnect(SQLite(), ":memory:")
  .init_sql_warehouse(con_sql)
  
  session$onSessionEnded(function() {
    try(dbDisconnect(con_sql), silent = TRUE)
  })
  
  output$analysis_tabs_ui <- renderUI({
    tabs <- list(
      tabPanel(
        "Data preview",
        h4("Counts"),
        DTOutput("counts_preview"),
        h4("Metadata"),
        DTOutput("meta_preview")
      ),
      
      tabPanel(
        "PCA & outliers",
        h4("PCA before outlier removal"),
        plotlyOutput("pca_before"),
        fluidRow(
          column(
            width = 4,
            selectInput(
              "pca_before_format",
              "Download format",
              choices = c("PNG" = "png", "TIFF" = "tiff", "SVG" = "svg", "JPEG" = "jpeg"),
              selected = "png"
            ),
            downloadButton("download_pca_before", "Download PCA (before)")
          )
        ),
        verbatimTextOutput("outlier_log"),
        h4("PCA after outlier removal"),
        plotlyOutput("pca_after"),
        fluidRow(
          column(
            width = 4,
            selectInput(
              "pca_after_format",
              "Download format",
              choices = c("PNG" = "png", "TIFF" = "tiff", "SVG" = "svg", "JPEG" = "jpeg"),
              selected = "png"
            ),
            downloadButton("download_pca_after", "Download PCA (after)")
          )
        )
      ),
      
      tabPanel(
        "DE results",
        verbatimTextOutput("deg_summary"),
        DTOutput("deg_table"),
        br(),
        fluidRow(
          column(
            width = 6,
            h5("DEGs"),
            radioButtons(
              "deg_format",
              "Format",
              choices = c("CSV" = "csv", "Tab-separated TXT" = "txt"),
              selected = "csv",
              inline = TRUE
            ),
            downloadButton("download_degs_any", "Download DEGs")
          ),
          column(
            width = 6,
            h5("Normalized counts"),
            radioButtons(
              "norm_format",
              "Format",
              choices = c("CSV" = "csv", "Tab-separated TXT" = "txt"),
              selected = "csv",
              inline = TRUE
            ),
            downloadButton("download_norm_any", "Download Norm Counts")
          )
        ),
        br(),
        fluidRow(
          column(
            width = 6,
            h5("FULL results (all genes)"),
            downloadButton("download_full_results", "Download FULL results")
          ),
          column(
            width = 6,
            h5("VST matrix (all genes)"),
            downloadButton("download_vst", "Download VST")
          )
        ),
        br(),
        fluidRow(
          column(
            width = 6,
            h5("genes_tested.txt"),
            downloadButton("download_genes_tested", "Download genes tested")
          ),
          column(
            width = 6,
            h5("samples_used.txt"),
            downloadButton("download_samples_used", "Download samples used")
          )
        ),
        
        br(),
        h4("Heat-map: selected DEGs"),
        selectInput(
          "heatmap_deg_n",
          "Number of top DEGs to display",
          choices = c(
            "10" = 10,
            "25" = 25,
            "50" = 50,
            "100" = 100,
            "250" = 250,
            "500" = 500,
            "All" = "all"
          ),
          selected = 100
        ),
        
        imageOutput("heatmap_degs", height = "900px"),
        br(),
        fluidRow(
          column(
            width = 4,
            selectInput(
              "heat_degs_format",
              "Download format",
              choices = c("PNG" = "png", "TIFF" = "tiff", "SVG" = "svg", "JPEG" = "jpeg"),
              selected = "png"
            ),
            downloadButton("download_heat_degs", "Download heatmap (DEGs)")
          )
        )
      ),
      
      tabPanel(
        "Gene expression",
        tabsetPanel(
          tabPanel(
            "Top gene expression",
            sliderInput("nTopGenes", "Total DEGs to consider", min = 5, max = 100, step = 5, value = 20),
            uiOutput("topRangeUI"),
            radioButtons(
              "top_label_type",
              "Gene labels",
              choices = c("Ensembl ID" = "ens", "Gene symbol" = "sym"),
              selected = "ens",
              inline = TRUE
            ),
            plotlyOutput("topGenePlots"),
            br(),
            fluidRow(
              column(
                width = 4,
                selectInput(
                  "top_expr_format",
                  "Download format",
                  choices = c("PNG" = "png", "TIFF" = "tiff", "SVG" = "svg", "JPEG" = "jpeg"),
                  selected = "png"
                ),
                downloadButton("download_top_expr", "Download plot")
              )
            )
          ),
          
          tabPanel(
            "Average expression",
            sliderInput("nAvgGenes", "Total DEGs to consider", min = 5, max = 100, step = 5, value = 20),
            uiOutput("avgRangeUI"),
            radioButtons(
              "avg_label_type",
              "Gene labels",
              choices = c("Ensembl ID" = "ens", "Gene symbol" = "sym"),
              selected = "ens",
              inline = TRUE
            ),
            plotlyOutput("avgExprPlot"),
            br(),
            fluidRow(
              column(
                width = 4,
                selectInput(
                  "avg_expr_format",
                  "Download format",
                  choices = c("PNG" = "png", "TIFF" = "tiff", "SVG" = "svg", "JPEG" = "jpeg"),
                  selected = "png"
                ),
                downloadButton("download_avg_expr", "Download plot")
              )
            )
          ),
          
          tabPanel(
            "Volcano plot",
            radioButtons(
              "volcano_label_type",
              "Gene labels",
              choices = c("Ensembl ID" = "ens", "Gene symbol" = "sym"),
              selected = "sym",
              inline = TRUE
            ),
            plotlyOutput("volcanoPlot"),
            br(),
            fluidRow(
              column(
                width = 4,
                selectInput(
                  "volcano_format",
                  "Download format",
                  choices = c("PNG" = "png", "TIFF" = "tiff", "SVG" = "svg", "JPEG" = "jpeg"),
                  selected = "png"
                ),
                downloadButton("download_volcano", "Download plot")
              )
            )
          )
          
        )
      ),
      
      tabPanel(
        "Results",
        h4("Results (session only)"),
        p("This table shows all DEG results archived in your in-memory SQL database during the current session only (closing the app clears it)."),
        DTOutput("warehouse_table"),
        hr(),
        h4("Results Comparison"),
        p("Dataset comparison: compare Study A and Study B using DEGs stored in the Results table."),
        fluidRow(
          column(4, uiOutput("select_meta_A")),
          column(4, uiOutput("select_meta_B")),
          column(
            4,
            selectInput(
              "meta_question",
              "Discovery Question",
              choices = c(
                "Shared DEGs (A ∩ B)" = "shared",
                "Shared, same direction (Concordant)" = "concordant",
                "Shared, opposite direction (Discordant)" = "discordant",
                "Unique to Study A (A \\ B)" = "uniqueA",
                "Unique to Study B (B \\ A)" = "uniqueB",
                "Largest effect-size shift (|FC_A - FC_B| ≥ 1)" = "shift"
              )
            )
          )
        ),
        actionButton("runMetaDiscovery", "Run Discovery Query", class = "btn btn-success"),
        hr(),
        DTOutput("meta_discovery_table")
      )
    )
    
    if (isTRUE(input$enable_pathway_analysis)) {
      tabs <- c(
        tabs,
        list(
          tabPanel(
            "Pathway analysis",
            h4("Pathway analysis"),
            p("Pathway analysis results for the currently selected dataset are shown below."),
            p("Input DEG threshold: padj < 0.01. Pathway result filter: padj < 0.05."),
            verbatimTextOutput("pathway_summary"),
            DTOutput("pathway_table"),
            br(),
            fluidRow(
              column(
                width = 6,
                downloadButton("download_pathway_results", "Download pathway results")
              ),
              column(
                width = 6,
                downloadButton("download_pathway_input_genes", "Download pathway input genes")
              )
            )
          )
        )
      )
    }
    
    tabs <- c(
      tabs,
      list(
        tabPanel(
          "Summary",
          h4("Analysis summary"),
          verbatimTextOutput("analysis_summary"),
          br(),
          downloadButton("download_report", "Download summary (.txt)"),
          downloadButton("download_sessioninfo", "Download session info (.txt)")
        )
      )
    )
    
    do.call(tabsetPanel, tabs)
  })
  
  .pca_plot_gg <- function(dds, dc, condition1, condition2, title = NULL) {
    vsd <- vst(dds, blind = TRUE) |> assay()
    md <- as.data.frame(colData(dds))
    md$sample <- rownames(md)
    
    ctrl_label <- as.character(condition1)
    case_label <- as.character(condition2)
    
    md$group_legend <- ifelse(
      as.character(md[[dc]]) == condition1,
      ctrl_label,
      case_label
    )
    md$group_legend <- factor(md$group_legend, levels = c(ctrl_label, case_label))
    
    pca_obj <- pca(vsd, metadata = md, removeVar = 0.1)
    scores <- as.data.frame(pca_obj$rotated[, c("PC1", "PC2"), drop = FALSE])
    scores$sample <- rownames(scores)
    scores$group_legend <- md$group_legend[match(scores$sample, md$sample)]
    
    var <- pca_obj$variance
    var1 <- if (!is.null(names(var)) && "PC1" %in% names(var)) as.numeric(var["PC1"]) else as.numeric(var[1])
    var2 <- if (!is.null(names(var)) && "PC2" %in% names(var)) as.numeric(var["PC2"]) else as.numeric(var[2])
    
    xlab <- paste0("PC1 (", round(var1, 1), "%)")
    ylab <- paste0("PC2 (", round(var2, 1), "%)")
    
    scores$tooltip <- paste0(
      "Sample: ", scores$sample,
      "<br>Group: ", scores$group_legend,
      "<br>PC1: ", signif(scores$PC1, 4),
      "<br>PC2: ", signif(scores$PC2, 4)
    )
    
    ggplot(scores, aes(PC1, PC2, color = group_legend, text = tooltip)) +
      geom_point(size = 3, alpha = 0.9) +
      geom_text(
        aes(label = sample),
        vjust = -0.8,
        size = 3,
        check_overlap = TRUE,
        show.legend = FALSE,
        color = "black"
      ) +
      labs(x = xlab, y = ylab, color = "Group", title = title) +
      coord_fixed() +
      theme_minimal() +
      theme(
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black")
      )
  }
  
  .volcano_plot_gg <- function(dat) {
    dat$Significance <- factor(
      dat$Significance,
      levels = c("Down", "Not significant", "Up")
    )
    
    volcano_colors <- c(
      "Down" = "#00BFC4",
      "Not significant" = "grey70",
      "Up" = "#F8766D"
    )
    
    ggplot(
      dat,
      aes(
        x = log2FoldChange,
        y = neglog10_padj,
        color = Significance,
        text = tooltip
      )
    ) +
      geom_point(alpha = 0.7, na.rm = TRUE) +
      scale_color_manual(values = volcano_colors) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed") +
      labs(
        x = "log2 fold change",
        y = "-log10 adjusted p-value",
        color = "Status"
      ) +
      theme_minimal()
  }  
  output$dataset_file_inputs <- renderUI({
    n <- input$n_datasets
    tagList(
      lapply(seq_len(n), function(i) {
        std_cond <- sprintf("!input.use_time_dose_%d", i)
        td_cond  <- sprintf("input.use_time_dose_%d", i)
        
        wellPanel(
          tags$h4(paste0("Dataset ", i)),
          textInput(paste0("dataset_name_", i), "Dataset name", value = paste0("dataset_", i)),
          selectInput(
            paste0("organism_", i),
            "Organism (for gene symbols)",
            choices = c("Human" = "human", "Mouse" = "mouse", "Rat" = "rat"),
            selected = "human"
          ),
          checkboxInput(
            paste0("use_example_", i),
            "Use example dataset",
            value = FALSE
          ),
          conditionalPanel(
            condition = sprintf("input.use_example_%d", i),
            tags$small("")
          ),
          conditionalPanel(
            condition = sprintf("!input.use_example_%d", i),
            fileInput(paste0("counts_file_", i), "Counts (genes × samples)", accept = c(".csv", ".txt")),
            fileInput(paste0("meta_file_", i), "Metadata (samples × attributes)", accept = c(".csv", ".txt"))
          ),
          checkboxInput(
            paste0("use_time_dose_", i),
            "Enable timepoint-specific contrasts (Dose within selected timepoints)",
            value = FALSE
          ),
          conditionalPanel(
            condition = std_cond,
            textInput(paste0("design_col_", i), "Design column", value = "treatment"),
            textInput(paste0("condition1_", i), "Control label", value = "control"),
            textInput(paste0("condition2_", i), "Exposure label", value = "exposure")
          ),
          conditionalPanel(
            condition = td_cond,
            textInput(paste0("time_col_", i), "Time column", value = "time"),
            textInput(paste0("dose_col_", i), "Dose column", value = "dose"),
            uiOutput(paste0("time_dose_selectors_", i))
          )
        )
      })
    )
  })
  
  output$dataset_selector_ui <- renderUI({
    n <- input$n_datasets
    if (n <= 1) return(NULL)
    
    labels <- vapply(seq_len(n), function(i) {
      nm <- input[[paste0("dataset_name_", i)]]
      if (!is.null(nm) && nzchar(nm)) {
        paste0("Dataset ", i, ": ", nm)
      } else {
        paste0("Dataset ", i)
      }
    }, character(1))
    
    choices <- setNames(as.character(seq_len(n)), labels)
    
    selectInput(
      "active_dataset",
      "Dataset to display",
      choices = choices,
      selected = "1"
    )
  })
  
  active_idx <- reactive({
    if (!is.null(input$active_dataset)) as.integer(input$active_dataset) else 1L
  })
  
  observeEvent(input$n_datasets, {
    if (!is.null(input$active_dataset)) {
      if (as.integer(input$active_dataset) > input$n_datasets) {
        updateSelectInput(session, "active_dataset", selected = "1")
      }
    }
  })
  
  for (i in 1:5) {
    local({
      ii <- i
      output[[paste0("time_dose_selectors_", ii)]] <- renderUI({
        if (is.null(input$n_datasets) || input$n_datasets < ii) return(NULL)
        if (!isTRUE(input[[paste0("use_time_dose_", ii)]])) return(NULL)
        
        meta_file <- .get_meta_source(ii)
        req(meta_file)
        meta0 <- .read_meta(meta_file)
        
        time_col <- input[[paste0("time_col_", ii)]]
        dose_col <- input[[paste0("dose_col_", ii)]]
        
        validate(
          need(!is.null(time_col) && nzchar(time_col), "Time column name is required."),
          need(!is.null(dose_col) && nzchar(dose_col), "Dose column name is required."),
          need(time_col %in% colnames(meta0), paste0("Time column not found in metadata: ", time_col)),
          need(dose_col %in% colnames(meta0), paste0("Dose column not found in metadata: ", dose_col))
        )
        
        time_vals <- sort(unique(as.character(meta0[[time_col]])))
        dose_vals <- sort(unique(as.character(meta0[[dose_col]])))
        
        validate(
          need(length(time_vals) >= 1, "No time values detected in the selected Time column."),
          need(length(dose_vals) >= 2, "Dose column must contain at least 2 unique values (Control and Case).")
        )
        
        default_ctrl <- if ("control" %in% dose_vals) "control" else dose_vals[1]
        remaining <- setdiff(dose_vals, default_ctrl)
        default_case <- if ("high" %in% dose_vals) "high" else if (length(remaining) >= 1) remaining[1] else dose_vals[1]
        
        tagList(
          selectizeInput(
            paste0("time_vals_", ii),
            "Timepoints to include",
            choices = time_vals,
            selected = time_vals,
            multiple = TRUE,
            options = list(
              plugins = list("remove_button"),
              placeholder = "Select one or more timepoints"
            )
          ),
          selectInput(
            paste0("control_dose_", ii),
            "Control dose label",
            choices = dose_vals,
            selected = default_ctrl
          ),
          selectInput(
            paste0("case_dose_", ii),
            "Case dose label",
            choices = dose_vals,
            selected = default_case
          ),
          tags$small("This filters samples to the selected timepoint(s) and keeps ONLY the selected Control/Case dose groups before running DESeq2 (Case vs Control).")
        )
      })
    })
  }
  
  active_counts_data <- reactive({
    i <- active_idx()
    fi <- .get_counts_source(i)
    req(fi)
    .read_counts(fi)
  })
  
  active_meta_data <- reactive({
    i <- active_idx()
    fi <- .get_meta_source(i)
    req(fi)
    .read_meta(fi)
  })
  
  output$counts_preview <- renderDT({
    req(active_counts_data())
    datatable(active_counts_data(), options = list(scrollX = TRUE, pageLength = 5))
  })
  
  output$meta_preview <- renderDT({
    req(active_meta_data())
    datatable(active_meta_data(), options = list(scrollX = TRUE, pageLength = 5))
  })
  
  .run_single_dataset <- function(
    cnts0,
    meta0,
    dataset_name,
    organism_used,
    design_col,
    condition1,
    condition2,
    subset_col = NULL,
    subset_values = NULL,
    contrast_desc = "",
    filter_desc = "",
    pathway_input_fdr = PATHWAY_INPUT_FDR,
    run_pathway_analysis = FALSE
  ) {
    validate(
      need(
        tryCatch({ .get_orgdb(organism_used); TRUE }, error = function(e) FALSE),
        paste0(
          "Annotation package missing for selected Organism: ", organism_used, ". ",
          "Install the corresponding org.*.eg.db package or switch Organism."
        )
      )
    )
    
    common <- intersect(colnames(cnts0), rownames(meta0))
    validate(
      need(
        length(common) >= 2,
        paste0("[", dataset_name, "] No overlapping sample IDs between counts (columns) and metadata (row names).")
      )
    )
    
    cnts <- cnts0[, common, drop = FALSE]
    meta <- meta0[common, , drop = FALSE]
    
    cnts[is.na(cnts)] <- 0
    
    if (!is.null(subset_col) && nzchar(subset_col)) {
      validate(
        need(
          subset_col %in% colnames(meta),
          paste0("[", dataset_name, "] Subset column not found in metadata: ", subset_col)
        )
      )
      if (!is.null(subset_values) && length(subset_values) > 0) {
        keep_rows <- as.character(meta[[subset_col]]) %in% as.character(subset_values)
        meta <- meta[keep_rows, , drop = FALSE]
        cnts <- cnts[, rownames(meta), drop = FALSE]
      }
    }
    
    validate(
      need(nrow(meta) >= 2, paste0("[", dataset_name, "] Not enough samples after filtering/subsetting.")),
      need(design_col %in% colnames(meta), paste0("[", dataset_name, "] Design column not found in metadata: ", design_col))
    )
    
    meta[[design_col]] <- as.character(meta[[design_col]])
    keep_levels <- meta[[design_col]] %in% c(condition1, condition2)
    meta <- meta[keep_levels, , drop = FALSE]
    cnts <- cnts[, rownames(meta), drop = FALSE]
    
    validate(
      need(nrow(meta) >= 2, paste0("[", dataset_name, "] Not enough samples after keeping only contrast groups.")),
      need(
        all(c(condition1, condition2) %in% unique(meta[[design_col]])),
        paste0(
          "[", dataset_name, "] Both condition labels must be present in metadata column '",
          design_col, "' after filtering. Missing: ",
          paste(setdiff(c(condition1, condition2), unique(meta[[design_col]])), collapse = ", ")
        )
      )
    )
    
    meta[[design_col]] <- factor(meta[[design_col]], levels = c(condition1, condition2))
    
    keep_samp <- colSums(cnts) > MIN_COVERAGE
    cnts <- cnts[, keep_samp, drop = FALSE]
    meta <- meta[keep_samp, , drop = FALSE]
    
    validate(
      need(nrow(meta) >= 2, paste0("[", dataset_name, "] Not enough samples after coverage filtering (MIN_COVERAGE)."))
    )
    
    design_f <- as.formula(paste0("~", design_col))
    dds_init <- DESeqDataSetFromMatrix(round(cnts), meta, design_f)
    dds_init <- DESeq(dds_init, quiet = TRUE)
    
    dds_cur <- dds_init
    cnts_cur <- cnts
    meta_cur <- meta
    removed <- character()
    conds <- c(condition1, condition2)
    
    repeat {
      outs <- check_outliers_vst(dds_cur, as.data.frame(colData(dds_cur)), design_col, conds, VAR_THRESHOLD)
      if (!length(outs)) break
      
      prefixes <- unique(sub("_.*", "", outs))
      mask <- grepl(paste0("^(", paste(prefixes, collapse = "|"), ")_"), colnames(cnts_cur))
      if (!any(mask)) break
      
      removed <- c(removed, colnames(cnts_cur)[mask])
      cnts_cur <- cnts_cur[, !mask, drop = FALSE]
      meta_cur <- meta_cur[!mask, , drop = FALSE]
      
      validate(
        need(nrow(meta_cur) >= 2, paste0("[", dataset_name, "] Not enough samples after outlier removal."))
      )
      
      dds_cur <- DESeqDataSetFromMatrix(round(cnts_cur), meta_cur, design_f)
      dds_cur <- DESeq(dds_cur, quiet = TRUE)
    }
    
    norm_counts <- counts(dds_cur, normalized = TRUE)
    CPMnorm <- cpm(norm_counts)
    
    samp_tab <- table(meta_cur[[design_col]])
    
    keep_genes <- rownames(norm_counts)[
      sapply(rownames(norm_counts), function(g) {
        any(sapply(names(samp_tab), function(gr) {
          cols <- which(meta_cur[[design_col]] == gr)
          if (length(cols) == 0) return(FALSE)
          sum(CPMnorm[g, cols] >= MIN_CPM) >= 0.75 * length(cols)
        }))
      })
    ]
    
    validate(
      need(length(keep_genes) >= 2, paste0("[", dataset_name, "] Too few genes after CPM filtering."))
    )
    
    res_all <- results(
      dds_cur[keep_genes, ],
      contrast = c(design_col, condition2, condition1),
      alpha = FDR_CUTOFF,
      independentFiltering = FALSE,
      cooksCutoff = FALSE,
      pAdjustMethod = "fdr"
    )
    
    res_df <- as.data.frame(res_all)
    
    max_sig_fdr <- max(FDR_CUTOFF, pathway_input_fdr, na.rm = TRUE)
    sig_df_for_filters <- subset(res_df, !is.na(padj) & padj < max_sig_fdr)
    
    if (nrow(sig_df_for_filters) > 0) {
      Filt <- data.frame(
        low = rep(1, nrow(sig_df_for_filters)),
        quant = rep(0, nrow(sig_df_for_filters)),
        spike = rep(1, nrow(sig_df_for_filters)),
        row.names = rownames(sig_df_for_filters),
        stringsAsFactors = FALSE
      )
      
      for (g in rownames(sig_df_for_filters)) {
        gr1 <- which(meta_cur[[design_col]] == condition1)
        gr2 <- which(meta_cur[[design_col]] == condition2)
        if (length(gr1) == 0 || length(gr2) == 0) next
        
        m1 <- median(norm_counts[g, gr1]); q2 <- quantile(norm_counts[g, gr2], .75)
        m2 <- median(norm_counts[g, gr2]); q1 <- quantile(norm_counts[g, gr1], .75)
        
        Filt[g, "quant"] <- as.integer((m1 > q2) | (m2 > q1))
        
        bad_spike <- sapply(list(gr1, gr2), function(cols) {
          v <- norm_counts[g, cols]
          if (max(v) == 0) FALSE else (max(v) / sum(v)) >= 1.4 * (length(cols))^(-0.66)
        })
        
        Filt[g, "spike"] <- as.integer(sum(bad_spike) <= 1)
      }
      
      keep_final <- rownames(Filt)[rowSums(Filt) == 3]
      filtered_sig_df <- sig_df_for_filters[rownames(sig_df_for_filters) %in% keep_final, , drop = FALSE]
      degs_final <- filtered_sig_df[!is.na(filtered_sig_df$padj) & filtered_sig_df$padj < FDR_CUTOFF, , drop = FALSE]
      pathway_degs_final <- filtered_sig_df[!is.na(filtered_sig_df$padj) & filtered_sig_df$padj < pathway_input_fdr, , drop = FALSE]
    } else {
      degs_final <- sig_df_for_filters
      pathway_degs_final <- sig_df_for_filters
    }
    
    all_sym <- map_to_symbol(rownames(norm_counts), organism_used)
    norm_hgnc <- data.frame(
      GeneID = rownames(norm_counts),
      GeneSymbol = unname(all_sym),
      norm_counts,
      row.names = NULL,
      check.names = FALSE
    )
    
    if (nrow(res_df) > 0) {
      res_sym <- map_to_symbol(rownames(res_df), organism_used)
      all_res_hgnc <- data.frame(
        GeneID = rownames(res_df),
        GeneSymbol = unname(res_sym),
        res_df,
        row.names = NULL,
        check.names = FALSE
      )
    } else {
      all_res_hgnc <- data.frame()
    }
    
    if (nrow(degs_final) > 0) {
      deg_sym <- map_to_symbol(rownames(degs_final), organism_used)
      degs_hgnc <- data.frame(
        GeneID = rownames(degs_final),
        GeneSymbol = unname(deg_sym),
        degs_final,
        row.names = NULL,
        check.names = FALSE
      )
    } else {
      degs_hgnc <- data.frame()
    }
    
    if (nrow(pathway_degs_final) > 0) {
      pathway_deg_sym <- map_to_symbol(rownames(pathway_degs_final), organism_used)
      pathway_degs_hgnc <- data.frame(
        GeneID = rownames(pathway_degs_final),
        GeneSymbol = unname(pathway_deg_sym),
        pathway_degs_final,
        row.names = NULL,
        check.names = FALSE
      )
    } else {
      pathway_degs_hgnc <- data.frame()
    }
    
    if (isTRUE(run_pathway_analysis)) {
      pathway_obj <- .run_pathway_enrichr(
        pathway_degs_hgnc = pathway_degs_hgnc,
        organism = organism_used,
        input_fdr = pathway_input_fdr,
        result_fdr = PATHWAY_RESULT_FDR,
        min_genes = MIN_PATHWAY_GENES
      )
    } else {
      pathway_obj <- list(
        input_genes = character(0),
        results = data.frame(),
        libraries_used = character(0)
      )
    }
    
    list(
      dataset_name = dataset_name,
      organism_used = organism_used,
      design_col = design_col,
      condition1 = condition1,
      condition2 = condition2,
      contrast_desc = contrast_desc,
      filter_desc = filter_desc,
      dds_init = dds_init,
      dds_final = dds_cur,
      norm = norm_counts,
      norm_hgnc = norm_hgnc,
      all_res = res_df,
      all_res_hgnc = all_res_hgnc,
      degs = degs_final,
      degs_hgnc = degs_hgnc,
      pathway_degs_hgnc = pathway_degs_hgnc,
      pathway_input_fdr = pathway_input_fdr,
      pathway_input_genes = pathway_obj$input_genes,
      pathway_results = pathway_obj$results,
      pathway_libraries = pathway_obj$libraries_used,
      pathway_enabled = isTRUE(run_pathway_analysis),
      removed = unique(removed),
      n_samples_used = ncol(norm_counts)
    )
  }
  
  results_all <- eventReactive(input$runAnalysis, {
    run_start <- Sys.time()
    run_stamp <- .make_run_stamp(run_start)
    pathway_input_fdr <- PATHWAY_INPUT_FDR
    run_pathway_analysis <- isTRUE(input$enable_pathway_analysis)
    
    withProgress(message = "Running R-ODAF analysis", value = 0, {
      n <- input$n_datasets
      res_list <- vector("list", n)
      
      for (i in seq_len(n)) {
        ds_start <- Sys.time()
        incProgress((i - 1) / max(1, n), detail = paste0("Dataset ", i, ": reading inputs..."))
        
        dataset_name <- input[[paste0("dataset_name_", i)]]
        if (is.null(dataset_name) || !nzchar(dataset_name)) dataset_name <- paste0("dataset_", i)
        
        organism_used <- input[[paste0("organism_", i)]]
        counts_file <- .get_counts_source(i)
        meta_file   <- .get_meta_source(i)
        
        req(counts_file, meta_file, organism_used)
        
        cnts0 <- .read_counts(counts_file)
        meta0 <- .read_meta(meta_file)
        
        use_td <- isTRUE(input[[paste0("use_time_dose_", i)]])
        design_col_run <- input[[paste0("design_col_", i)]]
        condition1_run <- input[[paste0("condition1_", i)]]
        condition2_run <- input[[paste0("condition2_", i)]]
        
        subset_col <- NULL
        subset_vals <- NULL
        contrast_desc <- paste0(condition2_run, " vs ", condition1_run)
        filter_desc <- "None"
        contrast_tag <- ""
        
        if (use_td) {
          time_col   <- input[[paste0("time_col_", i)]]
          dose_col   <- input[[paste0("dose_col_", i)]]
          time_vals  <- input[[paste0("time_vals_", i)]]
          ctrl_dose  <- input[[paste0("control_dose_", i)]]
          case_dose  <- input[[paste0("case_dose_", i)]]
          
          req(time_col, dose_col, ctrl_dose, case_dose)
          
          if (is.null(time_vals) || length(time_vals) == 0) {
            validate(
              need(time_col %in% colnames(meta0), paste0("[", dataset_name, "] Time column not in metadata: ", time_col))
            )
            time_vals <- sort(unique(as.character(meta0[[time_col]])))
          }
          
          design_col_run <- dose_col
          condition1_run <- ctrl_dose
          condition2_run <- case_dose
          subset_col <- time_col
          subset_vals <- time_vals
          contrast_desc <- paste0(dose_col, ": ", case_dose, " vs ", ctrl_dose)
          filter_desc <- paste0(time_col, " ∈ {", paste(time_vals, collapse = ", "), "}")
          contrast_tag <- paste0("TP_", paste(time_vals, collapse = "+"), "__", case_dose, "_vs_", ctrl_dose)
        }
        
        req(design_col_run, condition1_run, condition2_run)
        
        incProgress((i - 1) / max(1, n) + (0.5 / max(1, n)), detail = paste0("Dataset ", i, ": running pipeline..."))
        
        one <- .run_single_dataset(
          cnts0 = cnts0,
          meta0 = meta0,
          dataset_name = dataset_name,
          organism_used = organism_used,
          design_col = design_col_run,
          condition1 = condition1_run,
          condition2 = condition2_run,
          subset_col = subset_col,
          subset_values = subset_vals,
          contrast_desc = contrast_desc,
          filter_desc = filter_desc,
          pathway_input_fdr = pathway_input_fdr,
          run_pathway_analysis = run_pathway_analysis
        )
        
        tag_base <- paste0("DS", i, "_", .sanitize_label(dataset_name))
        tag_ct <- if (nzchar(contrast_tag)) paste0("__", .sanitize_label(contrast_tag)) else ""
        one$dataset_tag <- paste0(tag_base, tag_ct)
        
        ds_end <- Sys.time()
        one$dataset_start <- ds_start
        one$dataset_end <- ds_end
        one$dataset_duration_sec <- as.numeric(difftime(ds_end, ds_start, units = "secs"))
        
        res_list[[i]] <- one
      }
      
      run_end <- Sys.time()
      run_duration_sec <- as.numeric(difftime(run_end, run_start, units = "secs"))
      
      .init_sql_warehouse(con_sql)
      
      for (j in seq_len(n)) {
        res_list[[j]]$run_start <- run_start
        res_list[[j]]$run_end <- run_end
        res_list[[j]]$run_stamp <- run_stamp
        res_list[[j]]$run_duration_sec <- run_duration_sec
        
        if (nrow(res_list[[j]]$degs_hgnc) > 0) {
          sql_df <- res_list[[j]]$degs_hgnc
          sql_df$dataset_tag <- res_list[[j]]$dataset_tag
          sql_df$organism <- res_list[[j]]$organism_used
          sql_df$run_stamp <- run_stamp
          sql_df$timestamp <- as.character(run_start)
          sql_df$GeneSymbol_Key <- .normalize_symbol(sql_df$GeneSymbol)
          
          if (dbExistsTable(con_sql, "master_degs")) {
            dbExecute(con_sql, sprintf("DELETE FROM master_degs WHERE dataset_tag = '%s'", sql_df$dataset_tag[1]))
          }
          
          dbWriteTable(con_sql, "master_degs", sql_df, append = TRUE, row.names = FALSE)
        }
      }
      
      incProgress(1, detail = "Done. Saved to Results (session only, in-memory).")
      showNotification("Results successfully archived to Results (SESSION-ONLY, in-memory).", type = "message")
      res_list
    })
  })
  
  observeEvent(input$runAnalysis, {
    results_all()
  })
  
  active_res <- reactive({
    req(results_all())
    i <- active_idx()
    req(i <= length(results_all()))
    results_all()[[i]]
  })
  
  active_design_col <- reactive({
    req(active_res())
    active_res()$design_col
  })
  
  pathway_table_data <- reactive({
    req(input$enable_pathway_analysis, active_res())
    tbl <- active_res()$pathway_results
    if (is.null(tbl) || nrow(tbl) == 0) {
      return(data.frame())
    }
    tbl$dataset_name <- active_res()$dataset_name
    tbl$dataset_tag <- active_res()$dataset_tag
    tbl$organism <- active_res()$organism_used
    tbl$input_fdr <- active_res()$pathway_input_fdr
    tbl
  })
  
  pathway_input_gene_table <- reactive({
    req(input$enable_pathway_analysis, active_res())
    if (is.null(active_res()$pathway_input_genes) || length(active_res()$pathway_input_genes) == 0) {
      return(data.frame())
    }
    data.frame(
      dataset_name = active_res()$dataset_name,
      dataset_tag = active_res()$dataset_tag,
      organism = active_res()$organism_used,
      input_fdr = active_res()$pathway_input_fdr,
      GeneSymbol = active_res()$pathway_input_genes,
      stringsAsFactors = FALSE
    )
  })
  
  heatmap_deg_matrix <- reactive({
    req(active_res())
    
    dg <- active_res()$degs
    if (nrow(dg) < 2) return(NULL)
    
    deg_order <- rownames(dg)[order(dg$padj, na.last = NA)]
    n_choice <- if (is.null(input$heatmap_deg_n)) "100" else as.character(input$heatmap_deg_n)
    
    genes_to_plot <- if (n_choice == "all") {
      deg_order
    } else {
      head(deg_order, as.integer(n_choice))
    }
    
    active_res()$norm[genes_to_plot, , drop = FALSE]
  })
  
  output$topRangeUI <- renderUI({
    sliderInput(
      "topRange",
      "Gene range",
      min = 1,
      max = input$nTopGenes,
      value = c(1, min(WINDOW_SIZE, input$nTopGenes)),
      step = 1
    )
  })
  
  output$avgRangeUI <- renderUI({
    sliderInput(
      "avgRange",
      "Gene range",
      min = 1,
      max = input$nAvgGenes,
      value = c(1, min(WINDOW_SIZE, input$nAvgGenes)),
      step = 1
    )
  })
  
  observeEvent(input$nTopGenes, {
    rng <- input$topRange
    if (is.null(rng)) rng <- c(1, WINDOW_SIZE)
    rng <- pmin(pmax(rng, 1), input$nTopGenes)
    updateSliderInput(session, "topRange", max = input$nTopGenes, value = rng)
  })
  
  observeEvent(input$nAvgGenes, {
    rng <- input$avgRange
    if (is.null(rng)) rng <- c(1, WINDOW_SIZE)
    rng <- pmin(pmax(rng, 1), input$nAvgGenes)
    updateSliderInput(session, "avgRange", max = input$nAvgGenes, value = rng)
  })
  
  ranked_genes <- reactive({
    req(active_res(), nrow(active_res()$degs) > 0)
    rownames(active_res()$degs)[order(active_res()$degs$padj)]
  })
  
  topGenesWin <- reactive({
    req(ranked_genes(), input$topRange)
    g <- head(ranked_genes(), input$nTopGenes)
    idx <- seq(input$topRange[1], min(input$topRange[2], length(g)))
    g[idx]
  })
  
  topAvgGenesWin <- reactive({
    req(ranked_genes(), input$avgRange)
    g <- head(ranked_genes(), input$nAvgGenes)
    idx <- seq(input$avgRange[1], min(input$avgRange[2], length(g)))
    g[idx]
  })
  
  symbol_map <- reactive({
    req(active_res())
    nh <- active_res()$norm_hgnc
    setNames(nh$GeneSymbol, nh$GeneID)
  })
  
  topExprDF <- reactive({
    req(topGenesWin(), active_res())
    mat <- active_res()$norm[topGenesWin(), , drop = FALSE]
    gene_ids <- rownames(mat)
    sym_map <- symbol_map()
    
    labels <- gene_ids
    if (!is.null(input$top_label_type) && input$top_label_type == "sym") {
      tmp <- sym_map[gene_ids]
      labels <- ifelse(is.na(tmp) | tmp == "", gene_ids, tmp)
    }
    
    df <- as.data.frame(mat)
    df$gene <- labels
    long <- pivot_longer(df, -gene, names_to = "sample", values_to = "expression")
    
    meta <- as.data.frame(colData(active_res()$dds_final))
    meta$sample <- rownames(meta)
    
    left_join(long, meta, by = "sample")
  })
  
  output$topGenePlots <- renderPlotly({
    req(topExprDF(), active_design_col(), active_res())
    dc <- active_design_col()
    
    top_colors <- setNames(
      c("#00BFC4", "#F8766D"),
      c(active_res()$condition1, active_res()$condition2)
    )
    
    p <- ggplot(topExprDF(), aes_string(x = dc, y = "expression", color = dc, text = "sample")) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      facet_wrap(
        ~ gene,
        ncol = WINDOW_SIZE,
        scales = "free_y",
        labeller = labeller(gene = function(x) str_wrap(x, width = 10))
      ) +
      scale_color_manual(values = top_colors) +
      scale_x_discrete(labels = NULL) +
      labs(x = dc, y = "Normalized count") +
      theme_minimal() +
      theme(
        strip.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(1.1, "lines")
      )
    
    ggplotly(p, tooltip = c("text", "gene", "expression"))
  })
  
  avgExprDF <- reactive({
    req(topAvgGenesWin(), active_res(), active_design_col())
    dc <- active_design_col()
    mat <- active_res()$norm[topAvgGenesWin(), , drop = FALSE]
    gene_ids <- rownames(mat)
    sym_map <- symbol_map()
    
    labels <- gene_ids
    if (!is.null(input$avg_label_type) && input$avg_label_type == "sym") {
      tmp <- sym_map[gene_ids]
      labels <- ifelse(is.na(tmp) | tmp == "", gene_ids, tmp)
    }
    
    label_vec <- setNames(labels, gene_ids)
    meta <- as.data.frame(colData(active_res()$dds_final))
    
    do.call(
      rbind,
      lapply(gene_ids, function(g) {
        mu <- tapply(mat[g, ], meta[[dc]], mean)
        data.frame(
          gene = label_vec[g],
          condition = names(mu),
          mean_expression = as.numeric(mu),
          stringsAsFactors = FALSE
        )
      })
    )
  })
  
  output$avgExprPlot <- renderPlotly({
    req(avgExprDF(), active_res())
    
    avg_dat <- avgExprDF()
    avg_dat$condition <- factor(
      avg_dat$condition,
      levels = c(active_res()$condition1, active_res()$condition2)
    )
    
    avg_colors <- setNames(
      c("#00BFC4", "#F8766D"),
      c(active_res()$condition1, active_res()$condition2)
    )
    
    p <- ggplot(avg_dat, aes_string(x = "gene", y = "mean_expression", fill = "condition")) +
      geom_col(position = position_dodge()) +
      scale_fill_manual(values = avg_colors) +
      labs(x = "Gene", y = "Mean normalized count") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.title.x = element_text(margin = margin(t = 10))
      )
    
    ggplotly(p, tooltip = c("gene", "condition", "mean_expression"))
  })
  
  volcano_df <- reactive({
    req(active_res())
    dat <- active_res()$all_res_hgnc
    if (is.null(dat) || nrow(dat) == 0) {
      return(data.frame())
    }
    
    label_type <- if (!is.null(input$volcano_label_type) && input$volcano_label_type == "sym") "sym" else "ens"
    
    dat <- dat %>%
      mutate(
        padj_safe = case_when(
          is.na(padj) ~ NA_real_,
          padj <= 0 ~ .Machine$double.xmin,
          TRUE ~ padj
        ),
        neglog10_padj = -log10(padj_safe),
        Significance = case_when(
          !is.na(padj) & padj < FDR_CUTOFF & !is.na(log2FoldChange) & log2FoldChange > 0 ~ "Up",
          !is.na(padj) & padj < FDR_CUTOFF & !is.na(log2FoldChange) & log2FoldChange < 0 ~ "Down",
          TRUE ~ "Not significant"
        ),
        plot_label = case_when(
          label_type == "sym" & !is.na(GeneSymbol) & GeneSymbol != "" ~ GeneSymbol,
          TRUE ~ GeneID
        ),
        tooltip = paste0(
          "GeneID: ", GeneID,
          ifelse(is.na(GeneSymbol) | GeneSymbol == "", "", paste0("<br>GeneSymbol: ", GeneSymbol)),
          "<br>log2FC: ", signif(log2FoldChange, 4),
          "<br>padj: ", signif(padj, 4),
          "<br>Status: ", Significance
        )
      )
    dat
  })
  

  output$volcanoPlot <- renderPlotly({
    req(volcano_df())
    dat <- volcano_df()
    validate(
      need(nrow(dat) > 0, "No results available for volcano plot.")
    )
    
    p <- .volcano_plot_gg(dat = dat)
    ggplotly(p, tooltip = "text")
  })
  
 
  output$pca_before <- renderPlotly({
    req(active_res(), active_design_col())
    dc <- active_design_col()
    g <- .pca_plot_gg(
      dds = active_res()$dds_init,
      dc = dc,
      condition1 = active_res()$condition1,
      condition2 = active_res()$condition2,
      title = "PCA before outlier removal"
    )
    
    ggplotly(g, tooltip = "text") %>%
      layout(legend = list(title = list(text = "Group"), font = list(color = "black")))
  })
  
  output$outlier_log <- renderPrint({
    req(active_res())
    if (length(active_res()$removed) == 0) {
      cat("No outliers removed.")
    } else {
      cat("Removed samples:", paste(active_res()$removed, collapse = ", "))
    }
  })
  
  output$pca_after <- renderPlotly({
    req(active_res(), active_design_col())
    dc <- active_design_col()
    g <- .pca_plot_gg(
      dds = active_res()$dds_final,
      dc = dc,
      condition1 = active_res()$condition1,
      condition2 = active_res()$condition2,
      title = "PCA after outlier removal"
    )
    
    ggplotly(g, tooltip = "text") %>%
      layout(legend = list(title = list(text = "Group"), font = list(color = "black")))
  })
  
  output$deg_summary <- renderPrint({
    req(active_res())
    cat("Dataset:", active_res()$dataset_name, "\n")
    cat("Dataset tag:", active_res()$dataset_tag, "\n")
    cat("Organism used for gene symbols:", active_res()$organism_used, "\n")
    cat("Design column:", active_res()$design_col, "\n")
    if (!is.null(active_res()$contrast_desc) && nzchar(active_res()$contrast_desc)) {
      cat("Contrast:", active_res()$contrast_desc, "\n")
    } else {
      cat("Conditions:", active_res()$condition1, "vs", active_res()$condition2, "\n")
    }
    if (.has_real_filter(active_res()$filter_desc)) {
      cat("Filter:", active_res()$filter_desc, "\n")
    }
    cat("Samples used:", active_res()$n_samples_used, "\n")
    cat("Total genes tested:", nrow(active_res()$all_res), "\n")
    cat("Final DEGs (FDR <", FDR_CUTOFF, "):", nrow(active_res()$degs), "\n")
  })
  
  output$pathway_summary <- renderPrint({
    req(input$enable_pathway_analysis, active_res())
    cat("Dataset:", active_res()$dataset_name, "\n")
    cat("Dataset tag:", active_res()$dataset_tag, "\n")
    cat("Input DEG threshold:", active_res()$pathway_input_fdr, "\n")
    cat("Pathway result threshold:", PATHWAY_RESULT_FDR, "\n")
    cat("Input genes (total):", length(active_res()$pathway_input_genes), "\n")
    cat("Libraries used:", if (length(active_res()$pathway_libraries) > 0) paste(active_res()$pathway_libraries, collapse = ", ") else "None", "\n")
    cat("Significant pathway terms:", if (!is.null(active_res()$pathway_results)) nrow(active_res()$pathway_results) else 0, "\n")
  })
  
  output$deg_table <- renderDT({
    req(active_res())
    datatable(active_res()$degs_hgnc, options = list(scrollX = TRUE, pageLength = 5))
  })
  
  output$pathway_table <- renderDT({
    req(input$enable_pathway_analysis)
    dat <- pathway_table_data()
    datatable(dat, options = list(scrollX = TRUE, pageLength = 10), filter = "top")
  })
  
  output$warehouse_table <- renderDT({
    input$runAnalysis
    if (dbExistsTable(con_sql, "master_degs")) {
      df <- dbReadTable(con_sql, "master_degs")
      datatable(df, options = list(scrollX = TRUE, pageLength = 10), filter = "top")
    } else {
      NULL
    }
  })
  
  output$select_meta_A <- renderUI({
    input$runAnalysis
    if (!dbExistsTable(con_sql, "master_degs")) return(NULL)
    choices <- dbGetQuery(con_sql, "SELECT DISTINCT dataset_tag FROM master_degs")[, 1]
    if (length(choices) < 1) return(NULL)
    selectInput("metaA", "Study A", choices = choices)
  })
  
  output$select_meta_B <- renderUI({
    input$runAnalysis
    if (!dbExistsTable(con_sql, "master_degs")) return(NULL)
    choices <- dbGetQuery(con_sql, "SELECT DISTINCT dataset_tag FROM master_degs")[, 1]
    if (length(choices) < 1) return(NULL)
    selectInput("metaB", "Study B", choices = choices)
  })
  
  output$meta_discovery_table <- renderDT({
    req(input$runMetaDiscovery, input$metaA, input$metaB)
    
    cte <- sprintf("
      WITH
      a_raw AS (
        SELECT GeneSymbol, GeneSymbol_Key, log2FoldChange, padj
        FROM master_degs
        WHERE dataset_tag = '%s'
          AND GeneSymbol_Key IS NOT NULL
          AND GeneSymbol_Key <> ''
      ),
      b_raw AS (
        SELECT GeneSymbol, GeneSymbol_Key, log2FoldChange, padj
        FROM master_degs
        WHERE dataset_tag = '%s'
          AND GeneSymbol_Key IS NOT NULL
          AND GeneSymbol_Key <> ''
      ),
      a AS (
        SELECT
          GeneSymbol_Key AS GeneKey,
          MIN(GeneSymbol) AS GeneSymbol_A,
          AVG(log2FoldChange) AS FC_A,
          MIN(padj) AS FDR_A
        FROM a_raw
        GROUP BY GeneSymbol_Key
      ),
      b AS (
        SELECT
          GeneSymbol_Key AS GeneKey,
          MIN(GeneSymbol) AS GeneSymbol_B,
          AVG(log2FoldChange) AS FC_B,
          MIN(padj) AS FDR_B
        FROM b_raw
        GROUP BY GeneSymbol_Key
      )
    ", input$metaA, input$metaB)
    
    final_query <- switch(
      input$meta_question,
      
      "shared" = paste0(cte, "
        SELECT
          a.GeneKey,
          a.GeneSymbol_A,
          b.GeneSymbol_B,
          a.FC_A,
          b.FC_B,
          a.FDR_A,
          b.FDR_B,
          (a.FC_A - b.FC_B) AS DeltaFC,
          ABS(a.FC_A - b.FC_B) AS AbsDeltaFC
        FROM a
        JOIN b ON a.GeneKey = b.GeneKey
        ORDER BY AbsDeltaFC DESC
      "),
      
      "concordant" = paste0(cte, "
        SELECT
          a.GeneKey,
          a.GeneSymbol_A,
          b.GeneSymbol_B,
          a.FC_A,
          b.FC_B,
          a.FDR_A,
          b.FDR_B,
          CASE
            WHEN a.FC_A > 0 AND b.FC_B > 0 THEN 'Up'
            WHEN a.FC_A < 0 AND b.FC_B < 0 THEN 'Down'
            ELSE 'Mixed'
          END AS Direction,
          (a.FC_A - b.FC_B) AS DeltaFC,
          ABS(a.FC_A - b.FC_B) AS AbsDeltaFC
        FROM a
        JOIN b ON a.GeneKey = b.GeneKey
        WHERE (a.FC_A > 0 AND b.FC_B > 0)
           OR (a.FC_A < 0 AND b.FC_B < 0)
        ORDER BY (ABS(a.FC_A) + ABS(b.FC_B)) DESC
      "),
      
      "discordant" = paste0(cte, "
        SELECT
          a.GeneKey,
          a.GeneSymbol_A,
          b.GeneSymbol_B,
          a.FC_A,
          b.FC_B,
          a.FDR_A,
          b.FDR_B,
          (a.FC_A - b.FC_B) AS DeltaFC,
          ABS(a.FC_A - b.FC_B) AS AbsDeltaFC
        FROM a
        JOIN b ON a.GeneKey = b.GeneKey
        WHERE (a.FC_A > 0 AND b.FC_B < 0)
           OR (a.FC_A < 0 AND b.FC_B > 0)
        ORDER BY AbsDeltaFC DESC
      "),
      
      "uniqueA" = paste0(cte, "
        SELECT
          a.GeneKey,
          a.GeneSymbol_A,
          a.FC_A,
          a.FDR_A
        FROM a
        LEFT JOIN b ON a.GeneKey = b.GeneKey
        WHERE b.GeneKey IS NULL
        ORDER BY ABS(a.FC_A) DESC
      "),
      
      "uniqueB" = paste0(cte, "
        SELECT
          b.GeneKey,
          b.GeneSymbol_B,
          b.FC_B,
          b.FDR_B
        FROM b
        LEFT JOIN a ON b.GeneKey = a.GeneKey
        WHERE a.GeneKey IS NULL
        ORDER BY ABS(b.FC_B) DESC
      "),
      
      "shift" = paste0(cte, "
        SELECT
          a.GeneKey,
          a.GeneSymbol_A,
          b.GeneSymbol_B,
          a.FC_A,
          b.FC_B,
          a.FDR_A,
          b.FDR_B,
          (a.FC_A - b.FC_B) AS DeltaFC,
          ABS(a.FC_A - b.FC_B) AS AbsDeltaFC
        FROM a
        JOIN b ON a.GeneKey = b.GeneKey
        WHERE ABS(a.FC_A - b.FC_B) >= 1
        ORDER BY AbsDeltaFC DESC
      ")
    )
    
    res <- dbGetQuery(con_sql, final_query)
    
    datatable(
      res,
      options = list(pageLength = 10, scrollX = TRUE),
      caption = paste("Query Result for:", input$meta_question, " | A:", input$metaA, " | B:", input$metaB)
    )
  }, server = TRUE)
  
  output$heatmap_degs <- renderImage({
    mat <- heatmap_deg_matrix()
    
    if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 2) {
      outfile <- tempfile(fileext = ".png")
      png(outfile, width = 1200, height = 500, res = 150)
      plot.new()
      text(0.5, 0.5, "Not enough DEGs for heat-map.")
      dev.off()
      
      return(list(
        src = outfile,
        contentType = "image/png",
        alt = "No heatmap available"
      ))
    }
    
    n_choice <- if (is.null(input$heatmap_deg_n)) "100" else as.character(input$heatmap_deg_n)
    main_title <- if (n_choice == "all") {
      "Heat-map: all DEGs"
    } else {
      paste0("Heat-map: top ", n_choice, " DEGs")
    }
    
    outfile <- tempfile(fileext = ".png")
    .render_heatmap_png(mat, main_title, outfile, res = 150)
    
    list(
      src = outfile,
      contentType = "image/png",
      alt = main_title
    )
  }, deleteFile = TRUE)
  
  output$download_degs_any <- downloadHandler(
    filename = function() {
      req(active_res())
      ext <- if (is.null(input$deg_format) || input$deg_format == "csv") ".csv" else ".txt"
      paste0(active_res()$dataset_tag, "_DEGs_FDR", FDR_CUTOFF, "_", active_res()$run_stamp, ext)
    },
    content = function(file) {
      req(active_res())
      if (is.null(input$deg_format) || input$deg_format == "csv") {
        write.csv(active_res()$degs_hgnc, file, row.names = FALSE)
      } else {
        write.table(active_res()$degs_hgnc, file, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    }
  )
  
  output$download_norm_any <- downloadHandler(
    filename = function() {
      req(active_res())
      ext <- if (is.null(input$norm_format) || input$norm_format == "csv") ".csv" else ".txt"
      paste0(active_res()$dataset_tag, "_norm_counts_", active_res()$run_stamp, ext)
    },
    content = function(file) {
      req(active_res())
      if (is.null(input$norm_format) || input$norm_format == "csv") {
        write.csv(active_res()$norm_hgnc, file, row.names = FALSE)
      } else {
        write.table(active_res()$norm_hgnc, file, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    }
  )
  
  output$download_full_results <- downloadHandler(
    filename = function() {
      req(active_res())
      paste0(active_res()$dataset_tag, "_FULL_DESeq2_results_FDR", FDR_CUTOFF, "_", active_res()$run_stamp, ".txt")
    },
    content = function(file) {
      req(active_res())
      write.table(active_res()$all_res_hgnc, file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )
  
  output$download_vst <- downloadHandler(
    filename = function() {
      req(active_res())
      paste0(active_res()$dataset_tag, "_VST_matrix_", active_res()$run_stamp, ".txt")
    },
    content = function(file) {
      req(active_res())
      vst_mat <- vst(active_res()$dds_final, blind = TRUE) |> assay()
      gene_ids <- rownames(vst_mat)
      gene_sym <- map_to_symbol(gene_ids, active_res()$organism_used)
      
      vst_df <- data.frame(
        GeneID = gene_ids,
        GeneSymbol = unname(gene_sym),
        vst_mat,
        row.names = NULL,
        check.names = FALSE
      )
      
      write.table(vst_df, file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )
  
  output$download_genes_tested <- downloadHandler(
    filename = function() {
      req(active_res())
      paste0(active_res()$dataset_tag, "_genes_tested_", active_res()$run_stamp, ".txt")
    },
    content = function(file) {
      req(active_res())
      genes <- character()
      if (!is.null(active_res()$all_res_hgnc) && nrow(active_res()$all_res_hgnc) > 0) {
        genes <- as.character(active_res()$all_res_hgnc$GeneID)
      } else if (!is.null(active_res()$all_res) && nrow(active_res()$all_res) > 0) {
        genes <- as.character(rownames(active_res()$all_res))
      }
      genes <- unique(genes[!is.na(genes) & genes != ""])
      writeLines(genes, con = file)
    }
  )
  
  output$download_samples_used <- downloadHandler(
    filename = function() {
      req(active_res())
      paste0(active_res()$dataset_tag, "_samples_used_", active_res()$run_stamp, ".txt")
    },
    content = function(file) {
      req(active_res())
      samples <- as.character(colnames(active_res()$norm))
      samples <- unique(samples[!is.na(samples) & samples != ""])
      writeLines(samples, con = file)
    }
  )
  
  output$download_pathway_results <- downloadHandler(
    filename = function() {
      req(input$enable_pathway_analysis, active_res())
      paste0(active_res()$dataset_tag, "_PathwayResults_inputFDR_", active_res()$pathway_input_fdr, "_", active_res()$run_stamp, ".csv")
    },
    content = function(file) {
      req(input$enable_pathway_analysis)
      write.csv(pathway_table_data(), file, row.names = FALSE)
    }
  )
  
  output$download_pathway_input_genes <- downloadHandler(
    filename = function() {
      req(input$enable_pathway_analysis, active_res())
      paste0(active_res()$dataset_tag, "_PathwayInputGenes_inputFDR_", active_res()$pathway_input_fdr, "_", active_res()$run_stamp, ".csv")
    },
    content = function(file) {
      req(input$enable_pathway_analysis)
      write.csv(pathway_input_gene_table(), file, row.names = FALSE)
    }
  )
  
  output$analysis_summary <- renderPrint({
    req(active_res())
    cat(.build_summary_text(active_res()))
  })
  
  output$download_report <- downloadHandler(
    filename = function() {
      req(active_res())
      paste0(active_res()$dataset_tag, "_R-ODAF_summary_", active_res()$run_stamp, ".txt")
    },
    content = function(file) {
      req(active_res())
      writeLines(.build_summary_text(active_res()), con = file)
    }
  )
  
  output$download_sessioninfo <- downloadHandler(
    filename = function() {
      req(active_res())
      paste0("R-ODAF_sessionInfo_", active_res()$run_stamp, ".txt")
    },
    content = function(file) {
      writeLines(capture.output(utils::sessionInfo()), con = file)
    }
  )
  
  output$download_pca_before <- downloadHandler(
    filename = function() {
      req(active_res())
      paste0(active_res()$dataset_tag, "_PCA_before_", active_res()$run_stamp, ".", input$pca_before_format)
    },
    content = function(file) {
      req(active_res(), active_design_col())
      dc <- active_design_col()
      g <- .pca_plot_gg(
        dds = active_res()$dds_init,
        dc = dc,
        condition1 = active_res()$condition1,
        condition2 = active_res()$condition2,
        title = "PCA before outlier removal"
      )
      ggsave(file, plot = g, dpi = 300, width = 6, height = 5, units = "in")
    }
  )
  
  output$download_pca_after <- downloadHandler(
    filename = function() {
      req(active_res())
      paste0(active_res()$dataset_tag, "_PCA_after_", active_res()$run_stamp, ".", input$pca_after_format)
    },
    content = function(file) {
      req(active_res(), active_design_col())
      dc <- active_design_col()
      g <- .pca_plot_gg(
        dds = active_res()$dds_final,
        dc = dc,
        condition1 = active_res()$condition1,
        condition2 = active_res()$condition2,
        title = "PCA after outlier removal"
      )
      ggsave(file, plot = g, dpi = 300, width = 6, height = 5, units = "in")
    }
  )
  
  output$download_top_expr <- downloadHandler(
    filename = function() {
      req(active_res())
      paste0(active_res()$dataset_tag, "_TopGeneExpression_", active_res()$run_stamp, ".", input$top_expr_format)
    },
    content = function(file) {
      req(topExprDF(), active_design_col(), active_res())
      dc <- active_design_col()
      
      top_colors <- setNames(
        c("#00BFC4", "#F8766D"),
        c(active_res()$condition1, active_res()$condition2)
      )
      
      p <- ggplot(
        topExprDF(),
        aes_string(x = dc, y = "expression", color = dc, text = "sample")
      ) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.6) +
        facet_wrap(
          ~ gene,
          ncol = WINDOW_SIZE,
          scales = "free_y",
          labeller = labeller(gene = function(x) str_wrap(x, width = 10))
        ) +
        scale_color_manual(values = top_colors) +
        scale_x_discrete(labels = NULL) +
        labs(x = dc, y = "Normalized count") +
        theme_minimal() +
        theme(
          strip.text = element_text(size = 8),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing.x = unit(1.1, "lines")
        )      
      ggsave(file, plot = p, dpi = 300, width = 7, height = 5, units = "in")
    }
  )
  
  output$download_avg_expr <- downloadHandler(
    filename = function() {
      req(active_res())
      paste0(active_res()$dataset_tag, "_AvgExpression_", active_res()$run_stamp, ".", input$avg_expr_format)
    },
    content = function(file) {
      req(avgExprDF(), active_res())
      
      avg_dat <- avgExprDF()
      avg_dat$condition <- factor(
        avg_dat$condition,
        levels = c(active_res()$condition1, active_res()$condition2)
      )
      
      avg_colors <- setNames(
        c("#00BFC4", "#F8766D"),
        c(active_res()$condition1, active_res()$condition2)
      )
      
      p <- ggplot(avg_dat, aes_string(x = "gene", y = "mean_expression", fill = "condition")) +
        geom_col(position = position_dodge()) +
        scale_fill_manual(values = avg_colors) +
        labs(x = "Gene", y = "Mean normalized count") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          axis.title.x = element_text(margin = margin(t = 10))
        )
      
      ggsave(file, plot = p, dpi = 300, width = 7, height = 5, units = "in")
    }
    
  )
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      req(active_res())
      paste0(active_res()$dataset_tag, "_VolcanoPlot_", active_res()$run_stamp, ".", input$volcano_format)
    },
    content = function(file) {
      req(volcano_df(), active_res())
      dat <- volcano_df()
      p <- .volcano_plot_gg(dat = dat)
      ggsave(file, plot = p, dpi = 300, width = 7, height = 5, units = "in")
    }
  )
  
  output$download_heat_degs <- downloadHandler(
    filename = function() {
      req(active_res())
      n_choice <- if (is.null(input$heatmap_deg_n)) "100" else as.character(input$heatmap_deg_n)
      paste0(
        active_res()$dataset_tag,
        "_Heatmap_",
        if (n_choice == "all") "allDEGs" else paste0("top", n_choice, "DEGs"),
        "_",
        active_res()$run_stamp,
        ".",
        input$heat_degs_format
      )
    },
    content = function(file) {
      mat <- heatmap_deg_matrix()
      n_choice <- if (is.null(input$heatmap_deg_n)) "100" else as.character(input$heatmap_deg_n)
      main_title <- if (n_choice == "all") {
        "Heat-map: all DEGs"
      } else {
        paste0("Heat-map: top ", n_choice, " DEGs")
      }
      
      if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 2) {
        png(file, width = 1200, height = 500, res = 150)
        plot.new()
        text(0.5, 0.5, "Not enough DEGs for heat-map.")
        dev.off()
        return(invisible(NULL))
      }
      
      ext <- input$heat_degs_format
      
      if (ext == "png") {
        .render_heatmap_png(mat, main_title, file, res = 300)
      } else {
        if (ext == "svg") {
          svg(file, width = 10, height = 8)
        } else if (ext == "tiff") {
          tiff(file, width = 10, height = 8, units = "in", res = 300, compression = "lzw")
        } else if (ext == "jpeg") {
          jpeg(file, width = 10, height = 8, units = "in", res = 300, quality = 100)
        } else {
          png(file, width = 10, height = 8, units = "in", res = 300)
        }
        
        on.exit(dev.off(), add = TRUE)
        
        pheatmap(
          mat,
          scale = "row",
          show_rownames = FALSE,
          main = main_title,
          border_color = NA,
          fontsize_col = 10
        )
      }
    }
  )
}

shinyApp(ui, server)

