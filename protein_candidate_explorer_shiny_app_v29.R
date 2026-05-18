
# Protein Candidate Explorer Shiny App (v29)
# Modes:
# 1) Gene-list mode: existing annotation/enrichment workflow for pre-filtered gene/protein lists
# 2) Raw-data mode: raw abundance table -> QC/cleaning -> differential analysis -> volcano -> enrichment -> ion annotation

required_pkgs <- c(
  "shiny", "DT", "readxl", "readr", "dplyr", "stringr", "ggplot2",
  "AnnotationDbi", "org.Hs.eg.db", "GO.db", "clusterProfiler", "enrichplot", "openxlsx", "devEMF"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkgs, collapse = ", "),
    "\nPlease install them first. Example:\n",
    "install.packages(c('shiny','DT','readxl','readr','dplyr','stringr','ggplot2','openxlsx','devEMF'))\n",
    "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
    "BiocManager::install(c('AnnotationDbi','org.Hs.eg.db','GO.db','clusterProfiler','enrichplot'))"
  )
}

library(shiny)
library(DT)
library(readxl)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(clusterProfiler)
library(enrichplot)
library(openxlsx)
library(devEMF)

read_any_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    df <- readxl::read_xlsx(path)
  } else if (ext %in% c("csv")) {
    df <- readr::read_csv(path, show_col_types = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    df <- readr::read_tsv(path, show_col_types = FALSE)
  } else {
    stop("Unsupported file type: ", path)
  }
  colnames(df) <- make.unique(colnames(df))
  as.data.frame(df)
}

extract_genename <- function(df, desc_col = NULL) {
  colnames(df) <- make.unique(colnames(df))
  if (is.null(desc_col) || identical(desc_col, "")) {
    cand <- c("Description", "description", "Protein Description", "protein_description")
    hit <- cand[cand %in% colnames(df)]
    if (length(hit) > 0) {
      desc_col <- hit[1]
    } else {
      idx <- which(vapply(df, function(x) any(grepl("GN=", as.character(x), fixed = TRUE), na.rm = TRUE), logical(1)))
      if (length(idx) == 0) stop("Cannot find a column containing 'GN=' (Description-like column).")
      desc_col <- colnames(df)[idx[1]]
    }
  }
  if (!desc_col %in% colnames(df)) stop("Selected Description column not found: ", desc_col)

  df$GeneName <- stringr::str_extract(as.character(df[[desc_col]]), "GN=[^ ;]+")
  df$GeneName <- stringr::str_remove(df$GeneName, "^GN=")
  df
}

normalize_ion_name <- function(x) {
  x <- tolower(trimws(x))
  dplyr::case_when(
    x %in% c("calcium", "ca", "钙") ~ "calcium",
    x %in% c("magnesium", "mg", "镁") ~ "magnesium",
    x %in% c("manganese", "mn", "锰") ~ "manganese",
    x %in% c("zinc", "zn", "锌") ~ "zinc",
    TRUE ~ NA_character_
  )
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(as.character(x[!is.na(x) & x != ""]))
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = sep)
}

make_pattern <- function(x) paste0("(", paste(unique(x), collapse = "|"), ")")

empty_enrich_df <- function() {
  data.frame(
    ID = character(), Description = character(), GeneRatio = character(), BgRatio = character(),
    pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(), geneID = character(), Count = integer(),
    stringsAsFactors = FALSE
  )
}

empty_export_df <- function(message = "No results") {
  data.frame(Message = message, stringsAsFactors = FALSE)
}

add_direction_column <- function(df, direction) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  out <- as.data.frame(df)
  out$Direction <- direction
  out
}

combine_up_down_df <- function(up_df, down_df) {
  parts <- Filter(Negate(is.null), list(
    add_direction_column(up_df, "Up"),
    add_direction_column(down_df, "Down")
  ))
  if (length(parts) == 0) return(empty_export_df())
  dplyr::bind_rows(parts)
}

summarise_ion_go_terms <- function(go_hits_df) {
  if (is.null(go_hits_df) || nrow(go_hits_df) == 0) {
    return(data.frame(
      GO = character(),
      GO_TERM = character(),
      GO_ONTO = character(),
      protein_n = integer(),
      proteins = character(),
      ion = character(),
      stringsAsFactors = FALSE
    ))
  }
  req_cols <- c("GO", "GO_TERM", "GO_ONTOLOGY", "GeneName", "ion")
  missing_cols <- setdiff(req_cols, colnames(go_hits_df))
  if (length(missing_cols) > 0) {
    return(data.frame(
      GO = character(),
      GO_TERM = character(),
      GO_ONTO = character(),
      protein_n = integer(),
      proteins = character(),
      ion = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- go_hits_df %>%
    dplyr::group_by(GO, GO_TERM, GO_ONTOLOGY, ion) %>%
    dplyr::summarise(
      protein_n = dplyr::n_distinct(GeneName),
      proteins = collapse_unique(GeneName),
      .groups = "drop"
    ) %>%
    dplyr::rename(GO_ONTO = GO_ONTOLOGY) %>%
    dplyr::arrange(ion, dplyr::desc(protein_n), GO_TERM)
  as.data.frame(out)
}


write_named_table <- function(wb, sheet, title, df, start_row = 1) {
  openxlsx::writeData(wb, sheet, x = data.frame(Section = title, stringsAsFactors = FALSE), startRow = start_row, colNames = FALSE)
  out_df <- if (is.null(df) || nrow(df) == 0) empty_export_df() else as.data.frame(df)
  openxlsx::writeData(wb, sheet, x = out_df, startRow = start_row + 1, withFilter = TRUE)
  invisible(start_row + nrow(out_df) + 3)
}

safe_parse_numeric <- function(x) {
  if (is.numeric(x)) {
    y <- as.numeric(x)
  } else {
    ch <- trimws(as.character(x))
    y <- suppressWarnings(as.numeric(ch))
    need_fallback <- is.na(y) & !is.na(ch) & ch != ""
    if (any(need_fallback)) {
      y[need_fallback] <- suppressWarnings(
        readr::parse_number(ch[need_fallback], na = c("", "NA", "NaN", "NULL", "null", "N/A"))
      )
    }
  }
  y[is.infinite(y)] <- NA_real_
  # Keep raw zeros in the parsed data. Mode A may use 0 in mean abundance / FC,
  # while Modes B/C still treat 0 as undetected in counts, cleaning, and tests.
  y
}


detect_confidence_column <- function(df) {
  nms <- colnames(df)
  nms_lc <- tolower(nms)
  exact_hits <- nms[grepl("protein fdr confidence", nms_lc) & grepl("combined", nms_lc)]
  if (length(exact_hits) > 0) return(exact_hits[1])
  broad_hits <- nms[grepl("confidence", nms_lc)]
  if (length(broad_hits) > 0) return(broad_hits[1])
  NULL
}

apply_confidence_filter <- function(df, confidence_col = NULL, confidence_filter = c("high", "high_medium", "all")) {
  confidence_filter <- match.arg(confidence_filter)
  if (is.null(confidence_col) || !confidence_col %in% colnames(df) || identical(confidence_filter, "all")) {
    return(list(df = df, confidence_col = confidence_col, filter_applied = FALSE, levels_kept = "All"))
  }
  vals <- trimws(as.character(df[[confidence_col]]))
  vals_uc <- toupper(vals)
  keep_levels <- switch(
    confidence_filter,
    high = c("HIGH"),
    high_medium = c("HIGH", "MEDIUM"),
    all = unique(vals_uc)
  )
  out <- df[!is.na(vals_uc) & vals_uc %in% keep_levels, , drop = FALSE]
  list(
    df = out,
    confidence_col = confidence_col,
    filter_applied = TRUE,
    levels_kept = paste(keep_levels, collapse = " + ")
  )
}

detect_description_column <- function(df, preferred = NULL) {
  if (!is.null(preferred) && nzchar(preferred) && preferred %in% colnames(df)) return(preferred)
  cand <- c("Description", "description", "Protein Description", "protein_description")
  hit <- cand[cand %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  idx <- which(vapply(df, function(x) any(grepl("GN=", as.character(x), fixed = TRUE), na.rm = TRUE), logical(1)))
  if (length(idx) > 0) return(colnames(df)[idx[1]])
  NULL
}

get_contaminant_rules <- function(contaminant_filter = c("common", "keratin", "none")) {
  contaminant_filter <- match.arg(contaminant_filter)
  if (identical(contaminant_filter, "none")) return(character())

  keratin_rules <- c(
    "keratin_gene" = "^KRT(?:[0-9A-Z].*)?$",
    "krtap_gene" = "^KRTAP(?:[0-9A-Z].*)?$",
    "keratin_desc" = "\\bkeratin\\b"
  )

  common_extra <- c(
    "trypsin_or_lysc" = "\\btrypsin\\b|trypsinogen|lys-?c",
    "albumin_bsa_casein" = "^ALB$|bovine serum albumin|\\bbsa\\b|\\bcasein\\b|serum albumin",
    "hemoglobin_gene" = "^HBA[0-9A-Z]*$|^HBB[0-9A-Z]*$|^HBD[0-9A-Z]*$|^HBE[0-9A-Z]*$|^HBG[0-9A-Z]*$|^HBM[0-9A-Z]*$|^HBQ[0-9A-Z]*$|^HBZ[0-9A-Z]*$",
    "hemoglobin_desc" = "hemoglobin",
    "immunoglobulin_gene" = "^IGH|^IGK|^IGL",
    "immunoglobulin_desc" = "immunoglobulin"
  )

  if (identical(contaminant_filter, "keratin")) return(keratin_rules)
  c(keratin_rules, common_extra)
}

apply_contaminant_filter <- function(df, gene_col = "GeneName", desc_col = NULL, contaminant_filter = c("common", "keratin", "none")) {
  contaminant_filter <- match.arg(contaminant_filter)
  if (identical(contaminant_filter, "none")) {
    return(list(
      df = df,
      removed_df = data.frame(),
      filter_applied = FALSE,
      filter_label = "None"
    ))
  }

  rules <- get_contaminant_rules(contaminant_filter)
  if (length(rules) == 0) {
    return(list(
      df = df,
      removed_df = data.frame(),
      filter_applied = FALSE,
      filter_label = "None"
    ))
  }

  gene_txt <- if (gene_col %in% colnames(df)) as.character(df[[gene_col]]) else rep("", nrow(df))
  desc_txt <- if (!is.null(desc_col) && desc_col %in% colnames(df)) as.character(df[[desc_col]]) else rep("", nrow(df))

  matched_rule <- rep(NA_character_, nrow(df))
  for (nm in names(rules)) {
    hit <- is.na(matched_rule) & (
      grepl(rules[[nm]], gene_txt, ignore.case = TRUE, perl = TRUE) |
      grepl(rules[[nm]], desc_txt, ignore.case = TRUE, perl = TRUE)
    )
    matched_rule[hit] <- nm
  }

  removed_df <- df[!is.na(matched_rule), , drop = FALSE]
  if (nrow(removed_df) > 0) {
    removed_df$matched_contaminant_rule <- matched_rule[!is.na(matched_rule)]
  }
  kept_df <- df[is.na(matched_rule), , drop = FALSE]

  list(
    df = kept_df,
    removed_df = removed_df,
    filter_applied = TRUE,
    filter_label = if (identical(contaminant_filter, "common")) "Common contaminants (293T default)" else "Keratin only"
  )
}

detect_group_columns <- function(df) {
  test_cols <- colnames(df)[grepl("^test\\d+$", colnames(df), ignore.case = TRUE)]
  ctrl_cols <- colnames(df)[grepl("^ctrl\\d+$", colnames(df), ignore.case = TRUE)]
  list(test_cols = test_cols, ctrl_cols = ctrl_cols)
}

safe_bitr <- function(genes, to_types = c("ENTREZID")) {
  genes <- unique(as.character(genes[!is.na(genes) & genes != ""]))
  if (length(genes) == 0) return(data.frame())
  out <- tryCatch(
    suppressMessages(
      clusterProfiler::bitr(
        genes,
        fromType = "SYMBOL",
        toType = to_types,
        OrgDb = org.Hs.eg.db
      )
    ),
    error = function(e) data.frame()
  )
  if (nrow(out) == 0) return(data.frame())
  colnames(out) <- make.unique(colnames(out))
  as.data.frame(out)
}

safe_dotplot <- function(enrich_obj, title = "Enrichment", show_n = 12, wrap_width = 35, y_text_size = 9) {
  if (is.null(enrich_obj)) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = paste0(title, "\nNo enriched terms"), size = 6) +
        theme_void()
    )
  }
  enrichplot::dotplot(enrich_obj, showCategory = show_n) +
    ggtitle(title) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = wrap_width)) +
    theme(axis.text.y = element_text(size = y_text_size))
}

safe_enrich_go <- function(entrez, universe = NULL, ont = "BP", p_cut = 0.05, q_cut = 0.2) {
  entrez <- unique(as.character(entrez[!is.na(entrez) & entrez != ""]))
  universe <- unique(as.character(universe[!is.na(universe) & universe != ""]))
  if (length(entrez) == 0) return(list(obj = NULL, df = empty_enrich_df(), error = NULL))

  args <- list(
    gene = entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = p_cut,
    qvalueCutoff = q_cut,
    readable = TRUE
  )
  if (length(universe) > 0) args$universe <- universe

  out <- tryCatch(do.call(enrichGO, args), error = function(e) e)
  if (inherits(out, "error")) return(list(obj = NULL, df = empty_enrich_df(), error = conditionMessage(out)))
  df <- as.data.frame(out)
  if (nrow(df) == 0) return(list(obj = NULL, df = empty_enrich_df(), error = NULL))
  list(obj = out, df = df, error = NULL)
}

safe_enrich_kegg <- function(entrez, universe = NULL, p_cut = 0.05) {
  entrez <- unique(as.character(entrez[!is.na(entrez) & entrez != ""]))
  universe <- unique(as.character(universe[!is.na(universe) & universe != ""]))
  if (length(entrez) == 0) return(list(obj = NULL, df = empty_enrich_df(), error = NULL))

  args <- list(
    gene = entrez,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = p_cut
  )
  if (length(universe) > 0) args$universe <- universe

  out <- tryCatch(do.call(enrichKEGG, args), error = function(e) e)
  if (inherits(out, "error")) return(list(obj = NULL, df = empty_enrich_df(), error = conditionMessage(out)))
  df <- as.data.frame(out)
  if (nrow(df) == 0) return(list(obj = NULL, df = empty_enrich_df(), error = NULL))
  list(obj = out, df = df, error = NULL)
}


save_plot_emf <- function(plot_obj, file, width = 8, height = 6, bg = "white") {
  devEMF::emf(file = file, width = width, height = height, bg = bg)
  print(plot_obj)
  grDevices::dev.off()
}

save_plot_pdf <- function(plot_obj, file, width = 8, height = 6, bg = "white") {
  grDevices::pdf(file = file, width = width, height = height, bg = bg, useDingbats = FALSE)
  print(plot_obj)
  grDevices::dev.off()
}

save_plot_file <- function(plot_obj, file, format = c("png", "jpg", "emf", "pdf"), width = 8, height = 6, dpi = 300, bg = "white") {
  format <- match.arg(format)
  if (format == "emf") {
    save_plot_emf(plot_obj, file, width = width, height = height, bg = bg)
  } else if (format == "pdf") {
    save_plot_pdf(plot_obj, file, width = width, height = height, bg = bg)
  } else {
    ggplot2::ggsave(filename = file, plot = plot_obj, width = width, height = height, dpi = dpi, bg = bg, device = format)
  }
}

run_classic_enrichment <- function(symbol_genes, universe_symbols = NULL, p_cut = 0.05, q_cut = 0.2, show_n = 12) {
  gene_map <- safe_bitr(symbol_genes, to_types = c("ENTREZID"))
  universe_map <- safe_bitr(universe_symbols, to_types = c("ENTREZID"))

  entrez <- unique(as.character(gene_map$ENTREZID))
  universe_entrez <- unique(as.character(universe_map$ENTREZID))

  bp_res <- safe_enrich_go(entrez, universe = universe_entrez, ont = "BP", p_cut = p_cut, q_cut = q_cut)
  mf_res <- safe_enrich_go(entrez, universe = universe_entrez, ont = "MF", p_cut = p_cut, q_cut = q_cut)
  cc_res <- safe_enrich_go(entrez, universe = universe_entrez, ont = "CC", p_cut = p_cut, q_cut = q_cut)
  kegg_res <- safe_enrich_kegg(entrez, universe = universe_entrez, p_cut = p_cut)

  list(
    gene_map = gene_map,
    universe_map = universe_map,
    entrez_count = length(entrez),
    universe_entrez_count = length(universe_entrez),
    bp = bp_res,
    mf = mf_res,
    cc = cc_res,
    kegg = kegg_res,
    show_n = show_n,
    p_cut = p_cut,
    q_cut = q_cut
  )
}

get_go_info_table <- function() {
  all_go_info <- AnnotationDbi::toTable(GO.db::GOTERM)
  colnames(all_go_info) <- make.unique(colnames(all_go_info))
  name_map <- tolower(colnames(all_go_info))
  col_go <- colnames(all_go_info)[match("go_id", name_map)]
  col_term <- colnames(all_go_info)[match("term", name_map)]
  col_ont <- colnames(all_go_info)[match("ontology", name_map)]
  if (any(is.na(c(col_go, col_term, col_ont)))) stop("Cannot find GO term columns in GO.db::GOTERM table.")

  go_info <- data.frame(
    GO = as.character(all_go_info[[col_go]]),
    GO_TERM = as.character(all_go_info[[col_term]]),
    GO_ONTOLOGY = as.character(all_go_info[[col_ont]]),
    stringsAsFactors = FALSE
  )
  go_info <- go_info[!is.na(go_info$GO_TERM) & !is.na(go_info$GO_ONTOLOGY), , drop = FALSE]
  unique(go_info)
}

get_ion_patterns <- function(ion) {
  ion <- normalize_ion_name(ion)
  if (is.na(ion)) stop("Unsupported ion.")

  if (ion == "calcium") {
    strict_patterns <- c(
      "response to calcium ion",
      "cellular response to calcium ion",
      "calcium ion transport",
      "calcium ion transmembrane transport",
      "calcium ion homeostasis",
      "regulation of cytosolic calcium ion concentration",
      "cytosolic calcium ion concentration",
      "calcium-mediated signaling",
      "calcium-dependent cell-cell adhesion"
    )
    broad_patterns <- c(
      strict_patterns,
      "calcium ion binding",
      "calcium binding",
      "calmodulin",
      "calcium-dependent",
      "voltage-gated calcium",
      "store-operated calcium",
      "sarcoplasmic reticulum calcium",
      "endoplasmic reticulum calcium",
      "calcineurin"
    )
  } else if (ion == "magnesium") {
    strict_patterns <- c(
      "response to magnesium ion",
      "cellular response to magnesium ion",
      "magnesium ion transport",
      "magnesium ion transmembrane transport",
      "magnesium ion homeostasis",
      "regulation of magnesium ion transport",
      "regulation of magnesium ion homeostasis"
    )
    broad_patterns <- c(strict_patterns, "magnesium ion binding", "magnesium binding", "magnesium-dependent")
  } else if (ion == "manganese") {
    strict_patterns <- c(
      "response to manganese ion",
      "cellular response to manganese ion",
      "manganese ion transport",
      "manganese ion transmembrane transport",
      "manganese ion homeostasis",
      "regulation of manganese ion transport",
      "regulation of manganese ion homeostasis"
    )
    broad_patterns <- c(strict_patterns, "manganese ion binding", "manganese binding", "manganese-dependent")
  } else if (ion == "zinc") {
    strict_patterns <- c(
      "response to zinc ion",
      "cellular response to zinc ion",
      "zinc ion transport",
      "zinc ion transmembrane transport",
      "zinc ion homeostasis",
      "regulation of zinc ion transport",
      "regulation of zinc ion homeostasis"
    )
    broad_patterns <- c(strict_patterns, "zinc ion binding", "zinc binding", "zinc-dependent")
  }

  list(
    ion = ion,
    strict_patterns = unique(strict_patterns),
    broad_patterns = unique(broad_patterns)
  )
}

annotate_genes_dataframe <- function(df_in, gene_col = "GeneName", target_ions = c("calcium", "magnesium", "manganese", "zinc")) {
  if (!gene_col %in% colnames(df_in)) stop("Gene column not found: ", gene_col)
  normalized_ions <- unique(stats::na.omit(vapply(target_ions, normalize_ion_name, character(1))))
  genes <- unique(as.character(df_in[[gene_col]][!is.na(df_in[[gene_col]]) & df_in[[gene_col]] != ""]))

  if (length(genes) == 0 || length(normalized_ions) == 0) {
    return(list(
      annotated_df = df_in,
      protein_summary = data.frame(),
      go_hits = data.frame(),
      ion_term_summary = data.frame(),
      overview = data.frame(),
      gene_map_go = data.frame(),
      gene_map_entrez = data.frame(),
      selected_ions = normalized_ions
    ))
  }

  gene_map_go <- safe_bitr(genes, to_types = c("ENTREZID", "GO"))
  gene_map_entrez <- safe_bitr(genes, to_types = c("ENTREZID"))
  if (nrow(gene_map_go) == 0) {
    annotated_df <- df_in
    for (ion in normalized_ions) {
      annotated_df[[paste0(ion, "_strict")]] <- FALSE
      annotated_df[[paste0(ion, "_broad")]] <- FALSE
    }
    return(list(
      annotated_df = annotated_df,
      protein_summary = data.frame(),
      go_hits = data.frame(),
      ion_term_summary = data.frame(),
      overview = data.frame(
        ion = normalized_ions,
        strict_protein_n = 0,
        broad_protein_n = 0,
        go_hit_n = 0,
        stringsAsFactors = FALSE
      ),
      gene_map_go = gene_map_go,
      gene_map_entrez = gene_map_entrez,
      selected_ions = normalized_ions
    ))
  }

  go_info <- get_go_info_table()

  gene2go <- gene_map_go[, c("SYMBOL", "ENTREZID", "GO"), drop = FALSE]
  colnames(gene2go) <- c("GeneName", "ENTREZID", "GO")
  gene2go$GO <- as.character(gene2go$GO)

  merged_go <- merge(gene2go, go_info, by = "GO", all.x = TRUE)
  colnames(merged_go) <- make.unique(colnames(merged_go))
  merged_go$GO_TERM_lc <- tolower(merged_go$GO_TERM)

  annotated_df <- df_in
  combined_protein_summary <- list()
  combined_hit_table <- list()
  overview <- list()

  for (ion in normalized_ions) {
    pat <- get_ion_patterns(ion)
    strict_regex <- make_pattern(pat$strict_patterns)
    broad_regex <- make_pattern(pat$broad_patterns)

    ion_go <- merged_go
    ion_go$strict_response <- grepl(strict_regex, ion_go$GO_TERM_lc)
    ion_go$broad_related <- grepl(broad_regex, ion_go$GO_TERM_lc)

    ion_hits_all <- ion_go[ion_go$broad_related %in% TRUE, , drop = FALSE]

    if (nrow(ion_hits_all) == 0) {
      protein_summary <- data.frame(
        GeneName = character(),
        hit_type = character(),
        matched_GO_terms = character(),
        matched_GO_ids = character(),
        matched_ontology = character(),
        stringsAsFactors = FALSE
      )
    } else {
      ion_hits_all$hit_type <- ifelse(ion_hits_all$strict_response, "strict", "broad")

      protein_summary <- ion_hits_all %>%
        group_by(GeneName) %>%
        summarise(
          hit_type = ifelse(any(strict_response), "strict", "broad"),
          matched_GO_terms = collapse_unique(GO_TERM),
          matched_GO_ids = collapse_unique(GO),
          matched_ontology = collapse_unique(GO_ONTOLOGY),
          .groups = "drop"
        ) %>%
        arrange(hit_type, GeneName)
    }

    strict_genes <- unique(ion_go$GeneName[ion_go$strict_response %in% TRUE])
    broad_genes <- unique(ion_go$GeneName[ion_go$broad_related %in% TRUE])

    annotated_df[[paste0(ion, "_strict")]] <- annotated_df[[gene_col]] %in% strict_genes
    annotated_df[[paste0(ion, "_broad")]] <- annotated_df[[gene_col]] %in% broad_genes

    ion_hits_all$ion <- rep(ion, nrow(ion_hits_all))
    protein_summary$ion <- rep(ion, nrow(protein_summary))

    combined_hit_table[[ion]] <- ion_hits_all
    combined_protein_summary[[ion]] <- protein_summary
    overview[[ion]] <- data.frame(
      ion = ion,
      strict_protein_n = length(strict_genes),
      broad_protein_n = length(broad_genes),
      go_hit_n = nrow(ion_hits_all),
      stringsAsFactors = FALSE
    )
  }

  combined_go_hits <- dplyr::bind_rows(combined_hit_table)
  list(
    annotated_df = annotated_df,
    protein_summary = dplyr::bind_rows(combined_protein_summary),
    go_hits = combined_go_hits,
    ion_term_summary = summarise_ion_go_terms(combined_go_hits),
    overview = dplyr::bind_rows(overview),
    gene_map_go = gene_map_go,
    gene_map_entrez = gene_map_entrez,
    selected_ions = normalized_ions
  )
}

resolve_duplicate_genes <- function(df, gene_col = "GeneName", abundance_cols = character()) {
  df$.row_id_original <- seq_len(nrow(df))
  if (!gene_col %in% colnames(df)) stop("Gene column not found: ", gene_col)

  if (length(abundance_cols) == 0) {
    dedup_df <- df %>% filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "")
    dup_summary <- data.frame()
    return(list(df = dedup_df, duplicate_summary = dup_summary))
  }

  score <- apply(df[, abundance_cols, drop = FALSE], 1, function(x) {
    x <- as.numeric(x)
    if (all(is.na(x))) return(-Inf)
    mean(x, na.rm = TRUE)
  })
  df$.duplicate_score <- score

  duplicated_genes <- names(which(table(df[[gene_col]]) > 1))
  dup_summary <- data.frame()

  if (length(duplicated_genes) > 0) {
    dup_df <- df[df[[gene_col]] %in% duplicated_genes & !is.na(df[[gene_col]]) & df[[gene_col]] != "", , drop = FALSE]
    dup_summary <- dup_df %>%
      group_by(.data[[gene_col]]) %>%
      summarise(
        duplicate_n = dplyr::n(),
        kept_row_id = .row_id_original[which.max(.duplicate_score)],
        kept_score = max(.duplicate_score, na.rm = TRUE),
        candidate_row_ids = paste(.row_id_original, collapse = "; "),
        candidate_scores = paste(round(.duplicate_score, 4), collapse = "; "),
        .groups = "drop"
      )
    colnames(dup_summary)[1] <- gene_col
  }

  dedup_df <- df %>%
    filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") %>%
    group_by(.data[[gene_col]]) %>%
    slice_max(order_by = .duplicate_score, n = 1, with_ties = FALSE) %>%
    ungroup()

  dedup_df$.duplicate_score <- NULL
  list(df = as.data.frame(dedup_df), duplicate_summary = as.data.frame(dup_summary))
}

clean_group_values <- function(vals, sample_names, outlier_dev_log2 = 1, remaining_range_log2 = 0.5) {
  vals <- as.numeric(vals)
  names(vals) <- sample_names
  removed_samples <- character()
  removal_reason <- NA_character_

  valid_idx <- which(!is.na(vals) & vals > 0)
  if (length(valid_idx) >= 3) {
    log_vals <- log2(vals[valid_idx])
    med_val <- stats::median(log_vals)

    candidate_flag <- logical(length(valid_idx))
    candidate_dev <- rep(NA_real_, length(valid_idx))

    for (i in seq_along(valid_idx)) {
      others <- log_vals[-i]
      if (length(others) >= 2) {
        dev_i <- abs(log_vals[i] - med_val)
        remaining_range <- diff(range(others))
        if (is.finite(dev_i) && is.finite(remaining_range) &&
            dev_i > outlier_dev_log2 && remaining_range <= remaining_range_log2) {
          candidate_flag[i] <- TRUE
          candidate_dev[i] <- dev_i
        }
      }
    }

    if (sum(candidate_flag, na.rm = TRUE) == 1) {
      drop_pos <- valid_idx[which(candidate_flag)]
      removed_samples <- names(vals)[drop_pos]
      vals[drop_pos] <- NA_real_
      removal_reason <- "single_outlier_removed"
    }
  }

  list(values = vals, removed_samples = removed_samples, removal_reason = removal_reason)
}

collapse_named_values <- function(vals) {
  if (length(vals) == 0) return(NA_character_)
  paste(paste0(names(vals), "=", ifelse(is.na(vals), "NA", signif(vals, 6))), collapse = "; ")
}

calc_group_ratio <- function(vals) {
  vals <- vals[!is.na(vals) & vals > 0]
  if (length(vals) < 2) return(NA_real_)
  max(vals) / min(vals)
}




compare_single_gene <- function(row_df, test_cols, ctrl_cols,
                                total_test_reps = length(test_cols),
                                total_ctrl_reps = length(ctrl_cols),
                                outlier_dev_log2 = 1,
                                remaining_range_log2 = 0.5,
                                test_pair_fc = 5,
                                test_inconsistent_fc = 5,
                                fc_cut = 2,
                                p_cut = 0.05,
                                modes_to_compute = c("A", "B", "C")) {
  test_raw <- as.numeric(row_df[1, test_cols, drop = TRUE])
  ctrl_raw <- as.numeric(row_df[1, ctrl_cols, drop = TRUE])
  names(test_raw) <- test_cols
  names(ctrl_raw) <- ctrl_cols
  modes_to_compute <- unique(as.character(modes_to_compute))
  compute_A <- "A" %in% modes_to_compute
  compute_B <- any(c("B", "C") %in% modes_to_compute)
  compute_C <- "C" %in% modes_to_compute

  # Mode A (lenient): missing values are explicitly imputed to 0 before downstream handling.
  test_raw_A <- test_raw
  ctrl_raw_A <- ctrl_raw
  test_raw_A[is.na(test_raw_A)] <- 0
  ctrl_raw_A[is.na(ctrl_raw_A)] <- 0

  n_test_raw <- sum(!is.na(test_raw_A) & test_raw_A > 0)
  n_ctrl_raw <- sum(!is.na(ctrl_raw_A) & ctrl_raw_A > 0)

  empty_metrics <- function() {
    list(
      n_test_clean = NA_real_, n_ctrl_clean = NA_real_,
      n_test_signal = NA_real_, n_ctrl_signal = NA_real_,
      n_test_nonmissing = NA_real_, n_ctrl_nonmissing = NA_real_,
      test_present = FALSE, ctrl_present = FALSE,
      mean_test_abundance = NA_real_, mean_ctrl_abundance = NA_real_,
      mean_test_log2 = NA_real_, mean_ctrl_log2 = NA_real_,
      FC = NA_real_, log2FC = NA_real_,
      pvalue = NA_real_, neglog10p = NA_real_,
      test_ratio_clean = NA_real_, ctrl_ratio_clean = NA_real_,
      replicate_flag = NA_character_, ttest_note = NA_character_
    )
  }

  make_no_clean <- function(vals) list(values = vals, removed_samples = character(0))
  calc_metrics_mode_a <- function(test_vals, ctrl_vals) {
    # Mode A is intentionally independent from Modes B/C:
    # - keep raw 0 values
    # - Test: 0 does NOT count as a valid replicate
    # - Ctrl: 0 DOES count as a valid observation and must enter the mean
    n_test_signal <- sum(!is.na(test_vals) & test_vals > 0)
    n_ctrl_signal <- sum(!is.na(ctrl_vals) & ctrl_vals > 0)
    n_test_nonmissing <- sum(!is.na(test_vals))
    n_ctrl_nonmissing <- sum(!is.na(ctrl_vals))

    n_test_clean <- n_test_signal
    n_ctrl_clean <- n_ctrl_nonmissing

    test_present <- n_test_signal > 0
    ctrl_present <- n_ctrl_signal > 0

    mean_test_abundance <- if (n_test_signal > 0 && test_present) {
      mean(test_vals[!is.na(test_vals) & test_vals > 0])
    } else {
      NA_real_
    }
    mean_ctrl_abundance <- if (n_ctrl_nonmissing > 0) {
      mean(ctrl_vals[!is.na(ctrl_vals)])
    } else {
      NA_real_
    }

    FC <- if (is.finite(mean_test_abundance) && is.finite(mean_ctrl_abundance) &&
              mean_test_abundance > 0 && mean_ctrl_abundance > 0) {
      mean_test_abundance / mean_ctrl_abundance
    } else {
      NA_real_
    }

    log2FC <- if (is.finite(FC) && FC > 0) log2(FC) else NA_real_
    mean_test_log2 <- if (is.finite(mean_test_abundance) && mean_test_abundance > 0) log2(mean_test_abundance) else NA_real_
    mean_ctrl_log2 <- if (is.finite(mean_ctrl_abundance) && mean_ctrl_abundance > 0) log2(mean_ctrl_abundance) else NA_real_

    pvalue <- NA_real_
    ttest_note <- NA_character_
    if (n_test_signal >= 2 && n_ctrl_signal >= 2) {
      tt <- tryCatch(
        t.test(
          x = log2(test_vals[!is.na(test_vals) & test_vals > 0]),
          y = log2(ctrl_vals[!is.na(ctrl_vals) & ctrl_vals > 0]),
          var.equal = FALSE
        ),
        error = function(e) e
      )
      if (!inherits(tt, "error")) {
        pvalue <- tt$p.value
      } else {
        ttest_note <- conditionMessage(tt)
      }
    }

    list(
      n_test_clean = n_test_clean,
      n_ctrl_clean = n_ctrl_clean,
      n_test_signal = n_test_signal,
      n_ctrl_signal = n_ctrl_signal,
      n_test_nonmissing = n_test_nonmissing,
      n_ctrl_nonmissing = n_ctrl_nonmissing,
      test_present = test_present,
      ctrl_present = ctrl_present,
      mean_test_abundance = mean_test_abundance,
      mean_ctrl_abundance = mean_ctrl_abundance,
      mean_test_log2 = mean_test_log2,
      mean_ctrl_log2 = mean_ctrl_log2,
      FC = FC,
      log2FC = log2FC,
      pvalue = pvalue,
      neglog10p = if (!is.na(pvalue) && pvalue > 0) -log10(pvalue) else NA_real_,
      test_ratio_clean = calc_group_ratio(test_vals),
      ctrl_ratio_clean = calc_group_ratio(ctrl_vals),
      replicate_flag = if ((n_test_signal < 3 || n_ctrl_signal < 3)) "low_replicate" else "normal",
      ttest_note = ttest_note
    )
  }

  calc_metrics <- function(test_vals, ctrl_vals,
                               include_zero_in_mean = FALSE,
                               count_zero_as_valid_test = FALSE,
                               count_zero_as_valid_ctrl = FALSE) {
    n_test_signal <- sum(!is.na(test_vals) & test_vals > 0)
    n_ctrl_signal <- sum(!is.na(ctrl_vals) & ctrl_vals > 0)
    n_test_nonmissing <- sum(!is.na(test_vals))
    n_ctrl_nonmissing <- sum(!is.na(ctrl_vals))

    n_test_clean <- if (isTRUE(count_zero_as_valid_test)) n_test_nonmissing else n_test_signal
    n_ctrl_clean <- if (isTRUE(count_zero_as_valid_ctrl)) n_ctrl_nonmissing else n_ctrl_signal

    test_present <- n_test_signal > 0
    ctrl_present <- n_ctrl_signal > 0

    test_mean_idx <- if (isTRUE(include_zero_in_mean)) {
      !is.na(test_vals)
    } else {
      !is.na(test_vals) & test_vals > 0
    }
    ctrl_mean_idx <- if (isTRUE(include_zero_in_mean)) {
      !is.na(ctrl_vals)
    } else {
      !is.na(ctrl_vals) & ctrl_vals > 0
    }

    mean_test_abundance <- if (any(test_mean_idx) && (test_present || isTRUE(count_zero_as_valid_test))) {
      mean(test_vals[test_mean_idx])
    } else {
      NA_real_
    }
    mean_ctrl_abundance <- if (any(ctrl_mean_idx) && (ctrl_present || isTRUE(count_zero_as_valid_ctrl))) {
      mean(ctrl_vals[ctrl_mean_idx])
    } else {
      NA_real_
    }

    FC <- if (is.finite(mean_test_abundance) && is.finite(mean_ctrl_abundance) && mean_test_abundance > 0 && mean_ctrl_abundance > 0) {
      mean_test_abundance / mean_ctrl_abundance
    } else {
      NA_real_
    }
    log2FC <- if (is.finite(FC) && FC > 0) log2(FC) else NA_real_
    mean_test_log2 <- if (is.finite(mean_test_abundance) && mean_test_abundance > 0) log2(mean_test_abundance) else NA_real_
    mean_ctrl_log2 <- if (is.finite(mean_ctrl_abundance) && mean_ctrl_abundance > 0) log2(mean_ctrl_abundance) else NA_real_
    pvalue <- NA_real_
    ttest_note <- NA_character_
    if (n_test_signal >= 2 && n_ctrl_signal >= 2) {
      tt <- tryCatch(
        t.test(
          x = log2(test_vals[!is.na(test_vals) & test_vals > 0]),
          y = log2(ctrl_vals[!is.na(ctrl_vals) & ctrl_vals > 0]),
          var.equal = FALSE
        ),
        error = function(e) e
      )
      if (!inherits(tt, "error")) {
        pvalue <- tt$p.value
      } else {
        ttest_note <- conditionMessage(tt)
      }
    }
    list(
      n_test_clean = n_test_clean,
      n_ctrl_clean = n_ctrl_clean,
      n_test_signal = n_test_signal,
      n_ctrl_signal = n_ctrl_signal,
      n_test_nonmissing = n_test_nonmissing,
      n_ctrl_nonmissing = n_ctrl_nonmissing,
      test_present = test_present,
      ctrl_present = ctrl_present,
      mean_test_abundance = mean_test_abundance,
      mean_ctrl_abundance = mean_ctrl_abundance,
      mean_test_log2 = mean_test_log2,
      mean_ctrl_log2 = mean_ctrl_log2,
      FC = FC,
      log2FC = log2FC,
      pvalue = pvalue,
      neglog10p = if (!is.na(pvalue) && pvalue > 0) -log10(pvalue) else NA_real_,
      test_ratio_clean = calc_group_ratio(test_vals),
      ctrl_ratio_clean = calc_group_ratio(ctrl_vals),
      replicate_flag = if ((n_test_signal < 3 || n_ctrl_signal < 3)) "low_replicate" else "normal",
      ttest_note = ttest_note
    )
  }

  # Mode A: first impute missing values to 0; if original positive-signal n < 3, do not perform outlier removal
  test_clean_A <- if (isTRUE(compute_A)) {
    if (n_test_raw < 3) make_no_clean(test_raw_A) else clean_group_values(test_raw_A, test_cols, outlier_dev_log2, remaining_range_log2)
  } else {
    make_no_clean(test_raw_A)
  }
  ctrl_clean_A <- make_no_clean(ctrl_raw_A)

  # Modes B/C: keep balanced/strict cleaning behavior on the original NA-preserving vectors
  test_clean_B <- if (isTRUE(compute_B)) clean_group_values(test_raw, test_cols, outlier_dev_log2, remaining_range_log2) else make_no_clean(test_raw)
  ctrl_clean_B <- if (isTRUE(compute_B)) clean_group_values(ctrl_raw, ctrl_cols, outlier_dev_log2, remaining_range_log2) else make_no_clean(ctrl_raw)

  test_vals_A <- test_clean_A$values
  ctrl_vals_A <- ctrl_clean_A$values
  test_vals_B <- test_clean_B$values
  ctrl_vals_B <- ctrl_clean_B$values

  # Mode A uses raw zeros in mean abundance / FC to match the user's lenient manual workflow.
  metA <- if (isTRUE(compute_A)) calc_metrics_mode_a(test_vals_A, ctrl_vals_A) else empty_metrics()
  metB <- if (isTRUE(compute_B)) calc_metrics(test_vals_B, ctrl_vals_B, include_zero_in_mean = FALSE) else empty_metrics()

  filter_reason_A <- NA_character_
  filter_reason_B <- NA_character_

  # Mode A: lenient
  if (isTRUE(compute_A) && total_test_reps == 2) {
    if (n_test_raw > 0 && metA$n_test_clean < 2) {
      filter_reason_A <- "test_pair_incomplete"
    } else if (metA$n_test_clean == 2 && is.finite(metA$test_ratio_clean) && metA$test_ratio_clean > test_pair_fc) {
      filter_reason_A <- "test_pair_difference_gt_threshold"
    }
  } else if (isTRUE(compute_A) && total_test_reps >= 3) {
    if (metA$n_test_clean < 2) {
      filter_reason_A <- "test_insufficient_reps_after_cleaning"
    } else if (metA$n_test_clean < 3) {
      if (is.finite(metA$test_ratio_clean) && metA$test_ratio_clean > test_inconsistent_fc) {
        filter_reason_A <- "test_after_outlier_n_lt3_difference_gt_threshold"
      }
    }
  }

  # Mode B/C: balanced/strict
  if (isTRUE(compute_B) && total_test_reps == 2 && n_test_raw > 0 && metB$n_test_clean < 2) {
    filter_reason_B <- "test_pair_incomplete"
  } else if (total_test_reps == 2 && metB$n_test_clean == 2 && is.finite(metB$test_ratio_clean) && metB$test_ratio_clean > test_pair_fc) {
    filter_reason_B <- "test_pair_difference_gt_threshold"
  } else if (total_test_reps >= 3 && metB$n_test_clean >= 3 && is.finite(metB$test_ratio_clean) && metB$test_ratio_clean > test_inconsistent_fc) {
    filter_reason_B <- "test_group_inconsistent_after_cleaning"
  }

  sufficient_test_presence_A <- if (total_test_reps == 2) metA$n_test_clean == 2 else metA$n_test_clean >= 2
  sufficient_ctrl_presence_A <- if (total_ctrl_reps == 2) metA$n_ctrl_clean == 2 else metA$n_ctrl_clean >= 2
  sufficient_test_presence_B <- if (total_test_reps == 2) metB$n_test_clean == 2 else metB$n_test_clean >= 2
  sufficient_ctrl_presence_B <- if (total_ctrl_reps == 2) metB$n_ctrl_clean == 2 else metB$n_ctrl_clean >= 2

  fc_threshold_log2 <- log2(fc_cut)

  mode_a_status <- if (isTRUE(compute_A)) "ns" else NA_character_
  mode_b_status <- if (isTRUE(compute_B)) "ns" else NA_character_
  mode_c_status <- if (isTRUE(compute_C)) "ns" else NA_character_
  mode_a_subtype <- if (isTRUE(compute_A)) "ns" else NA_character_
  mode_b_subtype <- if (isTRUE(compute_B)) "ns" else NA_character_
  mode_c_subtype <- if (isTRUE(compute_C)) "ns" else NA_character_
  filter_reason <- NA_character_

  # Mode A assignment
  if (isTRUE(compute_A) && !is.na(filter_reason_A)) {
    mode_a_status <- "filtered_out"
    mode_a_subtype <- filter_reason_A
  } else if (metA$test_present && !metA$ctrl_present) {
    if (sufficient_test_presence_A) {
      mode_a_status <- "presence_up"
      mode_a_subtype <- "presence_up"
    } else {
      mode_a_status <- "filtered_out"
      mode_a_subtype <- "presence_up_insufficient_test_reps"
    }
  } else if (!metA$test_present && metA$ctrl_present) {
    if (sufficient_ctrl_presence_A) {
      mode_a_status <- "presence_down"
      mode_a_subtype <- "presence_down"
    } else {
      mode_a_status <- "filtered_out"
      mode_a_subtype <- "presence_down_insufficient_ctrl_reps"
    }
  } else if (metA$test_present && metA$ctrl_present) {
    if (is.finite(metA$log2FC) && metA$log2FC >= fc_threshold_log2) {
      mode_a_status <- "up"
      mode_a_subtype <- "fc_up"
    } else if (is.finite(metA$log2FC) && metA$log2FC <= -fc_threshold_log2) {
      mode_a_status <- "down"
      mode_a_subtype <- "fc_down"
    }
  } else {
    mode_a_subtype <- "all_missing"
  }

  # Mode B assignment
  if (isTRUE(compute_B) && !is.na(filter_reason_B)) {
    mode_b_status <- "filtered_out"
    mode_b_subtype <- filter_reason_B
  } else if (metB$test_present && !metB$ctrl_present) {
    if (sufficient_test_presence_B) {
      mode_b_status <- "presence_up"
      mode_b_subtype <- "presence_up"
    } else {
      mode_b_status <- "filtered_out"
      mode_b_subtype <- "presence_up_insufficient_test_reps"
    }
  } else if (!metB$test_present && metB$ctrl_present) {
    if (sufficient_ctrl_presence_B) {
      mode_b_status <- "presence_down"
      mode_b_subtype <- "presence_down"
    } else {
      mode_b_status <- "filtered_out"
      mode_b_subtype <- "presence_down_insufficient_ctrl_reps"
    }
  } else if (metB$test_present && metB$ctrl_present) {
    if (is.finite(metB$log2FC) && metB$log2FC >= fc_threshold_log2) {
      mode_b_status <- "up"
      mode_b_subtype <- "fc_up"
    } else if (is.finite(metB$log2FC) && metB$log2FC <= -fc_threshold_log2) {
      mode_b_status <- "down"
      mode_b_subtype <- "fc_down"
    }
  } else {
    mode_b_subtype <- "all_missing"
  }

  # Mode C assignment
  if (isTRUE(compute_C) && !is.na(filter_reason_B)) {
    mode_c_status <- "filtered_out"
    mode_c_subtype <- filter_reason_B
  } else if (metB$test_present && !metB$ctrl_present) {
    if (sufficient_test_presence_B) {
      mode_c_status <- "presence_up"
      mode_c_subtype <- "presence_up"
    } else {
      mode_c_status <- "filtered_out"
      mode_c_subtype <- "presence_up_insufficient_test_reps"
    }
  } else if (!metB$test_present && metB$ctrl_present) {
    if (sufficient_ctrl_presence_B) {
      mode_c_status <- "presence_down"
      mode_c_subtype <- "presence_down"
    } else {
      mode_c_status <- "filtered_out"
      mode_c_subtype <- "presence_down_insufficient_ctrl_reps"
    }
  } else if (metB$test_present && metB$ctrl_present) {
    if (!is.na(metB$pvalue) && metB$pvalue <= p_cut && is.finite(metB$log2FC) && metB$log2FC >= fc_threshold_log2) {
      mode_c_status <- "up"
      mode_c_subtype <- "fc_up"
    } else if (!is.na(metB$pvalue) && metB$pvalue <= p_cut && is.finite(metB$log2FC) && metB$log2FC <= -fc_threshold_log2) {
      mode_c_status <- "down"
      mode_c_subtype <- "fc_down"
    }
  } else {
    mode_c_subtype <- "all_missing"
  }

  filter_reason <- if (identical(mode_a_status, "filtered_out")) {
    mode_a_subtype
  } else if (identical(mode_b_status, "filtered_out")) {
    mode_b_subtype
  } else if (identical(mode_c_status, "filtered_out")) {
    mode_c_subtype
  } else {
    NA_character_
  }

  out <- data.frame(
    n_test_raw = n_test_raw,
    n_ctrl_raw = n_ctrl_raw,
    mode_a_n_test_clean = metA$n_test_clean,
    mode_a_n_ctrl_clean = metA$n_ctrl_clean,
    mode_b_n_test_clean = metB$n_test_clean,
    mode_b_n_ctrl_clean = metB$n_ctrl_clean,
    mode_a_mean_test_log2 = metA$mean_test_log2,
    mode_a_mean_ctrl_log2 = metA$mean_ctrl_log2,
    mode_a_mean_test_abundance = metA$mean_test_abundance,
    mode_a_mean_ctrl_abundance = metA$mean_ctrl_abundance,
    mode_a_log2FC = metA$log2FC,
    mode_a_FC = metA$FC,
    mode_a_pvalue = metA$pvalue,
    mode_a_neglog10p = metA$neglog10p,
    mode_a_replicate_flag = metA$replicate_flag,
    mode_a_test_ratio_clean = metA$test_ratio_clean,
    mode_a_ctrl_ratio_clean = metA$ctrl_ratio_clean,
    mode_b_mean_test_log2 = metB$mean_test_log2,
    mode_b_mean_ctrl_log2 = metB$mean_ctrl_log2,
    mode_b_mean_test_abundance = metB$mean_test_abundance,
    mode_b_mean_ctrl_abundance = metB$mean_ctrl_abundance,
    mode_b_log2FC = metB$log2FC,
    mode_b_FC = metB$FC,
    mode_b_pvalue = metB$pvalue,
    mode_b_neglog10p = metB$neglog10p,
    mode_b_replicate_flag = metB$replicate_flag,
    mode_b_test_ratio_clean = metB$test_ratio_clean,
    mode_b_ctrl_ratio_clean = metB$ctrl_ratio_clean,
    # legacy/general columns default to mode B metrics before active-mode overwrite
    n_test_clean = metB$n_test_clean,
    n_ctrl_clean = metB$n_ctrl_clean,
    mean_test_log2 = metB$mean_test_log2,
    mean_ctrl_log2 = metB$mean_ctrl_log2,
    mean_test_abundance = metB$mean_test_abundance,
    mean_ctrl_abundance = metB$mean_ctrl_abundance,
    log2FC = metB$log2FC,
    FC = metB$FC,
    pvalue = metB$pvalue,
    neglog10p = metB$neglog10p,
    replicate_flag = metB$replicate_flag,
    mode_a_status = mode_a_status,
    mode_b_status = mode_b_status,
    mode_c_status = mode_c_status,
    mode_a_subtype = mode_a_subtype,
    mode_b_subtype = mode_b_subtype,
    mode_c_subtype = mode_c_subtype,
    filter_reason = filter_reason,
    mode_a_removed_test_samples = if (length(test_clean_A$removed_samples) > 0) paste(test_clean_A$removed_samples, collapse = "; ") else NA_character_,
    mode_a_removed_ctrl_samples = if (length(ctrl_clean_A$removed_samples) > 0) paste(ctrl_clean_A$removed_samples, collapse = "; ") else NA_character_,
    removed_test_samples = if (length(test_clean_B$removed_samples) > 0) paste(test_clean_B$removed_samples, collapse = "; ") else NA_character_,
    removed_ctrl_samples = if (length(ctrl_clean_B$removed_samples) > 0) paste(ctrl_clean_B$removed_samples, collapse = "; ") else NA_character_,
    mode_a_test_values_clean = collapse_named_values(test_vals_A),
    mode_a_ctrl_values_clean = collapse_named_values(ctrl_vals_A),
    test_values_clean = collapse_named_values(test_vals_B),
    ctrl_values_clean = collapse_named_values(ctrl_vals_B),
    test_values_raw = collapse_named_values(test_raw),
    ctrl_values_raw = collapse_named_values(ctrl_raw),
    ttest_note = metB$ttest_note,
    stringsAsFactors = FALSE
  )

  cleaned_cols <- data.frame(
    as.list(setNames(as.list(test_vals_A), paste0("ModeA_Clean_", test_cols))),
    as.list(setNames(as.list(ctrl_vals_A), paste0("ModeA_Clean_", ctrl_cols))),
    as.list(setNames(as.list(test_vals_B), paste0("Clean_", test_cols))),
    as.list(setNames(as.list(ctrl_vals_B), paste0("Clean_", ctrl_cols))),
    stringsAsFactors = FALSE
  )

  cbind(row_df, cleaned_cols, out)
}


impute_left_tail_matrix <- function(df_num, width = 0.3, downshift = 1.8, seed = 123) {
  mat <- as.data.frame(df_num)
  for (nm in colnames(mat)) mat[[nm]] <- as.numeric(mat[[nm]])

  if (!is.null(seed) && !is.na(seed)) set.seed(seed)

  global_obs <- unlist(mat, use.names = FALSE)
  global_obs <- global_obs[!is.na(global_obs) & global_obs > 0]
  if (length(global_obs) == 0) global_obs <- c(1, 2, 4)
  global_log <- log2(global_obs)
  global_mu <- mean(global_log, na.rm = TRUE)
  global_sd <- stats::sd(global_log, na.rm = TRUE)
  if (!is.finite(global_sd) || global_sd <= 0) global_sd <- 1

  for (nm in colnames(mat)) {
    x <- as.numeric(mat[[nm]])
    miss <- which(is.na(x) | x <= 0)
    obs <- x[!is.na(x) & x > 0]

    if (length(miss) > 0) {
      if (length(obs) >= 5) {
        log_obs <- log2(obs)
        mu <- mean(log_obs, na.rm = TRUE)
        sdv <- stats::sd(log_obs, na.rm = TRUE)
        if (!is.finite(sdv) || sdv <= 0) sdv <- global_sd
      } else {
        mu <- global_mu
        sdv <- global_sd
      }
      if (!is.finite(sdv) || sdv <= 0) sdv <- 1
      imp_log <- stats::rnorm(length(miss), mean = mu - downshift * sdv, sd = width * sdv)
      x[miss] <- 2^imp_log
    }
    mat[[nm]] <- x
  }
  mat
}




impute_fixed_min_matrix <- function(df_num, min_value = NA_real_) {
  mat <- as.data.frame(df_num)
  for (nm in colnames(mat)) mat[[nm]] <- as.numeric(mat[[nm]])

  if (!is.finite(min_value) || is.na(min_value) || min_value <= 0) {
    global_obs <- unlist(mat, use.names = FALSE)
    global_obs <- global_obs[!is.na(global_obs) & global_obs > 0]
    if (length(global_obs) == 0) {
      min_value <- NA_real_
    } else {
      min_value <- min(global_obs, na.rm = TRUE)
    }
  }

  if (!is.finite(min_value) || is.na(min_value) || min_value <= 0) {
    return(mat)
  }

  for (nm in colnames(mat)) {
    x <- as.numeric(mat[[nm]])
    miss <- which(is.na(x) | x <= 0)
    if (length(miss) > 0) x[miss] <- min_value
    mat[[nm]] <- x
  }
  mat
}

calc_exploratory_metrics <- function(row_df, test_cols, ctrl_cols, min_nonzero = NA_real_) {
  test_vals <- as.numeric(row_df[1, test_cols, drop = TRUE])
  ctrl_vals <- as.numeric(row_df[1, ctrl_cols, drop = TRUE])

  if (!is.finite(min_nonzero) || is.na(min_nonzero) || min_nonzero <= 0) {
    return(data.frame(
      explore_mean_test = NA_real_,
      explore_mean_ctrl = NA_real_,
      explore_log2FC = NA_real_,
      explore_FC = NA_real_,
      explore_pvalue = NA_real_,
      explore_neglog10p = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  test_adj <- ifelse(is.na(test_vals) | test_vals == 0, min_nonzero, test_vals)
  ctrl_adj <- ifelse(is.na(ctrl_vals) | ctrl_vals == 0, min_nonzero, ctrl_vals)

  mean_test <- mean(test_adj, na.rm = TRUE)
  mean_ctrl <- mean(ctrl_adj, na.rm = TRUE)
  explore_FC <- if (is.finite(mean_test) && is.finite(mean_ctrl) && mean_test > 0 && mean_ctrl > 0) mean_test / mean_ctrl else NA_real_
  explore_log2FC <- if (is.finite(explore_FC) && explore_FC > 0) log2(explore_FC) else NA_real_

  explore_pvalue <- NA_real_
  tt <- tryCatch(t.test(test_adj, ctrl_adj, var.equal = TRUE), error = function(e) e)
  if (!inherits(tt, "error")) explore_pvalue <- tt$p.value
  explore_neglog10p <- if (!is.na(explore_pvalue) && explore_pvalue > 0) -log10(explore_pvalue) else NA_real_

  data.frame(
    explore_mean_test = mean_test,
    explore_mean_ctrl = mean_ctrl,
    explore_log2FC = explore_log2FC,
    explore_FC = explore_FC,
    explore_pvalue = explore_pvalue,
    explore_neglog10p = explore_neglog10p,
    stringsAsFactors = FALSE
  )
}

calc_visual_metrics <- function(row_df, test_cols, ctrl_cols) {
  test_vals <- as.numeric(row_df[1, test_cols, drop = TRUE])
  ctrl_vals <- as.numeric(row_df[1, ctrl_cols, drop = TRUE])

  mean_test_abundance <- mean(test_vals, na.rm = TRUE)
  mean_ctrl_abundance <- mean(ctrl_vals, na.rm = TRUE)
  mean_test_log2 <- log2(mean_test_abundance)
  mean_ctrl_log2 <- log2(mean_ctrl_abundance)
  viz_FC <- mean_test_abundance / mean_ctrl_abundance
  viz_log2FC <- log2(viz_FC)

  viz_pvalue <- NA_real_
  if (length(test_vals) >= 2 && length(ctrl_vals) >= 2) {
    tt <- tryCatch(
      t.test(log2(test_vals), log2(ctrl_vals), var.equal = FALSE),
      error = function(e) e
    )
    if (!inherits(tt, "error")) viz_pvalue <- tt$p.value
  }
  viz_neglog10p <- if (!is.na(viz_pvalue) && viz_pvalue > 0) -log10(viz_pvalue) else NA_real_

  data.frame(
    viz_mean_test_log2 = mean_test_log2,
    viz_mean_ctrl_log2 = mean_ctrl_log2,
    viz_log2FC = viz_log2FC,
    viz_FC = viz_FC,
    viz_pvalue = viz_pvalue,
    viz_neglog10p = viz_neglog10p,
    stringsAsFactors = FALSE
  )
}


make_volcano_plot <- function(df, fc_cut = 2, p_cut = 0.05,
                              col_up = "#F39B7F", col_down = "#4DBBD5", col_ns = "#A6A6A6",
                              col_line = "#000000",
                              axis_text_size = 18, axis_title_size = 24, legend_text_size = 14,
                              pt_alpha = 0.7, pt_size = 3,
                              presence_display = "count",
                              use_visual_imputation = TRUE,
                              active_mode_label = "A",
                              label_top_genes = FALSE,
                              label_top_n = 15,
                              display_style = c("candidate", "exploratory")) {
  display_style <- match.arg(display_style)
  df_plot <- df

  if (display_style == "exploratory" && all(c("explore_log2FC", "explore_neglog10p") %in% colnames(df_plot))) {
    df_plot$plot_x <- df_plot$explore_log2FC
    df_plot$plot_y <- df_plot$explore_neglog10p
  } else {
    if ("active_log2FC" %in% colnames(df_plot)) {
      df_plot$plot_x <- df_plot$active_log2FC
    } else {
      df_plot$plot_x <- df_plot$log2FC
    }
    if ("active_neglog10p" %in% colnames(df_plot)) {
      df_plot$plot_y <- df_plot$active_neglog10p
    } else {
      df_plot$plot_y <- df_plot$neglog10p
    }
  }

  df_plot$point_type <- ifelse(df_plot$status %in% c("presence_up", "presence_down"), "presence", "measured")
  df_plot$volcano_group <- dplyr::case_when(
    df_plot$status %in% c("up", "presence_up") ~ "up",
    df_plot$status %in% c("down", "presence_down") ~ "down",
    TRUE ~ "ns"
  )

  if (display_style == "candidate" && presence_display %in% c("plot_fixed", "plot_random") && isTRUE(use_visual_imputation)) {
    idx_presence <- which(df_plot$status %in% c("presence_up", "presence_down"))
    if (length(idx_presence) > 0) {
      if (identical(presence_display, "plot_fixed") && all(c("viz_log2FC", "viz_neglog10p") %in% colnames(df_plot))) {
        has_viz <- !is.na(df_plot$viz_log2FC[idx_presence]) & !is.na(df_plot$viz_neglog10p[idx_presence])
        if (any(has_viz)) {
          df_plot$plot_x[idx_presence[has_viz]] <- df_plot$viz_log2FC[idx_presence[has_viz]]
          df_plot$plot_y[idx_presence[has_viz]] <- df_plot$viz_neglog10p[idx_presence[has_viz]]
        }
      }
      if (identical(presence_display, "plot_random") && all(c("viz_rand_log2FC", "viz_rand_neglog10p") %in% colnames(df_plot))) {
        has_viz <- !is.na(df_plot$viz_rand_log2FC[idx_presence]) & !is.na(df_plot$viz_rand_neglog10p[idx_presence])
        if (any(has_viz)) {
          df_plot$plot_x[idx_presence[has_viz]] <- df_plot$viz_rand_log2FC[idx_presence[has_viz]]
          df_plot$plot_y[idx_presence[has_viz]] <- df_plot$viz_rand_neglog10p[idx_presence[has_viz]]
        }
      }
    }
  }

  if (display_style == "exploratory") {
    df_plot <- df_plot %>% filter(status %in% c("up", "down", "ns", "presence_up", "presence_down"))
  } else if (presence_display %in% c("plot_fixed", "plot_random")) {
    df_plot <- df_plot %>% filter(status %in% c("up", "down", "ns", "presence_up", "presence_down"))
  } else {
    df_plot <- df_plot %>% filter(status %in% c("up", "down", "ns"))
  }

  df_plot <- df_plot %>% filter(is.finite(plot_x), !is.na(plot_y))

  if (nrow(df_plot) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Volcano plot\nNo finite points to display", size = 6) +
        theme_void()
    )
  }

  label_df <- data.frame()
  if (isTRUE(label_top_genes) && "GeneName" %in% colnames(df_plot)) {
    label_up <- df_plot %>%
      filter(volcano_group == "up", !is.na(GeneName), GeneName != "") %>%
      arrange(desc(plot_y), desc(abs(plot_x))) %>%
      slice_head(n = label_top_n)
    label_down <- df_plot %>%
      filter(volcano_group == "down", !is.na(GeneName), GeneName != "") %>%
      arrange(desc(plot_y), desc(abs(plot_x))) %>%
      slice_head(n = label_top_n)
    label_df <- bind_rows(label_up, label_down)
  }

  fc_threshold_log2 <- log2(fc_cut)
  subtitle_txt <- paste0(
    "Mode: ", active_mode_label,
    " | style: ", ifelse(display_style == "exploratory", "Exploratory", "Candidate-aware"),
    " | presence_up: ", sum(df$status == "presence_up", na.rm = TRUE),
    " | presence_down: ", sum(df$status == "presence_down", na.rm = TRUE),
    " | filtered_out: ", sum(df$status == "filtered_out", na.rm = TRUE)
  )

  p <- ggplot(df_plot, aes(x = plot_x, y = plot_y, color = volcano_group, shape = point_type)) +
    geom_point(alpha = pt_alpha, size = pt_size) +
    geom_vline(xintercept = c(-fc_threshold_log2, fc_threshold_log2), color = col_line, linetype = "dashed") +
    geom_hline(yintercept = -log10(p_cut), color = col_line, linetype = "dashed") +
    scale_color_manual(values = c(up = col_up, down = col_down, ns = col_ns)) +
    scale_shape_manual(values = c(measured = 16, presence = 17)) +
    labs(
      title = "Volcano plot",
      subtitle = subtitle_txt,
      x = "log2 Fold Change (Test vs Ctrl)",
      y = "-log10(p value)",
      color = NULL,
      shape = NULL
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = axis_text_size),
      axis.title = element_text(size = axis_title_size),
      legend.text = element_text(size = legend_text_size),
      plot.title = element_text(face = "bold", size = axis_title_size),
      plot.subtitle = element_text(size = legend_text_size)
    )

  if (nrow(label_df) > 0) {
    p <- p + geom_text(
      data = label_df,
      aes(label = GeneName),
      check_overlap = TRUE,
      vjust = -0.5,
      size = 3.5,
      show.legend = FALSE
    )
  }

  p
}

prepare_mode_view_df <- function(df, mode_label = c("A", "B", "C")) {
  mode_label <- match.arg(mode_label)
  out <- as.data.frame(df)

  if (!"GeneName" %in% colnames(out)) return(out)

  if (mode_label == "A") {
    out$status <- out$mode_a_status
    out$subtype <- out$mode_a_subtype
    out$active_log2FC <- out$mode_a_log2FC
    out$active_FC <- out$mode_a_FC
    out$active_pvalue <- out$mode_a_pvalue
    out$active_neglog10p <- out$mode_a_neglog10p
    out$active_mean_test_abundance <- out$mode_a_mean_test_abundance
    out$active_mean_ctrl_abundance <- out$mode_a_mean_ctrl_abundance
    out$active_replicate_flag <- out$mode_a_replicate_flag
  } else if (mode_label == "B") {
    out$status <- out$mode_b_status
    out$subtype <- out$mode_b_subtype
    out$active_log2FC <- out$mode_b_log2FC
    out$active_FC <- out$mode_b_FC
    out$active_pvalue <- out$mode_b_pvalue
    out$active_neglog10p <- out$mode_b_neglog10p
    out$active_mean_test_abundance <- out$mode_b_mean_test_abundance
    out$active_mean_ctrl_abundance <- out$mode_b_mean_ctrl_abundance
    out$active_replicate_flag <- out$mode_b_replicate_flag
  } else {
    out$status <- out$mode_c_status
    out$subtype <- out$mode_c_subtype
    out$active_log2FC <- out$mode_b_log2FC
    out$active_FC <- out$mode_b_FC
    out$active_pvalue <- out$mode_b_pvalue
    out$active_neglog10p <- out$mode_b_neglog10p
    out$active_mean_test_abundance <- out$mode_b_mean_test_abundance
    out$active_mean_ctrl_abundance <- out$mode_b_mean_ctrl_abundance
    out$active_replicate_flag <- out$mode_b_replicate_flag
  }

  out
}


run_gene_list_mode <- function(input_path, desc_col = NULL,
                               target_ions = c("calcium", "magnesium", "manganese", "zinc"),
                               p_cut = 0.05, q_cut = 0.2, show_n = 12) {
  df_raw <- read_any_table(input_path)
  df <- extract_genename(df_raw, desc_col = desc_col)

  annot_res <- annotate_genes_dataframe(df, gene_col = "GeneName", target_ions = target_ions)

  genes <- unique(df$GeneName[!is.na(df$GeneName) & df$GeneName != ""])
  enrich_res <- run_classic_enrichment(
    symbol_genes = genes,
    universe_symbols = NULL,
    p_cut = p_cut,
    q_cut = q_cut,
    show_n = show_n
  )

  list(
    input_preview = utils::head(df_raw, 20),
    annotated_original = annot_res$annotated_df,
    protein_summary = annot_res$protein_summary,
    go_hits = annot_res$go_hits,
    ion_term_summary = annot_res$ion_term_summary,
    overview = annot_res$overview,
    gene_map = annot_res$gene_map_go,
    gene_map_entrez = annot_res$gene_map_entrez,
    selected_ions = annot_res$selected_ions,
    gene_count = length(genes),
    entrez_count = enrich_res$entrez_count,
    enrich = enrich_res
  )
}

run_raw_mode <- function(input_path, desc_col = NULL,
                         target_ions = c("calcium", "magnesium", "manganese", "zinc"),
                         analysis_mode = c("A", "B", "C"),
                         precompute_all_modes = FALSE,
                         confidence_filter = c("high", "high_medium", "all"),
                         contaminant_filter = c("common", "keratin", "none"),
                         deduplicate_by_gene = TRUE,
                         outlier_dev_log2 = 1,
                         remaining_range_log2 = 0.5,
                         test_pair_fc = 5,
                         test_inconsistent_fc = 5,
                         fc_cut = 2,
                         p_cut = 0.05,
                         enrich_p_cut = 0.05,
                         enrich_q_cut = 0.2,
                         enrich_strategy = c("exploratory", "context"),
                         enrich_scope = c("both", "up", "down"),
                         show_n = 12,
                         impute_for_visual = TRUE,
                         impute_width = 0.3,
                         impute_downshift = 1.8,
                         impute_seed = 123) {

  analysis_mode <- match.arg(analysis_mode)
  confidence_filter <- match.arg(confidence_filter)
  contaminant_filter <- match.arg(contaminant_filter)
  enrich_strategy <- match.arg(enrich_strategy)
  enrich_scope <- match.arg(enrich_scope)
  df_raw <- read_any_table(input_path)
  df <- extract_genename(df_raw, desc_col = desc_col)
  description_col <- detect_description_column(df, preferred = desc_col)
  confidence_col <- detect_confidence_column(df)
  conf_res <- apply_confidence_filter(df, confidence_col = confidence_col, confidence_filter = confidence_filter)
  df <- conf_res$df
  contam_res <- apply_contaminant_filter(df, gene_col = "GeneName", desc_col = description_col, contaminant_filter = contaminant_filter)
  df <- contam_res$df
  groups <- detect_group_columns(df)
  test_cols <- groups$test_cols
  ctrl_cols <- groups$ctrl_cols

  if (length(test_cols) < 2 || length(ctrl_cols) < 2) {
    stop("Raw-data mode requires at least 2 Test columns (Test1...) and 2 Ctrl columns (Ctrl1...).")
  }

  abundance_cols <- c(test_cols, ctrl_cols)
  for (col in abundance_cols) df[[col]] <- safe_parse_numeric(df[[col]])

  dup_res <- resolve_duplicate_genes(df, gene_col = "GeneName", abundance_cols = abundance_cols)
  duplicate_summary <- dup_res$duplicate_summary
  if (isTRUE(deduplicate_by_gene)) {
    df_dedup <- dup_res$df
  } else {
    df$.row_id_original <- seq_len(nrow(df))
    df_dedup <- df %>% filter(!is.na(GeneName), GeneName != "")
  }

  modes_to_compute <- if (isTRUE(precompute_all_modes)) {
    c("A", "B", "C")
  } else if (identical(analysis_mode, "A")) {
    c("A")
  } else if (identical(analysis_mode, "B")) {
    c("B")
  } else {
    c("B", "C")
  }

  analysis_list <- lapply(seq_len(nrow(df_dedup)), function(i) {
    compare_single_gene(
      row_df = df_dedup[i, , drop = FALSE],
      test_cols = test_cols,
      ctrl_cols = ctrl_cols,
      total_test_reps = length(test_cols),
      total_ctrl_reps = length(ctrl_cols),
      outlier_dev_log2 = outlier_dev_log2,
      remaining_range_log2 = remaining_range_log2,
      test_pair_fc = test_pair_fc,
      test_inconsistent_fc = test_inconsistent_fc,
      fc_cut = fc_cut,
      p_cut = p_cut,
      modes_to_compute = modes_to_compute
    )
  })
  diff_df <- dplyr::bind_rows(analysis_list)

  all_vals <- unlist(df_dedup[, abundance_cols, drop = FALSE], use.names = FALSE)
  all_vals <- suppressWarnings(as.numeric(all_vals))
  min_nonzero <- min(all_vals[is.finite(all_vals) & !is.na(all_vals) & all_vals > 0], na.rm = TRUE)
  if (!is.finite(min_nonzero)) min_nonzero <- NA_real_

  if (isTRUE(impute_for_visual)) {
    impute_df_fixed <- df_dedup
    impute_df_fixed[, abundance_cols] <- impute_fixed_min_matrix(
      impute_df_fixed[, abundance_cols, drop = FALSE],
      min_value = min_nonzero
    )
    viz_list <- lapply(seq_len(nrow(impute_df_fixed)), function(i) {
      calc_visual_metrics(impute_df_fixed[i, , drop = FALSE], test_cols = test_cols, ctrl_cols = ctrl_cols)
    })
    viz_df <- dplyr::bind_rows(viz_list)

    impute_df_rand <- df_dedup
    impute_df_rand[, abundance_cols] <- impute_left_tail_matrix(
      impute_df_rand[, abundance_cols, drop = FALSE],
      width = impute_width,
      downshift = impute_downshift,
      seed = impute_seed
    )
    viz_rand_list <- lapply(seq_len(nrow(impute_df_rand)), function(i) {
      calc_visual_metrics(impute_df_rand[i, , drop = FALSE], test_cols = test_cols, ctrl_cols = ctrl_cols)
    })
    viz_rand_df <- dplyr::bind_rows(viz_rand_list)
    colnames(viz_rand_df) <- c("viz_rand_mean_test_log2", "viz_rand_mean_ctrl_log2", "viz_rand_log2FC", "viz_rand_FC", "viz_rand_pvalue", "viz_rand_neglog10p")
  } else {
    viz_df <- data.frame(
      viz_mean_test_log2 = rep(NA_real_, nrow(diff_df)),
      viz_mean_ctrl_log2 = rep(NA_real_, nrow(diff_df)),
      viz_log2FC = rep(NA_real_, nrow(diff_df)),
      viz_FC = rep(NA_real_, nrow(diff_df)),
      viz_pvalue = rep(NA_real_, nrow(diff_df)),
      viz_neglog10p = rep(NA_real_, nrow(diff_df))
    )
    viz_rand_df <- data.frame(
      viz_rand_mean_test_log2 = rep(NA_real_, nrow(diff_df)),
      viz_rand_mean_ctrl_log2 = rep(NA_real_, nrow(diff_df)),
      viz_rand_log2FC = rep(NA_real_, nrow(diff_df)),
      viz_rand_FC = rep(NA_real_, nrow(diff_df)),
      viz_rand_pvalue = rep(NA_real_, nrow(diff_df)),
      viz_rand_neglog10p = rep(NA_real_, nrow(diff_df))
    )
  }
  diff_df <- cbind(diff_df, viz_df, viz_rand_df)

  exploratory_available <- length(test_cols) == 2 && length(ctrl_cols) == 2
  if (isTRUE(exploratory_available)) {
    exploratory_list <- lapply(seq_len(nrow(df_dedup)), function(i) {
      calc_exploratory_metrics(df_dedup[i, , drop = FALSE], test_cols = test_cols, ctrl_cols = ctrl_cols, min_nonzero = min_nonzero)
    })
    exploratory_df <- dplyr::bind_rows(exploratory_list)
  } else {
    exploratory_df <- data.frame(
      explore_mean_test = rep(NA_real_, nrow(diff_df)),
      explore_mean_ctrl = rep(NA_real_, nrow(diff_df)),
      explore_log2FC = rep(NA_real_, nrow(diff_df)),
      explore_FC = rep(NA_real_, nrow(diff_df)),
      explore_pvalue = rep(NA_real_, nrow(diff_df)),
      explore_neglog10p = rep(NA_real_, nrow(diff_df)),
      stringsAsFactors = FALSE
    )
  }
  diff_df <- cbind(diff_df, exploratory_df)

  diff_df$active_mode <- analysis_mode
  diff_df$active_status <- dplyr::case_when(
    analysis_mode == "A" ~ diff_df$mode_a_status,
    analysis_mode == "B" ~ diff_df$mode_b_status,
    TRUE ~ diff_df$mode_c_status
  )
  diff_df$active_subtype <- dplyr::case_when(
    analysis_mode == "A" ~ diff_df$mode_a_subtype,
    analysis_mode == "B" ~ diff_df$mode_b_subtype,
    TRUE ~ diff_df$mode_c_subtype
  )
  diff_df$active_mean_test_abundance <- dplyr::case_when(
    analysis_mode == "A" ~ diff_df$mode_a_mean_test_abundance,
    TRUE ~ diff_df$mode_b_mean_test_abundance
  )
  diff_df$active_mean_ctrl_abundance <- dplyr::case_when(
    analysis_mode == "A" ~ diff_df$mode_a_mean_ctrl_abundance,
    TRUE ~ diff_df$mode_b_mean_ctrl_abundance
  )
  diff_df$active_mean_test_log2 <- dplyr::case_when(
    analysis_mode == "A" ~ diff_df$mode_a_mean_test_log2,
    TRUE ~ diff_df$mode_b_mean_test_log2
  )
  diff_df$active_mean_ctrl_log2 <- dplyr::case_when(
    analysis_mode == "A" ~ diff_df$mode_a_mean_ctrl_log2,
    TRUE ~ diff_df$mode_b_mean_ctrl_log2
  )
  diff_df$active_FC <- dplyr::case_when(
    analysis_mode == "A" ~ diff_df$mode_a_FC,
    TRUE ~ diff_df$mode_b_FC
  )
  diff_df$active_log2FC <- dplyr::case_when(
    analysis_mode == "A" ~ diff_df$mode_a_log2FC,
    TRUE ~ diff_df$mode_b_log2FC
  )
  diff_df$active_pvalue <- dplyr::case_when(
    analysis_mode == "A" ~ diff_df$mode_a_pvalue,
    TRUE ~ diff_df$mode_b_pvalue
  )
  diff_df$active_neglog10p <- dplyr::case_when(
    analysis_mode == "A" ~ diff_df$mode_a_neglog10p,
    TRUE ~ diff_df$mode_b_neglog10p
  )
  diff_df$active_replicate_flag <- dplyr::case_when(
    analysis_mode == "A" ~ diff_df$mode_a_replicate_flag,
    TRUE ~ diff_df$mode_b_replicate_flag
  )
  diff_df$status <- diff_df$active_status
  diff_df$subtype <- diff_df$active_subtype
  diff_df$mean_test_abundance <- diff_df$active_mean_test_abundance
  diff_df$mean_ctrl_abundance <- diff_df$active_mean_ctrl_abundance
  diff_df$mean_test_log2 <- diff_df$active_mean_test_log2
  diff_df$mean_ctrl_log2 <- diff_df$active_mean_ctrl_log2
  diff_df$FC <- diff_df$active_FC
  diff_df$log2FC <- diff_df$active_log2FC
  diff_df$pvalue <- diff_df$active_pvalue
  diff_df$neglog10p <- diff_df$active_neglog10p
  diff_df$replicate_flag <- diff_df$active_replicate_flag

  filtered_out_df <- diff_df %>% filter(status == "filtered_out") %>% arrange(GeneName)
  active_up_df <- diff_df %>% filter(status %in% c("up", "presence_up")) %>% arrange(desc(status), desc(FC), GeneName)
  active_down_df <- diff_df %>% filter(status %in% c("down", "presence_down")) %>% arrange(status, FC, GeneName)

  # Keep full-table ion flags in the annotated differential table for browsing,
  # but drive ion summary/hits/term-summary by the same scope used for enrichment.
  ion_res_all <- annotate_genes_dataframe(diff_df, gene_col = "GeneName", target_ions = target_ions)
  annotated_diff_df <- ion_res_all$annotated_df
  active_up_df <- annotated_diff_df %>% filter(status %in% c("up", "presence_up")) %>% arrange(desc(status), desc(FC), GeneName)
  active_down_df <- annotated_diff_df %>% filter(status %in% c("down", "presence_down")) %>% arrange(status, FC, GeneName)

  universe_genes <- unique(annotated_diff_df$GeneName[annotated_diff_df$status != "filtered_out"])
  enrich_universe <- if (identical(enrich_strategy, "context")) universe_genes else NULL
  up_genes <- unique(active_up_df$GeneName)
  down_genes <- unique(active_down_df$GeneName)
  if (identical(enrich_scope, "up")) {
    ion_res_scope <- annotate_genes_dataframe(active_up_df, gene_col = "GeneName", target_ions = target_ions)
    ion_res_scope$overview <- add_direction_column(ion_res_scope$overview, "Up")
    ion_res_scope$protein_summary <- add_direction_column(ion_res_scope$protein_summary, "Up")
    ion_res_scope$go_hits <- add_direction_column(ion_res_scope$go_hits, "Up")
    ion_res_scope$ion_term_summary <- add_direction_column(ion_res_scope$ion_term_summary, "Up")
  } else if (identical(enrich_scope, "down")) {
    ion_res_scope <- annotate_genes_dataframe(active_down_df, gene_col = "GeneName", target_ions = target_ions)
    ion_res_scope$overview <- add_direction_column(ion_res_scope$overview, "Down")
    ion_res_scope$protein_summary <- add_direction_column(ion_res_scope$protein_summary, "Down")
    ion_res_scope$go_hits <- add_direction_column(ion_res_scope$go_hits, "Down")
    ion_res_scope$ion_term_summary <- add_direction_column(ion_res_scope$ion_term_summary, "Down")
  } else {
    ion_res_up <- annotate_genes_dataframe(active_up_df, gene_col = "GeneName", target_ions = target_ions)
    ion_res_down <- annotate_genes_dataframe(active_down_df, gene_col = "GeneName", target_ions = target_ions)
    ion_res_scope <- list(
      overview = dplyr::bind_rows(
        add_direction_column(ion_res_up$overview, "Up"),
        add_direction_column(ion_res_down$overview, "Down")
      ),
      protein_summary = dplyr::bind_rows(
        add_direction_column(ion_res_up$protein_summary, "Up"),
        add_direction_column(ion_res_down$protein_summary, "Down")
      ),
      go_hits = dplyr::bind_rows(
        add_direction_column(ion_res_up$go_hits, "Up"),
        add_direction_column(ion_res_down$go_hits, "Down")
      ),
      ion_term_summary = dplyr::bind_rows(
        add_direction_column(ion_res_up$ion_term_summary, "Up"),
        add_direction_column(ion_res_down$ion_term_summary, "Down")
      ),
      selected_ions = unique(c(ion_res_up$selected_ions, ion_res_down$selected_ions))
    )
  }

  enrich_up <- if (enrich_scope %in% c("both", "up")) {
    run_classic_enrichment(
      symbol_genes = up_genes,
      universe_symbols = enrich_universe,
      p_cut = enrich_p_cut,
      q_cut = enrich_q_cut,
      show_n = show_n
    )
  } else {
    run_classic_enrichment(
      symbol_genes = character(),
      universe_symbols = enrich_universe,
      p_cut = enrich_p_cut,
      q_cut = enrich_q_cut,
      show_n = show_n
    )
  }
  enrich_down <- if (enrich_scope %in% c("both", "down")) {
    run_classic_enrichment(
      symbol_genes = down_genes,
      universe_symbols = enrich_universe,
      p_cut = enrich_p_cut,
      q_cut = enrich_q_cut,
      show_n = show_n
    )
  } else {
    run_classic_enrichment(
      symbol_genes = character(),
      universe_symbols = enrich_universe,
      p_cut = enrich_p_cut,
      q_cut = enrich_q_cut,
      show_n = show_n
    )
  }

  summary_df <- data.frame(
    metric = c(
      "input_rows", "analysis_rows", "deduplicate_by_gene", "duplicated_gene_sets_detected",
      "contaminant_filter", "contaminants_removed_n",
      "test_column_n", "ctrl_column_n",
      "mode_a_up_n", "mode_a_down_n", "mode_b_up_n", "mode_b_down_n", "mode_c_up_n", "mode_c_down_n",
      "current_presence_up_n", "current_presence_down_n", "current_ns_n", "current_filtered_out_n",
      "active_up_n", "active_down_n", "universe_gene_n", "confidence_filter", "precompute_all_modes", "enrich_strategy", "enrich_scope", "ion_annotation_scope"
    ),
    value = c(
      nrow(df_raw), nrow(df_dedup), if (isTRUE(deduplicate_by_gene)) "TRUE" else "FALSE", if (nrow(duplicate_summary) > 0) nrow(duplicate_summary) else 0,
      contam_res$filter_label, nrow(contam_res$removed_df),
      length(test_cols), length(ctrl_cols),
      sum(diff_df$mode_a_status == "up", na.rm = TRUE),
      sum(diff_df$mode_a_status == "down", na.rm = TRUE),
      sum(diff_df$mode_b_status == "up", na.rm = TRUE),
      sum(diff_df$mode_b_status == "down", na.rm = TRUE),
      sum(diff_df$mode_c_status == "up", na.rm = TRUE),
      sum(diff_df$mode_c_status == "down", na.rm = TRUE),
      sum(diff_df$status == "presence_up", na.rm = TRUE),
      sum(diff_df$status == "presence_down", na.rm = TRUE),
      sum(diff_df$status == "ns", na.rm = TRUE),
      sum(diff_df$status == "filtered_out", na.rm = TRUE),
      nrow(active_up_df), nrow(active_down_df), length(universe_genes), conf_res$levels_kept, if (isTRUE(precompute_all_modes)) "TRUE" else "FALSE", enrich_strategy, enrich_scope, enrich_scope
    ),
    stringsAsFactors = FALSE
  )

  list(
    input_preview = utils::head(df_raw, 20),
    raw_df = df_raw,
    dedup_df = df_dedup,
    contaminants_removed_df = contam_res$removed_df,
    contaminant_filter = contam_res$filter_label,
    duplicate_summary = duplicate_summary,
    diff_df = diff_df,
    annotated_diff_df = annotated_diff_df,
    filtered_out_df = filtered_out_df,
    active_up_df = active_up_df,
    active_down_df = active_down_df,
    ion_overview = ion_res_scope$overview,
    ion_term_summary = ion_res_scope$ion_term_summary,
    ion_protein_summary = ion_res_scope$protein_summary,
    ion_go_hits = ion_res_scope$go_hits,
    selected_ions = ion_res_scope$selected_ions,
    test_cols = test_cols,
    ctrl_cols = ctrl_cols,
    universe_genes = universe_genes,
    up_genes = up_genes,
    down_genes = down_genes,
    enrich_up = enrich_up,
    enrich_down = enrich_down,
    summary_df = summary_df,
    analysis_mode = analysis_mode,
    precompute_all_modes = precompute_all_modes,
    confidence_col = confidence_col,
    confidence_filter = conf_res$levels_kept,
    description_col = description_col
  )
}

ui <- navbarPage(
  "Protein Candidate Explorer App",
  tabPanel(
    "Gene-list mode",
    sidebarLayout(
      sidebarPanel(
        fileInput("file_gene", "上传已筛好的基因/蛋白表", accept = c(".xlsx", ".xls", ".csv", ".tsv", ".txt")),
        textInput("desc_col_gene", "Description 列名（可留空自动识别）", value = ""),
        numericInput("p_cut_gene", "enrichment pvalue cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput("q_cut_gene", "GO qvalue cutoff", value = 0.2, min = 0, max = 1, step = 0.01),
        numericInput("show_n_gene", "dotplot 显示条目数", value = 12, min = 5, max = 50, step = 1),
        selectInput(
          "plot_format_gene", "图形导出格式（PDF/EMF 更适合后续编辑）",
          choices = c("PNG" = "png", "JPG" = "jpg", "EMF" = "emf", "PDF" = "pdf"),
          selected = "png"
        ),
        tags$hr(),
        tags$strong("离子注释"),
        checkboxGroupInput(
          "ions_gene", "选择离子",
          choices = c("钙 calcium" = "calcium", "镁 magnesium" = "magnesium", "锰 manganese" = "manganese", "锌 zinc" = "zinc"),
          selected = c("calcium", "magnesium", "manganese", "zinc")
        ),
        actionButton("run_gene", "运行 Gene-list 分析", class = "btn-primary"),
        br(), br(),
        helpText("适用于已经筛好的基因/蛋白列表。程序会保留你当前已有的离子注释 + BP/MF/CC/KEGG enrichment 工作流。"),
        verbatimTextOutput("status_gene"),
        hr(),
        downloadButton("download_gene_annotated", "下载 annotated 原表"),
        br(), br(),
        downloadButton("download_gene_summary", "下载 ion protein summary"),
        br(), br(),
        downloadButton("download_gene_hits", "下载 ion GO hits"),
        br(), br(),
        downloadButton("download_gene_bp", "下载 GO BP 表格"),
        br(), br(),
        downloadButton("download_gene_bp_plot", "下载 GO BP 图"),
        br(), br(),
        downloadButton("download_gene_mf", "下载 GO MF 表格"),
        br(), br(),
        downloadButton("download_gene_mf_plot", "下载 GO MF 图"),
        br(), br(),
        downloadButton("download_gene_cc", "下载 GO CC 表格"),
        br(), br(),
        downloadButton("download_gene_cc_plot", "下载 GO CC 图"),
        br(), br(),
        downloadButton("download_gene_kegg", "下载 KEGG 表格"),
        br(), br(),
        downloadButton("download_gene_kegg_plot", "下载 KEGG 图")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("概览", br(), DTOutput("gene_overview_tbl")),
          tabPanel("输入预览", br(), DTOutput("gene_input_tbl")),
          tabPanel("Ion Protein summary", br(), DTOutput("gene_summary_tbl")),
          tabPanel("Ion GO hits", br(), DTOutput("gene_hits_tbl")),
          tabPanel("Annotated 原表", br(), DTOutput("gene_annotated_tbl")),
          tabPanel("GO BP", br(), plotOutput("gene_bp_plot", height = "520px"), br(), DTOutput("gene_bp_tbl")),
          tabPanel("GO MF", br(), plotOutput("gene_mf_plot", height = "520px"), br(), DTOutput("gene_mf_tbl")),
          tabPanel("GO CC", br(), plotOutput("gene_cc_plot", height = "520px"), br(), DTOutput("gene_cc_tbl")),
          tabPanel("KEGG", br(), plotOutput("gene_kegg_plot", height = "520px"), br(), DTOutput("gene_kegg_tbl"))
        )
      )
    )
  ),
  tabPanel(
    "Raw-data mode",
    sidebarLayout(
      sidebarPanel(
        fileInput("file_raw", "上传 raw abundance 表", accept = c(".xlsx", ".xls", ".csv", ".tsv", ".txt")),
        textInput("desc_col_raw", "Description 列名（可留空自动识别）", value = ""),
        tags$hr(),
        tags$strong("分析模式"),
        radioButtons(
          "analysis_mode_raw", "分析模式",
          choices = c("模式 A：宽松" = "A", "模式 B：平衡" = "B", "模式 C：严格" = "C"),
          selected = "A"
        ),
        selectInput(
          "confidence_filter_raw", "Protein FDR Confidence: Combined",
          choices = c("High only" = "high", "High + Medium" = "high_medium", "All" = "all"),
          selected = "high"
        ),
        selectInput(
          "contaminant_filter_raw", "Contaminant filter",
          choices = c("Common contaminants（293T 默认）" = "common", "Keratin only" = "keratin", "不过滤" = "none"),
          selected = "common"
        ),
        checkboxInput("precompute_all_modes_raw", "预计算全部模式（用于跨模式对比）", value = FALSE),
        checkboxInput("deduplicate_by_gene_raw", "按 GeneName 去重（保留整体丰度最高的一行）", value = FALSE),
        tags$hr(),
        tags$strong("差异分析阈值"),
        numericInput("fc_cut_raw", "Fold change cutoff", value = 2, min = 1, step = 0.1),
        numericInput("p_cut_raw", "Differential pvalue cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
        tags$hr(),
        tags$strong("重复/异常值规则"),
        numericInput("outlier_dev_log2_raw", "单点离群偏离阈值（log2）", value = 1, min = 0, step = 0.1),
        numericInput("remaining_range_log2_raw", "剩余重复接近阈值（log2 range）", value = 0.5, min = 0, step = 0.1),
        numericInput("test_pair_fc_raw", "两重复时 Test 组最大允许倍数差", value = 5, min = 1, step = 0.5),
        numericInput("test_inconsistent_fc_raw", "三重复及以上 Test 组最大允许倍数差", value = 5, min = 1, step = 0.5),
        tags$hr(),
        tags$strong("缺失值可视化处理"),
        checkboxInput("impute_visual_raw", "仅用于 volcano 的左尾随机填补", value = TRUE),
        numericInput("impute_width_raw", "Left-tail width", value = 0.3, min = 0.01, step = 0.05),
        numericInput("impute_downshift_raw", "Left-tail downshift", value = 1.8, min = 0.1, step = 0.1),
        numericInput("impute_seed_raw", "Random seed", value = 123, min = 1, step = 1),
        selectInput(
          "presence_display_raw", "Candidate-aware 中 presence 点显示",
          choices = c("隐藏" = "hide", "只统计数量" = "count", "使用固定最小值绘制（推荐）" = "plot_fixed", "使用左尾随机值绘制" = "plot_random"),
          selected = "plot_fixed"
        ),
        checkboxInput("label_top_genes_raw", "标注 up/down 前15个基因名", value = FALSE),
        numericInput("label_top_n_raw", "每个方向标注前 N 个基因", value = 15, min = 1, max = 100, step = 1),
        tags$hr(),
        tags$strong("Volcano plot settings"),
        selectInput(
          "volcano_view_raw", "Volcano 视图",
          choices = c("当前模式" = "current", "候选视图（模式 A）" = "A", "严格视图（模式 C）" = "C"),
          selected = "current"
        ),
        selectInput(
          "volcano_style_raw", "Volcano display style",
          choices = c("Candidate-aware" = "candidate", "Exploratory style" = "exploratory"),
          selected = "candidate"
        ),
        numericInput("vol_p_cut_raw", "Volcano 横线 pvalue cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
        textInput("col_up_raw", "上调点颜色", value = "#F39B7F"),
        textInput("col_down_raw", "下调点颜色", value = "#4DBBD5"),
        textInput("col_ns_raw", "NS 点颜色", value = "#A6A6A6"),
        textInput("col_line_raw", "阈值线颜色", value = "#000000"),
        numericInput("axis_text_size_raw", "坐标轴字号", value = 18, min = 8, step = 1),
        numericInput("axis_title_size_raw", "坐标标题字号", value = 24, min = 8, step = 1),
        numericInput("legend_text_size_raw", "图例字号", value = 14, min = 8, step = 1),
        numericInput("pt_alpha_raw", "点透明度", value = 0.7, min = 0, max = 1, step = 0.05),
        numericInput("pt_size_raw", "点大小", value = 3, min = 0.5, step = 0.5),
        tags$hr(),
        tags$strong("Enrichment settings"),
        selectInput(
          "enrich_strategy_raw", "Enrichment 策略",
          choices = c(
            "Exploratory enrichment（更适合候选解释）" = "exploratory",
            "Context-aware enrichment（更强调实验背景）" = "context"
          ),
          selected = "exploratory"
        ),
        selectInput(
          "enrich_scope_raw", "Enrichment 方向",
          choices = c("全部（Up + Down）" = "both", "只算 Up" = "up", "只算 Down" = "down"),
          selected = "both"
        ),
        numericInput("enrich_p_cut_raw", "Enrichment pvalue cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput("enrich_q_cut_raw", "GO qvalue cutoff", value = 0.2, min = 0, max = 1, step = 0.01),
        numericInput("show_n_raw", "dotplot 显示条目数", value = 12, min = 5, max = 50, step = 1),
        selectInput(
          "plot_format_raw", "图形导出格式（PDF/EMF 更适合后续编辑）",
          choices = c("PNG" = "png", "JPG" = "jpg", "EMF" = "emf", "PDF" = "pdf"),
          selected = "png"
        ),
        tags$hr(),
        tags$strong("离子注释"),
        checkboxGroupInput(
          "ions_raw", "选择离子",
          choices = c("钙 calcium" = "calcium", "镁 magnesium" = "magnesium", "锰 manganese" = "manganese", "锌 zinc" = "zinc"),
          selected = c("calcium", "magnesium", "manganese", "zinc")
        ),
        actionButton("run_raw", "运行 Raw-data 分析", class = "btn-danger"),
        br(), br(),
        helpText("Raw-data mode 要求丰度列命名为 Test1...TestN 和 Ctrl1...CtrlN。模式 A 会先把空白补 0，再按宽松规则做候选筛选；其中 Test 组的均值只用有效正值计算，Ctrl 组继续按已确认口径把 0 计入均值。Candidate-aware 是当前分析框架下的火山图；其中 presence 点可用固定最小值（推荐）或左尾随机值进行可视化。Exploratory style 只改变坐标画法，不改变 up/down/presence 的归类，并且当前仅支持 2×2 比较。默认会启用 293T 常见污染物过滤。模式 A 只清洗 Test 组，不清洗 Ctrl 组。Exploratory enrichment 更适合候选解释；Context-aware enrichment 会引入本次分析的实验背景。离子注释会跟随同一个 Up/Down 选择范围。默认只计算当前模式；若要跨模式切换火山图，请勾选“预计算全部模式”。"),
        verbatimTextOutput("status_raw"),
        hr(),
        downloadButton("download_raw_diff", "下载差异结果总表"),
        br(), br(),
        downloadButton("download_raw_filtered", "下载 filtered out 清单"),
        br(), br(),
        downloadButton("download_raw_contaminants", "下载 contaminants removed 清单"),
        br(), br(),
        downloadButton("download_raw_duplicates", "下载 duplicate 处理记录"),
        br(), br(),
        downloadButton("download_raw_annotated", "下载 annotated 差异结果"),
        br(), br(),
        downloadButton("download_raw_up", "下载当前模式 Up 列表"),
        br(), br(),
        downloadButton("download_raw_down", "下载当前模式 Down 列表"),
        br(), br(),
        downloadButton("download_raw_volcano", "下载火山图"),
        br(), br(),
        downloadButton("download_raw_up_bp_plot", "下载 Up BP 图"),
        br(), br(),
        downloadButton("download_raw_up_mf_plot", "下载 Up MF 图"),
        br(), br(),
        downloadButton("download_raw_up_cc_plot", "下载 Up CC 图"),
        br(), br(),
        downloadButton("download_raw_up_kegg_plot", "下载 Up KEGG 图"),
        br(), br(),
        downloadButton("download_raw_down_bp_plot", "下载 Down BP 图"),
        br(), br(),
        downloadButton("download_raw_down_mf_plot", "下载 Down MF 图"),
        br(), br(),
        downloadButton("download_raw_down_cc_plot", "下载 Down CC 图"),
        br(), br(),
        downloadButton("download_raw_down_kegg_plot", "下载 Down KEGG 图"),
        br(), br(),
        downloadButton("download_raw_go_xlsx", "下载 GO enrichment.xlsx"),
        br(), br(),
        downloadButton("download_raw_kegg_xlsx", "下载 KEGG enrichment.xlsx"),
        br(), br(),
        downloadButton("download_raw_ion_xlsx", "下载 Ion annotation.xlsx")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("概览", br(), DTOutput("raw_summary_tbl")),
          tabPanel("输入预览", br(), DTOutput("raw_input_tbl")),
          tabPanel("Duplicate 处理", br(), DTOutput("raw_dup_tbl")),
          tabPanel("Contaminants removed", br(), DTOutput("raw_contaminants_tbl")),
          tabPanel("差异结果总表", br(), DTOutput("raw_diff_tbl")),
          tabPanel("当前模式 Up", br(), DTOutput("raw_up_tbl")),
          tabPanel("当前模式 Down", br(), DTOutput("raw_down_tbl")),
          tabPanel("Filtered out", br(), DTOutput("raw_filtered_tbl")),
          tabPanel("Annotated 差异结果", br(), DTOutput("raw_annotated_tbl")),
          tabPanel(
            "Volcano plot",
            br(),
            plotOutput("raw_volcano_plot", height = "620px")
          ),
          tabPanel(
            "离子注释",
            br(),
            tabsetPanel(
              tabPanel("钙 calcium",
                       h4("GO term summary"), DTOutput("raw_ion_terms_calcium_tbl"), br(),
                       h4("Protein summary"), DTOutput("raw_ion_summary_calcium_tbl"), br(),
                       h4("GO hits"), DTOutput("raw_ion_hits_calcium_tbl")),
              tabPanel("镁 magnesium",
                       h4("GO term summary"), DTOutput("raw_ion_terms_magnesium_tbl"), br(),
                       h4("Protein summary"), DTOutput("raw_ion_summary_magnesium_tbl"), br(),
                       h4("GO hits"), DTOutput("raw_ion_hits_magnesium_tbl")),
              tabPanel("锰 manganese",
                       h4("GO term summary"), DTOutput("raw_ion_terms_manganese_tbl"), br(),
                       h4("Protein summary"), DTOutput("raw_ion_summary_manganese_tbl"), br(),
                       h4("GO hits"), DTOutput("raw_ion_hits_manganese_tbl")),
              tabPanel("锌 zinc",
                       h4("GO term summary"), DTOutput("raw_ion_terms_zinc_tbl"), br(),
                       h4("Protein summary"), DTOutput("raw_ion_summary_zinc_tbl"), br(),
                       h4("GO hits"), DTOutput("raw_ion_hits_zinc_tbl"))
            )
          ),
          tabPanel(
            "Enrichment: Up",
            br(),
            tabsetPanel(
              tabPanel("BP", br(), plotOutput("raw_up_bp_plot", height = "500px"), br(), DTOutput("raw_up_bp_tbl")),
              tabPanel("MF", br(), plotOutput("raw_up_mf_plot", height = "500px"), br(), DTOutput("raw_up_mf_tbl")),
              tabPanel("CC", br(), plotOutput("raw_up_cc_plot", height = "500px"), br(), DTOutput("raw_up_cc_tbl")),
              tabPanel("KEGG", br(), plotOutput("raw_up_kegg_plot", height = "500px"), br(), DTOutput("raw_up_kegg_tbl"))
            )
          ),
          tabPanel(
            "Enrichment: Down",
            br(),
            tabsetPanel(
              tabPanel("BP", br(), plotOutput("raw_down_bp_plot", height = "500px"), br(), DTOutput("raw_down_bp_tbl")),
              tabPanel("MF", br(), plotOutput("raw_down_mf_plot", height = "500px"), br(), DTOutput("raw_down_mf_tbl")),
              tabPanel("CC", br(), plotOutput("raw_down_cc_plot", height = "500px"), br(), DTOutput("raw_down_cc_tbl")),
              tabPanel("KEGG", br(), plotOutput("raw_down_kegg_plot", height = "500px"), br(), DTOutput("raw_down_kegg_tbl"))
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  dt_opts <- list(pageLength = 15, scrollX = TRUE)

  gene_result <- eventReactive(input$run_gene, {
    req(input$file_gene)
    run_gene_list_mode(
      input_path = input$file_gene$datapath,
      desc_col = if (nzchar(input$desc_col_gene)) input$desc_col_gene else NULL,
      target_ions = input$ions_gene,
      p_cut = input$p_cut_gene,
      q_cut = input$q_cut_gene,
      show_n = input$show_n_gene
    )
  }, ignoreNULL = TRUE)

  output$status_gene <- renderText({
    req(input$file_gene)
    if (input$run_gene < 1) {
      paste0("已选择文件：", input$file_gene$name, "\n等待运行。")
    } else {
      res <- gene_result()
      msg <- c(
        paste0("文件：", input$file_gene$name),
        paste0("离子：", paste(res$selected_ions, collapse = ", ")),
        paste0("提取到的唯一 GeneName 数：", res$gene_count),
        paste0("可映射 ENTREZID 数：", res$entrez_count),
        paste0("BP term 数：", nrow(res$enrich$bp$df)),
        paste0("MF term 数：", nrow(res$enrich$mf$df)),
        paste0("CC term 数：", nrow(res$enrich$cc$df)),
        paste0("KEGG term 数：", nrow(res$enrich$kegg$df)),
        "分析完成。"
      )
      if (!is.null(res$enrich$kegg$error) && nzchar(res$enrich$kegg$error)) {
        msg <- c(msg, paste0("KEGG 提示：", res$enrich$kegg$error))
      }
      paste(msg, collapse = "\n")
    }
  })

  output$gene_overview_tbl <- renderDT({ req(gene_result()); datatable(gene_result()$overview, options = dt_opts, filter = "top") })
  output$gene_input_tbl <- renderDT({ req(gene_result()); datatable(gene_result()$input_preview, options = dt_opts, filter = "top") })
  output$gene_summary_tbl <- renderDT({ req(gene_result()); datatable(gene_result()$protein_summary, options = dt_opts, filter = "top") })
  output$gene_hits_tbl <- renderDT({ req(gene_result()); datatable(gene_result()$go_hits, options = dt_opts, filter = "top") })
  output$gene_annotated_tbl <- renderDT({ req(gene_result()); datatable(gene_result()$annotated_original, options = dt_opts, filter = "top") })

  output$gene_bp_plot <- renderPlot({ req(gene_result()); safe_dotplot(gene_result()$enrich$bp$obj, "GO BP", gene_result()$enrich$show_n) })
  output$gene_mf_plot <- renderPlot({ req(gene_result()); safe_dotplot(gene_result()$enrich$mf$obj, "GO MF", gene_result()$enrich$show_n) })
  output$gene_cc_plot <- renderPlot({ req(gene_result()); safe_dotplot(gene_result()$enrich$cc$obj, "GO CC", gene_result()$enrich$show_n) })
  output$gene_kegg_plot <- renderPlot({ req(gene_result()); safe_dotplot(gene_result()$enrich$kegg$obj, "KEGG", gene_result()$enrich$show_n) })

  output$gene_bp_tbl <- renderDT({ req(gene_result()); datatable(gene_result()$enrich$bp$df, options = dt_opts, filter = "top") })
  output$gene_mf_tbl <- renderDT({ req(gene_result()); datatable(gene_result()$enrich$mf$df, options = dt_opts, filter = "top") })
  output$gene_cc_tbl <- renderDT({ req(gene_result()); datatable(gene_result()$enrich$cc$df, options = dt_opts, filter = "top") })
  output$gene_kegg_tbl <- renderDT({ req(gene_result()); datatable(gene_result()$enrich$kegg$df, options = dt_opts, filter = "top") })

  output$download_gene_annotated <- downloadHandler(
    filename = function() "gene_list_annotated.csv",
    content = function(file) write.csv(gene_result()$annotated_original, file, row.names = FALSE)
  )
  output$download_gene_summary <- downloadHandler(
    filename = function() "gene_list_ion_protein_summary.csv",
    content = function(file) write.csv(gene_result()$protein_summary, file, row.names = FALSE)
  )
  output$download_gene_hits <- downloadHandler(
    filename = function() "gene_list_ion_go_hits.csv",
    content = function(file) write.csv(gene_result()$go_hits, file, row.names = FALSE)
  )
  output$download_gene_bp <- downloadHandler(
    filename = function() "gene_list_GO_BP.csv",
    content = function(file) write.csv(gene_result()$enrich$bp$df, file, row.names = FALSE)
  )
  output$download_gene_mf <- downloadHandler(
    filename = function() "gene_list_GO_MF.csv",
    content = function(file) write.csv(gene_result()$enrich$mf$df, file, row.names = FALSE)
  )
  output$download_gene_cc <- downloadHandler(
    filename = function() "gene_list_GO_CC.csv",
    content = function(file) write.csv(gene_result()$enrich$cc$df, file, row.names = FALSE)
  )
  output$download_gene_kegg <- downloadHandler(
    filename = function() "gene_list_KEGG.csv",
    content = function(file) write.csv(gene_result()$enrich$kegg$df, file, row.names = FALSE)
  )
  output$download_gene_bp_plot <- downloadHandler(
    filename = function() {
      paste0("gene_list_GO_BP.", input$plot_format_gene)
    },
    content = function(file) {
      req(gene_result())
      plot_obj <- safe_dotplot(gene_result()$enrich$bp$obj, "GO BP", gene_result()$enrich$show_n)
      save_plot_file(plot_obj, file, format = input$plot_format_gene, width = 8, height = 6)
    }
  )
  output$download_gene_mf_plot <- downloadHandler(
    filename = function() {
      paste0("gene_list_GO_MF.", input$plot_format_gene)
    },
    content = function(file) {
      req(gene_result())
      plot_obj <- safe_dotplot(gene_result()$enrich$mf$obj, "GO MF", gene_result()$enrich$show_n)
      save_plot_file(plot_obj, file, format = input$plot_format_gene, width = 8, height = 6)
    }
  )
  output$download_gene_cc_plot <- downloadHandler(
    filename = function() {
      paste0("gene_list_GO_CC.", input$plot_format_gene)
    },
    content = function(file) {
      req(gene_result())
      plot_obj <- safe_dotplot(gene_result()$enrich$cc$obj, "GO CC", gene_result()$enrich$show_n)
      save_plot_file(plot_obj, file, format = input$plot_format_gene, width = 8, height = 6)
    }
  )
  output$download_gene_kegg_plot <- downloadHandler(
    filename = function() {
      paste0("gene_list_KEGG.", input$plot_format_gene)
    },
    content = function(file) {
      req(gene_result())
      plot_obj <- safe_dotplot(gene_result()$enrich$kegg$obj, "KEGG", gene_result()$enrich$show_n)
      save_plot_file(plot_obj, file, format = input$plot_format_gene, width = 8, height = 6)
    }
  )

  raw_result <- eventReactive(input$run_raw, {
    req(input$file_raw)
    run_raw_mode(
      input_path = input$file_raw$datapath,
      desc_col = if (nzchar(input$desc_col_raw)) input$desc_col_raw else NULL,
      target_ions = input$ions_raw,
      analysis_mode = input$analysis_mode_raw,
      precompute_all_modes = isTRUE(input$precompute_all_modes_raw),
      confidence_filter = input$confidence_filter_raw,
      contaminant_filter = input$contaminant_filter_raw,
      deduplicate_by_gene = isTRUE(input$deduplicate_by_gene_raw),
      outlier_dev_log2 = input$outlier_dev_log2_raw,
      remaining_range_log2 = input$remaining_range_log2_raw,
      test_pair_fc = input$test_pair_fc_raw,
      test_inconsistent_fc = input$test_inconsistent_fc_raw,
      fc_cut = input$fc_cut_raw,
      p_cut = input$p_cut_raw,
      enrich_p_cut = input$enrich_p_cut_raw,
      enrich_q_cut = input$enrich_q_cut_raw,
      enrich_strategy = input$enrich_strategy_raw,
      enrich_scope = input$enrich_scope_raw,
      show_n = input$show_n_raw,
      impute_for_visual = isTRUE(input$impute_visual_raw),
      impute_width = input$impute_width_raw,
      impute_downshift = input$impute_downshift_raw,
      impute_seed = input$impute_seed_raw
    )
  }, ignoreNULL = TRUE)

  output$status_raw <- renderText({
    req(input$file_raw)
    if (input$run_raw < 1) {
      paste0("已选择文件：", input$file_raw$name, "\n等待运行。")
    } else {
      res <- raw_result()
      msg <- c(
        paste0("文件：", input$file_raw$name),
        paste0("当前模式：", dplyr::case_when(
          res$analysis_mode == "A" ~ "模式 A：宽松",
          res$analysis_mode == "B" ~ "模式 B：平衡",
          TRUE ~ "模式 C：严格"
        )),
        paste0("Protein FDR Confidence 过滤：", res$confidence_filter),
        paste0("Contaminant 过滤：", res$contaminant_filter, "（移除 ", nrow(res$contaminants_removed_df), " 条）"),
        paste0("预计算全部模式：", if (isTRUE(res$precompute_all_modes)) "是" else "否"),
        paste0("检测到 Test 列：", paste(res$test_cols, collapse = ", ")),
        paste0("检测到 Ctrl 列：", paste(res$ctrl_cols, collapse = ", ")),
        paste0(if (isTRUE(res$deduplicate_by_gene)) "按 GeneName 去重后条目数：" else "参与分析条目数（未按 GeneName 去重）：", nrow(res$dedup_df)),
        paste0("模式 A up/down：", sum(res$diff_df$mode_a_status == "up", na.rm = TRUE), " / ", sum(res$diff_df$mode_a_status == "down", na.rm = TRUE)),
        paste0("模式 B up/down：", sum(res$diff_df$mode_b_status == "up", na.rm = TRUE), " / ", sum(res$diff_df$mode_b_status == "down", na.rm = TRUE)),
        paste0("模式 C up/down：", sum(res$diff_df$mode_c_status == "up", na.rm = TRUE), " / ", sum(res$diff_df$mode_c_status == "down", na.rm = TRUE)),
        paste0("当前模式 Presence up/down：", sum(res$diff_df$status == "presence_up", na.rm = TRUE), " / ", sum(res$diff_df$status == "presence_down", na.rm = TRUE)),
        paste0("当前模式 Up/Down：", nrow(res$active_up_df), " / ", nrow(res$active_down_df)),
        paste0("Filtered out：", sum(res$diff_df$status == "filtered_out", na.rm = TRUE)),
        paste0("用于 enrichment 的 up/down 基因数：", length(res$up_genes), " / ", length(res$down_genes)),
        paste0("Exploratory style 可用：", if (isTRUE(res$exploratory_available)) "是（2×2）" else "否（当前不是 2×2）"),
        "分析完成。"
      )
      if (!is.null(res$enrich_up$kegg$error) && nzchar(res$enrich_up$kegg$error)) {
        msg <- c(msg, paste0("Up KEGG 提示：", res$enrich_up$kegg$error))
      }
      if (!is.null(res$enrich_down$kegg$error) && nzchar(res$enrich_down$kegg$error)) {
        msg <- c(msg, paste0("Down KEGG 提示：", res$enrich_down$kegg$error))
      }
      paste(msg, collapse = "\n")
    }
  })

  output$raw_summary_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$summary_df, options = dt_opts, filter = "top") })
  output$raw_input_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$input_preview, options = dt_opts, filter = "top") })
  output$raw_dup_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$duplicate_summary, options = dt_opts, filter = "top") })
  output$raw_contaminants_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$contaminants_removed_df, options = dt_opts, filter = "top") })
  output$raw_diff_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$diff_df, options = dt_opts, filter = "top") })
  output$raw_up_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$active_up_df, options = dt_opts, filter = "top") })
  output$raw_down_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$active_down_df, options = dt_opts, filter = "top") })
  output$raw_filtered_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$filtered_out_df, options = dt_opts, filter = "top") })
  output$raw_annotated_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$annotated_diff_df, options = dt_opts, filter = "top") })

  get_ion_subset <- function(df, ion) {
    req(raw_result())
    if (is.null(df) || nrow(df) == 0 || !"ion" %in% colnames(df)) return(df[0, , drop = FALSE])
    out <- df[df$ion == ion, , drop = FALSE]
    out
  }

  output$raw_ion_terms_calcium_tbl <- renderDT({ req(raw_result()); datatable(get_ion_subset(raw_result()$ion_term_summary, "calcium"), options = dt_opts, filter = "top") })
  output$raw_ion_terms_magnesium_tbl <- renderDT({ req(raw_result()); datatable(get_ion_subset(raw_result()$ion_term_summary, "magnesium"), options = dt_opts, filter = "top") })
  output$raw_ion_terms_manganese_tbl <- renderDT({ req(raw_result()); datatable(get_ion_subset(raw_result()$ion_term_summary, "manganese"), options = dt_opts, filter = "top") })
  output$raw_ion_terms_zinc_tbl <- renderDT({ req(raw_result()); datatable(get_ion_subset(raw_result()$ion_term_summary, "zinc"), options = dt_opts, filter = "top") })

  output$raw_ion_summary_calcium_tbl <- renderDT({ req(raw_result()); datatable(get_ion_subset(raw_result()$ion_protein_summary, "calcium"), options = dt_opts, filter = "top") })
  output$raw_ion_summary_magnesium_tbl <- renderDT({ req(raw_result()); datatable(get_ion_subset(raw_result()$ion_protein_summary, "magnesium"), options = dt_opts, filter = "top") })
  output$raw_ion_summary_manganese_tbl <- renderDT({ req(raw_result()); datatable(get_ion_subset(raw_result()$ion_protein_summary, "manganese"), options = dt_opts, filter = "top") })
  output$raw_ion_summary_zinc_tbl <- renderDT({ req(raw_result()); datatable(get_ion_subset(raw_result()$ion_protein_summary, "zinc"), options = dt_opts, filter = "top") })

  output$raw_ion_hits_calcium_tbl <- renderDT({ req(raw_result()); datatable(get_ion_subset(raw_result()$ion_go_hits, "calcium"), options = dt_opts, filter = "top") })
  output$raw_ion_hits_magnesium_tbl <- renderDT({ req(raw_result()); datatable(get_ion_subset(raw_result()$ion_go_hits, "magnesium"), options = dt_opts, filter = "top") })
  output$raw_ion_hits_manganese_tbl <- renderDT({ req(raw_result()); datatable(get_ion_subset(raw_result()$ion_go_hits, "manganese"), options = dt_opts, filter = "top") })
  output$raw_ion_hits_zinc_tbl <- renderDT({ req(raw_result()); datatable(get_ion_subset(raw_result()$ion_go_hits, "zinc"), options = dt_opts, filter = "top") })

  output$raw_volcano_plot <- renderPlot({
    req(raw_result())
    mode_to_show <- if (identical(input$volcano_view_raw, "current")) raw_result()$analysis_mode else input$volcano_view_raw
    if (!isTRUE(raw_result()$precompute_all_modes) && !identical(input$volcano_view_raw, "current") && !identical(mode_to_show, raw_result()$analysis_mode)) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "该视图需要勾选“预计算全部模式”", size = 6) +
          theme_void()
      )
    }
    if (identical(input$volcano_style_raw, "exploratory") && !isTRUE(raw_result()$exploratory_available)) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "Exploratory style 当前仅支持 2×2 比较", size = 6) +
          theme_void()
      )
    }
    plot_df <- prepare_mode_view_df(raw_result()$annotated_diff_df, mode_to_show)
    make_volcano_plot(
      plot_df,
      fc_cut = input$fc_cut_raw,
      p_cut = input$vol_p_cut_raw,
      col_up = input$col_up_raw,
      col_down = input$col_down_raw,
      col_ns = input$col_ns_raw,
      col_line = input$col_line_raw,
      axis_text_size = input$axis_text_size_raw,
      axis_title_size = input$axis_title_size_raw,
      legend_text_size = input$legend_text_size_raw,
      pt_alpha = input$pt_alpha_raw,
      pt_size = input$pt_size_raw,
      presence_display = input$presence_display_raw,
      use_visual_imputation = isTRUE(input$impute_visual_raw),
      active_mode_label = mode_to_show,
      label_top_genes = isTRUE(input$label_top_genes_raw),
      label_top_n = input$label_top_n_raw,
      display_style = input$volcano_style_raw
    )
  })

  output$raw_up_bp_plot <- renderPlot({ req(raw_result()); safe_dotplot(raw_result()$enrich_up$bp$obj, "Up: GO BP", raw_result()$enrich_up$show_n) })
  output$raw_up_mf_plot <- renderPlot({ req(raw_result()); safe_dotplot(raw_result()$enrich_up$mf$obj, "Up: GO MF", raw_result()$enrich_up$show_n) })
  output$raw_up_cc_plot <- renderPlot({ req(raw_result()); safe_dotplot(raw_result()$enrich_up$cc$obj, "Up: GO CC", raw_result()$enrich_up$show_n) })
  output$raw_up_kegg_plot <- renderPlot({ req(raw_result()); safe_dotplot(raw_result()$enrich_up$kegg$obj, "Up: KEGG", raw_result()$enrich_up$show_n) })

  output$raw_down_bp_plot <- renderPlot({ req(raw_result()); safe_dotplot(raw_result()$enrich_down$bp$obj, "Down: GO BP", raw_result()$enrich_down$show_n) })
  output$raw_down_mf_plot <- renderPlot({ req(raw_result()); safe_dotplot(raw_result()$enrich_down$mf$obj, "Down: GO MF", raw_result()$enrich_down$show_n) })
  output$raw_down_cc_plot <- renderPlot({ req(raw_result()); safe_dotplot(raw_result()$enrich_down$cc$obj, "Down: GO CC", raw_result()$enrich_down$show_n) })
  output$raw_down_kegg_plot <- renderPlot({ req(raw_result()); safe_dotplot(raw_result()$enrich_down$kegg$obj, "Down: KEGG", raw_result()$enrich_down$show_n) })

  output$raw_up_bp_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$enrich_up$bp$df, options = dt_opts, filter = "top") })
  output$raw_up_mf_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$enrich_up$mf$df, options = dt_opts, filter = "top") })
  output$raw_up_cc_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$enrich_up$cc$df, options = dt_opts, filter = "top") })
  output$raw_up_kegg_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$enrich_up$kegg$df, options = dt_opts, filter = "top") })

  output$raw_down_bp_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$enrich_down$bp$df, options = dt_opts, filter = "top") })
  output$raw_down_mf_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$enrich_down$mf$df, options = dt_opts, filter = "top") })
  output$raw_down_cc_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$enrich_down$cc$df, options = dt_opts, filter = "top") })
  output$raw_down_kegg_tbl <- renderDT({ req(raw_result()); datatable(raw_result()$enrich_down$kegg$df, options = dt_opts, filter = "top") })

  output$download_raw_diff <- downloadHandler(
    filename = function() "raw_differential_results_all.csv",
    content = function(file) write.csv(raw_result()$diff_df, file, row.names = FALSE)
  )
  output$download_raw_filtered <- downloadHandler(
    filename = function() "raw_filtered_out_genes.csv",
    content = function(file) write.csv(raw_result()$filtered_out_df, file, row.names = FALSE)
  )
  output$download_raw_contaminants <- downloadHandler(
    filename = function() "raw_contaminants_removed.csv",
    content = function(file) write.csv(raw_result()$contaminants_removed_df, file, row.names = FALSE)
  )
  output$download_raw_duplicates <- downloadHandler(
    filename = function() "raw_duplicate_resolution.csv",
    content = function(file) write.csv(raw_result()$duplicate_summary, file, row.names = FALSE)
  )
  output$download_raw_annotated <- downloadHandler(
    filename = function() "raw_annotated_differential_results.csv",
    content = function(file) write.csv(raw_result()$annotated_diff_df, file, row.names = FALSE)
  )
  output$download_raw_up <- downloadHandler(
    filename = function() paste0("raw_mode_", raw_result()$analysis_mode, "_up_list.csv"),
    content = function(file) write.csv(raw_result()$active_up_df, file, row.names = FALSE)
  )
  output$download_raw_down <- downloadHandler(
    filename = function() paste0("raw_mode_", raw_result()$analysis_mode, "_down_list.csv"),
    content = function(file) write.csv(raw_result()$active_down_df, file, row.names = FALSE)
  )
  output$download_raw_volcano <- downloadHandler(
    filename = function() {
      paste0(
        "volcano_mode_",
        if (identical(input$volcano_view_raw, "current")) raw_result()$analysis_mode else input$volcano_view_raw,
        "_",
        input$volcano_style_raw,
        ".",
        input$plot_format_raw
      )
    },
    content = function(file) {
      req(raw_result())
      mode_to_show <- if (identical(input$volcano_view_raw, "current")) raw_result()$analysis_mode else input$volcano_view_raw
      if (identical(input$volcano_style_raw, "exploratory") && !isTRUE(raw_result()$exploratory_available)) {
        stop("Exploratory style 当前仅支持 2×2 比较")
      }
      plot_df <- prepare_mode_view_df(raw_result()$annotated_diff_df, mode_to_show)
      plot_obj <- make_volcano_plot(
        plot_df,
        fc_cut = input$fc_cut_raw,
        p_cut = input$vol_p_cut_raw,
        col_up = input$col_up_raw,
        col_down = input$col_down_raw,
        col_ns = input$col_ns_raw,
        col_line = input$col_line_raw,
        axis_text_size = input$axis_text_size_raw,
        axis_title_size = input$axis_title_size_raw,
        legend_text_size = input$legend_text_size_raw,
        pt_alpha = input$pt_alpha_raw,
        pt_size = input$pt_size_raw,
        presence_display = input$presence_display_raw,
        use_visual_imputation = isTRUE(input$impute_visual_raw),
        active_mode_label = mode_to_show,
        label_top_genes = isTRUE(input$label_top_genes_raw),
        label_top_n = input$label_top_n_raw,
        display_style = input$volcano_style_raw
      )
      save_plot_file(plot_obj, file, format = input$plot_format_raw, width = 10, height = 8)
    }
  )
  output$download_raw_up_bp_plot <- downloadHandler(
    filename = function() {
      paste0("raw_up_GO_BP.", input$plot_format_raw)
    },
    content = function(file) {
      req(raw_result())
      plot_obj <- safe_dotplot(raw_result()$enrich_up$bp$obj, "Up: GO BP", raw_result()$enrich_up$show_n)
      save_plot_file(plot_obj, file, format = input$plot_format_raw, width = 8, height = 6)
    }
  )
  output$download_raw_up_mf_plot <- downloadHandler(
    filename = function() {
      paste0("raw_up_GO_MF.", input$plot_format_raw)
    },
    content = function(file) {
      req(raw_result())
      plot_obj <- safe_dotplot(raw_result()$enrich_up$mf$obj, "Up: GO MF", raw_result()$enrich_up$show_n)
      save_plot_file(plot_obj, file, format = input$plot_format_raw, width = 8, height = 6)
    }
  )
  output$download_raw_up_cc_plot <- downloadHandler(
    filename = function() {
      paste0("raw_up_GO_CC.", input$plot_format_raw)
    },
    content = function(file) {
      req(raw_result())
      plot_obj <- safe_dotplot(raw_result()$enrich_up$cc$obj, "Up: GO CC", raw_result()$enrich_up$show_n)
      save_plot_file(plot_obj, file, format = input$plot_format_raw, width = 8, height = 6)
    }
  )
  output$download_raw_up_kegg_plot <- downloadHandler(
    filename = function() {
      paste0("raw_up_KEGG.", input$plot_format_raw)
    },
    content = function(file) {
      req(raw_result())
      plot_obj <- safe_dotplot(raw_result()$enrich_up$kegg$obj, "Up: KEGG", raw_result()$enrich_up$show_n)
      save_plot_file(plot_obj, file, format = input$plot_format_raw, width = 8, height = 6)
    }
  )
  output$download_raw_down_bp_plot <- downloadHandler(
    filename = function() {
      paste0("raw_down_GO_BP.", input$plot_format_raw)
    },
    content = function(file) {
      req(raw_result())
      plot_obj <- safe_dotplot(raw_result()$enrich_down$bp$obj, "Down: GO BP", raw_result()$enrich_down$show_n)
      save_plot_file(plot_obj, file, format = input$plot_format_raw, width = 8, height = 6)
    }
  )
  output$download_raw_down_mf_plot <- downloadHandler(
    filename = function() {
      paste0("raw_down_GO_MF.", input$plot_format_raw)
    },
    content = function(file) {
      req(raw_result())
      plot_obj <- safe_dotplot(raw_result()$enrich_down$mf$obj, "Down: GO MF", raw_result()$enrich_down$show_n)
      save_plot_file(plot_obj, file, format = input$plot_format_raw, width = 8, height = 6)
    }
  )
  output$download_raw_down_cc_plot <- downloadHandler(
    filename = function() {
      paste0("raw_down_GO_CC.", input$plot_format_raw)
    },
    content = function(file) {
      req(raw_result())
      plot_obj <- safe_dotplot(raw_result()$enrich_down$cc$obj, "Down: GO CC", raw_result()$enrich_down$show_n)
      save_plot_file(plot_obj, file, format = input$plot_format_raw, width = 8, height = 6)
    }
  )
  output$download_raw_down_kegg_plot <- downloadHandler(
    filename = function() {
      paste0("raw_down_KEGG.", input$plot_format_raw)
    },
    content = function(file) {
      req(raw_result())
      plot_obj <- safe_dotplot(raw_result()$enrich_down$kegg$obj, "Down: KEGG", raw_result()$enrich_down$show_n)
      save_plot_file(plot_obj, file, format = input$plot_format_raw, width = 8, height = 6)
    }
  )

  output$download_raw_go_xlsx <- downloadHandler(
    filename = function() paste0("GO_enrichment_mode", raw_result()$analysis_mode, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content = function(file) {
      req(raw_result())
      wb <- openxlsx::createWorkbook()
      go_sheets <- list(
        BP = combine_up_down_df(raw_result()$enrich_up$bp$df, raw_result()$enrich_down$bp$df),
        CC = combine_up_down_df(raw_result()$enrich_up$cc$df, raw_result()$enrich_down$cc$df),
        MF = combine_up_down_df(raw_result()$enrich_up$mf$df, raw_result()$enrich_down$mf$df)
      )
      for (nm in names(go_sheets)) {
        openxlsx::addWorksheet(wb, nm)
        openxlsx::writeData(wb, nm, go_sheets[[nm]], withFilter = TRUE)
      }
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )

  output$download_raw_kegg_xlsx <- downloadHandler(
    filename = function() paste0("KEGG_enrichment_mode", raw_result()$analysis_mode, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content = function(file) {
      req(raw_result())
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "KEGG")
      kegg_df <- combine_up_down_df(raw_result()$enrich_up$kegg$df, raw_result()$enrich_down$kegg$df)
      openxlsx::writeData(wb, "KEGG", kegg_df, withFilter = TRUE)
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )

  output$download_raw_ion_xlsx <- downloadHandler(
    filename = function() paste0("Ion_annotation_mode", raw_result()$analysis_mode, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content = function(file) {
      req(raw_result())
      wb <- openxlsx::createWorkbook()
      ion_map <- c(calcium = "calcium", magnesium = "magnesium", manganese = "manganese", zinc = "zinc")
      for (ion in names(ion_map)) {
        openxlsx::addWorksheet(wb, ion)
        ts <- raw_result()$ion_term_summary
        ts <- if (!is.null(ts) && nrow(ts) > 0 && "ion" %in% colnames(ts)) ts[ts$ion == ion, , drop = FALSE] else ts[0, , drop = FALSE]
        ps <- raw_result()$ion_protein_summary
        ps <- if (!is.null(ps) && nrow(ps) > 0 && "ion" %in% colnames(ps)) ps[ps$ion == ion, , drop = FALSE] else ps[0, , drop = FALSE]
        gh <- raw_result()$ion_go_hits
        gh <- if (!is.null(gh) && nrow(gh) > 0 && "ion" %in% colnames(gh)) gh[gh$ion == ion, , drop = FALSE] else gh[0, , drop = FALSE]
        next_row <- write_named_table(wb, ion, "GO term summary", ts, start_row = 1)
        next_row <- write_named_table(wb, ion, "Protein summary", ps, start_row = next_row)
        write_named_table(wb, ion, "GO hits", gh, start_row = next_row)
      }
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui = ui, server = server)
