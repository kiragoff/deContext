parse_gbff <- function(gbff_path, de_csv_path) {
  
  message("Step 1: Reading and parsing GenBank flat-file natively...")
  
  lines <- readLines(gbff_path)
  feat_start <- which(grepl("^FEATURES", lines))[1]
  origin_start <- which(grepl("^ORIGIN", lines))[1]
  
  if (is.na(feat_start) || is.na(origin_start)) {
    stop("Invalid GenBank format: Missing FEATURES or ORIGIN block.")
  }
  
  feat_lines <- lines[feat_start:(origin_start - 1)]
  locus_line <- grep("^LOCUS", lines, value = TRUE)[1]
  molecule_name <- strsplit(trimws(locus_line), "\\s+")[[1]][2]
  
  is_feature_row <- grepl("^\\s{5}[A-Za-z]", feat_lines)
  feat_indices <- which(is_feature_row)
  
  features_list <- list()
  
  for (i in seq_along(feat_indices)) {
    start_idx <- feat_indices[i]
    end_idx <- if (i < length(feat_indices)) feat_indices[i + 1] - 1 else length(feat_lines)
    
    block <- feat_lines[start_idx:end_idx]
    first_line_tokens <- strsplit(trimws(block[1]), "\\s+")[[1]]
    
    feat_type <- first_line_tokens[1]
    location_str <- if (length(first_line_tokens) >= 2) first_line_tokens[2] else ""
    
    is_complement <- grepl("complement", location_str)
    strand <- if (is_complement) "-" else "+"
    
    clean_coords <- gsub("[a-zA-Z\\(\\)]", "", location_str)
    coords <- as.numeric(unlist(strsplit(clean_coords, "\\.\\.|,|\\(|\\)")))
    coords <- coords[!is.na(coords)]
    
    if (length(coords) >= 2) {
      gene_start <- min(coords)
      gene_end <- max(coords)
      block_text <- paste(block, collapse = " ")
      
      extract_qual <- function(tag, text) {
        pattern <- sprintf('.*?%s="([^"]+)".*', tag)
        val <- sub(pattern, '\\1', text)
        if (val == text) return(NA_character_)
        return(val)
      }
      
      locus_tag <- extract_qual("locus_tag", block_text)
      gene_name <- extract_qual("gene", block_text)
      product   <- extract_qual("product", block_text)
      operon_id <- extract_qual("operon", block_text)
      
      if (!is.na(locus_tag) || feat_type == "CDS") {
        features_list[[length(features_list) + 1]] <- tibble(
          seqid = ifelse(is.na(molecule_name), "chromosome", molecule_name),
          type = feat_type,
          start = gene_start,
          end = gene_end,
          strand = strand,
          locus_tag = locus_tag,
          gene_name = gene_name,
          product = product,
          operon_id = as.character(operon_id)
        )
      }
    }
  }
  
  gb_df <- bind_rows(features_list) %>%
    filter(type == "CDS", !is.na(locus_tag)) %>%
    group_by(locus_tag) %>%
    summarise(
      seqid = first(seqid),
      start = min(start, na.rm = TRUE),
      end = max(end, na.rm = TRUE),
      strand = first(strand),
      gene_name = first(na.omit(gene_name)),
      product = first(na.omit(product)),
      operon_id = first(na.omit(operon_id)),
      .groups = "drop"
    )
  
  message("Step 2: Reading DESeq2 CSV under strict formatting rules...")
  deseq_data <- read_csv(de_csv_path, show_col_types = FALSE)
  
  # Ensure all column names are strictly unique to prevent duplicate naming errors
  colnames(deseq_data) <- make.unique(colnames(deseq_data))
  
  # Clean potential empty import columns
  names(deseq_data)[names(deseq_data) == ""] <- "row_index_unnamed"
  
  # Standardize gene identifier column name explicitly
  gene_col_candidates <- c("locus_tag", "Gene", "gene", "locus", "ID")
  matched_gene_col <- intersect(gene_col_candidates, colnames(deseq_data))
  if (length(matched_gene_col) > 0) {
    deseq_data <- rename(deseq_data, locus_tag = !!sym(matched_gene_col[1]))
  } else {
    stop("Error: Could not find a gene/locus_tag column in the DESeq2 CSV.")
  }
  
  # Strict Pivot: Expects columns formatted strictly as metric_timepoint using names_sep
  metric_cols <- setdiff(colnames(deseq_data), c("locus_tag", "operon", "row_index_unnamed", "Gene", "product"))
  
  deseq_long <- deseq_data %>%
    pivot_longer(
      cols = matches("_(0[1-9]h|[1-9][0-9]h|[0-9]+)$"),
      names_to = c(".value", "timepoint"),
      names_sep = "_"
    )
  
  # Normalize standard field names safely without creating duplicates
  curr_cols <- colnames(deseq_long)
  lfc_candidate <- curr_cols[grepl("lfc|log2fc|fold|fc", curr_cols, ignore.case = TRUE)][1]
  padj_candidate <- curr_cols[grepl("padj|p_adj|qval|pval", curr_cols, ignore.case = TRUE)][1]
  
  if (!is.na(lfc_candidate) && lfc_candidate != "log2FC" && !("log2FC" %in% curr_cols)) {
    deseq_long <- rename(deseq_long, log2FC = !!sym(lfc_candidate))
  }
  if (!is.na(padj_candidate) && padj_candidate != "padj" && !("padj" %in% curr_cols)) {
    deseq_long <- rename(deseq_long, padj = !!sym(padj_candidate))
  }
  
  if ("operon" %in% colnames(deseq_long)) {
    deseq_long <- deseq_long %>% mutate(operon_id = as.character(operon))
  } else {
    deseq_long$operon_id <- NA_character_
  }
  
  # Ensure deseq_long only contains unique columns needed for the join
  deseq_clean <- deseq_long %>%
    select(any_of(c("locus_tag", "timepoint", "log2FC", "padj", "operon_id"))) %>%
    distinct(locus_tag, timepoint, .keep_all = TRUE)
  
  genomic_df <- gb_df %>%
    left_join(deseq_clean, by = "locus_tag", suffix = c("", "_csv"))
  
  if ("operon_id_csv" %in% colnames(genomic_df)) {
    genomic_df <- genomic_df %>%
      mutate(operon_id = coalesce(operon_id, as.character(operon_id_csv))) %>%
      select(-operon_id_csv)
  }
  
  return(genomic_df)
}

extract_operon_neighborhood <- function(genomic_df, type = "operon", values = NULL, flanks = 3, fallback_method = "gene_count") {
  
  # Ensure type is valid
  type <- match.arg(type, choices = c("operon", "locus_tag", "ko"))
  
  # 1. Identify target gene(s) or operon(s) based on query type
  if (type == "operon") {
    target_rows <- genomic_df %>% filter(operon_id %in% values)
  } else if (type == "locus_tag") {
    target_rows <- genomic_df %>% filter(locus_tag %in% values)
  } else if (type == "ko") {
    if (!"ko" %in% colnames(genomic_df)) stop("Error: 'ko' column not found in genomic_df.")
    target_rows <- genomic_df %>% filter(ko %in% values)
  }
  
  if (nrow(target_rows) == 0) {
    stop("Error: No matching features found for the given query values.")
  }
  
  # 2. Get unique genomic sequence IDs and coordinate bounds for neighborhood extraction
  target_seqid <- unique(target_rows$seqid)[1]
  min_coord <- min(target_rows$start, na.rm = TRUE)
  max_coord <- max(target_rows$end, na.rm = TRUE)
  
  # 3. Extract the local neighborhood on the same chromosome/contig
  # Pull a wider window first, then slice by operon flanks
  chromosome_df <- genomic_df %>% 
    filter(seqid == target_seqid) %>% 
    arrange(start)
  
  # Find unique operons ordered along the chromosome
  unique_operons <- unique(na.omit(chromosome_df$operon_id))
  
  # Target operons check
  target_operons <- unique(na.omit(target_rows$operon_id))
  
  if (length(target_operons) == 0 || all(target_operons == "NA")) {
    
    if (fallback_method == "gene_count") {
      # --- NEW: Index-based gene count fallback ---
      # Find the row positions of our target genes in the sorted chromosome dataframe
      target_indices <- which(chromosome_df$locus_tag %in% target_rows$locus_tag)
      
      if (length(target_indices) == 0) {
        stop("Error: Could not locate target locus tags in the chromosome dataframe.")
      }
      
      # Expand by 'flanks' number of rows, ensuring we don't go past the ends of the chromosome
      min_idx <- max(1, min(target_indices) - flanks)
      max_idx <- min(nrow(chromosome_df), max(target_indices) + flanks)
      
      filtered_data <- chromosome_df[min_idx:max_idx, ]
      
    } else {
      # Old coordinate fallback (kept temporarily until we add architecture)
      filtered_data <- chromosome_df %>%
        filter(end >= (min_coord - (flanks * 5000)) & start <= (max_coord + (flanks * 5000)))
    }
    
  } else {
    # Existing operon-based slicing logic...
    target_indices <- which(unique_operons %in% target_operons)
    min_idx <- max(1, min(target_indices) - flanks)
    max_idx <- min(length(unique_operons), max(target_indices) + flanks)
    
    selected_operons <- unique_operons[min_idx:max_idx]
    
    filtered_data <- chromosome_df %>%
      filter(operon_id %in% selected_operons | (is.na(operon_id) & start >= min_coord & end <= max_coord))
  }
  
  return(filtered_data)
}
