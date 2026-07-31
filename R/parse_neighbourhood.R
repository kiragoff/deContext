#' Module 2: Advanced Query Engine with Operon-Level Flanking & NA Safe-Guards

# NEW: function for autodetecting operon boundaries

detect_operon_boundaries <- function(master_data, query_locus, max_gap = 100) {
  df_sorted <- master_data %>%
    arrange(seqid, start)
  
  # Safely find the exact row index for the query locus or gene name
  match_pos <- which(df_sorted$locus_tag == query_locus | df_sorted$gene_name == query_locus)
  
  if (length(match_pos) == 0) {
    stop("Query gene '", query_locus, "' not found in master data.")
  }
  
  # Take the first match to guarantee a single scalar integer
  target_row_global <- match_pos[1]
  target_seqid <- df_sorted$seqid[target_row_global]
  target_strand <- df_sorted$strand[target_row_global]
  
  replicon_df <- df_sorted %>%
    filter(seqid == target_seqid, strand == target_strand) %>%
    mutate(gene_idx = row_number())
  
  target_locus_tag <- df_sorted$locus_tag[target_row_global]
  
  target_idx <- replicon_df %>%
    filter(locus_tag == target_locus_tag) %>%
    pull(gene_idx)
  
  # Ensure target_idx is strictly a single scalar value
  target_idx <- target_idx[1]
  
  up_idx <- target_idx
  while (up_idx > 1) {
    current_start <- replicon_df$start[up_idx]
    prev_end <- replicon_df$end[up_idx - 1]
    gap <- current_start - prev_end
    
    if (length(gap) == 0 || is.na(gap) || gap > max_gap) break
    up_idx <- up_idx - 1
  }
  
  down_idx <- target_idx
  max_nrow <- nrow(replicon_df)
  while (down_idx < max_nrow) {
    current_end <- replicon_df$end[down_idx]
    next_start <- replicon_df$start[down_idx + 1]
    gap <- next_start - current_end
    
    if (length(gap) == 0 || is.na(gap) || gap > max_gap) break
    down_idx <- down_idx + 1
  }
  
  detected_locus_tags <- replicon_df$locus_tag[up_idx:down_idx]
  return(detected_locus_tags)
}
 
 
parse_neighbourhood <- function(master_data, query_type, query_values, flank_operons = 3, max_gap = 100) {
  
  target_locus_tags <- c()
  
  if (query_type == "operon_auto") {
    target_locus_tags <- c()
    # Keep track of the original query input before it gets overwritten
    original_query <- query_values[1] 
    
    for (val in query_values) {
      auto_locus <- detect_operon_boundaries(master_data, query_locus = val, max_gap = max_gap)
      target_locus_tags <- unique(c(target_locus_tags, auto_locus))
    }
    query_type <- "locus"
    query_values <- target_locus_tags
    
    was_auto_detected <- TRUE
  }
  # 1. Base filtering based on query type
  if (query_type == "operon") {
    base_filtered <- master_data %>% filter(as.character(operon_id) %in% as.character(query_values))
  } else if (query_type %in% c("locus_tag", "locus")) {
    base_filtered <- master_data %>% filter(locus_tag %in% query_values)
  } else if (query_type == "ko" && "KO" %in% colnames(master_data)) {
    base_filtered <- master_data %>% filter(KO %in% query_values)
  } else {
    stop("Invalid query_type specified or required column missing from dataset.")
  }
  
  if (nrow(base_filtered) == 0) {
    warning("No matches found for the given query values.")
    return(base_filtered)
  }
  
  # 2. If operon flanking is requested
  if (flank_operons > 0) {
    message(sprintf("Expanding window by %d complete operon(s) upstream and downstream...", flank_operons))
    
    operon_map <- master_data %>%
      arrange(seqid, start) %>%
      mutate(
        raw_op = if_else(is.na(operon_id) | operon_id == "", paste0("singleton_", row_number()), as.character(operon_id)),
        op_block = if_else(raw_op != lag(raw_op, default = first(raw_op)), 1, 0),
        operon_group_id = cumsum(op_block)
      ) %>%
      group_by(seqid, operon_group_id) %>%
      mutate(temp_unique_op = first(raw_op)) %>% # Keep track of the block identifier
      slice(1) %>% 
      ungroup() %>%
      group_by(seqid) %>%
      mutate(operon_index = row_number()) %>%
      ungroup()
    
    # Map target rows back to their operon block tokens
    target_tokens <- master_data %>%
      arrange(seqid, start) %>%
      mutate(
        raw_op = if_else(is.na(operon_id) | operon_id == "", paste0("singleton_", row_number()), as.character(operon_id)),
        op_block = if_else(raw_op != lag(raw_op, default = first(raw_op)), 1, 0),
        operon_group_id = cumsum(op_block)
      ) %>%
      inner_join(base_filtered %>% select(locus_tag), by = "locus_tag") %>%
      pull(operon_group_id) %>%
      unique()
    
    target_op_indices <- operon_map %>%
      filter(operon_group_id %in% target_tokens) %>%
      pull(operon_index)
    
    expanded_op_indices <- unlist(map(target_op_indices, ~ (.x - flank_operons):(.x + flank_operons)))
    expanded_op_indices <- unique(expanded_op_indices)
    
    valid_groups <- operon_map %>%
      filter(operon_index %in% expanded_op_indices) %>%
      pull(operon_group_id)
    
    final_filtered <- master_data %>%
      arrange(seqid, start) %>%
      mutate(
        raw_op = if_else(is.na(operon_id) | operon_id == "", paste0("singleton_", row_number()), as.character(operon_id)),
        op_block = if_else(raw_op != lag(raw_op, default = first(raw_op)), 1, 0),
        operon_group_id = cumsum(op_block)
      ) %>%
      filter(operon_group_id %in% valid_groups) %>%
      select(-raw_op, -op_block, -operon_group_id)
    
  } else {
    final_filtered <- base_filtered
  }
  
  # Determine if this was an auto-detected run
  mode_label <- if (exists("was_auto_detected") && was_auto_detected) {
    "Auto-Detected Operon"
  } else if (any(grepl("operon_auto", match.call()))) {
    "Auto-Detected Operon"
  } else {
    "Pre-Annotated Operon"
  }
  
  attr(final_filtered, "query_mode") <- if (exists("was_auto_detected") && was_auto_detected) "Auto-Detected Operon" else "Pre-Annotated"
  attr(final_filtered, "query_target") <- if (exists("original_query")) original_query else query_values[1]
  attr(final_filtered, "flank_count") <- flank_operons
  
  return(final_filtered)
}