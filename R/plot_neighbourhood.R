#' Module 3: Operon Visualization Engine with Automated Dynamic Titling & Saving
#' these are fallback options if the use doesn't provide anything different
plot_neighbourhood <- function(filtered_data, 
                                     p_thresh = 0.05, 
                                     annotation_col = "product", 
                                     style = c("none", "above", "table"),
                                     annot_subset = NULL,
                                     clean_short_genes = TRUE, 
                                     deseq_path = NULL,          # <-- Changed from DESEQ_PATH
                                     save_plot = TRUE,
                                     output_dir = "operon_results", # <-- Changed from out_dir
                                     width = 10,
                                     height = 6,
                                     dpi = 300,
                                     colour_limits = NULL,
                                     custom_palette = c("deepskyblue3", "#ffffbf", "darkred")) {
  if (nrow(filtered_data) == 0) {
    stop("Oops: Cannot plot empty dataset. Check your query parameters.")
  }
  
  style <- match.arg(style)
  
  if (style != "none") {
    if (is.null(deseq_path)) {
      stop("Oops: You must provide a `deseq_path` when using style = '", style, "'.")
    }
    if (file.exists(deseq_path)) {
      master_csv <- read.csv(deseq_path, stringsAsFactors = FALSE)
      
      matched_col <- names(master_csv)[tolower(names(master_csv)) == tolower(annotation_col)]
      if (length(matched_col) > 0) {
        annotation_col <- matched_col[1]
        
        annotation_lookup <- master_csv %>%
          select(locus_tag, all_of(annotation_col)) %>%
          distinct(locus_tag, .keep_all = TRUE)
        
        if (annotation_col %in% names(filtered_data)) {
          filtered_data[[annotation_col]] <- NULL
        }
        
        filtered_data <- filtered_data %>%
          left_join(annotation_lookup, by = "locus_tag")
      } else {
        stop("Oops: Annotation column '", annotation_col, "' not found in CSV file: ", deseq_path)
      }
    } else {
      stop("Oops: DESEQ_PATH file not found at: ", deseq_path)
    }
  }
  
  plot_df <- filtered_data %>%
    mutate(
      gene_label = if_else(!is.na(gene_name) & gene_name != "", gene_name, locus_tag),
      sig_star = if_else(!is.na(padj) & padj < p_thresh, "*", "")
    )
  
  timepoints_sorted <- sort(unique(plot_df$timepoint))
  plot_df$timepoint <- factor(plot_df$timepoint, levels = timepoints_sorted)
  
  # --- AUTOMATED DYNAMIC TITLE, SUBTITLE, CAPTION & FILENAME TAG ---
  flank_used <- attr(filtered_data, "flank_count")
  if (is.null(flank_used)) flank_used <- 0
  
  mode_attr <- attr(filtered_data, "query_mode")
  is_auto_detected <- identical(mode_attr, "Auto-Detected Operon")
  
  # Build the caption string using your p_thresh argument
  plot_caption <- paste0("Significance threshold: padj < ", p_thresh)
  
  if (is_auto_detected) {
    query_gene_tag <- attr(filtered_data, "query_target")
    if (is.null(query_gene_tag)) {
      query_gene_tag <- plot_df$locus_tag[1]
    }
    plot_title <- paste0("Auto-Detected Operon Neighborhood: ", query_gene_tag)
    plot_subtitle <- paste0("Flanking Window: ", flank_used, " operon(s)")
    file_tag <- paste0("auto_", gsub("[^[:alnum:]]", "_", query_gene_tag))
  } else {
    operon_tag <- if ("operon_id" %in% names(plot_df) && any(!is.na(plot_df$operon_id))) {
      as.character(na.omit(plot_df$operon_id)[1])
    } else {
      as.character(plot_df$locus_tag[1])
    }
    plot_title <- paste0("Predefined Operon Neighborhood: ", operon_tag)
    plot_subtitle <- paste0("Flanking Window: ", flank_used, " operon(s)")
    file_tag <- paste0("operon_", operon_tag)
  }
  
  # Base gggenomes / gggenes plot
  p <- ggplot(
    plot_df, 
    aes(
      xmin = start, 
      xmax = end, 
      y = timepoint,         
      fill = log2FC,         
      forward = (strand == "+" | strand == "forward" | strand == 1),
      gene = gene_label
    )
  ) +
    geom_gene_arrow(
      arrowhead_height = unit(6, "mm"), 
      arrow_body_height = unit(4, "mm")        
    )
  
  # Conditionally add internal labels with smart color-flipping for background contrast
  if (!clean_short_genes) {
    p <- p + geom_gene_label(
      aes(
        label = gene_label,
        color = if_else(abs(log2FC) >= (max(abs(log2FC), na.rm = TRUE) / 2), "white", "black")
      ), 
      align = "centre", 
      fontface = "italic",
      min.size = 3
    ) +
      scale_color_identity()
  }
  
  # Significance stars sit on the clean background, so they remain standard black/bold
  p <- p + geom_text(
    aes(
      x = (start + end) / 2, 
      y = timepoint, 
      label = sig_star
    ),
    vjust = if (clean_short_genes) 0.5 else 1.8, 
    size = 4, 
    fontface = "bold",
    color = "black", 
    show.legend = FALSE
  ) +
    scale_y_discrete(limits = rev(timepoints_sorted)) +
    scale_fill_gradient2(
      low = custom_palette[1],  
      mid = custom_palette[2],  
      high = custom_palette[3],  
      midpoint = 0,
      limits = colour_limits,   
      labels = scales::label_number(accuracy = 0.1), 
      name = "Log2 FC"
    ) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(face = "bold", size = 10)
    ) +
    labs(
      x = "Genomic Coordinates (bp)",
      title = plot_title,
      subtitle = plot_subtitle,
      caption = plot_caption
    )
  
  # Handling display styles and clean external labels independently
  if (style == "above") {
    annot_data <- plot_df %>%
      group_by(locus_tag, start, end) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(
        mid_x = (start + end) / 2,
        # Default to gene_label, but include annotation_col if it's distinct from gene names
        label_text = gene_label 
      )
    
    if (!is.null(annot_subset)) {
      annot_data <- annot_data %>% 
        mutate(label_text = ifelse(locus_tag %in% annot_subset, label_text, ""))
    }
    
    p_final <- p + 
      geom_text_repel(
        data = annot_data,
        aes(x = mid_x, y = timepoints_sorted[length(timepoints_sorted)], label = label_text, fill = NULL),
        inherit.aes = FALSE,
        nudge_y = 0.5,
        segment.color = "lightgrey",
        segment.size = 0.5,
        direction = "y",
        force = 2,
        size = 3,
        box.padding = 0.5,
        max.overlaps = Inf,
        angle = 0,
        hjust = 0.5
      )
    
  } else if (style == "table") {
    table_data <- plot_df %>%
      arrange(start) %>%
      select(locus_tag, start, end, operon_id, annotation = all_of(annotation_col)) %>%
      distinct(locus_tag, .keep_all = TRUE) %>%
      mutate(row_id = row_number())
    
    p_table <- ggplot(table_data, aes(y = factor(row_id))) +
      geom_text(aes(x = 0.7, label = locus_tag), hjust = 0, fontface = "bold", size = 3.5, angle = 0) +
      geom_text(aes(x = 2, label = ifelse(is.na(operon_id) | operon_id == "", "-", operon_id)), hjust = 0, size = 3.5, angle = 0) +
      geom_text(aes(x = 3, label = ifelse(is.na(annotation), "", annotation)), hjust = 0, size = 3, color = "grey30", angle = 0) +
      xlim(0.5, 5) + 
      scale_y_discrete(limits = rev(levels(factor(table_data$row_id)))) +
      theme_void() +
      theme(
        plot.margin = margin(t = 5, r = 5, b = 10, l = 10),
        plot.subtitle = element_text(face = "bold", size = 9, color = "black")
      ) +
      labs(subtitle = paste0("Metadata Summary: Locus Tag | Operon ID | ", annotation_col))
    
    # If clean_short_genes is TRUE alongside table style, add external gene labels to the top gene panel track
    if (clean_short_genes) {
      annot_data <- plot_df %>%
        group_by(locus_tag, start, end) %>%
        slice(1) %>%
        ungroup() %>%
        mutate(
          mid_x = (start + end) / 2,
          label_text = gene_label
        )
      
      p_genes_annotated <- p + 
        geom_text_repel(
          data = annot_data,
          aes(x = mid_x, y = timepoints_sorted[length(timepoints_sorted)], label = label_text, fill = NULL),
          inherit.aes = FALSE,
          nudge_y = 0.5,
          segment.color = "lightgrey",
          segment.size = 0.5,
          direction = "y",
          force = 2,
          size = 3,
          box.padding = 0.5,
          max.overlaps = Inf,
          angle = 0,
          hjust = 0.5
        )
    } else {
      p_genes_annotated <- p
    }
    
    p_final <- plot_grid(p_genes_annotated + theme(legend.position = "bottom"), 
                         p_table, 
                         ncol = 1, 
                         rel_heights = c(1, 0.7))
    
  } else if (clean_short_genes) {
    annot_data <- plot_df %>%
      group_by(locus_tag, start, end) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(
        mid_x = (start + end) / 2,
        label_text = gene_label
      )
    
    p_final <- p + 
      geom_text_repel(
        data = annot_data,
        aes(x = mid_x, y = timepoints_sorted[length(timepoints_sorted)], label = label_text, fill = NULL),
        inherit.aes = FALSE,
        nudge_y = 0.5,
        segment.color = "lightgrey",
        segment.size = 0.5,
        direction = "y",
        force = 2,
        size = 3,
        box.padding = 0.5,
        max.overlaps = Inf,
        angle = 0,
        hjust = 0.5
      )
  } else {
    p_final <- p
  }
  
  # Automated Dynamic Saving with custom output directory support
  if (save_plot) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    tp_string <- paste(sort(unique(as.character(plot_df$timepoint))), collapse = "-")
    clean_suffix <- if (clean_short_genes) "_clean" else ""
    
    filename <- paste0(file_tag, "_flank_", flank_used, "_tps_", tp_string, "_style_", style, clean_suffix, ".png")
    file_path <- file.path(output_dir, filename)
    
    ggsave(filename = file_path, plot = p_final, width = width, height = height, dpi = dpi)
    message("Plot successfully saved to: ", file_path)
  }
  
  return(p_final)
}