rm(list=ls())
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(grid)

# ===== celltypes compares and conditions that needed to be interated =====
compares <- c("AC_to_MG", "MG_to_AC", 'AC_to_OL', 'OL_to_AC', 
              'MG_to_OL', 'OL_to_MG')
conditions <- c("WT_Veh","WT_LPS","5XFAD_Veh","5XFAD_LPS")
comparisions <- c('5XFAD_LPS_vs_5XFAD_Veh', 
                  '5XFAD_LPS_vs_WT_LPS',
                  '5XFAD_Veh_vs_WT_Veh', 
                  'WT_LPS_vs_WT_Veh')

# ===== specify input directory and output directory =====
input_root  <- "/Users/tianmi/Desktop/Plots/Cross_Talk"
output_root <- "/Users/tianmi/Desktop/Plots/Cross_Talk"  

# ===== specify input directory and output directory of Male=====
input_root_M  <- "/Users/tianmi/Desktop/Plots/Cross_Talk/Male"
output_root_M <- "/Users/tianmi/Desktop/Plots/Cross_Talk/Male" 

# ===== specify input directory and output directory of Female=====
input_root_F  <- "/Users/tianmi/Desktop/Plots/Cross_Talk/Female"
output_root_F <- "/Users/tianmi/Desktop/Plots/Cross_Talk/Female" 

# ===== mapping for celltype names=====
cell_label <- list(
  AC  = list(title = "Astrocytes",        legend = "Astrocyte"),
  MG  = list(title = "Microglia",         legend = "Microglia"),
  OL  = list(title = "Oligodendrocytes",  legend = "Oligodendrocytes")
)
get_label <- function(code, which=c("title","legend")){
  which <- match.arg(which)
  if (code %in% names(cell_label)) cell_label[[code]][[which]] else code
}

# # ==== specify the conditional files====
# files <- c(
#   "WT_Veh"     = "/Users/tianmi/Desktop/Plots/Cross_Talk/Geno_Treat/AC_to_MG/WT_Veh/AC_to_MG_filtered.csv",
#   "WT_LPS"     = "/Users/tianmi/Desktop/Plots/Cross_Talk/Geno_Treat/AC_to_MG/WT_LPS/AC_to_MG_filtered.csv",
#   "5XFAD_Veh"  = "/Users/tianmi/Desktop/Plots/Cross_Talk/Geno_Treat/AC_to_MG/5XFAD_Veh/AC_to_MG_filtered.csv",
#   "5XFAD_LPS"  = "/Users/tianmi/Desktop/Plots/Cross_Talk/Geno_Treat/AC_to_MG/5XFAD_LPS/AC_to_MG_filtered.csv"
# )

# ==== read in df and convert into chordDiagram ====
read_interactions <- function(csv_path) {
  x <- read.csv(csv_path, stringsAsFactors = FALSE)
  df <- data.frame(
    from = x$ligand_complex,
    to   = x$receptor_complex,
    regulatory_potential = x$lr_logfc
  )
  
  df
}

dfs <- lapply(files, read_interactions)

# ==== assign palette for regulatory legend ==== 
col_fun <- circlize::colorRamp2(c(0, 2), c("white", "red"))  

# ==== function for chord ====
plot_chord <- function(df, cond, sector_order, grid.col, col_fun,
                       from_title, to_title, 
                       from_legend, to_legend,
                       big_gap = 20, small_gap = 2,
                       split_duplicate_names = FALSE,
                       out_pdf = NULL) {
  df$color <- col_fun(df$regulatory_potential)
  
  circos.clear()
  
  if (isTRUE(split_duplicate_names)) {

    df$from_node <- paste0(df$from, " (L)")
    df$to_node   <- paste0(df$to,   " (R)")
    from_nodes <- unique(df$from_node)
    to_nodes   <- unique(df$to_node)
    sector_order2 <- c(from_nodes, to_nodes)
  } else {

    from_nodes <- unique(df$from)
    to_nodes   <- setdiff(unique(df$to), from_nodes)
    sector_order2 <- c(from_nodes, to_nodes)
  }
  
  if (isTRUE(split_duplicate_names)) {
    if (is.null(grid.col)) {
      grid.col <- c(
        setNames(rep("brown2",   length(from_nodes)), from_nodes),
        setNames(rep("darkgreen",length(to_nodes)),   to_nodes)
      )
    } else {

      missing <- setdiff(sector_order2, names(grid.col))
      if (length(missing) > 0) {
        fill_cols <- ifelse(missing %in% from_nodes, "brown2", "darkgreen")
        grid.col  <- c(grid.col, setNames(fill_cols, missing))
      }
    }
  } else {
    if (is.null(grid.col)) {
      grid.col <- c(
        setNames(rep("brown2",   length(from_nodes)), from_nodes),
        setNames(rep("darkgreen",length(to_nodes)),   to_nodes)
      )
    } else {
      missing <- setdiff(sector_order2, names(grid.col))
      if (length(missing) > 0) {
        fill_cols <- ifelse(missing %in% from_nodes, "brown2", "darkgreen")
        grid.col  <- c(grid.col, setNames(fill_cols, missing))
      }
    }
  }
  
  gaps <- c(
    rep(small_gap, max(length(from_nodes) - 1, 0)), big_gap,
    rep(small_gap, max(length(to_nodes)   - 1, 0)), big_gap
  )
  
  circos.par(canvas.xlim = c(-1.5, 1.5),
             canvas.ylim = c(-1.5, 1.5),
             track.margin = c(0.01, 0.01),
             start.degree = -90,
             "track.height" = 0.1)
  
  if (!is.null(out_pdf)) pdf(out_pdf, width = 8, height = 8)
  
  if (isTRUE(split_duplicate_names)) {
    chordDiagram(
      df[, c("from_node","to_node","regulatory_potential","color")],
      order = sector_order2,
      directional = 1,
      col = df$color,
      grid.col = grid.col,
      direction.type = c("diffHeight", "arrows"),
      link.arr.type = "big.arrow",
      link.zindex = rank(df$regulatory_potential, ties.method = "first"),
      annotationTrack = "grid",
      annotationTrackHeight = mm_h(2),
      preAllocateTracks = list(track.height = mm_h(3), track.margin = c(mm_h(1), 0))
    )
  } else {
    chordDiagram(
      df,
      order = sector_order2,
      directional = 1,
      col = df$color,
      grid.col = grid.col,
      direction.type = c("diffHeight", "arrows"),
      link.arr.type = "big.arrow",
      link.zindex = rank(df$regulatory_potential, ties.method = "first"),
      annotationTrack = "grid",
      annotationTrackHeight = mm_h(2),
      preAllocateTracks = list(track.height = mm_h(3), track.margin = c(mm_h(1), 0))
    )
  }
  
  
  highlight.sector(unique(df$from), track.index = 1, col = "brown2",  border = "black", lwd = 1.2)
  highlight.sector(unique(df$to),   track.index = 1, col = "darkgreen",border = "black", lwd = 1.2)
  
  circos.track(track.index = 2, ylim = c(0, 1.5), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + mm_y(12),
                CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, cex = 1.2, adj = 0)
  }, bg.border = NA)
  
  grid.text(paste0(from_title, " to ", to_title, " in Male Mice", "\n", cond),
            x = 0.5, y = 0.95, gp = gpar(fontsize = 18, fontface = "bold"))
  
  lgd_links <- Legend(at = c(0, 0.5, 1), col_fun = col_fun,
                      title_position = "topcenter", title = "LR LogFC")
  draw(lgd_links, x = unit(20, "mm"), y = unit(140, "mm"), just = c("left","bottom"))
  
  lgd_LR <- Legend(at = c(paste(from_legend, "ligands"), paste(to_legend, "receptors")),
                   legend_gp = gpar(fill = c("brown2","darkgreen")),
                   title = "LR pairs", title_position = "topcenter")
  draw(lgd_LR, x = unit(20, "mm"), y = unit(40, "mm"),  just = c("left","top"))
  
  if (!is.null(out_pdf)) dev.off()
}


# ===== interation for celltypes compare × condition =====
for (compare in compares) {
  # split parts, preparing for full celltype names
  parts <- strsplit(compare, "_to_", fixed = TRUE)[[1]]
  from_celltype <- parts[1]; to_celltype <- parts[2]
  from_title  <- get_label(from_celltype, "title")
  to_title    <- get_label(to_celltype,   "title")
  from_legend <- get_label(from_celltype, "legend")
  to_legend   <- get_label(to_celltype,   "legend")
  
  # ==== assemble the input file path under each condition ==== 
  files <- setNames(file.path(input_root, compare, comparisions, paste0(compare, "_filtered.csv")),
                    comparisions)
  
  # ==== read on ==== 
  dfs <- lapply(files, read_interactions)
  nonnull <- Filter(Negate(is.null), dfs)
  if (!length(nonnull)) { message("[SKIP] ", compare, "No data"); next }
  
  # ==== uniform the sectors and colors ====
  sectors_all <- unique(unlist(lapply(dfs, function(df) union(df$from, df$to))))
  base_colors <- brewer.pal(8, "Dark2")
  cols <- colorRampPalette(base_colors)(length(sectors_all))
  # assign specific colors for each sector
  grid.col <- setNames(as.character(cols), sectors_all)
  sector_order <- sectors_all
  
  # ==== assemble the output file path under each condition ==== 
  out_dir <- file.path(output_root, compare, "Chord_Plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # output PDFs
  mapply(function(df, cond){
    if (is.null(df)) return(NULL)
    out_pdf <- file.path(out_dir, paste0(cond, "_chord.pdf"))
    plot_chord(df, cond, sector_order, grid.col, col_fun,
               from_title, to_title, from_legend, to_legend,
               out_pdf = out_pdf)
  }, dfs, names(dfs))
  
  message("[OK] Done：", out_dir)
}


# ===== interation for celltypes compare × condition in Male=====
for (compare in compares) {
  # split parts, preparing for full celltype names
  parts <- strsplit(compare, "_to_", fixed = TRUE)[[1]]
  from_celltype <- parts[1]; to_celltype <- parts[2]
  from_title  <- get_label(from_celltype, "title")
  to_title    <- get_label(to_celltype,   "title")
  from_legend <- get_label(from_celltype, "legend")
  to_legend   <- get_label(to_celltype,   "legend")
  
  # ==== assemble the input file path under each condition ==== 
  files <- setNames(file.path(input_root_M, compare, comparisions, paste0(compare, "_filtered.csv")),
                    comparisions)
  
  # ==== read on ==== 
  dfs <- lapply(files, read_interactions)
  nonnull <- Filter(Negate(is.null), dfs)
  if (!length(nonnull)) { message("[SKIP] ", compare, " No data"); next }
  
  # ==== uniform the sectors and colors ====
  sectors_all <- unique(unlist(lapply(dfs, function(df) union(df$from, df$to))))
  base_colors <- brewer.pal(8, "Dark2")
  cols <- colorRampPalette(base_colors)(length(sectors_all))
  # assign specific colors for each sector
  grid.col <- setNames(as.character(cols), sectors_all)
  sector_order <- sectors_all
  
  # ==== assemble the output file path under each condition ==== 
  out_dir <- file.path(output_root_M, compare, "Chord_Plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # output PDFs
  mapply(function(df, cond){
    if (is.null(df)) return(NULL)
    out_pdf <- file.path(out_dir, paste0(cond, "_chord.pdf"))
    plot_chord(df, cond, sector_order, grid.col, col_fun,
               from_title, to_title, from_legend, to_legend,
               out_pdf = out_pdf)
  }, dfs, names(dfs))
  
  message("[OK] Done：", out_dir)
}


# ===== interation for celltypes compare × condition in Female=====
for (compare in compares) {
  # split parts, preparing for full celltype names
  parts <- strsplit(compare, "_to_", fixed = TRUE)[[1]]
  from_celltype <- parts[1]; to_celltype <- parts[2]
  from_title  <- get_label(from_celltype, "title")
  to_title    <- get_label(to_celltype,   "title")
  from_legend <- get_label(from_celltype, "legend")
  to_legend   <- get_label(to_celltype,   "legend")
  
  # ==== assemble the input file path under each condition ==== 
  files <- setNames(file.path(input_root_F, compare, comparisions, paste0(compare, "_filtered.csv")),
                    comparisions)
  
  # ==== read on ==== 
  dfs <- lapply(files, read_interactions)
  nonnull <- Filter(Negate(is.null), dfs)
  if (!length(nonnull)) { message("[SKIP] ", compare, " No data"); next }
  
  # ==== uniform the sectors and colors ====
  sectors_all <- unique(unlist(lapply(dfs, function(df) union(df$from, df$to))))
  base_colors <- brewer.pal(8, "Dark2")
  cols <- colorRampPalette(base_colors)(length(sectors_all))
  # assign specific colors for each sector
  grid.col <- setNames(as.character(cols), sectors_all)
  sector_order <- sectors_all
  
  # ==== assemble the output file path under each condition ==== 
  out_dir <- file.path(output_root_F, compare, "Chord_Plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # output PDFs
  mapply(function(df, cond){
    if (is.null(df)) return(NULL)
    out_pdf <- file.path(out_dir, paste0(cond, "_chord.pdf"))
    plot_chord(df, cond, sector_order, grid.col, col_fun,
               from_title, to_title, from_legend, to_legend,
               out_pdf = out_pdf)
  }, dfs, names(dfs))
  
  message("[OK] Done：", out_dir)
}



