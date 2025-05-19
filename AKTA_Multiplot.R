library(dplyr)
library(viridis)
library(readr)

plotChromMultiple <- function(filenames, 
                              wkDir, 
                              samples, 
                              plots = c('A280', 'A254'),
                              plotStartEnd = NULL, 
                              plotPostInj = TRUE,
                              axNames = c("Volume (mL)", "mAU", 'Conductivity (mS/cm)'),
                              plotFracs = FALSE,
                              rotateFracs = FALSE,
                              outputLocation = NULL,
                              overlay = TRUE) {
  
  #--------------------------------------------
  # Helper function to extract columns
  #--------------------------------------------
  get_column_data <- function(df, identifier, subtract_start = FALSE, subsetStart = 0, data_is_numeric = TRUE) {
    index <- which(sapply(df, function(x) dplyr::first(x) == identifier))
    if (length(index) == 0) {
      warning(paste("Identifier", identifier, "not found in data frame"))
      return(NULL)
    }
    # Row 1 is the column header (identifier), row 2 is typically units, row 3+ are data
    mL_column   <- df[[index]][3:nrow(df)]
    data_column <- df[[index + 1]][3:nrow(df)]
    
    # Convert the first column (mL) to numeric
    mL <- suppressWarnings(as.numeric(mL_column)) - ifelse(subtract_start, subsetStart, 0)
    
    # Convert the second column conditionally to numeric or keep as character
    if (data_is_numeric) {
      data <- suppressWarnings(as.numeric(data_column))
    } else {
      data <- data_column
    }
    
    return(data.frame(mL = mL, data = data))
  }
  
  #--------------------------------------------
  # Color palette for each file
  #--------------------------------------------
  file_colors <- viridis::viridis(length(filenames))
  
  #--------------------------------------------
  # If outputLocation is provided, set up EPS
  #--------------------------------------------
  if (!is.null(outputLocation)) {
    if (!dir.exists(outputLocation)) {
      dir.create(outputLocation, recursive = TRUE)
    }
    eps_file <- file.path(outputLocation, "chromatograms.eps")
    cairo_ps(file = eps_file, width = 6, height = 12, onefile = TRUE, fallback_resolution = 300)
  }
  
  #--------------------------------------------
  # First pass: read data from files, store,
  # and compute global ranges if overlay=TRUE
  #--------------------------------------------
  all_data <- vector("list", length(filenames))
  
  global_minX <- Inf
  global_maxX <- -Inf
  global_maxY <- -Inf
  
  for (i in seq_along(filenames)) {
    filename <- filenames[i]
    sample   <- samples[i]
    
    # Load the data
    df_old <- read_tsv(file.path(wkDir, filename),
                       locale = locale(encoding = 'UTF-16LE'),
                       col_types = cols(.default = 'c'))
    
    # Determine subsetStart if plotPostInj = TRUE
    subsetStart <- 0
    if (plotPostInj) {
      runLog <- get_column_data(df_old, 'Run Log', data_is_numeric = FALSE)
      if (!is.null(runLog)) {
        sample_app_index <- which(runLog$data == 'Sample Application') + 1
        if (length(sample_app_index) > 0) {
          subsetStart <- runLog$mL[sample_app_index]
        }
      }
    }
    
    # Extract the data for the specified plots
    plot_data_list <- list()
    for (plot_name in plots) {
      plot_data <- get_column_data(df_old, plot_name, subtract_start = plotPostInj, subsetStart = subsetStart)
      if (!is.null(plot_data)) {
        plot_data_list[[plot_name]] <- plot_data
      }
    }
    
    # If fraction info is needed:
    frac_data <- NULL
    if (plotFracs) {
      frac_data <- get_column_data(df_old, 'Fraction',
                                   subtract_start = plotPostInj,
                                   subsetStart    = subsetStart,
                                   data_is_numeric = FALSE)
    }
    
    # Determine local start/end
    if (!is.null(plotStartEnd)) {
      minX <- plotStartEnd[1]
      maxX <- plotStartEnd[2]
    } else {
      minX <- min(sapply(plot_data_list, function(d) min(d$mL, na.rm = TRUE)), na.rm = TRUE)
      maxX <- max(sapply(plot_data_list, function(d) max(d$mL, na.rm = TRUE)), na.rm = TRUE)
    }
    
    # Shift each curve so min in chosen x-range is zero
    local_maxY <- 0
    for (pn in names(plot_data_list)) {
      pd <- plot_data_list[[pn]]
      in_range <- pd$mL >= minX & pd$mL <= maxX
      data_in_range <- pd$data[in_range]
      min_data_value <- min(data_in_range, na.rm = TRUE)
      pd$data <- pd$data - min_data_value
      plot_data_list[[pn]] <- pd
      local_maxY <- max(local_maxY, max(pd$data[in_range], na.rm = TRUE))
    }
    
    # Store in all_data
    all_data[[i]] <- list(
      sample = sample,
      col_main = file_colors[i],
      plot_data_list = plot_data_list,
      frac_data = frac_data,
      minX = minX,
      maxX = maxX,
      maxY = local_maxY
    )
    
    # Update global ranges (for overlay=TRUE)
    global_minX <- min(global_minX, minX)
    global_maxX <- max(global_maxX, maxX)
    global_maxY <- max(global_maxY, local_maxY)
  }
  
  #--------------------------------------------
  # Plotting
  #--------------------------------------------
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  # -- overlay mode --
  if (overlay) {
    
    # Prepare global axes
    xAxTix <- pretty(c(global_minX, global_maxX))
    yAxTix <- pretty(c(0, global_maxY * 1.05))
    
    plot(NA,
         xlim = c(global_minX, global_maxX),
         ylim = c(0, max(yAxTix)),
         xlab = axNames[1],
         ylab = axNames[2],
         axes = FALSE,
         main = "Overlay of All Samples")
    axis(1, at = xAxTix, cex.axis = 1)
    axis(2, at = yAxTix, las = 1, cex.axis = 1)
    
    # Overlaid lines for each sample
    for (i in seq_along(all_data)) {
      ad <- all_data[[i]]
      sample       <- ad$sample
      col_main     <- ad$col_main
      plot_data_li <- ad$plot_data_list
      frac_data    <- ad$frac_data
      
      main_rgb <- col2rgb(col_main)/255
      lighter_rgb <- (main_rgb + c(1,1,1))/2
      col_260 <- rgb(lighter_rgb[1], lighter_rgb[2], lighter_rgb[3])
      
      lty_280 <- 1 # solid
      lty_260 <- 3 # dotted (example)
      
      # Plot lines
      if ('UV 1_280' %in% names(plot_data_li)) {
        pd_280 <- plot_data_li[['UV 1_280']]
        lines(pd_280$mL, pd_280$data, col = col_main, lwd = 3, lty = lty_280)
      }
      if ('UV 2_260' %in% names(plot_data_li)) {
        pd_260 <- plot_data_li[['UV 2_260']]
        lines(pd_260$mL, pd_260$data, col = col_260, lwd = 3, lty = lty_260)
      }
      
      # Fractions (if requested)
      if (plotFracs && !is.null(frac_data)) {
        in_range <- !is.na(frac_data$mL) & frac_data$mL >= ad$minX & frac_data$mL <= ad$maxX
        frac_positions <- frac_data$mL[in_range]
        frac_labels    <- frac_data$data[in_range]
        
        abline(v = frac_positions, lty = 3, col = 'gray40')
        text(x = frac_positions,
             y = rep(max(yAxTix) * 0.95, length(frac_positions)),
             labels = frac_labels,
             srt = if (rotateFracs) 45 else 0,
             pos = 4,
             offset = 0.3,
             cex = 0.6)
      }
    }
    
    # A simple combined legend
    legend("topright",
           legend = samples,
           col    = file_colors,
           lty    = 1,
           lwd    = 3,
           bty    = 'n',
           cex    = 1)
    
  } else {
    
    # -- non-overlay mode --
    for (i in seq_along(all_data)) {
      ad <- all_data[[i]]
      sample       <- ad$sample
      col_main     <- ad$col_main
      plot_data_li <- ad$plot_data_list
      frac_data    <- ad$frac_data
      minX         <- ad$minX
      maxX         <- ad$maxX
      maxY         <- ad$maxY
      
      yAxTix <- pretty(c(0, maxY * 1.05))
      xAxTix <- pretty(c(minX, maxX))
      
      plot(NA,
           xlim = c(minX, maxX),
           ylim = c(0, max(yAxTix)), 
           xlab = axNames[1],
           ylab = axNames[2],
           axes = FALSE,
           main = sample,
           cex.main = 1.2)
      axis(1, at = xAxTix, cex.axis = 1)
      axis(2, at = yAxTix, las = 1, cex.axis = 1)
      
      main_rgb <- col2rgb(col_main)/255
      lighter_rgb <- (main_rgb + c(1,1,1))/2
      col_260 <- rgb(lighter_rgb[1], lighter_rgb[2], lighter_rgb[3])
      
      lty_280 <- 1 # solid
      lty_260 <- 3 # dotted
      
      # Plot lines
      if ('UV 1_280' %in% names(plot_data_li)) {
        pd_280 <- plot_data_li[['UV 1_280']]
        lines(pd_280$mL, pd_280$data, col = col_main, lwd = 4, lty = lty_280)
      }
      if ('UV 2_260' %in% names(plot_data_li)) {
        pd_260 <- plot_data_li[['UV 2_260']]
        lines(pd_260$mL, pd_260$data, col = col_260, lwd = 4, lty = lty_260)
      }
      
      legend("topright",
             legend = c("A280", "A260"),
             col = c(col_main, col_260),
             lty = c(lty_280, lty_260),
             lwd = 4,
             bty = 'n',
             cex = 1,
             xpd = NA,
             inset = c(0, -0.15))
      
      # Plot fraction lines
      if (plotFracs && !is.null(frac_data)) {
        in_range <- !is.na(frac_data$mL) & frac_data$mL >= minX & frac_data$mL <= maxX
        frac_positions <- frac_data$mL[in_range]
        frac_labels    <- frac_data$data[in_range]
        
        abline(v = frac_positions, lty = 3, col = 'gray40')
        text(x      = frac_positions,
             y      = rep(max(yAxTix) * 0.95, length(frac_positions)), 
             labels = frac_labels,
             srt    = if (rotateFracs) 45 else 0,
             pos    = 4,
             offset = 0.3,
             cex    = 0.6)
      }
    }
  }
  
  #--------------------------------------------
  # Close device if opened
  #--------------------------------------------
  if (!is.null(outputLocation)) {
    dev.off()
  }
}


# Example usage (adjust paths and filenames as needed):
plotChromMultiple(
  filenames      = c("20250430_HisTrap_HP_2x5mL_AP3497_benz.csv","20250430_HisTrap_HP_2x5mL_AP3497_PEI.csv"),
  wkDir          = '/Users/andrew.grassetti/Library/CloudStorage/OneDrive-SharedLibraries-AeraTherapeutics/Protein Sciences - Documents/Protein Sciences Lab/Protein Sciences Instruments/AKTA_Pure_Systems/AKTA5/Andrew/Chromatogram_CSVs/20250502/',
  samples        = c("No PEI","PEI"),
  plots          = c('UV 1_280', 'UV 2_260'),
  plotStartEnd   = c(0, 1200),
  plotFracs      = FALSE,
  rotateFracs    = TRUE,
  plotPostInj    = FALSE
  # outputLocation = "..."
)
