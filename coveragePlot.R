#version August 21th, 2025
#Modesto Redrejo Rodríguez

#We use here nanopore reads (after Porechop & Nanofilt) mapped against diverse reference genomes with Minimap2


########################################################
#########    STEP 1: load packages         #############
########################################################
paquetes <- c("ggplot2","data.table","Rsubread", "dplyr","tidyverse","BioCircos","circlize")
unavailable <- setdiff(paquetes, rownames(installed.packages()))
if ("Rsubread" %in% unavailable) {   BiocManager::install("Rsubread")}

install.packages(unavailable)
lapply(paquetes, library, character.only = TRUE)
#before going forward check that the directory where the script file is located is the WD
#setwd("~/data/integron2/all_samples/") 

########################################################
############   STEP 2: coverage plots   ################
########################################################

#reads coverage per nt 
#function to parse the coverage tables in a given directory (sample_strain1_cov.tsv.gz) into windows of 1000 nt
process_coverage_files <- function(directory) {
  # List all .tsv.gz files in the directory
  coverage_files <- list.files(directory, pattern = "\\.tsv\\.gz$", full.names = TRUE)
  coverage_list <- list()
  
  # Extract sample names from filenames
  samples <- gsub("\\_cov.tsv\\.gz$", "", basename(coverage_files))
  #check for results directory or create it
  if (!dir.exists("results")) {
    # If it doesn't exist, create it
    dir.create("results")
  }
  for(i in seq_along(coverage_files)) {
    # Read the coverage file
    coverage_data <- fread(coverage_files[i])
    
    # Get unique contigs for this sample
    contigs <- unique(coverage_data$V1)
    
    # Initialize coverage table
    cov <- data.frame(matrix(NA, nrow = 1, ncol = 4))
    
    # Process each contig
    row_counter <- 1
    for(contig in contigs) {
      # Get coverage data for this contig
      contig_coverage <- coverage_data[coverage_data$V1 == contig,]
      
      # Create windows for this contig
      for (j in 1:ceiling(nrow(contig_coverage)/1000)) {
        start_idx <- (j-1)*1000 + 1
        end_idx <- min(j*1000, nrow(contig_coverage))
        
        cov[row_counter,] <- c(
          contig,
          start_idx,
          end_idx,
          mean(contig_coverage$V5[start_idx:end_idx])
        )
        row_counter <- row_counter + 1
      }
      
      # Fix the last window of this contig if needed
      last_window <- which(cov$V1 == contig)[length(which(cov$V1 == contig))]
      cov[last_window,3] <- nrow(contig_coverage)
      cov[last_window,4] <- mean(contig_coverage$V5[(((j-1)*1000)+1):nrow(contig_coverage)])
    }
    
    names(cov) <- c("chr","start","end","value1")
    cov[,2:4] <- apply(cov[,2:4], 2, as.numeric)
    coverage_list[[samples[i]]] <- cov
    
    # Write individual files if needed
    write.table(cov, file=paste0(directory,"/../results/cov_", samples[i], ".tsv"), 
                sep="\t", row.names=FALSE, quote=FALSE)
  }
  
  return(coverage_list)
}

# Call the function with the directory containing the mapping files
coverage_list <- process_coverage_files('mapped/')


#function to plot the circos and save them in pdf files
create_circos_plot <- function(cov_data, sample_name) {
  #open the PDF
  pdf(paste0("results/circos_", sample_name, ".pdf"), width=5, height=5)
  #clear circos
  circos.clear()
  col_text <- "grey30"
  circos.par("track.height"=0.8, gap.degree=5, cell.padding=c(0, 0, 0, 0))
  
  # convert data to numeric
  cov_data$start <- as.numeric(cov_data$start)
  cov_data$end <- as.numeric(cov_data$end)
  cov_data$value1 <- as.numeric(cov_data$value1)
  
  # Get unique contigs and calculate limits
  unique_contigs <- unique(cov_data$chr)
  
  # Calculate xlimits with a small margin
  xlimits <- matrix(0, nrow=length(unique_contigs), ncol=2)
  for(i in seq_along(unique_contigs)) {
    contig_data <- cov_data[cov_data$chr == unique_contigs[i],]
    xlimits[i,] <- c(min(contig_data$start), max(contig_data$end) * 1.01) # 1% extra margen
  }
  
  circos.initialize(factors=unique_contigs, xlim=xlimits)
  
  # Track for contig names with alternating background colors
  bg_colors <- paste0(c("#FF0000", "#00FF00", "#0000FF", "#FF00FF", "#00FFFF"), "40")

  circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
    chr=CELL_META$sector.index
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
    circos.text(mean(xlim), mean(ylim), chr, cex=0.8, font=2, col=col_text, 
                facing="bending.inside", niceFacing=TRUE)
  }, bg.col=rep(bg_colors, length.out=length(unique_contigs)), 
  bg.border=F, track.height=0.1)
  
  # Axis track with adaptive breaks
  max_length <- max(xlimits[,2])
  if(max_length > 5e5) {
    break_interval <- max_length/10
  } else {
    break_interval <- max_length/15
  }
  brk <- seq(0, max_length, by=break_interval)
  
  circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
    circos.axis(h="top", major.at=brk, labels=round(brk/10^6, 1), labels.cex=0.6, 
                col=col_text, labels.col=col_text, lwd=2, labels.facing="clockwise")
  }, bg.border=F)
  
  # Calculate ylim for coverage track
  coverage_values <- log10(cov_data$value1 + 1)
  ylim_range <- range(coverage_values[is.finite(coverage_values)])
  ylim_range <- c(floor(ylim_range[1]), ceiling(ylim_range[2]))
  
  # Coverage track (log10 scale)
  circos.track(factors=cov_data$chr, 
               x=cov_data$start, 
               y=coverage_values, 
               panel.fun=function(x, y) {
                 circos.lines(x, y, col="#1F78B4", lwd=1)
               }, 
               ylim=ylim_range,
               track.height=0.3, 
               bg.border=F)

  # Add sample name in the center
  text(0, 0, sample_name, 
       col = col_text, 
       cex = 1.2,           # font size
       font = 2)          # 2 = bold

  # Close PDF device
  dev.off()
}

# Create plots for all samples
for(sample_name in 1:length(coverage_list)) {
  if (dir.exists(paste0("results/",names(coverage_list)[sample_name]))) {
    cat(sprintf("Directory '%s' already exists.\n", names(coverage_list)[sample_name]))
  } else  {
    dir.create(paste0("results/",names(coverage_list)[sample_name]), recursive = TRUE, showWarnings = FALSE)
    }  
  tryCatch({
    create_circos_plot(coverage_list[[sample_name]], names(coverage_list)[sample_name])
  }, error = function(e) {
    message(sprintf("Error en muestra %s: %s", sample_name, e$message))
  })
}

