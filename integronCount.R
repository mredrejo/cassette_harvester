#version August 22th, 2025
#Modesto Redrejo Rodríguez

# This script uses the function featureCounts() from package Rsubread and the annotations table from 
# previous script annotations.R to quantify the mapped reads (BAM file created with minimap2) in each CDS and integron Cassette.

########################################################
#########    STEP 1: load packages         #############
########################################################
paquetes <- c("ggplot2","data.table","gggenes","Rsubread", "dplyr","tidyverse","ggpubr")
unavailable <- setdiff(paquetes, rownames(installed.packages()))
if ("Rsubread" %in% unavailable) {   BiocManager::install("Rsubread")}
unavailable <- setdiff(paquetes, rownames(installed.packages()))
install.packages(unavailable)
lapply(paquetes, library, character.only = TRUE)


###################################################################
######   STEP 1: count mapped reads with featureCounts   ##########
###################################################################

# Annotations files generated with script "annotations.R" should be in "results" folder
# Mapping BAM files should be in "mapped" folder
#You can add a list of samples and strains in vectors or read them from the directory with the following function

samples_by_strain  <- function(base_dir) {
  # 1. Construct full paths and check if directories exist
  annot_dir <- file.path(base_dir, "results")
  bam_dir <- file.path(base_dir, "mapped")
  
  if (!dir.exists(annot_dir) || !dir.exists(bam_dir)) {
    stop("The specified directory structure is incorrect.
             'results' or 'mapped' subdirectory is missing.")
  }
  
  # 2. List the files in each directory
  annotation_files <- list.files(path = annot_dir, pattern = "^annotations_.*\\.csv$", full.names = TRUE)
  bam_files <- list.files(path = bam_dir, pattern = ".*_mapped\\.bam$", full.names = FALSE)
  
  # 3. Create clean vectors of names for matching and output
  # The regex for BAM files is modified to capture the full sample_strain name
  extract_clean_name <- function(file_name) {
    if (grepl("annotations_.*\\.csv$", file_name)) {
      return(gsub("^.*annotations_(.*)\\.csv$", "\\1", file_name))
    }
    if (grepl("_mapped\\.bam$", file_name)) {
      return(gsub("_mapped.bam","", file_name))
    }
    return(NA_character_)
  }
  
  clean_annot_strains <- sapply(annotation_files, extract_clean_name, USE.NAMES = FALSE)
  clean_bam_names <- sapply(bam_files, extract_clean_name, USE.NAMES = FALSE)
  
  # 4. Extract just the strain from the clean BAM names for matching
  bam_strains <- gsub(".*_", "", clean_bam_names)
  
  # 5. Match the clean annotation strains with the BAM file strains
  matched_bam_log <- rep(FALSE, length(bam_strains))
  matched_list <- list()
  
  for (i in seq_along(clean_annot_strains)) {
    annot_strain <- clean_annot_strains[i]
    
    matching_indices <- which(bam_strains == annot_strain)
    
    if (length(matching_indices) > 0) {
      matched_bam_log[matching_indices] <- TRUE
      # Store the full clean BAM names (sample_strain)
      matched_list[[annot_strain]] <- clean_bam_names[matching_indices]
    } else {
      warning(paste("No matching BAM file found for annotation strain:", annot_strain))
      matched_list[[annot_strain]] <- NA_character_
    }
  }
  
  # 6. Check for any unmatched BAM files
  unmatched_bam_indices <- which(!matched_bam_log)
  if (length(unmatched_bam_indices) > 0) {
    unmatched_bam_files <- bam_files[unmatched_bam_indices]
    warning(paste("The following BAM files have no annotation match:", 
                  paste(unmatched_bam_files, collapse = ", ")))
  }
  
  # 7. Return the final list
  matched_vec <- unlist(matched_list, use.names = FALSE)
  names(matched_vec) <- rep(names(matched_list), times = lengths(matched_list))
  matched_list <- matched_vec
  return(matched_list)
}

matched_list <- samples_by_strain(".")


# The following function will parse samples and annotation files from a list
# As in the function samples_by_strain(), elements of the list must be the names of the samples in format sample_strain
# Names of the elements of the list must be the matched strain names
process_feature_counts <- function(matched_list, bam_dir="mapped", results_dir="results") {
  results_list <- list()
  #check and create directories
  
  for(i in 1:length(matched_list)) {
    if (!dir.exists(paste0(results_dir,"/", matched_list[i]))) {
      # If it doesn't exist, create it.
      # recursive = TRUE ensures that parent directories are also created if they don't exist.
      dir.create(paste0(results_dir,"/", matched_list[i]), recursive = TRUE)
      message("Results directory for sample created successfully!")
    } else {
      message("Results directory for sample already exists.")
    }
    #log
    sink(paste0(results_dir,"/", matched_list[i], "/", matched_list[i],".log"), append=TRUE, split=TRUE)
    paste(Sys.Date(),"\n")
    #START processing  
      # Get BAM file path and annotations for this sample
      bam_file <- file.path(bam_dir, paste0(matched_list[i], "_mapped.bam"))
      sample_annot <- fread(file.path(results_dir, paste0("annotations_",names(matched_list[i]), ".csv")))
      (paste("Processing", matched_list[i]))
      
      # Run featureCounts allowing overlapping of reads between cassettes and CDS within the cassettes
      counts <- featureCounts(
        files = bam_file,
        annot.ext = na.omit(sample_annot),
        isGTFAnnotationFile = FALSE,
        isPairedEnd = FALSE,
        countMultiMappingReads = TRUE,
        allowMultiOverlap = TRUE,
        isLongRead = TRUE,
        fraction = FALSE,
        ignoreDup = FALSE,
        reportReads = "CORE",
        nthreads = 12
      )
      
      # Create counts table
      table_counts <- data.frame(
        Feature_name = row.names(counts$counts),
        Counts = counts$counts[,1],
        sample_annot[match(row.names(counts$counts), 
                           sample_annot$GeneID), c("Chr", "Start", "End", "Strand")]
      )
      # Ensure numeric columns
      table_counts$Start <- as.numeric(table_counts$Start)
      table_counts$End <- as.numeric(table_counts$End)
      
      # Add type information
      table_counts$Type <- "cds"  # Default type
      table_counts$Type[grep("attC_", table_counts$Feature_name)] <- "attC"
      table_counts$Type[grep("Cassette", table_counts$Feature_name)] <- "cassette"
      
      # Calculate CPM
      map_stats <-  readLines(gsub("mapped.bam","map_stats.txt",bam_file), warn = FALSE)
      mapped_reads <- gsub("^([0-9]+).*", "\\1", map_stats[7])
      # Handle cases where the line is not found
      if (length(mapped_reads) == 0) {
        warning(paste("Mapped reads line not found in:", stats_file))
        return(NA_integer_)
      }
      table_counts$CPM <- as.numeric(table_counts$Counts) * 1000000 / as.numeric(mapped_reads)
      
      # Find superintegron coordinates based on attC positions and strand
      attC_positions <- table_counts[grep("attC_", table_counts$Feature_name), ]
      integron_name <- unique(gsub("attC_[0-9]+_(integron_[0-9]+)_.*", "\\1", attC_positions$Feature_name))
      contig_name <- unique(gsub("attC_[0-9]+_integron_[0-9]+_(.*)", "\\1", attC_positions$Feature_name))
      integron_stats <- data.frame(cbind(integron_name,contig_name))
      for (j in 1:nrow(integron_stats)){
        integron_stats$start[j] <- min(attC_positions$Start[intersect(grep(integron_stats$integron_name[j],attC_positions$Feature_name) , which(attC_positions$Chr==integron_stats$contig_name[j]))])
        integron_stats$end[j] <- max(attC_positions$End[intersect(grep(integron_stats$integron_name[j],attC_positions$Feature_name) , which(attC_positions$Chr==integron_stats$contig_name[j]))])
        kk <- table_counts[table_counts$Chr==integron_stats$contig_name[j] & table_counts$Type=="cassette" & table_counts$Start >= integron_stats$start[j] & table_counts$End <= integron_stats$end[j], ]
        integron_stats$cassettes[j] <- nrow(kk)
        integron_stats$CPM_mean[j] <- round(mean(kk$CPM),3)
        integron_stats$Ratio[j] <- round(sum(kk$Counts)*100/(sum(kk$Counts)+sum(table_counts$Chr==integron_stats$contig_name[j] & table_counts$Counts[table_counts$Start < integron_stats$start[j] & table_counts$End < integron_stats$end[j] ])),2)
      }
      
      #write stats_table
      write.table(integron_stats,paste0(results_dir,"/",matched_list[i],"/integrons_stats_",matched_list[i],".csv"),row.names = F,quote=F,sep=";")
      
      # Store results
      results_list[[i]] <- table_counts
      # Write full table
        write.table(table_counts, paste0(results_dir,"/",matched_list[i],"/",matched_list[i],".csv"),
                    quote = FALSE, sep = ";", row.names = FALSE)
        
        
        # Create plots
        # CDS counts plot
        counts_plot <- ggplot(data = na.omit(table_counts[table_counts$Type == "cds",])) +
          geom_point(aes(x = reorder(Feature_name, as.numeric(Start)),
                         y = as.numeric(CPM),
                         colour = ifelse(Chr %in% contig_name & 
                                           inrange(na.omit(table_counts[table_counts$Type == "cds",]$Start), integron_stats$start, integron_stats$end, incbounds=TRUE) &
                                           inrange(na.omit(table_counts[table_counts$Type == "cds",]$End), integron_stats$start, integron_stats$end, incbounds=TRUE),
                                         "Integron(s)", "Chromosome"))) +
          facet_grid(~Chr, scales = "free", space = "free_x") +
          theme_classic() +
          scale_x_discrete(expand = c(0.01, 0)) +
          # scale_y_continuous(limits = c(0, 100000)) +
          theme(axis.text.x = element_blank(),
                axis.text.y = element_text(size = 12, face = "bold")) +
          ylab("Normalized CPM") + xlab("CDS") + labs(color = "")
        
        ggsave(paste0(results_dir,"/",matched_list[i],"/counts_chromosomes.pdf"),
               plot = counts_plot, width = 8, height = 5)
        ggsave(paste0(results_dir,"/",matched_list[i],"/counts_chromosomes_log.pdf"),
               plot = counts_plot + scale_y_log10(), width = 8, height = 5)
        
        # Cassette counts plot
        cassettes_plot <- ggplot(data = na.omit(table_counts[table_counts$Type == "cassette",])) +
          geom_point(aes(x = reorder(Feature_name, as.numeric(Start)),
                         y = as.numeric(CPM),
                         colour = "integron"),
                     size = 3) +
          theme_classic() +
          scale_x_discrete(expand = c(0.01, 0),
                           breaks = table_counts[table_counts$Type == "cassette",]$Feature_name[rep(c(F,F,F,T))]) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),
                axis.text.y = element_text(size = 12, face = "bold"),
                legend.position = "none") +
          ylab("Normalized CPM") + xlab("Cassette") + labs(color = "")
        
        ggsave(paste0(results_dir,"/",matched_list[i],"/","counts_cassettes_integron.pdf"),
               plot = cassettes_plot, width = 6, height = 8)
        ggsave(paste0(results_dir,"/",matched_list[i],"/", "counts_cassettes_integron_log10.pdf"),
               plot = cassettes_plot + scale_y_log10(), width = 6, height = 8)
        
        #Correlation cassette counts vs length
        table_counts$length <- table_counts$End - table_counts$Start
        ggscatter(na.omit(table_counts[table_counts$Type == "cassette",]), x = "CPM", y = "length", color = "#CC0000", shape = 21, size = 2,
                  add = "reg.line", add.params = list(color = "#CC0000", size = 2), conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman") + 
          xlab("Cassette CPM (Normalized counts per M reads)") + ylab("Cassette length")
        ggsave(paste0(results_dir,"/",matched_list[i],"/", "cor_cassette_length_count.pdf"),
                  width = 5, height = 5)
        
        #CDS multimapped reads
        multimap <- fread(paste0(matched_list[i],"_mapped.bam.featureCounts"),sep="\t",header=FALSE)
        
        #number of cassettes per read
        multimap <- multimap[multimap$V3>0,]
        cassettes <- stringr::str_count(multimap$V4, "Cassette")
        stat_cassette <- as.data.frame(table(cassettes))
        stat_cassette$perc <- stat_cassette$Freq*100/sum(stat_cassette$Freq)
        reads_cas <- ggplot(data=stat_cassette)+geom_bar(aes(x=cassettes,y=Freq),stat="identity", col="black",fill="grey60") + 
          geom_text(aes(x=cassettes,y=Freq, label=paste0(round(perc,2),"%")),vjust = -1, nudge_y  = -.5)+
          ylim(0,max(stat_cassette$Freq)*1.1)+
          xlab("Number of cassettes") +ylab("Mapped reads") + theme_bw()+theme(text=element_text(size=12))
        ggsave(paste0(results_dir,"/",matched_list[i],"/", "reads_per_cassete.pdf"),
               width = 6, height = 5)
        #reads_cas <- ggplot(data=stat_cassette)+geom_bar(aes(x=cassettes,y=Freq),stat="identity", col="black",fill="grey60") + 
        #  geom_text(aes(x=cassettes,y=Freq, label=paste0(round(perc,2),"%")),vjust = -1, nudge_y  = -.5)+
        #  scale_y_log10() +
        #  xlab("Number of cassettes") +ylab("Mapped reads") + theme_bw()+theme(text=element_text(size=12))
        #ggsave(paste0(results_dir,"/",matched_list[i],"/", "reads_per_cassete_log.pdf"),
        #       width = 6, height = 5)
        
        #number of CDS per read
        genome_tag <- unique(sub("_.*","",sample_annot$GeneID))[which(unique(sub("_.*","",sample_annot$GeneID))!="attC" & unique(sub("_.*","",sample_annot$GeneID))!= "Cassette" & unique(sub("_.*","",sample_annot$GeneID))!= "attI")]
        CDS <- stringr::str_count(multimap$V4, paste0(genome_tag,"_"))
        stat_CDS <- as.data.frame(table(CDS))
        stat_CDS$perc <- stat_CDS$Freq*100/sum(stat_CDS$Freq)
        reads_CDS <- ggplot(data=stat_CDS)+geom_bar(aes(x=CDS,y=Freq),stat="identity", col="black",fill="grey60") + 
          geom_text(aes(x=CDS,y=Freq, label=paste0(round(perc,2),"%")),vjust = -1, nudge_y  = -.5)+
          ylim(0,max(stat_CDS$Freq)*1.1)+
          xlab("Number of CDS") +ylab("Mapped reads") + theme_bw()+theme(text=element_text(size=12))
        ggsave(paste0(results_dir,"/",matched_list[i],"/", "reads_per_CDS.pdf"),
               width = 6, height = 5)
        #reads_CDS <- ggplot(data=stat_CDS)+geom_bar(aes(x=CDS,y=Freq),stat="identity", col="black",fill="grey60") + 
        #  geom_text(aes(x=CDS,y=Freq, label=paste0(round(perc,2),"%")),vjust = -1, nudge_y  = -.5)+
        #  scale_y_log10() +
        #  xlab("Number of CDS") +ylab("Mapped reads") + theme_bw()+theme(text=element_text(size=12))
        #ggsave(paste0(results_dir,"/",matched_list[i],"/", "reads_per_CDS_log.pdf"),
        #       width = 6, height = 5)
        
        (paste("Successfully processed sample: %s", matched_list[i]))
        sink()
        }
}



# Use the function with the files obtained in samples_by_strain() function above
process_feature_counts(matched_list)

