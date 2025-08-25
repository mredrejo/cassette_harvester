#version August 20th, 2025
#Modesto Redrejo Rodríguez


########################################################
#########    STEP 1: load packages         #############
########################################################
paquetes <- c("ggplot2","data.table","Rsubread", "dplyr","tidyverse")
unavailable <- setdiff(paquetes, rownames(installed.packages()))
if ("Rsubread" %in% unavailable) {   BiocManager::install("Rsubread")}

install.packages(unavailable)
lapply(paquetes, library, character.only = TRUE)
#before going forward check that the directory where the script file is located is the WD
#setwd("~/data/integron2/all_samples/") 



########################################################
######   STEP 2: load & prepare annotations   ##########
########################################################
#This block parse Bakta and IntegronFinder output tables to generate a combined table containing full genome annotation
#A log file will be created for each parsed genome
#Cassettes are defined as the sequence between two attC motifs
#check the filenames for bakta table: it should be strain.tsv (not strain.fna.tsv or strain.fasta.tsv)
#There two functions that should be loaded: 
#(1) process_annotations()
#(2) process_multiple_annotations()
#
process_annotations <-  function(strain, bakta_file, integron_finder_file, manual_attc_file = NULL) {
  # Read Bakta annotations (keep only CDS)
  bakta <- fread(bakta_file)
  bakta <- bakta[bakta$Type == "cds", ]
  if (!dir.exists("result")) {
    # If it doesn't exist, create it.
    dir.create("results")
    message("Results directory for sample created successfully!")
  } else {
    message("Results directory for sample already exists.")
  }
  log_con <- file(paste0("results/",strain,"_combined_annotation.log"),open="a")
  cat("\n", date(), "\n \n",file=log_con)
  cat("Parsing Bakta and IntegronFinder annotation files",bakta_file," and ",integron_finder_file,"\n","\n",file=log_con)
  # Rename columns for consistency
  names(bakta)[c(1,3,4,6,7,8)] <- c("Chr", "Start", "End", "GeneID","Gene_name","Product")
  
  # Read integron_finder results
  integrons <- fread(integron_finder_file, skip = 2)
  # Keep only attC and integrase and rename columns
  integrons <- subset(integrons,annotation=="attC") #|annotation=="intI")
  integrons <- data.frame(integrons[,c(1:6,11)])
  colnames(integrons) <- c("ID_integron","Chr", "GeneID", "Start", "End", "Strand","type")
  integrons$Strand <- ifelse(integrons$Strand == 1, '+', '-')
  integrons$merge <- paste0(integrons$ID_integron,"_",integrons$Chr)
  # Add manual annotations if provided
  if (!is.null(manual_attc_file)) {
    manual <- readxl::read_excel(manual_attc_file, sheet = 1)
    names(manual) <- names(integrons)
    integrons <- rbind(integrons, manual)
  }
  
  # Process attC IDs to remove 
  integrons$GeneID[grep("att",integrons$GeneID)] <- gsub("_[[:digit:]]+", "", integrons$GeneID[grep("att",integrons$GeneID)])
  integrons$GeneID <- gsub("attc", "attC", integrons$GeneID)
  integrons <- integrons[order(integrons$merge,integrons$Start), ]
  
  # Verify chromosome matches
  common_chr <- intersect(unique(bakta$Chr), unique(integrons$Chr))
  if(length(common_chr) == 0) {
    stop("No matching chromosomes found between Bakta and Integron Finder results")
  }  
  # Print unique chromosomes/contigs in both files
  message("Contigs in Bakta annotations:")
  print(unique(bakta$Chr))
  message("Integron Finder results:")
  cat(length(unique(integrons$merge)), "integrons in", length(unique(integrons$Chr)), "contigs")
  cat(paste("Contigs in Bakta annotations: ",length(unique(bakta$Chr))), "\n", file = log_con, append = TRUE)
  cat("Integron results:", length(unique(integrons$merge)), "integrons in", length(unique(integrons$Chr)), "contigs \n \n", file = log_con, append = TRUE)
  
  
  
  # Filter to keep only matching chromosomes
  #bakta <- bakta[bakta$Chr %in% common_chr, ]
  integrons <- integrons[integrons$Chr %in% common_chr, ]
  
  # Number attCs by strand
  for(integron in unique(integrons$merge)) {
    chr_attcs <- integrons[integrons$merge == integron, ]
    if(nrow(chr_attcs) > 0) {
      strand <- chr_attcs$Strand[chr_attcs$GeneID!="attC"]
      count <- nrow(chr_attcs[chr_attcs$GeneID=="attC",])
      if("+" %in% strand) {
        chr_attcs$GeneID[chr_attcs$GeneID=="attC"] <- paste0("attC_", count:1)
      } else {
        chr_attcs$GeneID[chr_attcs$GeneID=="attC"] <- paste0("attC_", 1:count)
      }
      integrons[integrons$merge == integron, ] <- chr_attcs
    }
  }
 
  # Add cassette annotations as inter-attC segments
  for(integron in unique(integrons$merge)) {
    chr_attcs <- integrons[integrons$merge == integron, ]
    if(nrow(chr_attcs) > 1) {
      for(i in 1:(nrow(chr_attcs)-1)) {
        if(chr_attcs[i, "End"] < chr_attcs[(i+1), "Start"]) {
          new_cassette <- data.frame(
            ID_integron = integron,
            Chr = chr_attcs$Chr,
            GeneID = paste0("Cassette_", i),
            Start = as.numeric(chr_attcs[i, "End"]) + 1,
            End = as.numeric(chr_attcs[(i+1), "Start"]) - 1,
            Strand = chr_attcs[i+1, "Strand"],
            merge = integron,
            type = chr_attcs$type
          )
          integrons <- rbind(integrons[,1:8], new_cassette)
        }
      }
    }
  }
  integrons <- unique(integrons)
  #merge feature & integron names to avoid duplicate feature names in table
  integrons$GeneID <- paste0(integrons$GeneID,"_",integrons$merge)
  #rename ID_integron to Product
  names(integrons)[7] <- "Annotation"
  integrons$Annotation[grep("complete",integrons$Annotation)] <- "Complete Integron"

  # Combine annotations
  names(bakta)[8] <- "Annotation"
  final_annotations <- rbind(
    bakta[, c("GeneID", "Chr", "Start", "End", "Strand","Annotation")],
    integrons[, c("GeneID", "Chr", "Start", "End", "Strand","Annotation")]
  )
  
  # Sort by chromosome and start position
  final_annotations <- unique(final_annotations[order(final_annotations$Chr, final_annotations$Start), ])
  
  # Print summary
  cat("\nAnnotations table is in the folder results","\n\n",file = log_con, append = TRUE)

  cat("Annotation summary: \n", file = log_con, append = TRUE)
  cat("Total features in genome: ", nrow(final_annotations),"\n",file = log_con, append = TRUE)
  cat("Total CDS: ", sum(grepl("^(?!attC_|Cassette)", final_annotations$GeneID, perl=TRUE)) ,"\n",file = log_con, append = TRUE)
  cat("Total attCs: ", sum(grepl("^attC_", final_annotations$GeneID, perl=TRUE)),"\n",file = log_con, append = TRUE)
  cat("Total Cassettes: ", sum(grepl("^Cassette", final_annotations$GeneID, perl=TRUE)),"\n",file = log_con, append = TRUE)
  cat("Note that cassettes are defined as the sequence between two attC motifs or the integrase gene and the first attC.","\n",file = log_con, append = TRUE)
 
  close(log_con)
  return(final_annotations)
 
}
  
process_multiple_annotations <- function(strains, base_dir = "refs") {
  annotation_list <- list()
  
  for(strain in strains) {
    tryCatch({
      # Construct file paths based on sample name
      strain <- strain
      
      bakta_file <- file.path(base_dir, "annotated", 
                              strain, 
                              paste0(strain, ".tsv"))
      
      integron_finder_file <- file.path(base_dir, "integron_finder_results",
                                        paste0("Results_Integron_Finder_", strain),
                                        paste0(strain, ".integrons"))
      
      manual_attc_file <- file.path(base_dir, "integron_finder_results",
                                    paste0(strain, "_attC_manually.xlsx"))
      
      # Check if manual attC file exists
      if (!file.exists(manual_attc_file)) {
        manual_attc_file <- NULL
        message(sprintf("No manual attC file found for sample %s", strain))
      }
      
      # Process annotations for this sample
      annotation_list[[strain]] <- process_annotations(
        strain = strain,
        bakta_file = bakta_file,
        integron_finder_file = integron_finder_file,
        manual_attc_file = manual_attc_file
      )
      write.csv(annotation_list[[strain]],paste0("results/annotations_",strain,".csv"),row.names = FALSE)
     
      
    }, error = function(e) {
      message(sprintf("Error processing sample %s: %s", sample_name, e$message))
    })
  }
  
  return(annotation_list)
}

# Use the function with the samples from coverage files
# You can provide the names of the Fasta files in a vector or obtain them from the 'refs/fasta' folder as follows
genomes <- unique(sub("\\..*", "", list.files("refs/fasta", pattern="*.fna", full.names=F)))

process_multiple_annotations(genomes)

