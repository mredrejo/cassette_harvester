#version August 26th, 2024
#We use nanopore reads (after Porechop & Nanofilt) mapped against the whole N16961 genome


########################################################
#########    STEP 1: load packages         #############
########################################################
paquetes <- c("ggplot2","data.table","gggenes","Rsubread", "dplyr","tidyverse","ggpubr")
unavailable <- setdiff(paquetes, rownames(installed.packages()))
if ("Rsubread" %in% unavailable) {   BiocManager::install("Rsubread")}

install.packages(unavailable)
lapply(paquetes, library, character.only = TRUE)

########################################################
######   STEP 2: load & prepare annotations   ##########
########################################################
#N16961 GCF_003063785.1.fna annotations from Bakta (keep only cds annotations)
N16961 <- fread("../bakta_annotations/N16961_GCF_000006745.1/N16961.tsv", skip=5)
N16961 <- N16961[N16961$Type=="cds", ]
N16961$`#Sequence Id` <- ifelse(N16961$`#Sequence Id`=='contig_1',"NC_002505.1","NC_002506.1")

#integron_finder2 results
N16961_superIn <- fread("../integron_finder_results/Results_Integron_Finder_GCF_000006745.1_ASM674v1_genomic_e4/GCF_000006745.1_ASM674v1_genomic.integrons", skip=2)

#keep only attC and create annotations of cassettes as inter-vcr stretches
N16961_superIn <- as.data.frame(N16961_superIn[N16961_superIn$type_elt=="attC",2:6])
colnames(N16961_superIn) <- c("Chr","GeneID","Start","End","Strand")
N16961_superIn$Strand <- ifelse(N16961_superIn$Strand==1,'+','-')

#add more manual annotations from Filipa and fix attC numbers
manual <- readxl::read_excel("../integron_finder_results/attC_manually.xlsx", sheet=1)
names(manual) <- names(N16961_superIn)
N16961_superIn <- rbind(N16961_superIn,manual)
N16961_superIn$GeneID <- gsub("_[[:digit:]]+","",N16961_superIn$GeneID)
N16961_superIn$GeneID <- gsub("attc","attC",N16961_superIn$GeneID)
N16961_superIn <- N16961_superIn[order(N16961_superIn$Start),]
N16961_superIn$GeneID[N16961_superIn$GeneID=="attC"] <- paste0("attC","_",1:table(N16961_superIn$GeneID=="attC")["TRUE"])

#add "cassette" annotations as inter-attC segments
for (i in 1:(nrow(N16961_superIn)-1)){
  if (N16961_superIn[i,4] < N16961_superIn[(i+1),3]){
    N16961_superIn[(nrow(N16961_superIn)+i),] <- c("NC_002506.1", 
                                "Cassette",
                                as.numeric(N16961_superIn[i,4]) + 1, 
                                as.numeric(N16961_superIn[(i+1),3]) - 1,
                                N16961_superIn[i+1,5])
  }
}
N16961_superIn <- na.omit(N16961_superIn[order(N16961_superIn$Start),]   ) 
#add the last cassette
N16961_superIn[(nrow(N16961_superIn)+1),] <- c("NC_002506.1", 
                                               "Cassette",
                                               as.numeric(N16961_superIn[N16961_superIn$GeneID=="attC_178",4]) + 1, 
                                               435387,"+")
 
N16961_superIn$GeneID[N16961_superIn$GeneID=="Cassette"] <- paste0("Cassette_",1:nrow(N16961_superIn[N16961_superIn$GeneID=="Cassette",]))   

#rename cols to allow merge
names(N16961)[c(1,3,4,6)] <- c("Chr","Start","End","GeneID")
colnames(N16961_superIn) <- c("Chr","GeneID","Start","End","Strand")
N16961_superIn$Strand <- ifelse(N16961_superIn$Strand==1,'+','-')

#final table
N16961full <- rbind(N16961[,c(6,1,3,4,5)],N16961_superIn)


########################################################
######       STEP 3: count mapped reads       ##########
########################################################
#minimap -ax map-ont 
counts <- featureCounts(files=paste0("../data/mapped_GCF_000006745_reads_chopped_filtered_minimap.bam"),annot.ext=na.omit(N16961full),isGTFAnnotationFile=FALSE,isPairedEnd=FALSE,countMultiMappingReads=TRUE, allowMultiOverlap=TRUE,isLongRead=TRUE, fraction = FALSE, ignoreDup = FALSE,reportReads="CORE",nthreads = 12)
table_counts <- data.frame(cbind(row.names(counts$counts), counts$counts, N16961full$Chr, N16961full$Start,N16961full$End))
table_counts$CPM <- as.numeric(table_counts$mapped_GCF_000006745_reads_chopped_filtered_minimap.bam)*1000000/72108

#add annotations to CDS
table_counts_final <- merge(table_counts,N16961[N16961$Type=="cds",], by.x=c("V4","V5"), by.y=c("Start","End"), all.x=TRUE)
#reorder and save
table_counts_final[grep("attC_",table_counts_final$V1),8] <- "attC"
table_counts_final[grep("Cassette",table_counts_final$V1),8] <- "cassette"
table_counts_final$Strand[table_counts_final$Type=="attC"] <- "+"
table_counts_final <- table_counts_final[,c(5,3,1,2,9,8,4,6,11,12)]
names(table_counts_final) <- c("Chromosome","Feature_name","Start","Stop","Strand","Type","Counts","CPM","Gene","Bakta Annotation")
table_counts_final[,c(3,4,7,8)] <- apply(table_counts_final[,c(3,4,7,8)],2,as.numeric)

#write full table of counts 
write.table(table_counts_final,"../results_GCF_000006745_179cass/coverage_N16961_IF2_nofrac.csv", quote=FALSE,sep=";",row.names = FALSE)

#write table of counts within superintegron
write.table(subset(table_counts_final,Chromosome=="NC_002506.1" & Start>309750 & Stop<435388),"../results_GCF_000006745_179cass/coverage_superintegron_nofrac.csv", quote=FALSE,sep=";",row.names = FALSE)

#ratio CDS mapping
mean(subset(table_counts_final,Chromosome=="NC_002506.1" & Start>309750 & Stop<435388 & Type=="cds")$Counts) / 
  mean(subset(table_counts_final,Chromosome=="NC_002506.1" & Start<309750 & Stop>435388 & Type=="cds" | Chromosome=="NC_002505.1")$Counts)

#plot counts
counts <- ggplot(data=na.omit(table_counts_final[table_counts_final$Type=="cds",])) +
       geom_point(aes(x=reorder(Feature_name,as.numeric(Start)),y=as.numeric(Counts), colour=ifelse(as.numeric(Start)>309750 & as.numeric(Stop)<435387 & Chromosome=="NC_002506.1" ,"Superintegron","Chromosome"))) +
       facet_grid(~Chromosome,scales="free")+ theme_classic() +
       scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits = c(0, 2500), breaks = seq(0, 2500, 500))+
       theme(axis.text.x = element_blank(),axis.text.y=element_text(size=12,face="bold")) +
       ylab("Counts") + xlab("CDS") + labs(color="")
ggsave("../results_GCF_000006745_179cass/counts_chromosomes.pdf", plot=counts, width=8,height=5)
ggsave("../results_GCF_000006745_179cass/counts_chromosomes_log10.pdf", plot=counts+ scale_y_log10(), width=8,height=5)

#plot counts/Cassette
casettes <- ggplot(data=na.omit(table_counts_final[table_counts_final$Type=="cassette",c(2,3,7)])) +
  geom_point(aes(x=reorder(Feature_name,as.numeric(Start)),y=(as.numeric(Counts)), colour="Superintegron"),size=3) +
   theme_classic() +
  scale_x_discrete(expand = c(0.01, 0),breaks = table_counts_final[table_counts_final$Type=="cassette",]$Feature_name[rep(c(F,F,F,T))]) +
  theme(axis.text.x=element_text(angle=90,hjust = 1),axis.text.y=element_text(size=12,face="bold"), legend.position = "none") +
  ylab("Counts") + xlab("Cassette") + labs(color="")
ggsave("../results_GCF_000006745_179cass/counts_cassettes_integron.pdf", plot=casettes, width=6,height=5)
ggsave("../results_GCF_000006745_179cass/counts_cassettes_integron_log10.pdf", plot=casettes+ scale_y_log10(), width=6,height=5)


#plot density/cassettes
ggplot(data=table_counts_final[table_counts_final$Type=="cassette",], aes(x=log10(as.numeric(Counts)+1))) + geom_density()


########################################################
######     STEP 4: check multimapped reads    ##########
########################################################

#CDS multimapped reads
multimap <- fread("mapped_GCF_000006745_reads_chopped_filtered_minimap.bam.featureCounts",sep="\t",header=FALSE)
#number of cassettes per read
multimap <- multimap[multimap$V3>0,]
cassettes <- stringr::str_count(multimap$V4, "Cassette")
table(cassettes)
stat_cassette <- as.data.frame(table(cassettes))
ggplot(data=stat_cassette)+geom_bar(aes(x=cassettes,y=(Freq)),stat="identity", col="grey30") + xlab("Number of cassettes") +ylab("Mapped reads (log10)") + theme_bw()+theme(text=element_text(size=12))

CDS <- stringr::str_count(multimap$V4, "N16961_")
table(CDS)
stat_CDS <- as.data.frame(table(CDS))
ggplot(data=stat_CDS)+geom_bar(aes(x=CDS,y=(Freq)),stat="identity", col="grey30") + xlab("Number of CDS") +ylab("Mapped reads (log10)") + theme_bw()+theme(text=element_text(size=12))



#reads that map in >1 site
table(duplicated(multimap$V1))
which(duplicated(multimap$V1))

#write table
write.table(multimap[multimap$V2=="Assigned" & multimap$V3>0,],"../results_GCF_000006745_179cass/multimapped_N16916_final.csv", quote=FALSE,sep=";")

#check insertion points outside superintegron
alignments <- fread(input = "../data/mapped_GCF_000006745_reads_chopped_filtered_minimap.bed")

spurious <- alignments[alignments$V1=="NC_002505.1" | alignments$V1=="NC_002506.1" & alignments$V2 < 309750 | alignments$V1=="NC_002506.1" & alignments$V2 >435388, ]
spurious <- spurious[spurious$V5>10,] #filter by alignment quality

spurious_GTF <- data.frame(spurious$V1,"N16961","insertions",spurious$V2-100,spurious$V2+100,".",spurious$V6,".","bla")

write.table(spurious_GTF,"../results_GCF_000006745/spurious_insertion_sites_Q30_100.gtf",col.names = FALSE,row.names = FALSE,quote=FALSE,sep="\t")
writeLines(spurious$V4,"../results_GCF_000006745/spurious_reads.txt")

#how often the outside SI alignments correspond with duplicated alignments?
table(spurious$V4 %in% alignments[which(duplicated(alignments$V4)),V4])

#reads that map ONLY outside integron
onlyoutside <- spurious$V4[-which(spurious$V4 %in% alignments[which(duplicated(alignments$V4)),V4])]
writeLines(onlyoutside,"../results_GCF_000006745_179cass/spurious_reads_onlyoutside.txt")

spurious_only <- spurious[-which(spurious$V4 %in% alignments[which(duplicated(alignments$V4)),V4]),]
spurious_only_GTF <- data.frame(spurious_only$V1,"N16961","insertions",spurious_only$V2-100,spurious_only$V2+100,".",spurious_only$V6,".","bla")
write.table(spurious_only_GTF,"../results_GCF_000006745_179cass/spurious_ONLY_insertion_sites_Q30_100.gtf",col.names = FALSE,row.names = FALSE,quote=FALSE,sep="\t")

########################################################
###########           EXTRA Steps        ###############
########################################################

#check any reorder
multimap <- multimap[multimap$V2=="Assigned" & multimap$V3>1,]

col_number <- max(str_count(multimap$V4, ",") + 1)
kk <- multimap %>% separate(V4, 
                      into = paste0("feature", seq_len(col_number)), 
                      sep = ",")
synteny <- factor(table_counts_final$Feature_name[as.numeric(table_counts_final$Start)>309750 & as.numeric(table_counts_final$Stop)<435387  & table_counts_final$Chromosome=="Chr_2"])
captured <- c()
switches <- c()
for (i in 1:nrow(multimap)){
  captured[i] <- ordered(kk[i,4:15],levels=synteny)
  switches[i] <- is.ordered(captured)
}
table(switches)
#no reorderings

#correlation plots
N16961_superIn <- fread("../integron_finder_results/Results_Integron_Finder_GCF_000006745.1_ASM674v1_genomic_e4/GCF_000006745.1_ASM674v1_genomic.integrons", skip=2)
N16961_superIn <- N16961_superIn[N16961_superIn$type_elt=="attC",]
names(N16961_superIn)[c(3,4,5,7,8)] <- c("element","Start","Stop","evalue","Type")

kk <- table_counts_final[as.numeric(table_counts_final$Start)>309750 & as.numeric(table_counts_final$Stop)<435388 & table_counts_final$Chromosome=="NC_002506.1",]
kk <- kk[kk$Type=="cassette" | kk$Type=="attC",]
kkk <- merge(kk[,c(2,3,4,6,7)],N16961_superIn[,c(3,4,5,7,8)], by=c("Start","Stop"),all.x=TRUE)
kkk$length <- kkk$Stop - kkk$Start
#attC counts vs evalue
ggscatter(kkk[kkk$Type.x=="attC",], x = "Counts", y = "evalue", color = "#CC0000", shape = 21, size = 2,
          add = "reg.line", add.params = list(color = "#CC0000", size = 2), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman") + yscale("log10", .format = TRUE) + 
  xlab("AttC counts") + ylab("AttC e-value")
ggsave("../results_GCF_000006745_179cass/cor_attC_evalue_count.pdf",  width=5,height=5)

#cassette counts vs length
ggscatter(kkk[kkk$Type.x=="cassette",], x = "Counts", y = "length", color = "#CC0000", shape = 21, size = 2,
          add = "reg.line", add.params = list(color = "#CC0000", size = 2), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman") + 
  xlab("Cassette counts") + ylab("Cassette length")
ggsave("../results_GCF_000006745_179cass/cor_cassette_length_count.pdf",  width=5,height=5)



## plots with evalue fill gradient of previous attC
#change NA evalues in manual attC
kkk <- kkk[-nrow(kkk),]
eval <- log10(ifelse(is.na(kkk[kkk$Type.x=="attC",7]),10,kkk[kkk$Type.x=="attC",7]))
eval <- c(1,eval)
cassettes2 <- ggplot(data=kkk[kkk$Type.x=="cassette",]) +
  geom_point(aes(x=reorder(Feature_name,Start),y=as.numeric(Counts), colour="black", fill=eval),shape=21,size=3) +
  theme_bw() + scale_fill_gradient2() +
  scale_x_discrete(expand = c(0.01, 0),breaks = table_counts_final[table_counts_final$Type=="cassette",]$Feature_name[rep(c(F,F,F,T))]) +
  theme(axis.text.x=element_text(angle=90,hjust = 1),axis.text.y=element_text(size=12,face="bold"), legend.position = "none") +
  ylab("Counts") + xlab("") + labs(color="")

ggsave("../results_GCF_000006745_179cass/counts_cassettes_integron_evalues.pdf", plot=cassettes2, width=6,height=5)
ggsave("../results_GCF_000006745_179cass/counts_cassettes_integron_log10_evalues.pdf", plot=cassettes2+ scale_y_log10()+ylab("Counts"), width=6,height=5)

