# Cassette harvester

This repo contains the code for analyzing long-read sequencing of a new biotechnological tool for capturing Integron cassettes. This tool was developed by [MBA lab](https://ucm.es/mbalab) and used as a proof of concept with *Vibrio* Super-Integron (SI). Further details can be found in a new paper that is currently being prepared.

The software packages for the analysis were installed either by Conda, CRAN or Bioconductor. R code for generation of final table and plots is available in the script [superintegron_mapping_GCF_000006745_179curated.R](superintegron_mapping_GCF_000006745_179curated.R). 

Briefly, we used here the *Vibrio cholerae* O1 biovar El Tor str. N16961 genome assembly [GCF_000006745.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000006745.1/). The genome was re-annotated using [Bakta 1.9.2](https://github.com/oschwengers/bakta) (DB v.5.1) and [IntegronFinder 2.0.5](https://github.com/gem-pasteur/Integron_Finder). Two additional attC described in the literature were added manually. The final SI contains 179 annotated cassettes.

The nanopore raw reads were analyzed with [Nanoplot 1.42.0](https://github.com/wdecoster/NanoPlot) and processed with [Porechop 0.2.4](https://github.com/rrwick/Porechop), and [NanoFilt 2.8.0](https://github.com/wdecoster/nanofilt)  (`-q 10`). 

**QC output from final reads is available [here](https://mredrejo.github.io/cassette_harvester/results/nanoplot_trimmed/NanoPlot-report.html).**


Teads were then mapped against the reference genome using [Minimap 2.28-r1209](https://github.com/lh3/minimap2) (`-ax map-ont`).  Mapping coverage per nucleotide was analyzed with [Samtools 1.16](https://github.com/samtools/samtools) and plotted with the R package `circlize` to obtain a **circular [coverage plot](https://mredrejo.github.io/cassette_harvester/results/circos_cov.pdf)**. More detailed coverage plots were obtained with [BamDash 0.2.4](https://github.com/jonas-fuchs/BAMdash), and are available as interactive plots (note the different Y scale) : [**Chromosome 1**](https://www2.iib.uam.es/mredrejo_lab/plots_cov_nt/NC_002505.1_plot.html)  and [**Chromosome 2**](https://www2.iib.uam.es/mredrejo_lab/plots_cov_nt/NC_002506.1_plot.html).

Mapped reads were transformed in *Counts* per feature (CDS, attC or cassette) using the R package `Rsubread v. 2.18.0`. **The final [Counts/feature table](results/coverage_N16961_IF2_manual.xlsx) and plots are available in the [**results**](results) folder. This folder also contains some interactive plots, like the [counts/cassette](https://mredrejo.github.io/cassette_harvester/results/counts_cassettes_fill_evalue.html) (also in Log10 scale [here](https://mredrejo.github.io/cassette_harvester/results/counts_cassettes_log10_fill_evalue.html)).**
