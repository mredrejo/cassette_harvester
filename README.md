# Integron Capturer
This repo contains the code for the analysis of long reads sequencing of a new biotechnological tool to capture integron cassettes. This tool has been developed by the [MBA lab](https://ucm.es/mbalab) and used as a proof of concept with *Vibrio* super-integron (SI). Full details in a new paper currently under preparation.

Software packages for analysis were installed either from Conda, CRAN or Bioconductor. 


Here, we used the *Vibrio cholerae* O1 biovar El Tor str. N16961 genome assembly [GCF_000006745.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000006745.1/). The genome was reannotated with [Bakta 1.9.2](https://github.com/oschwengers/bakta) (DB v.5.1) and [IntegronFinder 2.0.5](https://github.com/gem-pasteur/Integron_Finder). Two more attC described in the literature were added manually. The final SI contains 179 cassettes.
Raw reads were analyzed with [Nanoplot 1.42.0](https://github.com/wdecoster/NanoPlot) and processed with [Porechop 0.2.4](https://github.com/rrwick/Porechop), and [NanoFilt 2.8.0](https://github.com/wdecoster/nanofilt)  (`-q 10`). 

**QC output from final reads is available [here](results/nanoplot_trimmed/NanoPlot-report.html).**


Then, reads were mapped against the reference genome using [Minimap 2.28-r1209](https://github.com/lh3/minimap2) (`-ax map-ont`). Mapping coverage per nucleotide was analyzed using [BamDash 0.2.4](https://github.com/jonas-fuchs/BAMdash) and is available as interactive plots (note the different Y scale) : [**Chromosome 1**](https://www2.iib.uam.es/mredrejo_lab/plots_cov_nt/NC_002505.1_plot.html)  and [**Chromosome 2**](https://www2.iib.uam.es/mredrejo_lab/plots_cov_nt/NC_002506.1_plot.html).

Mapped reads were transformed in *Counts* per feature (CDS, attC or cassette) using the R package `Rsubread v. 2.18.0`. Code is available in the script [superintegron_mapping_GCF_000006745_179curated.R](superintegron_mapping_GCF_000006745_179curated.R). 

The final [Counts/feature table](results/coverage_N16961_IF2_manual.xlsx) and plots are available in the folder [**results**](results). The folder also contains some interactive plots, like the [counts/cassette](results/counts_cassettes_fill_evalue.html) (also in Log10 scale [here](results/counts_cassettes_log10_fill_evalue.html)).
