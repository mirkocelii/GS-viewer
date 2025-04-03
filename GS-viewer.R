
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    cat("No arguments provided. \n\n  Usage: Rscript GS-viewer.R <command> [options] \n")
    cat("\n commands:\n")
    cat("    scaffold         \t  # Visualize RagTag Scaffold(s) \n")
    cat("    repeat-profile    \t # Identify and visualize Repeats in an assembly/scaffold \n")
    cat("    rDNA-finder       \t # Identify rDNA loci using a reference rDNA library \n")
    cat("    coverage-plot     \t # Plot Contig coverage from Hifiasm .gfa file \n")
    cat("    ploidy-level      \t # Estimates the ploidy of an assembly based on RagTag output \n")
    cat("    busco-duplication \t # Quantifies the duplication level based on BuSCo output \n")
    
    cat("\n Dependencies:\n")
    cat("   - R (data.table, seqinr) \n")
    cat("   - samtools \n")
    cat("   - bedtools \n")
    cat("   - seqtk \n")
    cat("   - blast \n")
    cat("   - repeatmasker \n")
    cat("   - ragtag \n")
    cat("   - trf \n")
    
    cat("\n options:\n")
    cat("   -v, version\n")
    cat("   -c, citation \n")


stop()


}

command <- args[1]

switch(command,
       "scaffold" = source("Ragtag_visualization_v10.r"),
       "repeat-profile" = source("Telomere_finder.R"),
       "rDNA-finder" = source("rDNA_finder.R"),
       "coverage-plot" = source("Coverage_plot.R"),
       "ploidy-level" = source("Ploidy_levels.R"),
       "busco-duplication" = source("Plot_Busco_dup.R"),
       stop(paste("Command not found:", command))
)
