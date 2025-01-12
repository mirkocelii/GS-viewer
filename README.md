a Tool for visualizing and validating Genome Scaffolds

Genome Scaffold Viewer (GS-viewer)  is R tool that can visualize and compare in a few minutes several scaffolds created by RagTag Scaffold

Determine the best scaffold according to:
- scaffold length
- number of gaps
- percentage of reference genome covered

It will also help for:
- Determining the optimal minimal contig length threshold
- Marking duplicated contigs
- Marking overlapping contigs that need to be manually edited
- Identifing contigs that could be used for gap-filling 
- generating repeat profile (including Telomere, centromeres, rDNA)
- rDNA identification
- Ploidy estimation
- Haplotype and Subgenome separation
- Coverage analysis
- Integrate Hic and Optical maps scaffold (in progress)

Modules:
 - scaffold
 - repeat-profile
 - rDNA-finder
 - coverage-plot
 - ploidy-level
 - busco-duplication

Dependencies:
- R
- R-data-table
- R-seqinr
- ragtag
- samtoolds
- bedtools
- seqtk
- blast
- repeatmasker
- trf

