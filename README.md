<img width="944" alt="image" src="https://github.com/user-attachments/assets/45042d17-4e32-4f0a-9200-c7375c60c83a" /># GS-viewer
a Tool for visualizing and validating Genome Scaffolds

Genome Scaffold Viewer (GS-viewer)  is R tool that can visualize and compare in a few minutes several scaffolds created by RagTag Scaffold

Determine the best scaffold according to:
      scaffold length
      number of gaps
      percentage of reference genome covered

- Determine the optimal minimal contig length threshold
- Mark duplicated contigs
- Mark overlapping contigs that need to be manually edited
- Identify contigs that could be used for gap-filling 
- Repeat Profile (including Telomere, centromeres, rDNA)
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

