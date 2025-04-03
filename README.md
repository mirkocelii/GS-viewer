
# GS-viewer
<p align="center">
  <img src="GSviewer logo.png" alt="GS-viewer logo" width="300"/>
</p>

GS-viewer is a modular toolkit designed for the visualization and analysis of genome assembly outputs. It encompasses various functionalities, including:

- **RagTag Scaffold Visualization**:
  - Helps to determine the best assembly and scaffolding parameters according to scaffold length, number of gaps and percentage of reference  covered.
  - Determine the optimal minimal contig length threshold.
  - Identify duplicated contigs in the assembly
  - Identify overlapping contigs that need to be manually edited
  - Identify contigs that could be used for gap-filling

- **rDNA Loci Identification**: 
  - Detect rRNA 45S, 18S, 5S and 5.8S using a reference rDNA library.

- **Repeat Profiling**:
  - Generate complete repeat profile by combining Repeat Masker and Tandem Repeat Finder outputs.
  - Identify Telomeric sequences
  - Identify rDNA loci (if added to the repeat library)
  - Plot gene density from an external .gff 
  - Plot costum data

- **Coverage Plotting**: 
  - Generate contig coverage plots extrapolated from Hifiasm `.gfa` files, in order to flag potential anomalies

- **Ploidy Level Estimation**: 
  - Assess ploidy based on RagTag output files.
  - Helpful in separating haplotypes and subgenomes

- **BUSCO Exact Duplication Analysis**:
  - Quantifies and plot the exact copy number of each duplicated gene from BUSCO output files.

## Getting Started

To utilize GS-viewer, ensure that the following dependencies are installed:

### R Packages:
- `data.table`
- `ggplot2`
- `seqinr`

### External Tools:
- `samtools`
- `bedtools`
- `seqtk`
- `blast`
- `RepeatMasker`
- `RagTag`
- `trf`

After installing the necessary packages and tools, you can execute GS-viewer commands as follows:

### USAGE:
```bash
Rscript GS-viewer.R <module> [options]
``` 
### GS-viewer Modules:
 - scaffold
 - repeat-profile
 - rDNA-finder
 - coverage-plot
 - ploidy-level
 - busco-duplication

### GS-viewer.R scaffold 
```bash
Identify rDNA sequences from a given genome using rDNA sequences from a close species

USAGE: Rscript GS-viewer.r scaffold --ref=<reference.fa> [ options ]

 positional arguments:
   --ref=              |  REQUIRED            [ .fasta or .fa file       ]  reference genome fasta file used with Ragtag. It must be the same for all scaffold loaded
   --grep=             |  no default          [ comma separated strings  ]  greps folders containing ragtag.scaffold.fasta file, linux regular expression allowed,
   --skip=             |  no default          [ comma separated strings  ]  skips folders containing these strings
   --path=             |  no default          [ linux path               ]  provide path of 1 or more folders, regular expression allowed, folders must be at the same level,--path will ignore --grep

 output options:
   --sort=             |  default TRUE         [ TRUE/FALSE              ]  sort samples names by character or numeric suffix/prefix
   --chr.filt=         |  default 0            [ integer                 ]  minimal Reference Chromosome size plotted. Numeric and Kb,Mb character accepted, e.g. 1Mb, 500Kb, 1000
   --organelles=       |  default FALSE        [ TRUE/FALSE              ]  include >organelles independetly of the chr.filt
   --chr.rename=       |  default FALSE        [ TRUE/FALSE              ]  replace GenBank identifier with Chromosome or Contig name from the fasta header
   --out.suffix=       |  no default           [ string                  ]  output file suffix
   --out.dir=          |  def. GS-viewer_Plots [ string                  ]  output directory name
   --rename=           |  default nothing      [ string                  ]  if only one scaffold is selected, it can be renamed
   --name.del=         |  default nothing      [ string                  ]  strings to be deleted from all names, comma separated
   --plot.h=           |  default 18           [ integer                 ]  Height in inches of PDF output files
   --plot.w=           |  default 35           [ integer                 ]  Width in inches of PDF output files
   --plot.sep=         |  default FALSE        [ TRUE/FALSE              ]  if TRUE, separate PDF for each sample Scaffold Length will be created
   --plot.mix=         |  default TRUE         [ TRUE/FALSE              ]  if TRUE, create also a PDF with both Coverage and Length Plots of all samples
   --short.name=       |  default TRUE         [ TRUE/FALSE              ]  if YES common substrings are removed from samples names
   --CTG.min=          |  default 0            [ integer                 ]  minimal Ctg alignment length to be plotted on the Reference coverage Plot
   --MIN.WIND=         |  default 5000         [ integer                 ]  minimal Ctg alignment length to be plotted on the Single Chr coverage Plot
   --colors=           |  default 11 colors    [ comma separated R colors]  default firebrick3,cornflowerblue,gold,blueviolet,darkturquoise,darkorange3,lightcoral,royalblue4,forestgreen,deepskyblue,mediumorchid1,then rainbow() will be used
   --ctg.gap=          |  default FALSE        [ TRUE/FALSE              ]  identify contigs of secondary scaffold overlapping with gaps in the main/best scaffold
   --main=             |  default nothing      [ string                  ]  string to indicate the main scaffold. If absent, the Scaffold with best compromise between low number of gaps and higher coverage of reference genome is selected
   --ctg.fa=           |  default TRUE         [ TRUE/FALSE              ]  create fasta files with unplaced contigs of each scaffold and with contigs that overlap with gaps on the main Scaffold

 REQUIRED FILES:
  - reference file
  - >=1 Ragtag_scaffold folders in the working directory

 Modules required:
  - R data.table package
  - R seqinr package  [ for Bionano data only ]
  - samtools
  - bedtools
  - seqtk [for Gap filling fasta files only]
```

### GS-viewer.R rDNA-finder 
```bash

 Identify rDNA sequences from a given genome using rDNA sequences from a close species

 USAGE: Rscript GS-viewer.R rDNA-finder --ref.rDNA=plant_rDNA.fasta. --genome=my_assembly.fasta  --out=my_species_rDNA.fasta

 Arguments:
   --ref.rDNA=            | rDNA fasta file of a close species ]
   --genome=              | genome
   --out=                 | output file name [ default = genome.rDNA.fasta ]

 blast parameters:
   --num_threads=         |  num_threads   [ default = 32      ]
   --perc_identity=       |  perc_identity [ default = 70      ]
   --evalue=              |  evalue        [ default = 0.000001]

 Modules required:
  - R
  - R packages: data.table
  - R packages: seqinr
  - blast
  - samtools
  - bedtools

you can find rDNA sequences on NCBI or SILVA [ https://www.arb-silva.de/ ]

 make sure you rDNA sequences have comparable values with the following:
     --> rDNA 18S  1800 bp
     --> rDNA 28S  3500 bp
     --> rDNA 5.8S  155 bp
     --> rDNA 5S    120 bp
```

### GS-viewer.R repeat-profile 
```bash

 Indentify and visualize Telomeres, Tandem repeats, rDNA and transposable elements in a genome

 USAGE: Rscript GS-viewer.R repeat-profile --fasta=<query_scaffold.fasta> --lib=<query_TE_lib.fasta> [ options ]

 NOTE: Before the scritp it will run Tandem Repeat Finder and Repeat Masker
  --> e.g. :
    # trf <query_scaffold.fasta>  2 3 5 80 10 30 2000 -l 6 -dat -h  :
    # RepeatMasker -qq -no_is --pa 64 -lib <query_TE_lib.fasta> <query_scaffold.fasta> -gff -xsmall   # Do not use -nolow

 mandatory arguments:
   --fasta=            |  scaffold of assembly .fasta or .fa file ]
   --lib=              |  TE library, possibily including rDNA

 add external rDNA-only repeat profile:
   --rDNA.fa =         |  rDNA library, if not included in TE library]
   --rDNA.out=         |  repeat masker out with rDNA library        ]

  other parameters:
   --sort.size=        |  sort by size                                                [ default = FALSE  ]
   --tel=              |  identify a specififc telomeric pattern                      [ default = absent ]
   --diploid=          |  group chromosomes by name ( e.g Chr1.1 and Chr1.2  )        [ default = absent ]
   --polyploid=        |  group chromosomes by name ( e.g Chr1.A and Chr1.B  )        [ default = absent ]
   --ref.name=         |  replace GenBank IDs with chromosome names from reference.fna[ default = absent ]
   --gene.gff=         |  gene prediction in gff3 format                              [ default = absent ]
   --add.dots=         |  custom.txt file --> Chr, pos,      , name, color, pch, size [ default = absent ]
   --add.segm=         |  custom.bed file --> Chr, start, end, name, color            [ default = absent ]
   --chr.filt=         |  Integer, minimal Chr size plotted                           [ default = 1Mb / 1000Kb / 1000000 ; numeric, Kb,Mb strings accepted  ]
   --sample.name=      |  Sample name to be plotted                                   [ default = main folder name  ]
   --plot.h=           |  Hight in inches of PDF output files                         [ default 18       ]
   --plot.w=           |  Width in inches of PDF output files                         [ default 35       ]
   --name.space=       |  space on the left for chromosome names [ inches ]           [ default 6        ]
   --leg.size=         |  increase (>0) or decrease (<0) legend size by               [ default 1        ]

  TRF parameters :
   --match=            |  matching weight                     [ default = 2    ]
   --mismatch=         |  mismatching penalty                 [ default = 3    ]
   --delta=            |  indel penalty                       [ default = 5    ]
   --pm=               |  match probability (whole number)    [ default = 80   ]
   --pi=               |  indel probability (whole number)    [ default = 10   ]
   --minscore=         |  minimum alignment score to report   [ default = 30   ]
   --maxperiod=        |  maximum period size to report       [ default = 2000 ]
   --l=                |  matching weight                     [ default = 6    ]
   --consensus.size=   |  consensus.size in bp                [ default = 7    ]

 Modules required:
  - R
  - R packages: data.table
  - trf
  - samtools
  - bedtools
  - Repeat Masker
  - RagTag 
```

### GS-viewer.R busco-duplication 
```bash

 Diplace busco scorere and exact duplication percentages

 USAGE: GS-viewer.R busco-duplication --grep=Busco_folder

 Arguments:
   --grep=             | grep folders containing mysummary.txt          [ default all folders having full_table.tsv inside ]   linux regular expression allowed,
   --skip=             | skip folders containing these strings          [ default absent  ]
   --main.dir=         | main directory containing busco folders        [ default current  ]
   --thr.plot=         | min. n° of dupl. genes to be writtin o barplot [ default 30      ]
   --out=              | output file suffix                             [ default Busco_plot ]
   --plot.h=           | Hight in inches of PDF output files            [ default 18      ]
   --plot.w=           | Width in inches of PDF output files            [ default 35      ]

 Modules required:
  - R packages: data.table

```

### GS-viewer.R coverage-plot 
```bash

 Plot contig coverage from HIFIASM .gfa file

 USAGE: Rscript GS-viewer.R coverage-plot --gfa=my_hifiasm_assembly.fasta.bp.p_ctg.noseq.gfa --agp=my_scaffold.agp

 Arguments:
   --gfa=              | *noseq.gfa from HIFIASM
   --agp=              | Scaffold AGP file                     [ optional, default absent ] if absentm, contigs will be plot one by one
   --ref=              | Reference FASTA file                  [ optional, default absent ] used for importing real Chr names rather than GenBank IDs 
   --chr.filt=         | default = 1Mb                         [ Integer                  ] minimal Chr size plotted [ default = 1Mb / 1000Kb / 1000000 ; numeric-only and Kb,Mb character accepted  ]
   --out=              | no default                            [ string                   ] output file suffix
   --wind=             | default 100Kb/100000)                 [ string                   ] window size [ default = 1000Kb / 1000000 ; numeric-only, Mb,Kb,bp character accepted  ]
   --plot.h=           | Height in inches of PDF output files  [ default 18               ]
   --plot.w=           | Width in inches of PDF output files   [ default 35               ]

Modules required:
  - R
  - R packages: data.table
  - samtools
  - bedtools

```

### GS-viewer.R ploidy-level 
```bash

 Plot ploidy level from a RagTag scaffold

 USAGE: Rscript GS-viewer.R ploidy-level --scaffold=ragtag.scaffold.fasta --ref=reference_genome.fasta

 Arguments:
   --scaffold=         |  ragtag.scaffold.fasta file
   --ref=              |  Reference FASTA file                           [ optional, default absent   used for importing real Chr names   ]
   --overlap=          |  Max contig overlap for beein at the same level [ 0-to-1, default 0.02                                           ]
   --ploidy=           |  Set ploidy level                               [ default automatic determined with 90-95% of reference covered  ]
   --chr.rename=       |  Replace GenBank identifier with Chr/Contig     [ default FALSE                                                  ]
   --chr.filt=         |  minimal Chr size plotted                       [ default = 1Mb/1000Kb/1000000; numeric,Kb,Mb character accepted ]
   --hap.fa=           |  output fasta file of each haplotype set draft  [ default absent   ]
   --plot.h=           |  Height in inches of PDF output files           [ default 18       ]
   --plot.w=           |  Width in inches of PDF output files            [ default 35       ]

Modules required:
  - R
  - R packages: data.table
  - samtools
  - bedtools
