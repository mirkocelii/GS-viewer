#!/bin/R

#R

args <- commandArgs(TRUE)


# HELP SECTION

   if("--help" %in% args) 
    {
        cat('\n')
        cat(
        "Visualization of Scaffold(s) created with Ragtag \n",
        " \n",
        "USAGE: Rscript GS-viewer.r scaffold --ref=<reference.fa> [ options ] \n",
        " \n",
        "positional arguments:\n",
        "  --ref=              |  REQUIRED            [ .fasta or .fa file       ]  reference genome fasta file used with Ragtag. It must be the same for all scaffold loaded","\n",
        "  --grep=             |  no default          [ comma separated strings  ]  grep folders containing ragtag.scaffold.fasta file, linux regular expression allowed,  ","\n",
        "  --skip=             |  no default          [ comma separated strings  ]  skip folders containing these strings ","\n",
        "  --path=             |  no default          [ linux path               ]  provide path of 1 or more folders, regular expression allowed, folders must be at the same level, using --path will ignore --grep","\n",
       " \n",
        "optional inputs:\n",
 #       "  --rm.ref=           |  default FALSE     [ RepeatMasker.gff file      ]  Repeat Masker profile of reference ","\n",
 #       "  --rm.scaf=          |  default FALSE     [ RepeatMasker.gff file      ]  Repeat Masker profile of scaffolds, must be in the same folder ","\n",
        "  --agp.bionano=      |  default FALSE        [ bionano map agp  file      ]  Bionano Map AGP file of scaffolds, must be in the same folder ","\n", 
        " \n",
        "output options:\n",
        "  --sort=             |  default TRUE         [ TRUE/FALSE              ]  sort samples names by character or numeric suffix/prefix ","\n",
        "  --organelles=       |  default FALSE        [ TRUE/FALSE              ]  include >organelles independetly of the chr.filt  ","\n",
        "  --chr.filt=         |  default 0            [ integer                 ]  minimal Reference Chromosome size plotted. Numeric and Kb,Mb character accepted, e.g. 1Mb, 500Kb, 1000             ","\n",
        "  --chr.rename=       |  default FALSE        [ TRUE/FALSE              ]  replace GenBank identifier with Chromosome or Contig name from the fasta header  ","\n",
   #     "  --not.rename=       |  default FALSE        [ TRUE/FALSE              ]  do not replace GenBank identifier with Chromosome or Contig name from the fasta header  ","\n",
        "  --out.suffix=       |  no default           [ string                  ]  output file suffix             ","\n",
        "  --out.dir=          |  def. GS-viewer_Plots [ string                  ]  output directory name          ","\n",
        "  --rename=           |  default nothing      [ string                  ]  if only one scaffold is selected, it can be renamed", "\n",
        "  --name.del=         |  default nothing      [ string                  ]  strings to be deleted from all names, comma separated", "\n",
        "  --rm.unpl.ctg=      |  default TRUE         [ TRUE/FALSE              ]  if TRUE unplaced contigs and ch0 will be removed from the scaffold","\n",
        "  --chr0_plot=        |  default FALSE        [ TRUE/FALSE              ]  if TRUE unplaced contigs will be plotted as a merged Chromosome, but not added in output file","\n",
        "  --plot.conf=        |  default TRUE         [ TRUE/FALSE              ]  if TRUE will add the location confidence to the barplot","\n",
        "  --thr.conf=         |  default 0.25         [ number from 0 to 1      ]  if plot.conf=TRUE, it's the threshold below which location confidence is plotted","\n",
        "  --plot.h=           |  default 18           [ integer                 ]  Height in inches of PDF output files","\n",
        "  --plot.w=           |  default 35           [ integer                 ]  Width in inches of PDF output files","\n",
        "  --plot.sep=         |  default FALSE        [ TRUE/FALSE              ]  if TRUE, separate PDF for each sample Scaffold Length will be created","\n",
        "  --plot.mix=         |  default TRUE         [ TRUE/FALSE              ]  if TRUE, create also a PDF with both Coverage and Length Plots of all samples","\n",
        "  --short.name=       |  default TRUE         [ TRUE/FALSE              ]  if YES common substrings are removed from samples names","\n",
        "  --CTG.min=          |  default 0            [ integer                 ]  minimal Ctg alignment length to be plotted on the Reference coverage Plot","\n",
        "  --MIN.WIND=         |  default 5000         [ integer                 ]  minimal Ctg alignment length to be plotted on the Single Chr coverage Plot","\n",
        "  --colors=           |  default 11 colors    [ comma separated R colors]  default = firebrick3,cornflowerblue,gold,blueviolet,darkturquoise,darkorange3,lightcoral,royalblue4,forestgreen,deepskyblue,mediumorchid1 if >11 samples, from 11th color rainbow() function will be used", "\n",
        "  --ctg.gap=          |  default FALSE        [ TRUE/FALSE              ]  identify contigs of secondary scaffold overlapping with gaps in the main/best scaffold", "\n",
        "  --main=             |  default nothing      [ string                  ]  string to indicate the main scaffold. If absent, the Scaffold with best compromise between low number of gaps and higher coverage of reference genome is selected", "\n",
        "  --ctg.fa=           |  default TRUE         [ TRUE/FALSE              ]  create fasta files with unplaced contigs of each scaffold and with contigs that overlap with gaps on the main Scaffold", "\n",
        "\n",
        "REQUIRED FILES: \n",
        " - reference file \n",
        " - >=1 Ragtag_scaffold folders in the working directory \n",
       "\n",
        "Modules required: \n",
        " - R data.table package \n",
        " - R seqinr package  [ for Bionano data only ]\n",
        " - samtools \n",
        " - bedtools \n",
        " - seqtk [for Gap filling fasta files only] ",
       "\n"
        )

        cat('\n')
        q(save="no")
    }

#
# BUG : if 1 fai has to be produced for the first time, the script will block because of missing different length of:
# Error in data.frame(CHR = chr.files, AGP = agp.files, CNF = cnf.files,  :
#   arguments imply differing number of rows: 5, 4
# Execution halted

cat('\n')
cat('\n')
cat("############################################################################\n")
cat("############################################################################\n")
cat("###                             GS-viewer.R                             ###\n")
cat("############################################################################\n")
cat("############################################################################\n")

START.time=Sys.time()



cat('\n')
cat("############################################################################\n")
cat('-------------------------------OPTIONS--------------------------------------\n')
cat("############################################################################\n")
if (length(args)>0 ) print(as.data.frame(args))
cat('\n')


# COMMANDLINE ARGS

    if (length(args)==0) ARG = data.frame("X1"="--option", "X2"="argument")
    if (length(args)>0 ) ARG = data.frame(do.call(rbind,strsplit(args,"=")))
    
cat('\n')
cat("############################################################################\n")
cat('-------------------------------OPTIONS check--------------------------------\n')
cat("############################################################################\n")
cat('\n')
if (length(args)==0) cat("--> all default \n" )
if (length(args)>0 ) print(ARG, justify = 3 ) 
    
	ref.file     = as.character(ARG[ ARG$X1=="--ref"         ,]$"X2")
	grep         = as.character(ARG[ ARG$X1=="--grep"        ,]$"X2")
	skip         = as.character(ARG[ ARG$X1=="--skip"        ,]$"X2")
	path         = as.character(ARG[ ARG$X1=="--path"        ,]$"X2")
	chr.filt     = as.character(ARG[ ARG$X1=="--chr.filt"    ,]$"X2")
	chr.rename   = as.character(ARG[ ARG$X1=="--chr.rename"  ,]$"X2")
#	not.rename   = as.character(ARG[ ARG$X1=="--not.rename"  ,]$"X2")
	out.suffix   = as.character(ARG[ ARG$X1=="--out.suffix"  ,]$"X2")
	rm.unpl_ctg  = as.character(ARG[ ARG$X1=="--rm.unpl.ctg" ,]$"X2")
	chr0_plot    = as.character(ARG[ ARG$X1=="--chr0_plot"   ,]$"X2")
	plot.conf    = as.character(ARG[ ARG$X1=="--plot.conf"   ,]$"X2")
	thr.conf     = as.character(ARG[ ARG$X1=="--thr.conf"    ,]$"X2")
	plot.h       = as.character(ARG[ ARG$X1=="--plot.h"      ,]$"X2")
	plot.w       = as.character(ARG[ ARG$X1=="--plot.w"      ,]$"X2")
	plot.sep     = as.character(ARG[ ARG$X1=="--plot.sep"    ,]$"X2")
	plot.mix     = as.character(ARG[ ARG$X1=="--plot.mix"    ,]$"X2")
	short.name   = as.character(ARG[ ARG$X1=="--short.name"  ,]$"X2")
	out.dir      = as.character(ARG[ ARG$X1=="--out.dir"     ,]$"X2")
	CTG.min      = as.character(ARG[ ARG$X1=="--CTG.min"     ,]$"X2")
    MIN.WIND     = as.character(ARG[ ARG$X1=="--MIN.WIND"    ,]$"X2")
    COL.LIST     = as.character(ARG[ ARG$X1=="--colors"      ,]$"X2")
    ctg.gap      = as.character(ARG[ ARG$X1=="--ctg.gap"     ,]$"X2")
    ctg.fa       = as.character(ARG[ ARG$X1=="--ctg.fa"      ,]$"X2")
    main         = as.character(ARG[ ARG$X1=="--main"        ,]$"X2")
# 	rm.ref       = as.character(ARG[ ARG$X1=="--rm.ref"      ,]$"X2")
# 	rm.scaf      = as.character(ARG[ ARG$X1=="--rm.scaf"     ,]$"X2")
	bionano      = as.character(ARG[ ARG$X1=="--agp.bionano" ,]$"X2")
	rename       = as.character(ARG[ ARG$X1=="--rename"      ,]$"X2")
	name.del     = as.character(ARG[ ARG$X1=="--name.del"    ,]$"X2")
	sort.sam     = as.character(ARG[ ARG$X1=="--sort"        ,]$"X2")
	organelles   = as.character(ARG[ ARG$X1=="--organelles"  ,]$"X2")

#  DEFAULT VALUES

#   if ( length(ref.file   )==0  ) ref.file = "AjwaAllChr.fasta.fai"
    if ( length(grep       )==0  ) grep = NA  
    if ( length(skip       )==0  ) skip = NA 
    if ( length(path       )==0  ) path = NA 
    if ( length(chr.filt   )==0  ) chr.filt = "0Mb"
    if ( length(chr.rename )==0  ) chr.rename = FALSE
#    if ( length(not.rename )==0  ) not.rename = FALSE
    if ( length(rm.unpl_ctg)==0  ) rm.unpl_ctg = TRUE
    if ( length(chr0_plot)==0    ) chr0_plot   = TRUE
    if ( length(plot.conf)==0    ) plot.conf   = FALSE
    if ( length(thr.conf)==0     ) thr.conf    = 0.25
    if ( length(plot.h)==0       ) plot.h      = 18
    if ( length(plot.w)==0       ) plot.w      = 35
    if ( length(plot.sep)==0     ) plot.sep    = FALSE
    if ( length(plot.mix)==0     ) plot.mix    = FALSE
    if ( length(short.name)==0   ) short.name  = TRUE
    if ( length(CTG.min)==0      ) CTG.min     = 0
    if ( length(MIN.WIND)==0     ) MIN.WIND    = 5000
    if ( length(COL.LIST)==0     ) COL.LIST    = c("firebrick3","cornflowerblue","gold","blueviolet","darkturquoise","darkorange3","lightcoral","royalblue4","forestgreen","deepskyblue","mediumorchid1") 
    if ( length(ctg.gap)==0      ) ctg.gap     = FALSE
    if ( length(ctg.fa)==0       ) ctg.fa      = TRUE
#     if ( length(rm.ref)==0       ) rm.ref      = NA
#     if ( length(rm.scaf)==0      ) rm.scaf     = NA
    if ( length(bionano)==0      ) bionano     = FALSE
    if ( length(rename)==0       ) rename      = FALSE
    if ( length(name.del )==0    ) name.del    = FALSE
    if ( length(sort.sam)==0     ) sort.sam    = TRUE
    if ( length(organelles)==0   ) organelles  = FALSE
  #  if ( length(main)==0         ) main        = FALSE # later in the script
  
    if ( length(out.suffix)==0  &  length(out.dir)>0   ) out.suffix  = gsub("GS-viewer_Plots.|Ragtag_Plot.","",out.dir)
    if ( length(out.suffix)>0   &  length(out.dir)==0  ) out.dir     = paste0("GS-viewer_Plots.",out.suffix)
    if ( length(out.suffix)==0  &  length(out.dir)==0  ) out.dir     = "GS-viewer_Plots"
    if ( length(out.suffix)==0  &  length(out.dir)==0  ) out.suffix  = NA 

    
    if( chr.rename %in% c("--chr.rename", "T","TRUE" )) chr.rename   = TRUE
#    if( not.rename %in% c( "--not.rename", "T","TRUE" )) not.rename   = TRUE else not.rename = FALSE
    if( plot.mix   %in% c( "--plot.mix", "T","TRUE" )   ) plot.mix     = TRUE
    if( plot.sep   %in% c( "--plot.sep", "T","TRUE" )   ) plot.sep     = TRUE
    if( plot.conf  %in% c( "--plot.conf", "T","TRUE" )  ) plot.conf    = TRUE
    if( chr0_plot  %in% c( "--chr0_plot", "T","TRUE" )  ) chr0_plot    = TRUE
    if(rm.unpl_ctg %in% c( "--rm.unpl_ctg", "T","TRUE" )) rm.unpl_ctg  = TRUE
    if( short.name %in% c( "--short.name", "T","TRUE" ) ) short.name   = TRUE
    if( ctg.gap    %in% c( "--ctg.gap", "T","TRUE" )    ) ctg.gap      = TRUE
    if( ctg.fa     %in% c( "--ctg.fa", "T","TRUE" )     ) ctg.fa       = TRUE

# as.numeric inputs
     plot.h   = as.numeric( plot.h )
     plot.w   = as.numeric( plot.w )
     CTG.min  = as.numeric( CTG.min )
     MIN.WIND = as.numeric( MIN.WIND )
     thr.conf = as.numeric( thr.conf )


# TEST PARAMETERS for interactive session
if( 2 ==1 )
{
    setwd("/ibex/scratch/projects/c2067/Rice_CIAT/RICE/Pool2/1158")
 #  setwd("/ibex/scratch/projects/c2067/celiim/Ziziphus_spina-christi")
 #   setwd("/ibex/scratch/projects/c2067/celiim/Ajwa_dp")
    setwd("/ibex/scratch/projects/c2067/celiim/Reziz_dp/")
  #  setwd("/ibex/scratch/projects/c2067/celiim")
   ref.file="/ibex/scratch/projects/c2067/celiim/references/AjwaAllChr.fasta"
    setwd("/ibex/scratch/projects/c2067/Rice_CIAT/RICE/Pool2/1160")
    setwd("/ibex/scratch/projects/c2067/celiim/")

   # setwd("/ibex/scratch/projects/c2067/celiim/Ajwa_dp")
   setwd("/ibex/scratch/projects/c2067/celiim/Ajwa_dp")
   ref.file= "/ibex/scratch/projects/c2067/celiim/references/Artificial_reference_nogaps2.fasta"
   ref.file= "/ibex/scratch/projects/c2067/celiim/references/Ajwa_220x_18seq.fasta"
   ref.file="/ibex/scratch/projects/c2067/celiim/references/Ajwa_220x_18chr_singleCTG18.fasta"

 # setwd("/ibex/scratch/projects/c2067/celiim/for_mirko") ; ref.file= "asclepias_syriaca_500k.fasta"
 
   
  # ref.file=("/home/celiim/Desktop/references/references.fasta")
  # setwd("/home/celiim/Desktop/DATEPALMS")
   # ref.file="/ibex/scratch/projects/c2067/Rice_CIAT/RICE/MH63RS3.fasta"
   
  # ref.file="Z_jujuba_Genome_ref.fasta"
  #  ref.file="/ibex/scratch/projects/c2067/Manjula/Reziz/references/Bahree_all_chromosomes.fa"
  # ref.file="/ibex/scratch/projects/c2067/Manjula/Reziz/references/AjwaAllChr.fasta"
   #ref.file="//ibex/scratch/projects/c2067/celiim/date_palm_artificial_ref_nogaps/Artificial_reference_nogaps.fasta"
   
 #  Rscript "$Ragtag_vis" --ref=$Arti --path=*_dp/rag*hifiasm.100*Arti*   --skip=.F  --out.dir=Ajwa.hifiasm_100Kb  --out.suffix=Ajwa.hifiasm_100Kb

#x="Ajwa"     ; cd $DP ; cd "$x"_dp ; sbatch -t 1:00:00 --mem=64GB -c 16 --export=Script=$Script,x=$x,Arti=$Arti -J "6.$x" --wrap 'Rscript "$Script" --ref="$Arti" --grep=manual3,clr,100Kb --skip=Art2,.F,Bah,manual2,fn --out.dir=Ragtag_Plot.Arti.manual_vs_all --out.suffix=Arti.manual_vs_all  --ctg.gap=TRUE --main=manual ' ;


#cd $DP ; sbatch -t 1:00:00 --mem=64GB -c 16 --export=Script=$Script,Ajw4a=$Ajw4a -J "6.$x" --wrap ' Rscript "$Script" --ref="$Ajw4a" --path=*_dp/ragtag_scaffold.*.hifiasm.100Kb.Ajw4a --skip=Bahree,Khalas --out.suffix=all_hifiasm.100Kb.Ajw4 ' ;


    grep           =  NA #"500" #"Ajwa.hifiasm.100Kb.Ajw3" #"100Kb.Art2.manual4,ref_ctgs,HIFIASM_MER"#"hifiasm.100Kb" #"hifiasm.100Kb" # ".100Kb" #"1Mb*Ajwa"  
  #  grep           =  "Ajwa.hifiasm.100Kb.Ajw3" ;     bionano        = TRUE
    skip           =  NA# "Art2,Arti" #".F,corr,less" #"merged,sm_,_hif"# ".F,.merged,sm_,_hif"
    path           =  NA; #"*ragtag_scaffold.*.hifiasm.100Kb.Ajw4a" # "*_dp/*ag*ag*aff*100*F*"
    rm.unpl_ctg    = TRUE                                    # optional, removes unplaced contigs or Chr0
    chr0_plot      = TRUE                                    # Visualize all unplaced contigs separated by GAPs
    plot.conf      = TRUE                                    # PLot confidence of contig placment on reference
    thr.conf       = 0.7                                     # PLot confidence of contig placment on reference
    plot.h         = 18
    plot.w         = 35
    plot.sep       = FALSE
    plot.mix       = FALSE
    short.name     = TRUE
    out.suffix     = "test"
    out.dir        = "test"
    CTG.min        = 0       
    MIN.WIND       = 5000
    COL.LIST       = c("firebrick3","cornflowerblue","gold","blueviolet","darkturquoise","darkorange3","lightcoral","royalblue4","forestgreen","deepskyblue","mediumorchid1") 
    ctg.gap        = TRUE
    ctg.fa         = TRUE
    main           = "220x"
    rm.ref         = NA
    rm.scaf        = NA
    chr.filt       = "5Mb"
    bionano        = FALSE
    rename         = FALSE
    organelles     = FALSE
    name.del       = FALSE
   }


if (  "data.table" %in% rownames(installed.packages()) == FALSE )  { cat('\n\n############################################################################\n######### Please install data.table package! \n    R\n    > install.packages("data.table") \n######### or\n    conda install r-data.table \n######### Then try again! \n############################################################################\n\n\n') }

library(data.table)
options(width=350)
options(scipen=999)

# chr.filt into numeric type
chr.filt = gsub(" ","",chr.filt) ; chr.filt = gsub("bp","",chr.filt) ; chr.filt = gsub("Kb|kb","000",chr.filt) ;  chr.filt = gsub("Mb|mb","000000",chr.filt) ;      chr.filt = as.numeric( chr.filt )



# Colors
COL.LIST = unlist(strsplit(COL.LIST,",")) 

# GREP AND SKIP processing

if ( is.na(grep)==TRUE  ) grep="*"

grep.split = unlist(strsplit(grep,","))
if( length(grep.split)==1) grep.words = grep
if( length(grep.split)>=2) grep.words = paste0("{",grep,"}")

SKIP=gsub(",","|", skip)
OUT=paste0(out.suffix,".")
if( OUT==".")    OUT=""
if( OUT=="NA.")  OUT=""

if( dir.exists(out.dir) ==FALSE ) dir.create(out.dir)


if( plot.sep == TRUE  ) { PLOT_TYPES = c("single", "multiple","separate") }  
if( plot.sep == FALSE ) { PLOT_TYPES = c("single", "multiple")            }  

cat('\n')
cat("############################################################################\n")
cat('------------------------ Summary of Parameters -----------------------------\n')
cat("############################################################################\n")
cat('\n')
    cat("   directory    = ", getwd(),"\n")
    cat("   ref.file     = ", ref.file,"\n")
    cat("   grep         = ", grep,"\n")
    cat("   skip         = ", skip,"\n")
    cat("   path         = ", path,"\n")
    cat("   out.suffix   = ", out.suffix,"\n")
    cat("   out.dir      = ", out.dir,"\n")
#   cat("   OUT          = ", OUT,"\n")
    cat("   chr.filt     = ", chr.filt,"\n")
    cat("   organelles   = ", organelles,"\n")
    cat("   chr.rename   = ", chr.rename,"\n")
#    cat("   not.rename   = ", not.rename,"\n")
    cat("   rm.unpl_ctg  = ", rm.unpl_ctg,"\n")
    cat("   chr0_plot    = ", chr0_plot,"\n")
    cat("   plot.conf    = ", plot.conf,"\n")
    cat("   thr.conf     = ", thr.conf,"\n")
    cat("   plot.h       = ", plot.h,"\n")
    cat("   plot.w       = ", plot.w,"\n")
    cat("   plot.sep     = ", plot.sep,"\n")
    cat("   plot.mix     = ", plot.mix,"\n")
    cat("   CTG.min      = ", CTG.min,"\n")
    cat("   MIN.WIND     = ", MIN.WIND,"\n")
    cat("   COL.LIST     = ", paste(COL.LIST,collapse=","),"\n")
    cat("   short.name   = ", short.name,"\n")
    cat("   ctg.gap      = ", ctg.gap,"\n")
    cat("   ctg.fa       = ", ctg.fa,"\n")
    cat("   main         = ", main,"\n")
    cat("   short.name   = ", short.name,"\n")
#     cat("   rm.ref       = ", rm.ref,"\n")
#     cat("   rm.scaf      = ", rm.scaf,"\n")
    cat("   rename       = ", rename,"\n")
    cat("   name.del     = ", name.del,"\n")
    cat("   agp.bionano  = ", bionano,"\n")
cat('\n')
cat('\n')
 

# select files with GREP or default

if ( is.na(path)==TRUE  ) 
{ 
#DIR_SELECTION = paste0("?ag?ag*?caffold*",grep.words ) 
DIR_SELECTION = paste0(grep.words ) 
if ( is.na(skip)==FALSE ) { DIR_SELECTION = grep( SKIP, DIR_SELECTION, value=TRUE, invert=TRUE) }

fas.files = system( paste0( "ls *",DIR_SELECTION,"*/ragtag.scaffold.fasta"          ) , inter=TRUE)
#chr.files = system( paste0( "ls *",DIR_SELECTION,"*/ragtag.scaffold.fasta.fai"      ) , inter=TRUE)
agp.files = system( paste0( "ls *",DIR_SELECTION,"*/ragtag.scaffold.agp"            ) , inter=TRUE)
cnf.files = system( paste0( "ls *",DIR_SELECTION,"*/ragtag.scaffold.confidence.txt" ) , inter=TRUE)
sta.files = system( paste0( "ls *",DIR_SELECTION,"*/ragtag.scaffold.stats"          ) , inter=TRUE)
asm.files = system( paste0( "ls *",DIR_SELECTION,"*/ragtag.scaffold.asm.paf"        ) , inter=TRUE) 
log.files = system( paste0( "ls *",DIR_SELECTION,"*/ragtag.scaffold.asm.paf.log"    ) , inter=TRUE) 
if( bionano== TRUE ) bio.files = system( paste0( "ls *",DIR_SELECTION,"*/*bionano*.agp" ) , inter=TRUE) 

}

# select files PATH

if ( is.na(path)==FALSE ) 
{
DIR_SELECTION = system( paste0( "ls -d ", path) , inter=TRUE)
if ( is.na(skip)==FALSE ) { DIR_SELECTION = grep( SKIP, DIR_SELECTION, value=TRUE, invert=TRUE) }

fas.files = as.character( sapply(  paste0( DIR_SELECTION,"/ragtag.scaffold.fasta"         ) , function(x)  system( paste0( "ls ", x  ) , inter=TRUE)))
#chr.files = as.character( sapply(  paste0( DIR_SELECTION,"/ragtag.scaffold.fasta.fai"     ) , function(x)  system( paste0( "ls ", x  ) , inter=TRUE)))
agp.files = as.character( sapply(  paste0( DIR_SELECTION,"/ragtag.scaffold.agp"           ) , function(x)  system( paste0( "ls ", x  ) , inter=TRUE)))
cnf.files = as.character( sapply(  paste0( DIR_SELECTION,"/ragtag.scaffold.confidence.txt") , function(x)  system( paste0( "ls ", x  ) , inter=TRUE)))
sta.files = as.character( sapply(  paste0( DIR_SELECTION,"/ragtag.scaffold.stats"         ) , function(x)  system( paste0( "ls ", x  ) , inter=TRUE)))
asm.files = as.character( sapply(  paste0( DIR_SELECTION,"/ragtag.scaffold.asm.paf"       ) , function(x)  system( paste0( "ls ", x  ) , inter=TRUE)))
log.files = as.character( sapply(  paste0( DIR_SELECTION,"/ragtag.scaffold.asm.paf.log"   ) , function(x)  system( paste0( "ls ", x  ) , inter=TRUE)))
if( bionano== TRUE ) bio.files = as.character( sapply(  paste0( DIR_SELECTION,"/*bionano*.agp"   ) , function(x)  system( paste0( "ls ", x  ) , inter=TRUE)))
}



fas.files = unique( fas.files)
#chr.files = unique( chr.files)
agp.files = unique( agp.files)
sta.files = unique( sta.files)
asm.files = unique( asm.files)
cnf.files = unique( cnf.files)
log.files = unique( log.files)
if( bionano== TRUE ) bio.files = unique( bio.files)


# skip 

if ( is.na(skip)==FALSE ) 
{
fas.files =  grep( SKIP, fas.files, value=TRUE, invert=TRUE) 
#chr.files =  grep( SKIP, chr.files, value=TRUE, invert=TRUE) 
agp.files =  grep( SKIP, agp.files, value=TRUE, invert=TRUE) 
cnf.files =  grep( SKIP, cnf.files, value=TRUE, invert=TRUE) 
sta.files =  grep( SKIP, sta.files, value=TRUE, invert=TRUE) 
asm.files =  grep( SKIP, asm.files, value=TRUE, invert=TRUE) 
log.files =  grep( SKIP, log.files, value=TRUE, invert=TRUE) 
if( bionano== TRUE )  bio.files = grep( SKIP, bio.files, value=TRUE, invert=TRUE) 
}




# remove unwanted folders

# file.check = data.frame( CHR = chr.files, AGP=agp.files , CNF = cnf.files, STA = sta.files, ASM = asm.files, LOG = log.files) 
# file.check = data.frame( CHR = chr.files, AGP=agp.files , CNF = cnf.files, ASM = asm.files, LOG = log.files) 
# do the faiindex later

file.check = data.frame( FAS = fas.files, AGP=agp.files , CNF = cnf.files, ASM = asm.files, LOG = log.files) 
if( bionano== TRUE ) file.check = data.frame( FAS = fas.files, AGP=agp.files , CNF = cnf.files, ASM = asm.files, LOG = log.files, BIO = bio.files) 
cat('\n')
cat("############################################################################\n")
cat(" ==> Files to be imported:\n")
cat("############################################################################\n")
cat('\n')
print(file.check)
cat('\n')


# cat('\n')
# cat("############################################################################\n")
# cat(" ==> Files Names and order:\n")
# cat("############################################################################\n")
# cat('\n')


# generate index if missing

# to.run = setdiff( fas.files ,gsub(".fai","",chr.files))
# if (length(to.run)>0) {
#     cat("\n ==> generating FAI index:\n") ; for ( fasta in to.run) {  cat(fasta,"\n") ;  system( paste0( "samtools faidx ",fasta )) } ;
#     if ( is.na(path)==FALSE )  chr.files = as.character( sapply(  paste0( DIR_SELECTION,"/ragtag.scaffold.fasta.fai"     ) , function(x)  system( paste0( "ls ", x  ) , inter=TRUE)))
#     if ( is.na(path)==TRUE  )  chr.files = system( paste0( "ls *",DIR_SELECTION,"*/ragtag.scaffold.fasta.fai"      ), inter=TRUE)
#     if ( is.na(skip)==FALSE )  chr.files = grep( SKIP, chr.files, value=TRUE, invert=TRUE)
# 
# }


# FILE CHECK AND FOLDER CHECK
folder.check = data.frame( RagTag_Folders = gsub("/ragtag.scaffold.fasta","",fas.files )) 
cat('\n')
cat("############################################################################\n")
cat(" ==> Folders :\n")
cat("############################################################################\n")
cat('\n')
print(folder.check)
cat('\n')

# Simplify names 
rm.str=""


SAMPLE        = data.frame( file=agp.files, name = gsub("/ragtag.scaffold.agp","", agp.files), stringsAsFactors=FALSE)
SAMPLE$"name" = gsub(".ag.ag..caffold.","", SAMPLE$"name")
SAMPLE$"x"    = gsub(",|_|-|/",".", SAMPLE$"name")
SAMP_2        = SAMPLE
SAMPLE        = split(SAMPLE$"x", SAMPLE$"name")
SAMPLE        = lapply(SAMPLE, function(x) unlist(strsplit(x,".", fixed=TRUE)))
SAMPLE        = lapply(SAMPLE, unique)
all.words     = unique(unlist(SAMPLE))


# if one word is repeat twice in only 1 or few names (e.g. same name for reference.suffix and sample), remove it before all the others
ALL.WORDS = list()
for (i in all.words) ALL.WORDS[[i]] = sapply(SAMPLE, function(x) { sum(x==i) })
ALL.WORDS = data.frame( do.call(cbind,ALL.WORDS))
any_rep = sapply(ALL.WORDS, function(x) any(x>1))
any_rep = any_rep [ any_rep == TRUE]
if (length(any_rep)>0)
{
    which_rep = lapply(SAMPLE, function(x) { which(x== names(any_rep))  })
    all.reps     = unique(unlist(which_rep))
    # remove the same shared position in all
    for (j in all.reps) if ( any ( sapply(which_rep, function(x) any(x==j))  < 1) ) all.reps = setdiff(all.reps, j)
    # remove position 
    
    SAMPLE        = lapply(SAMPLE, function(x) x[-all.reps])
    all.words     = unique(unlist(SAMPLE))
    
    
    }

for( word in all.words) { if(  all(sapply(SAMPLE, function(x) any(x==word) ))) { SAMPLE = lapply(SAMPLE, function(x) x[ x!= word ]); rm.str=paste(rm.str,word,sep="|.")     } }
if ( any(sapply(SAMPLE,function(x)length(x)==0)) ){  w=which(sapply(SAMPLE,function(x)length(x)==0)) ; xx = names(SAMPLE[w]) ;   xx = do.call(rbind,strsplit(xx,"/")) ; xx=xx[,ncol(xx)];  SAMPLE[[w]] = xx;  SS=SAMPLE[[w]] }  # if 1 short name == "", give the real name


SAMPLE  =  lapply(SAMPLE, paste, collapse=".")
SAMPLE  =  data.frame( long_name = names(SAMPLE), short_name = as.character(SAMPLE), stringsAsFactors=FALSE)

rm.str = paste0( rm.str,   gsub(".","",rm.str, fixed=TRUE),"|/")

if( short.name == TRUE  ) RM.STRING = rm.str
if( short.name == FALSE ) RM.STRING = ""

# rm ref name
www = which(SAMPLE$"long_name" ==SAMPLE$"short_name" ) ; 
SS1 = SAMPLE$"long_name"[www]
  x = unlist(strsplit(SS1,".",fixed=T))
  x = x[ -length(x)]; 
SS2 = paste(x, collapse=".") 

if (length(www)==1) { SAMPLE$"short_name"[www] = SS2 }


# for 1 sample only, use only file name
if (nrow(SAMPLE)==1) if( grepl("/",SAMPLE$"long_name")==TRUE)  SAMPLE$"short_name" = sapply( strsplit( gsub(".ag.ag.|.caffold.|","",SAMPLE$"long_name"),"/") , tail, 1)
# if 1 sample emply, use only file name
w.empty = which(SAMPLE$"short_name"=="")
if (length(w.empty)>0)  SAMPLE$"short_name"[w.empty] =  sapply( strsplit( gsub(".ag.ag.|.caffold.|","",SAMPLE$"long_name"[[w.empty] ]),"/") , tail, 1)



# invert long_name and name
xxx                = SAMP_2$"long_name"
SAMP_2$"long_name" = SAMP_2$"name"
SAMP_2$"name" = xxx

names(SAMP_2) = c("file","name","long_name")
SAMP_3 = merge(SAMP_2,SAMPLE, by="long_name",all.x=T)

# if only 1 name has not numeric part --> sort by number 
# sort by number if any
SAMPLE$"numeric.sort" = SAMPLE$"short_name"
SAMPLE$"numeric.sort" = gsub("Kb","000",SAMPLE$"numeric.sort")
SAMPLE$"numeric.sort" = gsub("Mb","000000",SAMPLE$"numeric.sort")
SAMPLE$"numeric.sort" = gsub("Gb","000000000",SAMPLE$"numeric.sort")
SAMPLE$"numeric.sort" = gsub("less.","1001",SAMPLE$"numeric.sort")
SAMPLE$"numeric.sort" = gsub("All|all","0",SAMPLE$"numeric.sort")
SAMPLE$"numeric.sort" = gsub(".F",".",SAMPLE$"numeric.sort")
#SAMPLE$"numeric.sort" = as.numeric(gsub("[a-z]|[A-Z]|[:punct:]|_|-|[.]","",SAMPLE$"numeric.sort"))
SAMPLE$"numeric.sort" = as.numeric(gsub("[a-z]|[A-Z]|[:punct:]|_|-|"   ,"",SAMPLE$"numeric.sort")) ## keep the . for decimal positions

if ( sort.sam == TRUE) 
{ if( sum(is.na(SAMPLE$"numeric.sort")) <=1 ) { SAMPLE = SAMPLE[ order(SAMPLE$"numeric.sort"),]  }
}



# if rename == something
if( nrow(SAMPLE)==1 & rename!=FALSE ) { SAMPLE$"NEW_short_name" = rename  }

# simplify names with

name.del=gsub(",","|", name.del)

if( name.del!=FALSE  ) { SAMPLE$"short_name" = gsub(name.del,"", SAMPLE$"short_name" )  }
if( name.del!=FALSE  ) { SAMP_3$"short_name" = gsub(name.del,"", SAMP_3$"short_name")  }



cat('\n')
cat("############################################################################\n")
cat(" ==> SAMPLE NAMES :\n")
cat("############################################################################\n")
cat('\n')


print(SAMPLE)
# print(SAMP_2)
# print(SAMP_3)
# 
if( nrow(SAMPLE)==1 & rename!="FALSE" ) { SAMPLE$"short_name" = rename  }
if( nrow(SAMPLE)==1 & rename!="FALSE" ) { SAMP_3$"short_name" = rename  }



SAMP_4 = SAMPLE
SAMP_4$"common" = "-"
for( rr in 1:nrow(SAMP_4)) SAMP_4$"common"[rr] = gsub(  SAMP_4$"short_name"[rr] , "|",  SAMP_4$"long_name"[rr])
SAMP_RM  = unique(SAMP_4$"common")
SAMP_RM  = gsub("-|_",".",SAMP_RM)
SAMP_RM  = gsub("/","|",SAMP_RM)
if (length(SAMP_RM)==1) RM.STRING = SAMP_RM



cat('\n')
cat("############################################################################\n")
cat(" ==> Reference:\n")
cat("############################################################################\n")
cat("\n")

# REFERENCE FAI FILE
ref.fai =paste0(ref.file,".fai")
if (file.exists(ref.fai)==FALSE ) { system( paste0( "samtools faidx ",ref.file )) }

ref=fread(ref.fai)[,1:2]
names(ref)=c("chr","ref.len")
#ref = ref[ rev(order(ref$"chr")) ,]
ref= as.data.frame(ref)
#revert order
#ref = ref[ order(rev(ref$"chr")) ,]
ref$"chr" = gsub("_RagTag","",ref$"chr")

# sort if the have progressive names

ref$"prefix" = substr(ref$"chr",1,2)
ref$"num" = as.numeric( gsub("CHR|Chr|LG|CM0|C|l|c|ptg|PTG|h1g|h2g|","",ref$"chr" ))
ref$"num2"= ref$"num"

cat("\n")
print(ref)
cat("\n")

# Sorting Reference 
cat("   --> Sorting Reference: ")

# All with numeric part (e.g. Canonical Chromosomes Chr1 Chr2 Chr3)
if ( sum(is.na(ref$"num")) ==0 ) { ref = ref[ rev(order(ref$"chr")) ,] ; cat(" Cannonical Chromosomes only --> Alphanumeric sorting: \n") }

# Modified Chromosome-like names --> already sorted and most of them (>70%) have same prefix
if ( sum(is.na(ref$"num")) >0 & sum(is.na(ref$"num")) == nrow(ref) & max(table(ref$prefix)/nrow(ref))> 0.7 ) { ref = ref[ rev(order(ref$"chr")) ,]  ; cat(" Non Cannonical Chromosome names --> Alphanumeric sorting: \n") }

# Some with no numeric part (etc, ChrC, ChrM, ChrUnkwn) --> sort the non-numeric by size from bigger to smaller
if ( sum(is.na(ref$"num")) >0 & sum(is.na(ref$"num"))  < nrow(ref) ) { wna = which(is.na(ref$"num")) ;  ref$"num2"[wna] = max( ref$"num2"[ -wna]) + 1 +  abs(ref$"ref.len"[wna] -  max( ref$"ref.len"[wna] )) / min(( ref$"ref.len"[wna] ))/1000 ; ref = ref[ rev(order(ref$"num2")) ,] ; cat(" Some Cannonical Chromosomes --> First Chromosomes, then organels and scaffolds/contigs \n")  }

# If contig names with same prefix, all numeric, and not already sorted FW or REV) --> sort by contig size
if ( length(unique(ref$"prefix"))==1 & sum(is.na(ref$"num"))==0 & any(order(ref$"ref.len") != 1:nrow(ref)) & any(order(ref$"ref.len") != nrow(ref):1  ) )  { ref = ref[ rev(order(ref$"num")) ,] ; ; cat(" Scaffolds or Contigs --> Sorted by length \n") }

cat("\n")

ref$"prefix" = NULL 
ref$"num" = NULL
ref$"num2" = NULL

#print(ref)


# ref repeat masker
# if ( is.na(rm.ref)==FALSE )
# { 
#     REF.RM = fread( paste0(ref.file,".gff") )
#     }

chr_prefix = unique(substring(ref$"chr",1,3))
ref.name = system( paste(' grep ">" ',ref.file), inter=TRUE) 
ref.name = gsub("\t"," ",ref.name)
ref.name = gsub("="," ",ref.name)                                 # some NCBI accessions have \t instead of spaces
# some NCBI accessions have \t instead of spaces
ref.df= data.frame( name=ref.name, stringsAsFactors=FALSE     )

if (all(grepl("LG",chr_prefix))==TRUE  ) chr_prefix="LG"


cat("   --> Genome fasta index",ref.fai,"\n")
cat("   --> Chromosome prefix:",chr_prefix,"\n")
cat("\n")
#print( head(ref[ rev(order(ref$"ref.len")), ],30))
#cat("\n")

#print(not.rename)

# if some prefix are not CHR or CTG, but you can find Chr LG etc in the fasta header
#if ( all(chr_prefix %in% c("Chr","chr","CTG")==FALSE)  & sum(grepl("Chr|chr|LG|chl|Chl|Mit|mit|CTG", ref.name))>1  & not.rename==FALSE )  
if ( all(chr_prefix %in% c("Chr","chr","CTG","LG")==FALSE)  & sum(grepl("Chr|chr|LG|chl|Chl|Mit|mit|CTG", ref.name))>1  )  
{
cat("   ==> Chromosome complete names:\n\n")
aa_chr = grep("LG|Chr"     ,ref.df$"name",ignore.case=TRUE)
aa_org = grep("mitoc|plast" ,ref.df$"name",ignore.case=TRUE)
aa_sca = grep("scaff|contig",ref.df$"name",ignore.case=TRUE)


aa_selection = c(aa_chr,aa_org, head(aa_sca,5))


ref.df2 = data.frame( name=ref.df[ aa_selection , ], stringsAsFactors=FALSE     )
ref.df3 = ref.df2
ref.df3$"chr" = NA

w_chr = grep("LG|Chr"       ,ref.df3$"name",ignore.case=TRUE) 
w_org = grep("mitoc|plast"  ,ref.df3$"name",ignore.case=TRUE)
w_sca = grep("scaff|contig" ,ref.df3$"name",ignore.case=TRUE)

if (length(w_chr)>0 ) ref.df3$"chr" [w_chr ] = " --> Main Chromosomes or Linkage Groups"
if (length(w_org)>0 ) ref.df3$"chr" [w_org ] = " --> Organelles"
if (length(w_sca)>0 ) ref.df3$"chr" [w_sca ] = " --> top Scaffold/Contigs"

#print(ref.df3)

ref.df3 = setDT(ref.df3)
ref.df3 = split(ref.df3 ,ref.df3$"chr")
ref.df3 = lapply(ref.df3, function(x) { x$chr=NULL; return(x)})

#if (length(aa_chr) >0  ) { print( ref.df2  )            } else print(head(ref.df,30))
if (length(aa_chr) >0  ) { print( ref.df3  )            } else print(head(ref.df,30))
if (length(aa_sca) >5  ) { cat("       --> total of",length(aa_sca),"Scaffold/Contigs in the reference\n") }


print(1)
print(head(ref.df3))

# remove common words
ref.words = strsplit(ref.name, " ")
ref.all.words = sort(unique(unlist(ref.words)))

# count words occurencies for each ">" line
 REF.ALL.WORDS = list()
 for (i in ref.all.words) REF.ALL.WORDS[[i]] = sapply(ref.words, function(x) { sum(x==i) })

# identiofy duplicated columns names
 REF.ALL.WORDS.sin = REF.ALL.WORDS [ which(sapply(REF.ALL.WORDS, function(x) ! all(x==1))) ] 
 REF.ALL.WORDS.dup = REF.ALL.WORDS [ which(sapply(REF.ALL.WORDS, function(x)   all(x>=1))) ] 
 ref.dup = names(REF.ALL.WORDS.dup)
 ref.dup = grep("Chr|chr|LG|CTG",ref.dup, invert=T, value = TRUE)
 ref.dup = paste(ref.dup, collapse= "|")
print(2)

ref.df$"simp_name" = ref.df$"name"
ref.df$"simp_name" = gsub(ref.dup, "", ref.df$"simp_name" )
ref.df$"simp_name" = gsub("       |     |   |  |,|;", "  ", ref.df$"simp_name" )

 ref.df$"simp_name" = gsub("whole|shotgun|sequence|,|;", "  ", ref.df$"simp_name" )
 ref.df$"simp_name" = gsub("\\s{2,}", "  ", ref.df$"simp_name")
 
 #remove space after works Chromosomes Chr LG
  ref.df$"simp_name" = gsub("mosome:  ", "mosome", ref.df$"simp_name")
  ref.df$"simp_name" = gsub("mosome: " , "mosome", ref.df$"simp_name")
  ref.df$"simp_name" = gsub("mosome  " , "mosome", ref.df$"simp_name")
  ref.df$"simp_name" = gsub("mosome "  , "mosome", ref.df$"simp_name")
  ref.df$"simp_name" = gsub("Chr ", "Chr", ref.df$"simp_name")
  ref.df$"simp_name" = gsub("chr ", "chr", ref.df$"simp_name")
  ref.df$"simp_name" = gsub("LG  ", "LG", ref.df$"simp_name")
  ref.df$"simp_name" = gsub("contig: " , "contig:", ref.df$"simp_name")

 
 
print(3)
print(head(ref.df))

ref.edit = suppressWarnings( data.frame(do.call(rbind, strsplit(ref.df$"simp_name"," " ))))
for ( ii in names(ref.edit)) if ( all( as.character(ref.edit[,ii])==" ") | all(as.character(ref.edit[,ii])=="")) ref.edit[ii] = NULL
print(4)
print(head(ref.edit))


CHR_PREFIX = paste(chr_prefix, collapse="|")

# colums with NCBI name and Chr name
w1.col = names( which(sapply(ref.edit,function(x) length(grep(CHR_PREFIX,x))) == nrow(ref)))
w2.col = names( which(sapply(ref.edit,function(x) length(grep("Chr|chr|LG",x))) > 0))
w3.col = names( which(sapply(ref.edit,function(x) length(grep(".caffold|contig",x))) > 0))

print(w1.col)
print(w2.col)
print(w3.col)


#if (length(w2.col)>1 ) { w2.col = w2.col[1] ; ref.edit = ref.edit[ , c(w1.col ,w2.col )] }
if (length(w2.col)>1 ) { w2.col = w2.col[1] ; ref.edit = ref.edit[ , c(w1.col ,w2.col )] }
print(w2.col)
print(5)

name.ref.edit = names(ref.edit)
name.ref.edit[ name.ref.edit == w1.col] = "chr"
name.ref.edit[ name.ref.edit == w2.col] = "chr.rename"
names(ref.edit) = name.ref.edit

print(head(ref.edit))

ref.edit$"chr"        = gsub(">","", ref.edit$"chr")
ref.edit$"chr.rename" = gsub("chromosome:  |Chromosome:  |chromosome: |Chromosome: |Chromosome |Chromosome|chromosome |chromosome","Chr", ref.edit$"chr.rename")
ref.edit$"chr.rename" = gsub("Chr: ","Chr", ref.edit$"chr.rename")
ref.edit$"chr.rename" = gsub(": ",":", ref.edit$"chr.rename")
ref.edit$"chr.rename" = gsub("contig:","", ref.edit$"chr.rename")
print(6)
print(head(ref.edit))
print(tail(ref.edit))

w_chl = grep("plastid|chloroplast"           ,ref.edit$"chr.rename",ignore.case = TRUE)
w_mit = grep("mitochondrion"                 ,ref.edit$"chr.rename",ignore.case = TRUE)
w_ctg = grep("ptg|h1tg|h2tg|h1tg|tig|contig" ,ref.edit$"chr.rename",ignore.case = TRUE) # put capital letter in order to have different names with the assembly

grep("plastid|chloroplast",ref.edit$"chr.rename",ignore.case = TRUE)
grep("mitochondrion"      ,ref.edit$"chr.rename",ignore.case = TRUE)
print(7)

if (length(w_chl) == 1 ) ref.edit$"chr.rename"[w_chl]="ChrC"
if (length(w_mit) == 1 ) ref.edit$"chr.rename"[w_mit]="ChrM"
if (length(w_ctg) == 1 ) ref.edit$"chr.rename"[w_ctg]=toupper(ref.edit$"chr.rename"[w_ctg])  #
print(head(ref.edit))

ref.edit$"nchar" = nchar(ref.edit$"chr.rename")
chr_lines = grep("Chr[1-9]", ref.edit$"chr.rename")
# nchar of Chr, to decide if adding a zero
min_char = min(ref.edit$"nchar"[chr_lines] )
max_char = max(ref.edit$"nchar"[chr_lines] )
print(min_char)
print(max_char)


if( min_char == max_char-1  ) ref.edit[ ref.edit$"nchar" == min_char & grepl("Chr[1-9]", ref.edit$"chr.rename"), ]$"chr.rename" = gsub("Chr","Chr0",ref.edit[ ref.edit$"nchar" == min_char & grepl("Chr[1-9]", ref.edit$"chr.rename"), ]$"chr.rename")
print(8)
print(head(ref.edit))
print(tail(ref.edit))

ref.edit$"nchar" = NULL
ref.edit.full = ref.edit 
#ref.edit = ref.edit [ , 1:2]
final_col=c("chr","chr.rename")
ref.edit = ref.edit [ , final_col]

###  Check Scaffold names, if all the same, take it from other columns
w_chl = grep("Chr"           ,ref.edit$"chr.rename",ignore.case = TRUE)
w_oth = grep("Chr"           ,ref.edit$"chr.rename",ignore.case = TRUE, invert=TRUE)

print(w3.col)
print(length(w_chl))
print(length(w_oth))

print( any(duplicated(ref.edit$"chr.rename"[w_oth])) )

   
if (any(duplicated(ref.edit$"chr.rename"[w_oth]))   ) 
{

if( length(w3.col) >1)  { scaff_edit = ref.edit.full[ w_oth , w3.col ] }
if( length(w3.col)==1)  { scaff_edit = data.frame( a=1:nrow( ref.edit.full)) ; scaff_edit[[w3.col]] = ref.edit.full[  , w3.col ]  ; scaff_edit = scaff_edit[ w_oth , ]  ; scaff_edit$"a" = NULL }
    print(8.5)
    print(head(scaff_edit))  
    print(tail(scaff_edit))  
    print(nrow(scaff_edit))  
    print( sapply(scaff_edit, function(x) length(unique(x)) ) ) 
    
    scaff_stat = sapply(scaff_edit, function(x) length(unique(x)) ) 
    scaff_stat = scaff_stat[ scaff_stat == nrow(scaff_edit)]
    
    print(8.6)
    print(scaff_stat)
    scaff_stat = scaff_stat[1]
    print(scaff_stat)
    scaff_col = names(scaff_stat)
   
    print(8.7)     
    print(length(ref.edit$"chr.rename"[w_oth]))
    print(8.8)     
    print(length(scaff_edit[,scaff_col] ))
    print(8.9)     
    print(head(scaff_edit))  
    print(tail(scaff_edit))  
print(head(ref.edit))
print(tail(ref.edit))    
    print(8.99)     
    print(length(scaff_col))
#    print(ref.edit$"chr.rename"[w_oth])
#    print(scaff_edit[,scaff_col] )
   
if( length(scaff_col)==1) {  ref.edit$"chr.rename"[w_oth] = as.character(scaff_edit[,scaff_col])  }
if( length(scaff_col)==0) {  ref.edit$"chr.rename"[w_oth] = paste(ref.edit$"chr.rename"[w_oth],1:nrow(ref.edit),sep="_") }  # no info, add progressive number
}



print(9)

print(head(ref.edit))
print(tail(ref.edit))

cat("\n\n")
  
cat('\n  ==> If you want to rename chromosomes as below  ===>> USE OPTION  --chr.rename !! \n\n')
print(ref.edit)
print(ref)
ref.temp = merge(ref.edit,ref,by="chr")
print(head(ref.temp,50))

# print(grep("ChrC|ChrM",ref.edit.full$"chr.rename", value=TRUE) )
# print(ref.edit[ grep("ChrC|ChrM",ref.edit$"chr.rename"), ])
        ref.edit[ grep("ChrC|ChrM",ref.edit$"chr.rename"), ]$"chr" -> organelle_ctgs
# print(organelle_ctgs)

} else { ref.temp=ref; print( head(ref[ rev(order(ref$"ref.len")), ],30)) }

cat("\n\n")

cat("   --> Chromosome length filter :",chr.filt,"\n")
cat("\n")
chr_to_rm = ref[ ref$"ref.len" < chr.filt  ,]$"chr"
chr_to_kp = ref[ ref$"ref.len" >= chr.filt ,]$"chr"

if (organelles == TRUE ){ 
    chr_to_rm =  grep("mitoc|plast"  ,chr_to_rm,ignore.case=TRUE,invert=TRUE, value=TRUE)
    chr_to_rm =  setdiff(chr_to_rm, organelle_ctgs)
    if ( length(organelle_ctgs) >0 ) { cat("   --> Organelles present\n\n") ; print(ref.edit[ grep("ChrC|ChrM",ref.edit$"chr.rename"), ]) ; cat("\n") }
    if ( length(organelle_ctgs)==0 ) { cat("   --> Organelles absent \n\n") }
    }     

ref = ref[ ! ref$"chr" %in% chr_to_rm ,]
ref.temp =  ref.temp[ ! ref.temp$"chr" %in% chr_to_rm ,]
#nrow(ref)
#print(ref)
print(ref.temp[ order(ref.temp$"chr"), ])
cat("\n")


if ( all(grepl("ChrC|ChrM",ref$"chr")==FALSE) )  { cat("   --> Organelles not selected \n") }
if ( any(grepl("ChrC"     ,ref$"chr")==TRUE ) )  { cat("   --> ChrC selected \n") }
if ( any(grepl("ChrM"     ,ref$"chr")==TRUE ) )  { cat("   --> ChrM selected \n") }


cat("\n")



if ( ctg.gap == TRUE & ctg.fa == FALSE)   cat('\n ==> WARNING ! Make sure you have bedtools installed/loaded ! \n')
if ( ctg.gap == TRUE & ctg.fa == TRUE )   cat('\n ==> WARNING ! Make sure you have bedtools and seqtk installed/loaded ! \n')

if ( ctg.gap == TRUE & ctg.fa == TRUE )   cat('\n ==> If you want to rename chromosomes, used --chr.rename! \n')
# 
# print(head(ref.df))
# print(head(ref.edit))
# 

# PRESS ENTER TO CONTINUE
cat('\n')

    pause = function()
    {
        if (interactive()) { invisible(readLines(file("stdin"), 1)) } else {
            cat("---------------------------------------------------------------------------\n")
            cat("--------Please check DIRECTORY,PATH,FILES PARAMETERS printed above---------\n")
            cat("----------------------- Press any key to continue -------------------------\n")
            cat("------------------------- Press Ctrl+Z to abort ---------------------------\n")
            cat("---------------------------------------------------------------------------\n")
            invisible(readLines(file("stdin"), 1))
        }
    }   
cat('\n')
pause()


# IMPORT FAI index files
# generate index if missing
# 
# if ( is.na(path)==TRUE  ) chr.files = system( paste0( "ls *",DIR_SELECTION,"*/ragtag.scaffold.fasta.fai"      ) , inter=TRUE)
# if ( is.na(path)==FALSE ) chr.files = as.character( sapply(  paste0( DIR_SELECTION,"/ragtag.scaffold.fasta.fai"     ) , function(x)  system( paste0( "ls ", x  ) , inter=TRUE)))
# if ( is.na(skip)==FALSE ) chr.files =  grep( SKIP, chr.files, value=TRUE, invert=TRUE) 
# chr.files = unique( chr.files)


# RUN FAIDX ALL THE TIME
# because if you overwrite the ragtag files, faidx will be relative to the old one
cat("\n ==> generating FAI index:\n") ; 
for ( fasta in fas.files) {  cat(fasta,"\n") ;  system( paste0( "samtools faidx ",fasta )) } ;

# RUN missing FAIDX 
# to.run = setdiff( fas.files ,gsub(".fai","",chr.files))
# if (length(to.run)>0) {
#     cat("\n ==> generating FAI index:\n") ; for ( fasta in to.run) {  cat(fasta,"\n") ;  system( paste0( "samtools faidx ",fasta )) } ;
#     if ( is.na(path)==FALSE )  chr.files = as.character( sapply(  paste0( DIR_SELECTION,"/ragtag.scaffold.fasta.fai"     ) , function(x)  system( paste0( "ls ", x  ) , inter=TRUE)))
#     if ( is.na(path)==TRUE  )  chr.files = system( paste0( "ls *",DIR_SELECTION,"*/ragtag.scaffold.fasta.fai"      ), inter=TRUE)
#     if ( is.na(skip)==FALSE )  chr.files = grep( SKIP, chr.files, value=TRUE, invert=TRUE)
# 
# }
fas.files = unique( fas.files)
chr.files = paste0 (fas.files,".fai") 


# FILE CHECK 

file.check = data.frame( CHR = chr.files, AGP=agp.files , CNF = cnf.files, ASM = asm.files, LOG = log.files) 
if( bionano== TRUE ) file.check = data.frame( CHR = chr.files, AGP=agp.files , CNF = cnf.files, ASM = asm.files, LOG = log.files, BIO = bio.files) 



cat('\n')
cat("############################################################################\n")
cat(" ==> Loading Data:\n")
cat("############################################################################\n")
cat('\n')

cat("==> Warning messages in fread() are expected at this step \n\n")


# load files, 
CHR= list()
AGP= list()
CNF= list()
STA= list()
ASM= list()
LOG= list()
# BIO= list()



for( i in chr.files ) {  j=gsub( paste0(".ag.ag..caffold.|/ragtag.scaffold.fasta.fai")     ,"",i) ; j=SAMPLE[ SAMPLE$"long_name" ==j,]$"short_name" ; if( substr(j,1,1) %in% c("_",".","-")) {j=substr(j,2,100)} ; j=gsub("_",".",j); j=gsub("/$","",j); CHR[[j]]=fread(i)}
for( i in agp.files ) {  j=gsub( paste0(".ag.ag..caffold.|/ragtag.scaffold.agp")           ,"",i) ; j=SAMPLE[ SAMPLE$"long_name" ==j,]$"short_name" ; if( substr(j,1,1) %in% c("_",".","-")) {j=substr(j,2,100)} ; j=gsub("_",".",j); j=gsub("/$","",j); AGP[[j]]=fread(i)}
for( i in cnf.files ) {  j=gsub( paste0(".ag.ag..caffold.|/ragtag.scaffold.confidence.txt"),"",i) ; j=SAMPLE[ SAMPLE$"long_name" ==j,]$"short_name" ; if( substr(j,1,1) %in% c("_",".","-")) {j=substr(j,2,100)} ; j=gsub("_",".",j); j=gsub("/$","",j); CNF[[j]]=fread(i)}
#for( i in sta.files) {  j=gsub( paste0(".ag.ag..caffold.|/ragtag.scaffold.stats")         ,"",i) ; j=SAMPLE[ SAMPLE$"long_name" ==j,]$"short_name" ; if( substr(j,1,1) %in% c("_",".","-")) {j=substr(j,2,100)} ; j=gsub("_",".",j); j=gsub("/$","",j); STA[[j]]=fread(i)}
for( i in asm.files ) {  j=gsub( paste0(".ag.ag..caffold.|/ragtag.scaffold.asm.paf")       ,"",i) ; j=SAMPLE[ SAMPLE$"long_name" ==j,]$"short_name" ; if( substr(j,1,1) %in% c("_",".","-")) {j=substr(j,2,100)} ; j=gsub("_",".",j); j=gsub("/$","",j); ASM[[j]]=suppressWarnings(  fread(i, fill=T) )}
for( i in log.files ) {  j=gsub( paste0(".ag.ag..caffold.|/ragtag.scaffold.asm.paf.log")   ,"",i) ; j=SAMPLE[ SAMPLE$"long_name" ==j,]$"short_name" ; if( substr(j,1,1) %in% c("_",".","-")) {j=substr(j,2,100)} ; j=gsub("_",".",j); j=gsub("/$","",j); LOG[[j]]=system( paste( "cat ", i ), inter=T)}

# if( bionano== TRUE ){
# for( i in bio.files ) {  j=gsub( paste0(".ag.ag..caffold.|/bionano.agp")                   ,"",i) ; j=SAMPLE[ SAMPLE$"long_name" ==j,]$"short_name" ; if( substr(j,1,1) %in% c("_",".","-")) {j=substr(j,2,100)} ; j=gsub("_",".",j); j=gsub("/$","",j); BIO[[j]]=system( paste( "cat ", i ), inter=T)} 
# }




########################################################################################################################################################


# Rename  needed
if (length(www)==1) if ( ! (length(CHR)==1 & all(SS1==names(CHR)) ) )  { names(CHR) = gsub(SS1,SS2,names(CHR)) }
if (length(www)==1) if ( ! (length(CHR)==1 & all(SS1==names(CHR)) ) )  { names(AGP) = gsub(SS1,SS2,names(AGP)) }
if (length(www)==1) if ( ! (length(CHR)==1 & all(SS1==names(CHR)) ) )  { names(CNF) = gsub(SS1,SS2,names(CNF)) }
if (length(www)==1) if ( ! (length(CHR)==1 & all(SS1==names(CHR)) ) )  { names(ASM) = gsub(SS1,SS2,names(ASM)) }
if (length(www)==1) if ( ! (length(CHR)==1 & all(SS1==names(CHR)) ) )  { names(LOG) = gsub(SS1,SS2,names(LOG)) }

# convert "." to "_"
# all(names(CHR) == SAMPLE$"short_name")
S_dif= setdiff(names(CHR),SAMPLE$"short_name")
S_mod= gsub("_",".",S_dif)

if (length(S_dif)>0) for ( ii in S_dif) { ii_mod= gsub("_",".",ii) ; SAMPLE$"short_name" = gsub(ii_mod,ii,SAMPLE$"short_name") ; SAMP_3$"short_name" = gsub(ii_mod,ii,SAMP_3$"short_name") ; cat(" ==> SAMPLE NAMES sprted by number :\n") ; print(SAMPLE)    }


# if no more than 1 name has not numeric part --> sort by number 
# sort by number if any
if( sum(is.na(SAMPLE$"numeric.sort")) <=1 ) { CHR = CHR [SAMPLE$"short_name"];  AGP = AGP [SAMPLE$"short_name"]; CNF = CNF [SAMPLE$"short_name"]; STA = STA [SAMPLE$"short_name"]; ASM = ASM [SAMPLE$"short_name"]; LOG = LOG [SAMPLE$"short_name"] }



#file.check = data.frame( CHR = names(CHR), AGP =names(AGP) , CNF = names(CNF), STA = names(STA), ASM = names(ASM)) 
file.check = data.frame( CHR = names(CHR), AGP =names(AGP) , CNF = names(CNF), ASM = names(ASM), LOG=names(LOG)) 
cat("\n")
cat("############################################################################\n")
cat(" ==> Files imported:\n")
cat("############################################################################\n")
cat("\n")
print(file.check)
cat("\n")
cat("\n")

########################################################################################################################################################


# agp.name == c("ref.chr","ref.start","ref.end","part_number","component_type","component_id/gap_length","component_beg/gap_type","component_end/linkage","orientation/linkage evidence")
agp.name = c("ref.chr","ref.start","ref.end","part_number","component_type","query","q.start","q.end","strand")
asm.name = c("query","q.len","q.start","q.end","strand","ref.chr","ref.len","ref.start","ref.end","l1","l2","V12","V13","V14","V15","V16","V17")
bio.name = c("superscaffold","sup.start","sup.end","part_number","component_type","query","q.start","q.end","strand")
rep.name = c()

# select assembly file from LOG
# LOG = lapply(LOG, function(x) tail( unlist( strsplit( grep("fasta",x, value=TRUE), " ")),1) )
# LOG = data.frame( "short_name" = names(LOG), "assembly_file" = as.character(LOG),stringsAsFactors=FALSE)
# LOG = merge(LOG, SAMP_3, by="short_name")

# select columns 1,2 and rename them
#lapply(CHR, head)
#print(RM.STRING)

CHR = lapply(CHR, function(x) x[,1:2])

########################################################################################################################################################


CHR = lapply(CHR, function(x) x[,1:2])
CHR = lapply(CHR, function(x) { names(x)=c("chr","len"); return(x)} )
AGP = lapply(AGP, function(x) { names(x)=agp.name ; return(x) } )


CHR = lapply(CHR, function(x) { x$"chr"= gsub("_RagTag","",x$"chr") ; return(x)}  )
AGP = lapply(AGP, function(x) { x$"ref.chr"= gsub("_RagTag","",x$"ref.chr") ; return(x)}  )

n_paf = unique(sapply(ASM,ncol))
if( all(n_paf==18)) { asm.name = c(asm.name,"V18")}



# rename AGP colums and add query.length
for ( cc in names(ASM)) { ASM.name = asm.name ; if ( ncol(ASM[[cc]]) > length(asm.name)  )  { ASM.name = c(ASM.name, paste0("V",17+(1:(ncol(ASM[[cc]])-17)) ) ) } ; names(ASM[[cc]])=ASM.name  }

# cat ("\n\n check ASM file \n\n")
# print(lapply(ASM, function(x) head(x$"ref.chr")))
# ASM = lapply(ASM, function(x) { x$"ref.chr"= gsub("_RagTag","",x$"ref.chr") ; return(x)}  )
# print(lapply(ASM, function(x) head(x$"ref.chr")))
# 
# 

# Rename  needed
if (length(www)==1) if ( ! (length(CHR)==1 & all(SS1==names(CHR)) ) )  { names(CHR) = gsub(SS1,SS2,names(CHR)) }
if (length(www)==1) if ( ! (length(CHR)==1 & all(SS1==names(CHR)) ) )  { names(AGP) = gsub(SS1,SS2,names(AGP)) }
if (length(www)==1) if ( ! (length(CHR)==1 & all(SS1==names(CHR)) ) )  { names(CNF) = gsub(SS1,SS2,names(CNF)) }
if (length(www)==1) if ( ! (length(CHR)==1 & all(SS1==names(CHR)) ) )  { names(ASM) = gsub(SS1,SS2,names(ASM)) }
if (length(www)==1) if ( ! (length(CHR)==1 & all(SS1==names(CHR)) ) )  { names(LOG) = gsub(SS1,SS2,names(LOG)) }

# convert "." to "_"
# all(names(CHR) == SAMPLE$"short_name")
S_dif= setdiff(names(CHR),SAMPLE$"short_name")
S_mod= gsub("_",".",S_dif)

if (length(S_dif)>0) for ( ii in S_dif) { ii_mod= gsub("_",".",ii) ; SAMPLE$"short_name" = gsub(ii_mod,ii,SAMPLE$"short_name") ; SAMP_3$"short_name" = gsub(ii_mod,ii,SAMP_3$"short_name") ; print(SAMPLE)    }



# if no more than 1 name has not numeric part --> sort by number 
# sort by number if any

if ( sort.sam == TRUE) 
{
 if( sum(is.na(SAMPLE$"numeric.sort")) <=1 ) 
   {
    CHR = CHR [SAMPLE$"short_name"];  
    AGP = AGP [SAMPLE$"short_name"]; 
    CNF = CNF [SAMPLE$"short_name"]; 
    STA = STA [SAMPLE$"short_name"]; 
    ASM = ASM [SAMPLE$"short_name"]; 
    LOG = LOG [SAMPLE$"short_name"] }
}

# CHR renaming 
if ( chr.rename == TRUE) 
{

cat("\n")
cat("############################################################################\n")
cat(" ==> replacing GenBank indentifier with Chromosome or Contig name from the fasta header  \n")
cat("############################################################################\n")
cat("\n")
    ref.edi2 = ref.edit
    names(ref.edi2) = c("ref.chr","chr.rename")
    
     
    ref = merge(ref,ref.edit, by="chr"    , all.x=T) ;    
    ref$"chr" =ref$"chr.rename";
    ref$"chr.rename"=NULL 
    
    # filter all files
    ref = ref[ ref$"ref.len" >= chr.filt ,]
    ref.edi2 = ref.edi2[ ref.edi2$"chr.rename" %in% ref$"chr", ]
    ref.edit = ref.edit[ ref.edit$"chr.rename" %in% ref$"chr", ]
    
    ref = ref[ rev(order(ref$"chr"))  , ]
    cat('\n    --> new Reference FAI   \n\n')
    print( ref) 


# rm chr to folter (still using old names)
    CHR = lapply(CHR, function(x) { x[ ! x$"chr" %in% chr_to_rm , ] }  )
    #CTG = lapply(CTG, function(x) { x[ ! x$"chr" %in% chr_to_rm , ] }  )
    ASM = lapply(ASM, function(x) { x[ ! x$"ref.chr" %in% chr_to_rm , ] }  )
    AGP = lapply(AGP, function(x) { x[ ! x$"ref.chr" %in% chr_to_rm , ] }  )
    CHR = lapply(CHR, function(x) { x=merge(x,ref.edit, by="chr"    , all.x=T) ; w_rep = which( ! is.na(x$"chr.rename")) ; x$"chr"     [w_rep] =x$"chr.rename"[w_rep] ; x$"chr.rename"=NULL ; return(x) } )
    AGP = lapply(AGP, function(x) { x=merge(x,ref.edi2, by="ref.chr", all.x=T) ; w_rep = which( ! is.na(x$"chr.rename")) ; x$"ref.chr" [w_rep] =x$"chr.rename"[w_rep] ; x$"chr.rename"=NULL ; return(x) } )
    ASM = lapply(ASM, function(x) { x=merge(x,ref.edi2, by="ref.chr", all.x=T) ; w_rep = which( ! is.na(x$"chr.rename")) ; x$"ref.chr" [w_rep] =x$"chr.rename"[w_rep] ; x$"chr.rename"=NULL ; return(x) } )  #  BIO no need
# #       
cat("\n")
cat("\n")
      cat("   --> Conversion name check: \n")    
      cat("     -> ragtag FAI \n")    
      print(sapply(CHR, function(x) { paste( sort(unique((x$"chr"))), collapse= " ")     } ))
      cat("     -> ragtag AGP \n")    
      print(sapply(AGP, function(x) { paste( sort(unique((x$"ref.chr"))), collapse= " ")  } ))
      cat("     -> ragtag ASM \n")    
      print(sapply(ASM, function(x) { paste( sort(unique((x$"ref.chr"))), collapse= " ")  } ))
cat("\n")
      
      chr_prefix = unique(substring(ref$"chr",1,3))
      cat("   --> New Chromosome prefix:",chr_prefix,"\n")    
     }





#file.check = data.frame( CHR = names(CHR), AGP =names(AGP) , CNF = names(CNF), STA = names(STA), ASM = names(ASM)) 
file.check = data.frame( CHR = names(CHR), AGP =names(AGP) , CNF = names(CNF), ASM = names(ASM), LOG=names(LOG)) 
cat("\n")
cat("############################################################################\n")
cat(" ==> Files imported:\n")
cat("############################################################################\n")
cat("\n")
print(file.check)
cat("\n")
cat("\n")




# select assembly file from LOG
LOG = lapply(LOG, function(x) tail( unlist( strsplit( grep("fasta",x, value=TRUE), " ")),1) )
LOG = data.frame( "short_name" = names(LOG), "assembly_file" = as.character(LOG),stringsAsFactors=FALSE)
LOG = merge(LOG, SAMP_3, by="short_name")

# select columns 1,2 and rename them
#lapply(CHR, head)
#print(RM.STRING)

CHR = lapply(CHR, function(x) { x=as.data.frame(x) ; return(x)}  )
CHR = lapply(CHR, function(x) { x=x[ order(x$"chr"),]  ; return(x)}  )
CHR = lapply(CHR, function(x) { x=setDT(x) ; return(x)}  )

# print((ASM))
# print(length(asm.name))
# print(ncol(ASM[[1]]))
# print(2)
# 


#ASM = lapply(ASM, function(x) { names(x)=asm.name ; return(x) } )
ASM = lapply(ASM, function(x) { x$"ctg.len"=  x$"q.end"- x$"q.start"+1; x$"scaf.perc"= round( x$"ctg.len"/x$"q.len",2);  return(x) } )
ASM = lapply(ASM, function(x) { x$"ctg.ref"=  x$"ref.end"- x$"ref.start"+1;   return(x) } )

AGP = lapply(AGP, function(x) { w = which(x$"query" !="100"); x$"q.len"[w]= as.numeric(x$"q.end"[w])- as.numeric(x$"q.start"[w])+1; x$"q.len"[-w] =100; return(x) } )

AGP = lapply(AGP, function(x) { x=as.data.frame(x) ; return(x)}  )
AGP = lapply(AGP, function(x) { x=x[ order(x$"ref.chr"),]  ; return(x)}  )
AGP = lapply(AGP, function(x) { x=setDT(x) ; return(x)}  )

ASM = lapply(ASM, function(x) { x=as.data.frame(x) ; return(x)}  )
ASM = lapply(ASM, function(x) { x=x[ order(x$"ref.chr"),]  ; return(x)}  )
ASM = lapply(ASM, function(x) { x=setDT(x) ; return(x)}  )

if( bionano ==TRUE ) BIO = lapply(BIO, function(x) { names(x)=bio.name ; return(x) } )


#print(BIO)

# 
# # ref repeat masker
# if ( is.na(rm.scaf)==FALSE )
# { 
#     rep.files = paste( chr.files,".gff")
#     REP = list()
#     for( i in rep.files ) { j=gsub( paste0(".ag.ag..caffold.",RM.STRING,"|/ragtag.scaffold.fasta.gff") ,"",i) ; if(j=="") {j=SS} ; if( substr(j,1,1) %in% c("_",".","-")) {j=substr(j,2,100)} ; REP[[j]]=fread(i)}
#     REP = lapply(REP, function(x) { names(x)=asm.name ; return(x) } )
#     }
#     

    

cat("\n")
cat("############################################################################\n")
cat(" ==> assembly files from LOG files \n")
cat("############################################################################\n")
cat("\n")
print( LOG )
cat("\n")
cat("\n")

cat("############################################################################\n")
cat(" ==> FAI files first + last rows \n")
cat("############################################################################\n")
cat("\n")
print( cbind( data.frame(t(sapply(CHR, head,1))) , data.frame(t(sapply(CHR, tail,1)))) )
cat("\n")
cat("\n")

cat("############################################################################\n")
cat(" ==> AGP files first rows \n")
cat("############################################################################\n")
cat("\n")
print( data.frame(t(sapply(AGP, head,1))))
cat("\n")
cat("\n")


cat("############################################################################\n")
cat(" ==> PAF files first rows \n")
cat("############################################################################\n")
cat("\n")
print( data.frame(t(sapply(ASM, head,1))))
cat("\n\n")



# Separate Scaffold and unplaced contigs
CHR.ctg = CHR
CTG = lapply(CHR.ctg, function(x) x[ ! x$"chr" %in% ref$"chr" ,])
CTG = lapply(CTG    , function(x) x[  rev(order(x$len)      ) ,])
CHR = lapply(CHR.ctg, function(x) x[   x$"chr" %in% ref$"chr" ,])
#for( i in names(CHR)) { CHR.ctg[[i]] = rbind(CHR[[i]], CTG[[i]]) }

CTG.max = max(sapply(CTG, function(x) sum(x$"len")))


# optional, select only scaffold (and not unplaced contigs or Chr0)

if( rm.unpl_ctg == TRUE & length(unique(chr_prefix))==1 ) { cat("\n ==> Selected CHRs/Scaffolds only!  \n") ; CHR = lapply(CHR, function(x) x[ grep( chr_prefix , x$"chr"),]) }



if( length(chr_to_rm) >0 )
{
CHR = lapply(CHR, function(x) { x[ ! x$"chr" %in% chr_to_rm , ] }  )
#CTG = lapply(CTG, function(x) { x[ ! x$"chr" %in% chr_to_rm , ] }  )
ASM = lapply(ASM, function(x) { x[ ! x$"ref.chr" %in% chr_to_rm , ] }  )
AGP = lapply(AGP, function(x) { x[ ! x$"ref.chr" %in% chr_to_rm , ] }  )
    }


# adding extra colors
nnn = length(CHR) - length(COL.LIST)
if( nnn>0) { COL.LIST = c(COL.LIST,rainbow(nnn)) ;     cat("\n## -->  ADDING NEW COLORS --> ", paste(rainbow(nnn), collapse=" "),"\n") }


# FIND THE MAX CHR LENGTH for the plot
CCC = CHR 
for (i in names(CHR)) names(CCC[[i]]) = c("chr",i)
CC1 = ref
#CC1$"chr" = paste0(CC1$"chr","_RagTag")
for( i in names(CCC)) CC1 = merge(CC1, CCC[[i]],by="chr", all.x=TRUE)
CCmax = max(CC1[,-1], na.rm=T)

#  print(ref)
#  print(CHR)
#  print(CCC)
#  print(CC1)
#  print(CCmax)

STAT = list()
GAPS = list()


#choose size of multiplot pannel (good up to 20 samples)
N.FILES = length(CHR)
# N.FILES = N.FILES +1

n.row=ceiling(N.FILES/3)
n.col=ceiling(N.FILES/2)
if( N.FILES <= (n.row*n.col)-n.col ) n.row=n.row-1
if( N.FILES <= (n.row*n.col)-n.col ) n.row=n.row-1
if( N.FILES <= (n.row*n.col)-n.row ) n.col=n.col-1
if( N.FILES <= (n.row*n.col)-n.row ) n.col=n.col-1
if( N.FILES <= (n.row*n.col)-n.col ) n.row=n.row-1
if( N.FILES <= (n.row*n.col)-n.col ) n.row=n.row-1
if( N.FILES <= (n.row*n.col)-n.row ) n.col=n.col-1
if( N.FILES <= (n.row*n.col)-n.row ) n.col=n.col-1
if( N.FILES <= (n.row*n.col)-n.row ) n.col=n.col-1
if( N.FILES  > (n.row*n.col)       ) n.col=n.col+1


cat("\n")
cat("############################################################################\n")
cat("###                   1) Chromosome length Barplots                      ###\n")
cat("############################################################################\n")


cat("\n ==> Plot panel: ",N.FILES,"files ;",n.row,"rows",n.col,"columns \n")

# PLOT_TYPES = c("single", "multiple","separate")

for (plot_type in PLOT_TYPES)
{
if ( plot_type == "single"  ) { cat("\n ==> Plot type:   ",plot_type,"") ; pdf( paste0(out.dir,"/GS-viewer_01.Scaffolds_length_and_gaps.single.",OUT,"pdf"  ) ,width=plot.w, height=plot.h) ; par(mfcol=c(1    ,1    ), mar=c(2,10,4,1) , oma=c(1,4,1,1) ); LG.size=1   }
if ( plot_type == "multiple") { cat("\n ==> Plot type: "  ,plot_type,"") ; pdf( paste0(out.dir,"/GS-viewer_01.Scaffolds_length_and_gaps.multiple.",OUT,"pdf") ,width=plot.w, height=plot.h) ; par(mfrow=c(n.row,n.col), mar=c(3,11,5,3) , oma=c(1,4,1,1) ); LG.size=0.5 }

    for( i in names(CHR))
    {
        # i =  names(CHR)[3];
        # long_i = 
        if ( plot_type == "separate")  { cat("\n ==> Plot type: ",plot_type,"") ; pdf( paste0(out.dir,"/GS-viewer_01.Scaffolds__length_and_gaps.",gsub(" ","_",i),".",OUT,"pdf") ,width=plot.w, height=plot.h) ; par(mfcol=c(1    ,1    ), mar=c(5,10,5,3), oma=c(1,4,1,1) ) }
        cat(i,"")
        
        chr = data.frame(CHR[[i]]) ; # CHR length
        agp = data.frame(AGP[[i]]) ; # AGP
        cnf = data.frame(CNF[[i]]) ; # CONFIDENCE of contig placement on the reference

     #   print(agp)
        
        chr$"chr" = gsub("_RagTag","",chr$"chr")

        # add missing chromosome
        miss_chr = setdiff(ref$"chr", chr$"chr")
        if ( length(miss_chr)>0) { chr=rbind(chr, data.frame("chr" =miss_chr, len=0 , stringsAsFactors=FALSE)) }

        chr = chr[ rev(order(chr$"chr")) ,]   # invert order for the plot (top to bottom)
        chr = chr[ chr$"chr" != "Chr0"  , ] # rm Chr=0
        
        # add Chr0 = concatenated ctg unplaced
        if (chr0_plot == TRUE)
        {
            ctg = data.frame(CHR.ctg[[i]]) ; # CHR length
            ctg = ctg [ ! ctg$"chr" %in% ref$"chr",]  # rm Chr=0
            ctg = ctg [ rev(order(ctg$"len")), ]
            ctg$"gap.pos"= cumsum(as.numeric(ctg$"len"))
            ch0 = data.frame( "chr"="other.Ctgs", len = sum(ctg$"len"), stringsAsFactors=FALSE)
            if( CTG.max > CCmax ) CCmax = CTG.max
            chr = chr [ chr$"chr" %in% ref$"chr",]         
            chr = rbind(ch0,chr)
            }
            
        # create barplot and get coordinates of each chromosome in the plot 
        
        BLUE = rgb(0.68 , 0.85, 0.90,   alpha=0.75)
        GRAY = rgb(0.00 , 0.00, 0.00,   alpha=0.40)
        
    
        lll =  0.9+(0:(nrow(chr)-0)*1.4)
        bar=barplot(chr$"len", names.arg=chr$"chr", las=1, horiz=T, main=i, cex.main=4 , cex.axis=2, ,xlim=c(0, CCmax*1.26), ylim=c(0,max(lll)*0.99 ), cex.names=2.4, space=0.4, col=GRAY)
        b2 = data.frame(chr.pos = bar[,1])
        b2$"chr" = chr$"chr"
        rrr = merge(ref, b2,by="chr", all.x=TRUE) # add reference chr length
        ch2 = merge(chr, b2,by="chr", all.x=TRUE) # add reference chr length
        
      # AGP file
        agp=  agp [  agp$"ref.chr" %in% ref$"chr",] # rm GAP lines
        agp$"ref.chr" = as.factor( agp$"ref.chr")
        gap=  agp[ agp$"query"=="100" ,] # rm CTG lines
        agp=  agp[ agp$"query"!="100" ,] # rm GAP lines
        agp$"chr" = gsub("_RagTag","",agp$"ref.chr")
        agp$"ref.mid" = (agp$"ref.start" + agp$"ref.end")/2  # mid point of the contig on the reference 
        agp = merge(agp, cnf, by="query", all.x=TRUE)        # add confidence columns
        agp = merge(agp, b2, by="chr", all.x=TRUE)           # add coordinated on the plot
        agp = agp[ order(agp$"ref.chr"), ]
        
        # Simply contig names by remove stretched of 0
        simp=all(grepl("0000",agp$"query")) ;  if (simp==TRUE) { agp$"qq" = gsub("0000",".", agp$"query" )} ; #print(simp)
        simp=all(grepl("000" ,agp$"query")) ;  if (simp==TRUE) { agp$"qq" = gsub("000" ,".", agp$"query" )}  #print(simp)
        simp=all(grepl("00"  ,agp$"query")) ;  if (simp==TRUE) { agp$"qq" = gsub("00"  ,".", agp$"query" ) ; agp$"qq" = gsub("[.][.]"  ,".", agp$"qq" )}  #print(simp)
              
        # Reference chr length
      #  segments(rep(0, nrow(chr)),rrr$"chr.pos", rrr$"ref.len", rrr$"chr.pos", col="gray", lwd=5,lend=1)
     #   segments(rrr$"ref.len", rrr$"chr.pos"-0.5, rrr$"ref.len",rrr$"chr.pos"+0.5, col="blue", lwd=2,lend=1)
     for( iii in 1:nrow(rrr))   polygon( c(0 , rrr$"ref.len"[iii],rrr$"ref.len"[iii],0) ,c( rrr$"chr.pos"[iii]+0.6, rrr$"chr.pos"[iii]+0.6, rrr$"chr.pos"[iii]-0.6 ,rrr$"chr.pos"[iii]-0.6) , lwd=2 , col=BLUE , border="blue", lwd=1, lend=1)   
      #  segments(ch2$"len"    , ch2$"chr.pos"-0.5, ch2$"len",    ch2$"chr.pos"+0.5, col="black", lwd=1,lend=1)
        
        # add GAPS to the plot
      #  points(agp$"ref.end", agp$"chr.pos", type="p", cex=1 , pch=19,col="black")

        segments(agp$"ref.end", agp$"chr.pos"-0.5, agp$"ref.end",agp$"chr.pos"+0.5, col="black", lwd=2, lend=1)
      #  segments(rrr$"ref.len", rrr$"chr.pos"-0.5, rrr$"ref.len",rrr$"chr.pos"+0.5, col="blue", lwd=4,lend=2,lend=2)
        if (chr0_plot == TRUE) if( nrow(ctg)>0) { ctg$"chr.pos"= b2[ b2$"chr"=="other.Ctgs",]$"chr.pos"; segments(ctg$"gap.pos", ctg$"chr.pos"-0.5, ctg$"gap.pos",ctg$"chr.pos"+0.5, col="black", lwd=1)  }

        par(new = T,lwd = 2)
        bar=barplot(chr$"len", names.arg=NA, las=1, horiz=T, cex.axis=2, ,xlim=c(0, CCmax*1.26), ylim=c(0,max(lll)*0.99 ), cex.names=2.4, space=0.4, col=NA, border="black")
        
    # Confidence for contig association to refencen      
        ag1 = agp [ agp$"chr" !="Chr0", ]
        ag1$"grp.conf" = round(ag1$"grouping_confidence",2)
        ag1$"loc.conf" = round(ag1$"location_confidence",2)
        ag1.g = ag1 [ ag1$"grouping_confidence" < thr.conf, ]
        ag1.l = ag1 [ ag1$"location_confidence" < thr.conf, ]
        CEX.txt = 1
        if( plot.conf ==TRUE & nrow(ag1.g)>0) segments(ag1.g$"ref.start",ag1.g$"chr.pos"-0.05, ag1.g$"ref.end", ag1.g$"chr.pos"-0.05, col="darkgreen", cex=3, lwd=2)
        if( plot.conf ==TRUE & nrow(ag1.g)>0)      text(ag1.g$"ref.mid", ag1.g$"chr.pos"-0.25, ag1.g$"grp.conf", col="darkgreen", cex=CEX.txt)
        if( plot.conf ==TRUE & nrow(ag1.l)>0) segments(ag1.l$"ref.start",ag1.l$"chr.pos"+0.05, ag1.l$"ref.end", ag1.l$"chr.pos"+0.05, col="darkred", cex=3, lwd=2)
        if( plot.conf ==TRUE & nrow(ag1.l)>0)      text(ag1.l$"ref.mid", ag1.l$"chr.pos"+0.25, ag1.l$"loc.conf", col="darkred", cex=CEX.txt)
        if( plot.conf ==TRUE         )             points(agp$"ref.end", agp$"chr.pos", type="p", pch="|", cex=1.5 ,col="black")
               
        gan = sapply(  split(gap, gap$"ref.chr"), nrow)
        gan = data.frame( "ref.chr"= names(gan) , gaps=as.numeric(gan))
        
              
        
        gan$"chr"= gsub("_RagTag","",gan$"ref.chr")
        gap$"chr"= gsub("_RagTag","",gap$"ref.chr")
        gan = merge(gan, rrr,by="chr", all.x=TRUE)
        gan = merge(gan, chr,by="chr", all.x=TRUE)
        gan$"max"= apply(gan[,c("len","ref.len")],1,max)
        gan$"perc"= round(100*gan$"len"/gan$"ref.len",3)
        gan$"per1"= round(100*gan$"len"/gan$"ref.len",0)
        
        gap = merge(gap, b2,by="chr", all.x=TRUE)
        gap = gap[ gap$"chr" !="Chr0",]
               
        stat=gan[,c("chr","gaps","ref.len","len","perc")] 
        s.all = data.frame( t(colSums(stat[,-1])))
        s.all$"perc"= round(100*s.all$"len"/s.all$"ref.len",3)
        s.all = cbind( data.frame( chr="all"), s.all)
        
        stat = rbind(stat, s.all)
        stat$"sam" = i
        
        GPmax = max(c(chr$"len",rrr$"ref.len"), na.rm=T)
        GPma1 = CCmax/19
        GPma2 = GPma1*3.5
        CEX.bar = 2
        if ( plot_type %in% c("multiple")) CEX.bar = 2.75
        
        GXmax =max(lll)
        GXmin = 0

        STAT[[i]] = stat
        GAPS[[i]] = gap
        
        tot.gap = sum(gan$"gaps")
        tot.per = round(100*sum(gan$"len")/sum(gan$"ref.len"),1)
        
#         
#         print(GPmax)
#         print(GPma1)
#         print(gan$chr.pos)
#         print(gan$gaps)
#          print(gap)
        
      if( nrow(gan)>0) text(GPmax+GPma1, gan$"chr.pos", gan$"gaps", col="darkred", cex=CEX.bar ) else 
      if( nrow(gan)>0) text(GPmax+GPma2, gan$"chr.pos", gan$"per1", col="darkblue", cex=CEX.bar )       
                       text(GPmax+GPma1, GXmax, "Gaps", col="darkred", cex=CEX.bar )
                       text(GPmax+GPma2, GXmax, "%Ref.len", col="darkblue", cex=CEX.bar )
                       text(GPmax-GPma1, GXmin, "tot", col="black", cex=CEX.bar )
                       text(GPmax+GPma1, GXmin, tot.gap, col="darkred", cex=CEX.bar )
                       text(GPmax+GPma2, GXmin, tot.per, col="darkblue", cex=CEX.bar )
                       legend( x= CCmax * 0.8, y= max(lll) * 0.5, c("Reference", "Scaffold"), fill=c(BLUE,GRAY), cex=3*LG.size)

    if ( plot_type == "separate")  dev.off()  
        }
  if ( plot_type %in% c("single","multiple")) dev.off()
   }
        

# system ( "evince ../../../Ragtag_scaffolds_and_gaps.multiple.pdf ")

STAT.GAP = as.data.frame(t(sapply(STAT, function(x) sum(x$"gaps"))))
write.table( STAT.GAP, file = paste0(out.dir,"/Gaps_number.",OUT,"tsv"  ) , col.names=TRUE, row.names=F, quote=F, sep="\t")

cat("\n\n")

cat("############################################################################\n")
cat("###         Calculating Scaffold coverage on reference genome            ###\n")
cat("############################################################################\n")

##  chr comparisons
# COL.LIST=c("firebrick3","black","dodgerblue4","darkorange3","darkturquoise","blueviolet","forestgreen","lightcoral","gold")

col.data = data.frame ("assembly" = names(CHR), stringsAsFactors=FALSE)
col.data$"col"= COL.LIST[1:nrow(col.data)]
rownames(col.data) = col.data$"assembly"
based.asm =  names(CHR)

A.list= list()
B.list= list()
g.list= list()

for( ii in based.asm ) 1

for( ii in based.asm )
{
# ii =  "hifiasm.1Mb.ch0";
# ii =  "hicanu.1Mb.ch0";
cat("\n",ii,"")

chr = data.frame(CHR[[ii]]) ; 
agp = setDT(AGP[[ii]]) ; 
cnf = setDT(CNF[[ii]]) ; 
asm = setDT(ASM[[ii]]) ; 
chr = chr[ rev(order(chr$"chr")) ,]
chr$"chr" = gsub("_RagTag","",chr$"chr")

CCmax.max = max(CCmax, chr$"len")

# separate GAPs annotation from AGP
agp$"chr" = gsub("_RagTag","",agp$"ref.chr")
gap = agp[ agp$"query" =="100",]
agp = agp[ agp$"query" !="100",]

# select valid contigs on ASM from AGP 

agp2 = split( agp, agp$"chr")
asm2 = split( asm, asm$"ref.chr")
for (i in names(asm2)) { asm2[[i]] = asm2[[i]] [ asm2[[i]]$"query" %in% agp2[[i]]$"query",    ]  }

for (i in names(asm2)) { 
    cat(i,"")
    a = asm2[[i]]  
    g = agp2[[i]]  
    
    if (nrow(a)> 0 & length(g)> 0)
    {        
        # contigs on ref
        x = 0.02
        
        a$"ya" = x*2
        a$"yb" = a$"ya"
        a$"col"="blue"
        if( any(a$"strand"=="-")) a[ a$"strand"=="-",]$"yb"= a[ a$"strand"=="-",]$"yb" *-1
        if( any(a$"strand"=="-")) a[ a$"strand"=="-",]$"col"="darkred"
        if( any(a$"strand"=="+")) a[ a$"strand"=="+",]$"col"="blue"
        a$"strand.query" = paste(a$"query", a$"strand")
        as= split(a,a$"strand.query")
        n.a = length(as)
        col.a = rainbow(n.a)
        chr.max =a$"ref.len"[1]
        chr.txt = -chr.max/8
        n.ctg = length(unique(a$query))
        for (j in 1:length(as)) { as[[j]]$"ya"= as[[j]]$"ya"+j/100  }
        
        a= do.call(rbind,as)
        
        # scaffold
        g$"za"= -x/2
        g$"zb"= g$"za" -x*2
        g$"zm"=g$"zb"-g$"za"
        g$"zs"=g$"zb"-0.02
        g$"ref.med"=(g$"ref.start"+g$"ref.end")/2
        g$"ref.ctg_start"=g$"ref.med"-g$"q.len"/2
        g$"ref.ctg_end"  =g$"ref.med"+g$"q.len"/2
        
        g$"col"="blue"
        g$"QUERY"=g$"query"
        if( any(g$"strand"=="-")) g[ g$"strand"=="-",]$"col"="darkred"
        if( any(g$"strand"=="+")) g[ g$"strand"=="+",]$"col"="blue"
        if( any(g$"strand"=="-")) g[ g$"strand"=="-",]$"QUERY"= paste0(g[ g$"strand"=="-",]$"QUERY", ".rev")
        
        # invs 
         as2= lapply( split(a$"strand",a$"query"), unique )
         qg= split(a,a$"strand.query")
    
        # contig names
        col.aa = c("query","q.len","strand.query","ya","yb","col")
        aa = a[ ,..col.aa]
        aa = aa[ ! duplicated(aa),]
        all( g$"q.len" %in% aa$"q.len" )
        
        gg = g[, c(2,3,6)] 
        names(gg) = c("scaf.start","scaf.end","query")
        a = merge(a, gg, by="query",all.x=T)
        a$"scaf.pos"  =a$"scaf.start"+a$"q.start"-1
        a$"scaf.end"  =a$"scaf.start"+a$"q.end"-1
     
    
      # A = a[ a$"query" %in% "utg000010l",]
        A=a
        A = A[ A$"ctg.len" > MIN.WIND , ] #before it was 30000
        A = A[ order(A$"ref.start") , ]
        A$"COL" = gsub("darkred","firebrick3",A$"col") 
        A$"COL" = gsub("blue","cornflowerblue",A$"COL") 
        y1 = x 
        y2 = x/2 
               
           # invert Strand of  MINUS contigs
        
        min.contigs = g[ g$"strand"=="-", ]$"query"
      
        B = A
       # B = B[ B$"query" == "utg000024l" , ]
        B$"Q.START"  = B$"q.start" 
        B$"Q.END"    = B$"q.end" 
        B$"SCAF.POS" = B$"scaf.pos"
        B$"SCAF.END" = B$"scaf.end" 
        w.min= which(B$"query" %in% min.contigs)
        
        if ( length(min.contigs)>0 )
        {  
            B$"Q.START" [w.min]  = B$"q.len"[w.min] - B$"q.end"[w.min] +1
            B$"Q.END"   [w.min]  = B$"q.len"[w.min] - B$"q.start"[w.min] +1
            B$"SCAF.POS"[w.min] = B$"scaf.start"[w.min] + B$"Q.START"[w.min] -1
            B$"SCAF.END"[w.min] = B$"scaf.start"[w.min] + B$"Q.END"[w.min] -1
#            B$"COL"[w.min]  = gsub("firebrick3","pink",B$"COL"[w.min] ) 
#            B$"COL"[w.min]  = gsub("cornflowerblue","lightskyblue",B$"COL"[w.min] ) 
            }
        } else { A =B = g = data.frame( )}
    A.list [[i]] [[ii]] = A
    B.list [[i]] [[ii]] = B
    g.list [[i]] [[ii]] = g   
       }
   }

names(A.list) = gsub("_RagTag","",names(A.list))
names(B.list) = gsub("_RagTag","",names(B.list))
names(g.list) = gsub("_RagTag","",names(g.list))

cat("\n\n")
cat("############################################################################\n")
cat("###     2) Scaffold coverage of all asseblies on reference genome        ###\n")
cat("############################################################################\n")

##  CHANGE THE ORDER of plot ( high cov to low cov) ??

#   CTG.min = 000
#pdf( paste0(out.dir,"/GS-viewer_02.Scaffold_Coverage_on_Reference.",OUT,"min",CTG.min,"bp.pdf"),width=plot.w, height=plot.h )
pdf( paste0(out.dir,"/GS-viewer_02.Scaffold_Coverage_on_Reference.",OUT,"pdf"),width=plot.w, height=plot.h )
par( mar=c(2,12,4,1) )

# chr = chr[ chr$"chr" != "Chr0", ]
bar=barplot(ref$"ref.len", names.arg=ref$"chr", las=1, horiz=T, main="Scaffold coverage on Reference Genome",col=gray(0.95), cex.main=3, xlim=c(0,max(ref$"ref.len")*1.15),  cex.axis=2, cex.names=2)
b2 = data.frame(chr.pos = bar[,1])
b2$"chr" = ref$"chr"
rrr = merge(ref,b2, by="chr")

text( max(ref$"ref.len")/2, max(b2$"chr.pos")*1.06,  paste( "min plotted window =",CTG.min,"bp"), cex=0.5 )

sh = max(ref$"ref.len")*0.01

n.ass = length( A.list[[1]] )

shifts = rev(seq(-0.5,0.5,length.out=n.ass+1))

for( p in shifts ) segments(rep(0, nrow(rrr)),rrr$"chr.pos"+p, rrr$"ref.len", rrr$"chr.pos"+p, col="gray", cex=2)

for (i in names(A.list)) 1==1

A.coo.list=list()

for (i in names(A.list)) { 
#  i=  names(A.list)[1]
    A.list [[i]]  -> AA
    B.list [[i]]  -> BB
    g.list [[i]]  -> gg
        
    bb = b2[ b2$"chr" == i, ]
    bb$"chr.pos"
    
    CEX.txt = 1.2
    if(length(CHR)>7) CEX.txt = 0.6

    x = 0.02
    cat("\n",gsub(" ","_",i),"")
     
    for (k in length(AA):1) 
    {  
        A=AA[[k]] ;
        A = A[ A$"ctg.len" > CTG.min , ]
        ass=names(AA)[k]
        cat(ass,"")
        COL =  unique(col.data [ col.data$"assembly" == ass ,]$"col" )
        s1=shifts[k];
        s2=shifts[k+1];
        y1 = bb$"chr.pos" + s1
        y2 = bb$"chr.pos" + s2
        y3 = (y1+y2)/2
      #   cat("\n")
     #    print(head(A,3))
      #   print(bb)
#         print(y1)
#         print(y2)
#         
        if( nrow(A)>0)
        {
         
            A.coord = list()
            for (j in 1:nrow(A)) A.coord[[j]]  = c(A$"ref.start"[j]: A$"ref.end"[j])           
            A.coord = unique(do.call(c,A.coord))
            A.perc =round(100* length(A.coord)/A$"ref.len"[1])
            A.coo.list [[ass]] [[i]] = A.coord
                    
           # A = A[ A$"ctg.len" > CTG.min , ]
    
            for (j in 1:nrow(A)) { polygon( c(A$"ref.start"[j], A$"ref.end"[j], A$"ref.end"[j],A$"ref.start"[j]) ,c( y2 , y2, y1, y1) , lwd=2 , col=COL , border=NA)    }
            text( A$"ref.len"[1]+sh  , y3 ,  paste0(A.perc,"%"), cex=CEX.txt ,col=COL)
            }   else { A.coo.list [[ass]] [[i]] = data.frame()}       
        }
}


#COL.LIST=c("firebrick3","black","dodgerblue4","darkorange3","darkturquoise","blueviolet","forestgreen","lightcoral")
col.list = COL.LIST[ 1:length(CHR)]

colors= data.frame(assembly=names(A.list [[1]] ), col=col.list[1:length(A.list [[1]] )], row.names=names(A.list [[1]] ), stringsAsFactors=FALSE)

colors = col.data

chr.sum = sapply(A.coo.list,sapply,length)
chr.sum = as.data.frame(chr.sum)
gen.sum = sapply(chr.sum,sum)
gen.sum = data.frame( assembly=names(gen.sum), len=gen.sum, row.names=names(gen.sum), stringsAsFactors=FALSE)
colors  = merge( colors, gen.sum, by="assembly")
colors$"perc" = round( 100* colors$"len" / sum(ref$"ref.len") ,2)
colors$"tag" = paste0(colors$"assembly","        ",colors$"perc","%")
colors$"ta2" = paste0(colors$"assembly","                 ")
colors$"ta3" = paste0(colors$"perc","%")

rownames(colors) = colors$"assembly"
colors = colors[  names(A.list[[1]]) , ]
CEX.leg=2.5
legend("right", gsub(".ch0","",colors$"ta2"), cex=CEX.leg, fill=colors$col  )
#legend("bottomright", gsub(".ch0","",colors$"tag"), fill=colors$col  )
legend("right", gsub(".ch0","",colors$"ta3"),cex=CEX.leg, bty="n" )

 dev.off()       

 cat("\n")
#cat("\n") ; print(colors)
    
# 
# 

cat("############################################################################\n")
cat("### 3)  All Scaffold length and gaps in one Figure                       ###\n")
cat("############################################################################\n")

cat("\n  --> Scaffold length","\n")

CHR2 = ref
if ( any(grepl("RagTag",do.call(rbind,CHR)$"chr") )) CHR2$"chr" = paste0(CHR2$"chr","_RagTag")
for (i in names(CHR)) { x=CHR[[i]] ; names(x)=c("chr",i) ; CHR2 = merge(CHR2,x,by="chr", all.x=TRUE) }
# col names may have changed
CHR2 = as.data.frame(CHR2)
CHR2 [ is.na(CHR2) ] = 0
CHR2$"chr" = gsub("_RagTag","",CHR2$"chr")
rownames(CHR2) = CHR2$"chr" 
CHR2$"ref.len" = NULL
CHR2$"chr" = NULL
CHR2 = t(as.matrix(CHR2))

CHR2.colnames= colnames(CHR2)
CHR2.rownames= rownames(CHR2)

CEX.GAP = 1
if(length(chr.files)==1) CEX.GAP = 2
if(length(chr.files)==2) CEX.GAP = 1.5


if(nrow(CHR2) >1) { CHR2 = CHR2[ ,  gsub("_RagTag","",ref$"chr") ] }
if(nrow(CHR2)==1) { CHR2 = t(as.matrix(CHR2[ ,  gsub("_RagTag","",ref$"chr")])) ; rownames(CHR2)= CHR2.rownames }


if(nrow(CHR2)>1) CHR2 = CHR2[   rev(rownames(CHR2)) , ]
if(ncol(CHR2)>1) CHR2 = CHR2[ , rev(colnames(CHR2)) ]

if( is.null(dim(CHR2)) )  { CHR2= as.matrix(t(CHR2)) ;  dimnames(CHR2)[[1]] = CHR2.rownames } 

CCmax = max(CC1[,-1], na.rm=T)
# print(ref$"ref.len")
# print(CHR2)
# print(CCmax)


#   CTG.min = 000
#pdf( paste0(out.dir,"/GS-viewer_03.Scaffold_Length_vs_Reference.",OUT,"min",CTG.min,"bp.pdf"), width=plot.w, height=plot.h )
pdf( paste0(out.dir,"/GS-viewer_03.Scaffold_Length_vs_Reference.",OUT,"pdf"), width=plot.w, height=plot.h )
par( mar=c(2,12,4,1))
# chr = chr[ chr$"chr" != "Chr0", ]
bar=barplot(ref$"ref.len", names.arg=ref$"chr", las=1, horiz=T, main="Scaffold length vs Reference length",col=gray(0.95), cex.main=3, xlim=c(0, max(c(ref$"ref.len",CHR2,CCmax)*1.20 )) ,  cex.axis=2, cex.names=2)
b2 = data.frame(chr.pos = bar[,1])
b2$"chr" = ref$"chr"
rrr = merge(ref,b2, by="chr")
rrr = rrr[ order(rrr$"chr.pos"),]

sh = max(ref$"ref.len")*0.01

n.ass = length( A.list[[1]] )

shifts = rev(seq(-0.5,0.5,length.out=n.ass+1))

#for( p in shifts ) segments(rep(0, nrow(rrr)),rrr$"chr.pos"+p, rrr$"ref.len", rrr$"chr.pos"+p, col="gray", cex=2)


GAP2 = list()
for (i in names(A.list)) 1==1
for (i in names(A.list)) { 
#  i=  names(A.list)[4]
    AA = CHR2[,i]
    if( nrow(CHR2)==1 ) { names(AA) = rownames(CHR2)}
    AA = AA[ rev(names(AA)) ]
    bb = rrr[ rrr$"chr" == i, ]
    bb$"chr.pos"
         
    CEX.txt = 1.1
    if(length(CHR)>7) CEX.txt = 0.6

    x = 0.02
    cat("\n",gsub(" ","_",i),"")
     
    
    for (k in length(AA):1) 
    {  
        A=AA[[k]] ;
        ass=names(AA)[k]
        cat(ass,"")
        COL =  unique(col.data [ col.data$"assembly" == ass ,]$"col" )
        s1=shifts[k];
        s2=shifts[k+1];
        y1 = bb$"chr.pos" + s1
        y2 = bb$"chr.pos" + s2
        y3 = bb$"chr.pos" + (s1+s2)/2
        Y1 = bb$"chr.pos" - 0.5 -0.01
        Y2 = bb$"chr.pos" + 0.5 +0.01
        CL = bb$"ref.len" 
        #cat(y3,"")

        if (nrow(bb)>0) 
        {   polygon( c(0 , A, A,0) ,c( y2 , y2, y1, y1) , lwd=2 , col=COL , border="black", lwd=0.2)   
        
            gap = GAPS[[ ass]]
            if ( nrow(gap)> 0) gap$"y3" = y3 else gap$"y3" = gap$"chr"
            gap = gap[ gap$"chr" == i, ]
            GAP2[[i]][[ ass]] = gap
            
        if ( nrow(gap)> 0)  points(gap$"ref.end",  gap$"y3", type="p", pch=18, cex=CEX.GAP ,col="black")
        }
        
        if ( nrow(bb) > 0) text( max(c(CL,A))+sh  , y3 ,  paste0(nrow(gap)," gaps"), cex=CEX.txt ,col=COL,adj=0)
        
       # segments(rep(0, nrow(rrr)),rrr$"chr.pos"+p, rrr$"ref.len", rrr$"chr.pos"+p, col="gray", cex=2)
        if ( nrow(bb) == 0)  CL = rrr[ rrr$"chr" ==i, ]$"ref.len" 
        if ( nrow(bb) == 0)  Y1 = rrr[ rrr$"chr" ==i, ]$"chr.pos" - 0.5 -0.01
        if ( nrow(bb) == 0)  Y2 = rrr[ rrr$"chr" ==i, ]$"chr.pos" + 0.5 +0.01


        }          
        polygon( c(0 , CL,CL,0) ,c( Y2 , Y2, Y1, Y1) , lwd=1 , border="blue", lwd=2, lty=1 , col=adjustcolor("gray",alpha.f=0.5) )   

       
}
cat("\n")
cat("\n  --> Adding GAPs","\n")
 for (i in names(A.list)) { AA = CHR2[,i] ;     if( nrow(CHR2)==1 ) { names(AA) = rownames(CHR2)} ; AA = AA[ rev(names(AA))] ; cat("\n",i,"") ; for (k in length(AA):1) { ass=names(AA)[k] ; cat(ass,"") ; gap = GAP2[[i]][[ ass]] ; if ( nrow(gap)> 0) points(gap$"ref.end",  gap$"y3", type="p", pch=23, cex=1 ,col="black", bg="yellow", lwd=2 )}}

GAP.tot  = sapply(GAPS,nrow)
GAP.tot  = data.frame("assembly"= names(GAP.tot), "gaps"=as.numeric(GAP.tot), stringsAsFactors=TRUE)

colors2 = colors
colors2 = merge(colors2, GAP.tot, by="assembly", all.x=TRUE)
colors2$"n.gap" = paste0(colors2$"assembly","   ",colors2$"gap"," gaps")

colors3 = colors2[, c("assembly","col","gaps","n.gap")]

colors3$"bg" = colors3$"col"
#colors3$"col" = "black"
colors3$"pch" = 22

rownames(colors3) = colors3$"assembly"
#cat("\n"); print(data.frame(colors3))

colors3 = colors3[   (SAMPLE$"short_name") , ]
#cat("\n");print(data.frame(colors3))

#colors3 = colors3[  rev(SAMPLE$"short_name") , ]

colors3 = rbind(colors3, c(" ","white"," "," ","white",22))
colors3 = rbind(colors3, c("Reference","blue","na","Reference","gray",22))
colors3 = rbind(colors3, c("Gaps","black","na","Gaps","yellow",23))

rownames(colors3) = colors3$"assembly"
#colors3 = colors3[  names(A.list[[1]]) , ]

#cat("\n");print(data.frame(colors3))


colors3$"tag" = paste0(colors3$"assembly","        ",colors3$"gaps"," gaps")
colors3$"ta2" = paste0(colors3$"assembly","                 ")
colors3$"ta3" = paste0(colors3$"gaps"," gaps")
colors3[ !colors3$"assembly" %in% GAP.tot$"assembly", ]$"ta3"=""

w_1gap = which(colors3$"ta3" ==1)
if(length(w_1gap)>0) colors3$"ta3"[w_1gap] = gsub("gaps","gap",colors3$"ta3"[w_1gap])

legend("right", colors3$"ta2", cex=CEX.leg*0.8, pt.bg=colors3$"bg", col=colors3$"col" , pch=as.numeric(colors3$"pch") , bg=adjustcolor("white",alpha.f=0.5), pt.lwd =3 )
legend("right", colors3$"ta3", cex=CEX.leg*0.8, pt.lwd =3 , bty="n"  )


# 
# colors$"tag" = paste0(colors$"assembly","        ",colors$"perc","%")
# colors$"ta2" = paste0(colors$"assembly","                 ")
# colors$"ta3" = paste0(colors$"perc","%")
# 
# rownames(colors) = colors$"assembly"
# colors = colors[  names(A.list[[1]]) , ]
# CEX.leg=2.5
# legend("bottomright", gsub(".ch0","",colors$"ta2"), cex=CEX.leg, fill=colors$col  )
# #legend("bottomright", gsub(".ch0","",colors$"tag"), fill=colors$col  )
# legend("bottomright", gsub(".ch0","",colors$"ta3"),cex=CEX.leg, bty="n" )
# 
# 



 dev.off()       

cat("############################################################################\n")
cat("### 4)    Reference Coverage + Scaffold Length/GAPs in one Figure        ###\n")
cat("############################################################################\n")


#   CTG.min = 000
#pdf( paste0(out.dir,"/GS-viewer_04.Scaffold_Coverage_on_Reference_plus_Scaffold_Length_and_Gaps.",OUT,"min.",CTG.min,"bp.pdf"),width=plot.w, height=plot.h )
pdf( paste0(out.dir,"/GS-viewer_04.Scaffold_Coverage_on_Reference_plus_Scaffold_Length_and_Gaps.",OUT,"pdf"),width=plot.w, height=plot.h )
par( mfcol=c(1,2) ,  mar=c(2,11,4,1))
# chr = chr[ chr$"chr" != "Chr0", ]

CEX.names = 2

if ( nchar(max(ref$"chr"))>7 ) CEX.names = 1.5

cat("\n  --> Scaffold length","\n")

# CHR2 = ref
# CHR2$"chr" = paste0(CHR2$"chr","_RagTag")
# for (i in names(CHR)) { x=CHR[[i]] ; names(x)=c("chr",i) ; CHR2 = merge(CHR2,x,by="chr", all.x=TRUE) }
# CHR2 = as.data.frame(CHR2)
# CHR2 [ is.na(CHR2) ] = 0
# CHR2$"chr" = gsub("_RagTag","",CHR2$"chr")
# rownames(CHR2) = CHR2$"chr" 
# CHR2$"ref.len" = NULL
# CHR2$"chr" = NULL
# CHR2 = t(as.matrix(CHR2))
# 
# CHR2.colnames= colnames(CHR2)
# CHR2.rownames= rownames(CHR2)
# 
# if(nrow(CHR2)>1) CHR2 = CHR2[  rev(rownames(CHR2)) , ]
# if(ncol(CHR2)>1) CHR2 = CHR2[ , rev(colnames(CHR2)) ]
# 
# if( is.null(dim(CHR2)) )  { CHR2= as.matrix(t(CHR2)) ;  dimnames(CHR2)[[1]] = CHR2.rownames } 
# 

# 
# print(CHR2)
# print(max(CHR2))
# print(max(c(ref$"ref.len",CHR2)))

#   CTG.min = 000
# chr = chr[ chr$"chr" != "Chr0", ]
bar=barplot(ref$"ref.len", names.arg=ref$"chr", las=1, horiz=T, main="Scaffold length vs Reference length",col=gray(0.95), cex.main=3, xlim=c(0, max(c(ref$"ref.len",CHR2,CCmax)*1.20 )),  cex.names=CEX.names, cex.axis=CEX.names)
#bar=barplot(ref$"ref.len", names.arg=ref$"chr", las=1, horiz=T, main="Scaffold length vs Reference length",col=gray(0.95), cex.main=3, xlim=c(0, max(c(ref$"ref.len",CHR2,CCmax)*1.20 )) ,  cex.axis=2, cex.names=2)

b2 = data.frame(chr.pos = bar[,1])
b2$"chr" = ref$"chr"
rrr = merge(ref,b2, by="chr")

sh = max(ref$"ref.len")*0.01

n.ass = length( A.list[[1]] )

shifts = rev(seq(-0.5,0.5,length.out=n.ass+1))

#for( p in shifts ) segments(rep(0, nrow(rrr)),rrr$"chr.pos"+p, rrr$"ref.len", rrr$"chr.pos"+p, col="gray", cex=2)


GAP2 = list()
for (i in names(A.list)) 1==1
for (i in names(A.list)) { 
#  i=  names(A.list)[4]
    AA = CHR2[,i]
    if( nrow(CHR2)==1 ) { names(AA) = rownames(CHR2)}
    AA = AA[ rev(names(AA))]
    bb = rrr[ rrr$"chr" == i, ]
    bb$"chr.pos"
    
       
    CEX.txt = 1.1
    if(length(CHR)>7) CEX.txt = 0.6

    x = 0.02
    cat("\n",i,"")
     
    for (k in length(AA):1) 
    {  
        A=AA[[k]] ;
        ass=names(AA)[k]
        cat(ass,"")
        COL =  unique(col.data [ col.data$"assembly" == ass ,]$"col" )
        s1=shifts[k];
        s2=shifts[k+1];
        y1 = bb$"chr.pos" + s1
        y2 = bb$"chr.pos" + s2
        y3 = bb$"chr.pos" + (s1+s2)/2
        Y1 = bb$"chr.pos" - 0.5 -0.01
        Y2 = bb$"chr.pos" + 0.5 +0.01
        CL = bb$"ref.len" 
    #    cat(y3,"")
     if (nrow(bb)>0) 
        {   polygon( c(0 , A, A,0) ,c( y2 , y2, y1, y1) , lwd=2 , col=COL , border="black", lwd=0.2)   
        
            gap = GAPS[[ ass]]
            if ( nrow(gap)> 0) gap$"y3" = y3 else gap$"y3" = gap$"chr"
            gap = gap[ gap$"chr" == i, ]
            GAP2[[i]][[ ass]] = gap
            
        if ( nrow(gap)> 0)  points(gap$"ref.end",  gap$"y3", type="p", pch=18, cex=CEX.GAP ,col="black")
        }
        
        if ( nrow(bb) > 0) text( max(c(CL,A))+sh  , y3 ,  paste0(nrow(gap)," gaps"), cex=CEX.txt ,col=COL,adj=0)
        
       # segments(rep(0, nrow(rrr)),rrr$"chr.pos"+p, rrr$"ref.len", rrr$"chr.pos"+p, col="gray", cex=2)
        if ( nrow(bb) == 0)  CL = rrr[ rrr$"chr" ==i, ]$"ref.len" 
        if ( nrow(bb) == 0)  Y1 = rrr[ rrr$"chr" ==i, ]$"chr.pos" - 0.5 -0.01
        if ( nrow(bb) == 0)  Y2 = rrr[ rrr$"chr" ==i, ]$"chr.pos" + 0.5 +0.01


        }          
        polygon( c(0 , CL,CL,0) ,c( Y2 , Y2, Y1, Y1) , lwd=1 , border="blue", lwd=2, lty=1 , col=adjustcolor("gray",alpha.f=0.5) )          
}
cat("\n")
cat("\n  --> Adding GAPs","\n")
for (i in names(A.list)) { AA = CHR2[,i] ;  if( nrow(CHR2)==1 ) { names(AA) = rownames(CHR2)} ; AA = AA[ rev(names(AA))] ; cat("\n",i,"") ; for (k in length(AA):1) { ass=names(AA)[k] ; cat(ass,"") ; gap = GAP2[[i]][[ ass]] ; if ( nrow(gap)> 0) points(gap$"ref.end",  gap$"y3", type="p", pch=23, cex=1 ,col="black", bg="yellow", lwd=2 )}}

GAP.tot  = sapply(GAPS,nrow)
GAP.tot  = data.frame("assembly"= names(GAP.tot), "gaps"=as.numeric(GAP.tot), stringsAsFactors=TRUE)

colors2 = colors
colors2 = merge(colors2, GAP.tot, by="assembly", all.x=TRUE)
colors2$"n.gap" = paste0(colors2$"assembly","   ",colors2$"gap"," gaps")

colors3 = colors2[, c("assembly","col","gaps","n.gap")]

colors3$"bg" = colors3$"col"
#colors3$"col" = "black"
colors3$"pch" = 22

rownames(colors3) = colors3$"assembly"
#cat("\n"); print(data.frame(colors3))

colors3 = colors3[   SAMPLE$"short_name"  , ]
#cat("\n");print(data.frame(colors3))

#colors3 = colors3[  rev(SAMPLE$"short_name") , ]

colors3 = rbind(colors3, c(" ","white"," "," ","white",22))
colors3 = rbind(colors3, c("Reference","blue","na","Reference","gray",22))
colors3 = rbind(colors3, c("Gaps","black","na","Gaps","yellow",23))

rownames(colors3) = colors3$"assembly"
#colors3 = colors3[  names(A.list[[1]]) , ]

#cat("\n");print(data.frame(colors3))


colors3$"tag" = paste0(colors3$"assembly","        ",colors3$"gaps"," gaps")
colors3$"ta2" = paste0(colors3$"assembly","                 ")
colors3$"ta3" = paste0(colors3$"gaps"," gaps")
colors3[ !colors3$"assembly" %in% GAP.tot$"assembly", ]$"ta3"=""

w_1gap = which(colors3$"ta3" ==1)
if(length(w_1gap)>0) colors3$"ta3"[w_1gap] = gsub("gaps","gap",colors3$"ta3"[w_1gap])


cat("\n");print(data.frame(colors3))

legend("right", colors3$"ta2", cex=CEX.leg*0.8, pt.bg=colors3$"bg", col=colors3$"col" , pch=as.numeric(colors3$"pch") , bg=adjustcolor("white",alpha.f=0.5), pt.lwd =3 )
legend("right", colors3$"ta3", cex=CEX.leg*0.8, pt.lwd =3 , bty="n"  )


cat("\n  -->  Reference Coverage","\n")


bar=barplot(ref$"ref.len", names.arg=ref$"chr", las=1, horiz=T, main="Scaffold coverage on Reference Genome",col=gray(0.95), cex.main=3, xlim=c(0,max(ref$"ref.len")*1.20), cex.names=CEX.names, cex.axis=CEX.names)
b2 = data.frame(chr.pos = bar[,1])
b2$"chr" = ref$"chr"
rrr = merge(ref,b2, by="chr")

text( max(ref$"ref.len")/2, max(b2$"chr.pos")*1.06,  paste( "min plotted window =",CTG.min,"bp"), cex=0.5 )

sh = max(ref$"ref.len")*0.01

n.ass = length( A.list[[1]] )

shifts = rev(seq(-0.5,0.5,length.out=n.ass+1))

for( p in shifts ) segments(rep(0, nrow(rrr)),rrr$"chr.pos"+p, rrr$"ref.len", rrr$"chr.pos"+p, col="gray", cex=2)

for (i in names(A.list)) 1==1

A.coo.list=list()

for (i in names(A.list)) { 
#  i=  names(A.list)[1]
    A.list [[i]]  -> AA
    B.list [[i]]  -> BB
    g.list [[i]]  -> gg
        
    bb = b2[ b2$"chr" == i, ]
    bb$"chr.pos"
    
    CEX.txt = 1
    if(length(CHR)>7) CEX.txt = 0.6

    x = 0.02
    cat("\n",i,"")
     
    for (k in length(AA):1) 
    {  
        A=AA[[k]] ;
        ass=names(AA)[k]
        cat(ass,"")
        COL =  unique(col.data [ col.data$"assembly" == ass ,]$"col" )
        s1=shifts[k];
        s2=shifts[k+1];
        y1 = bb$"chr.pos" + s1
        y2 = bb$"chr.pos" + s2
        y3 = (y1+y2)/2
       # cat(k,"")
        if( nrow(A)>0)
        {
         
        A.coord = list()
        for (j in 1:nrow(A)) A.coord[[j]]  = c(A$"ref.start"[j]: A$"ref.end"[j])           
        A.coord = unique(do.call(c,A.coord))
        A.perc =round(100* length(A.coord)/A$"ref.len"[1])
        A.coo.list [[ass]] [[i]] = A.coord
                
        A = A[ A$"ctg.len" > CTG.min , ]

        for (j in 1:nrow(A)) { polygon( c(A$"ref.start"[j], A$"ref.end"[j], A$"ref.end"[j],A$"ref.start"[j]) ,c( y2 , y2, y1, y1) , lwd=2 , col=COL , border=NA)    }
        text( A$"ref.len"[1]+sh*2  , y3 ,  paste0(A.perc,"%"), cex=CEX.txt ,col=COL)
        }   else { A.coo.list [[ass]] [[i]] = data.frame()}       
        }
}


#COL.LIST=c("firebrick3","black","dodgerblue4","darkorange3","darkturquoise","blueviolet","forestgreen","lightcoral")
col.list = COL.LIST[ 1:length(CHR)]

colors= data.frame(assembly=names(A.list [[1]] ), col=col.list[1:length(A.list [[1]] )], row.names=names(A.list [[1]] ), stringsAsFactors=FALSE)

colors = col.data

chr.sum = sapply(A.coo.list,sapply,length)
chr.sum = as.data.frame(chr.sum)
gen.sum = sapply(chr.sum,sum)
gen.sum = data.frame( assembly=names(gen.sum), len=gen.sum, row.names=names(gen.sum), stringsAsFactors=FALSE)
colors  = merge( colors, gen.sum, by="assembly")
colors$"perc" = round( 100* colors$"len" / sum(ref$"ref.len") ,2)
colors$"tag" = paste0(colors$"assembly","        ",colors$"perc","%")
colors$"ta2" = paste0(colors$"assembly","                 ")
colors$"ta3" = paste0(colors$"perc","%")

colors$"bg" = colors$"col"
colors$"pch" = 22

rownames(colors) = colors$"assembly"
colors = colors[  names(A.list[[1]]) , ]
CEX.leg=2.5
legend("right", gsub(".ch0","",colors$"ta2"), cex=CEX.leg*0.8, pt.bg=colors$"bg", col=colors$"col",  pch=as.numeric(colors$"pch")  )
#legend("bottomright", gsub(".ch0","",colors$"tag"), fill=colors$col  )
legend("right", gsub(".ch0","",colors$"ta3"),cex=CEX.leg*0.8, bty="n" )


cat("\n")

 
 
 dev.off()       



cat("\n")
cat("############################################################################\n")
cat("### 5) Single Chromosome reports                                       ###\n")
cat("############################################################################\n")
cat("\n")

# MIN.WIND = 

for( ii in names(CHR))
{

# ii =  "hifi200k";

# ii =  names(CHR)[2];
# ii =  names(CHR)[1];

cat(ii,"")

chr = data.frame(CHR[[ii]]) ; 
agp = setDT(AGP[[ii]]) ; 
cnf = setDT(CNF[[ii]]) ; 
asm = setDT(ASM[[ii]]) ; 
chr = chr[ rev(order(chr$"chr")) ,]
chr$"chr" = gsub("_RagTag","",chr$"chr")

CCmax.max = max(CCmax, chr$"len")

# chr = chr[ chr$"chr" != "Chr0", ]
bar=barplot(chr$"len", names.arg=chr$"chr", las=1, horiz=T, main=i, xlim=c(0, CCmax.max*1.10),col=gray(0.95) )
b2 = data.frame(chr.pos = bar[,1])
b2$"chr" = chr$"chr"
rrr = merge(ref, b2,by="chr", all.x=TRUE)

# separate GAPs annotation from AGP
agp$"chr" = gsub("_RagTag","",agp$"ref.chr")
gap = agp[ agp$"query" =="100",]
agp = agp[ agp$"query" !="100",]

# ASM2 select in ASM the same contigs as AGP 
# ASM3 select in ASM the same contigs as AGP 


agp2 = split( agp, agp$"chr")
asm2 = split( asm, asm$"ref.chr")
asm3 = split( asm, asm$"ref.chr")

sapply(asm2,nrow)
for (i in names(asm2)) { asm2[[i]] = asm2[[i]] [  asm2[[i]]$"query" %in% agp2[[i]]$"query",    ]  }
for (i in names(asm3)) { asm3[[i]] = asm3[[i]] [ !asm3[[i]]$"query" %in% agp2[[i]]$"query",    ]  }
sapply(asm2,nrow)
sapply(asm3,nrow)


ASM3 = do.call(rbind,asm3)


INVERT_MINUS_CTG = TRUE
SINGLE_CHR_PLOT = FALSE

 if (INVERT_MINUS_CTG ==FALSE) pdf( paste0(out.dir,"/GS-viewer_05.Chromosome_details.",ii,".",OUT,"ALL_CHR.pdf"), width=plot.w, height=plot.h )
#if (INVERT_MINUS_CTG ==TRUE ) pdf( paste0(out.dir,"/GS-viewer_05.Chromosome_details.",ii,".",OUT,"ALL_CHR.oriented.pdf"),width=plot.w, height=plot.h )
 if (INVERT_MINUS_CTG ==TRUE ) pdf( paste0(out.dir,"/GS-viewer_05.Chromosome_details.",ii,".",OUT,"ALL_CHR.pdf"), width=plot.w, height=plot.h )

FINAL.AGP = list()
B.DATA   = list()
for (i in names(asm2))  { 
#   i = names(asm2)[1]
    cat(i,"")
    a = asm2[[i]]  
    g = agp2[[i]]
    a$"query" = as.character(a$"query")
    g$"query" = as.character(g$"query")

    if (nrow(a)> 0 & length(g)> 0)
    {
    
        # contigs on ref
        x = 0.02
        
        a$"ya" = x*2
        a$"yb" = a$"ya"
        a$"col"="blue4"
        if( any(a$"strand"=="-")) a[ a$"strand"=="-",]$"yb"= a[ a$"strand"=="-",]$"yb" *-1
        if( any(a$"strand"=="-")) a[ a$"strand"=="-",]$"col"="darkred"
        if( any(a$"strand"=="+")) a[ a$"strand"=="+",]$"col"="blue4"
        
        # sort by position of biggest fragment
        a$"ref.width" = a$"ref.end"- a$"ref.start"+1
        a$"q.left"    = a$"q.len"- a$"q.end"
        a$"strand.query" = paste(a$"query", a$"strand")
        a$"chr.ctg" = paste0(a$"ref.chr",":", a$"query")
    
        as= split(a,a$"query")
        ass= lapply(as, function(x) x[ which.max(x$"ctg.len"), ])
        ass= sapply(ass, function(x) x$"ref.start")   
        as = as [ names(sort(ass)) ]
        as = as [ g$"query"]
        # add cordinates for plot
        AS = list()
        for (j in 1:length(as)) { aa = as[[j]] ; aa=split(aa,aa$"strand.query") ; for (kk in names(aa)) AS[[kk]]=aa[[kk]]  }
        for (j in 1:length(AS)) { AS[[j]]$"ya"= AS[[j]]$"ya"+j/100  }
        a= do.call(rbind,AS)
        
        length(unique(a$query))
        length(unique(a$ya))
        length(unique(a$strand.query))
    
        n.a = length(as)
        col.a = rainbow(n.a)
        chr.max =a$"ref.len"[1]
        max.len  =max( chr.max, g$"ref.end" )
        MAX.LEN = max.len
        if( max(g$"ref.end") / chr.max < 1.20 ) { max.len = chr.max*1.25}
        chr.txt = -max.len/8
        n.ctg = length(unique(a$query))
        for (j in 1:length(as)) { as[[j]]$"ya"= as[[j]]$"ya"+j/100  }
        
        a= do.call(rbind,as)
    
        
        # scaffold
        g$"za"= -x/2
        g$"zb"= g$"za" -x
        g$"zm"=g$"zb"-g$"za"
     #  g$"zs"=g$"zb"-0.02
        g$"zs"=g$"zb" -x/2
        g$"ref.med"=(g$"ref.start"+g$"ref.end")/2
        g$"ref.ctg_start"=g$"ref.med"-g$"q.len"/2
        g$"ref.ctg_end"  =g$"ref.med"+g$"q.len"/2
        g$"ZS"=g$"zs" -(1:nrow(g))/100
        
        g$"col"="blue4"
        g$"QUERY"=g$"query"
        if( any(g$"strand"=="-")) g[ g$"strand"=="-",]$"col"="darkred"
        if( any(g$"strand"=="+")) g[ g$"strand"=="+",]$"col"="blue4"
        if( any(g$"strand"=="-")) g[ g$"strand"=="-",]$"QUERY"= paste0(g[ g$"strand"=="-",]$"QUERY", ".rev")

        # invs 
         as2= lapply( split(a$"strand",a$"query"), unique )
         qg= split(a,a$"strand.query")
       
        gg = g[, c(2,3,6)] 
        names(gg) = c("scaf.start","scaf.end","query")
        a = merge(a, gg, by="query",all.x=T)
        a$"scaf.pos"  =a$"scaf.start"+a$"q.start"-1
        a$"scaf.end"  =a$"scaf.start"+a$"q.end"-1
        
        gg3.col= c("query","strand","col")
        gg3 = g[,..gg3.col] 
        names(gg3) = c("query","agp.str","apg.col")
        a = merge(a,gg3, by="query", all.x=T)
        a = a[ order(a$"ref.chr",a$"query",a$"ref.start",a$"ref.end") , ]
        as= split(a,a$"query")
    
        # contig names
        col.aa = c("query","q.len","strand.query","ya","yb","col","agp.str","apg.col")
        aa = a[ ,..col.aa]
        aa = aa[ ! duplicated(aa),]
        all( g$"q.len" %in% aa$"q.len" )
            
        #CTG real start-ends
        a$"ref.st.ext" = NA
        a$"ref.e.ext"   = NA
        w.minu = which(a$"agp.str"=="-")
        w.plus = which(a$"agp.str"=="+")
        if( length(w.plus)>0 ) { a$"ref.st.ext"[w.plus]=a$"ref.start"[w.plus]-a$"q.start"[w.plus] ; a$"ref.e.ext"[w.plus]=a$"ref.end"[w.plus]+a$"q.left" [w.plus] }
        if( length(w.minu)>0 ) { a$"ref.st.ext"[w.minu]=a$"ref.start"[w.minu]-a$"q.left" [w.minu] ; a$"ref.e.ext"[w.minu]=a$"ref.end"[w.minu]+a$"q.start"[w.minu] }
        as= split(a,a$"query")
        
               
        # PLOT   
        b = do.call(rbind, as) ;  b=b[ b$"ref.width">MIN.WIND*1.0 ,] 
        Yax1 = max( b$"ya") + 2*x
        Yax2 = min( g$"ZS") - 2*x*1.5
        if (Yax1 <= 0.1) Yax1 = 0.1
        if (Yax2 >= -0.1) Yax2 = -0.1
        Yax2 = -Yax1
        ##   SINGLE_CHR_PLOT = TRUE
        if( SINGLE_CHR_PLOT == TRUE )  dir.create(paste0("../../../ragtag_scaffold.",ii, "/PDF/"))
        if( SINGLE_CHR_PLOT == TRUE )  PDF =  paste0("../../../ragtag_scaffold.",ii, "/PDF/Scaffold_vs_reference_",i,".pdf")
        if( SINGLE_CHR_PLOT == TRUE )  pdf( PDF, 30,16)
            # --> REF TRACK
        plot(c(0, chr.max ), c(x,x), type="l", lwd=5, xlim=c(-max.len/10,max.len), ylim=c(Yax2, Yax1), lend=1, cex.axis=2)
        abline(h=0, col="gray")
        
        g$"Q.LEN"  = as.character(g$"q.len")
        if( any(g$"q.len">=1000   ) )  g[ g$"q.len">=1000    ,]$"Q.LEN"= paste( as.character(round( g[ g$"q.len">=1000    ,]$"q.len"/1000   )),"Kb")
        if( any(g$"q.len">=1000000) )  g[ g$"q.len">=1000000 ,]$"Q.LEN"= paste( as.character(round( g[ g$"q.len">=1000000 ,]$"q.len"/1000000)),"Mb")
      #  print(head(data.frame(g)))
        
        # --> CONTIG ALIGN on REF
       #  MIN.WIND = 5000
        # for (j in 1:length(as)){ b= as[[j]];  b=b[ b$"ref.width">0000     ,] ; segments( b$"ref.start", b$"ya"-0.000, b$"ref.end", b$"ya"-0.000, lwd=1, col=b$"col") }
    #      for (j in 1:length(as)){ b= as[[j]];  b=b[ b$"ref.width">MIN.WIND*2.0 ,] ; segments( b$"ref.st.ext", b$"ya"-0.000, b$"ref.e.ext", b$"ya"-0.000, lwd=0.5, col="gray" ) ; points( b$"ref.st.ext", b$"ya", pch="|", cex=0.3,col="gray" ); points( b$"ref.e.ext", b$"ya", pch="|", cex=0.3,col="gray" ) }
    #      for (j in 1:length(bs)){ b= bs[[j]];  b=b[ b$"ref.width">BED.WIND*1.5 ,] ; segments( b$"ref.start" , b$"ya"-0.000, b$"ref.end"  , b$"ya"-0.000, lwd=5.0, col="gold") }
    #      for (j in 1:length(as)){ b= as[[j]];  b=b[ b$"ref.width">MIN.WIND*1.0 ,] ; segments( b$"ref.start" , b$"ya"-0.000, b$"ref.end"  , b$"ya"-0.000, lwd=2.0 , col=b$"col") }
         rs = as
         rs = lapply(rs, function(x) head(x[ rev( order(x$"ref.width")),]) )
     #   b = do.call(rbind, rs) ;  b=b[ b$"ref.width">MIN.WIND*2.0 ,] ; segments( b$"ref.st.ext", b$"ya"-0.000, b$"ref.e.ext", b$"ya"-0.000, lwd=0.5, col="gray" ) ; points( b$"ref.st.ext", b$"ya", pch="|", cex=0.3,col="gray" ); points( b$"ref.e.ext", b$"ya", pch="|", cex=0.3,col="gray" ) 
     #   b = do.call(rbind, bs) ;  b=b[ b$"ref.width">BED.WIND*1.5 ,] ; segments( b$"ref.start" , b$"ya"-0.000, b$"ref.end"  , b$"ya"-0.000, lwd=5.0, col="gold") 
         b = do.call(rbind, as) ;  b=b[ b$"ref.width">MIN.WIND*1.0 ,] ; segments( b$"ref.start" , b$"ya"-0.000, b$"ref.end"  , b$"ya"-0.000, lwd=2.0 , col=b$"col") 
         
         bbb.col=c("ref.chr","ref.start","ref.end","strand","query","q.len","ref.width")
         bbb = b[ , ..bbb.col]
         
        # check overlaps
        if( dir.exists(paste(out.dir,"/ctg_coordinates_on_ref",sep="")) == FALSE ) {  dir.create(paste(out.dir,"/ctg_coordinates_on_ref",sep=""))}
        bname = paste(out.dir,"/ctg_coordinates_on_ref/",ii,".",i,".ctg_coordinates_on_ref.txt",sep="")
        b.int = paste(out.dir,"/ctg_coordinates_on_ref/",ii,".",i,".ctg_coordinates_on_ref.intersect.txt",sep="")
        write.table( bbb , bname  ,col.names=FALSE, row.names=F, quote=F, sep="\t")
        system( paste0( " bedtools intersect -wao -a ",bname," -b ",bname, " > ",b.int) )
        system( paste0( "wc -l ",b.int ), inter=TRUE) -> bi.length
        bi.length = as.numeric(as.character(unlist(strsplit(bi.length," ")))) [1]

        bbb = split(bbb,bbb$"query")
#         cat("print bbb \n")
#         print(length(bbb))
        if( length(bbb)>0) { bbb.stat = data.frame("query"= names(bbb),"n.all.fragment"      =as.numeric(sapply(bbb,nrow))                                               ,"all.align.bp"= as.numeric(sapply(bbb,function(x) sum(x$"ref.width"))), stringsAsFactors=FALSE) }
        bbb.empty = data.frame("query"=0,"n.overlap.fragment"  =0, "overlapp.align.bp"= 0, stringsAsFactors=FALSE)
        bbb.empty = bbb.empty[ bbb.empty$"query">0, ]
        if( bi.length > 0 & is.na(bi.length)==FALSE)
        {        
            bi = fread(b.int, stringsAsFactors=FALSE)
            names(bi) = c(   paste0("A", bbb.col) , paste0("B", bbb.col), "overlap" )
            bi = bi[ bi$"Aquery" != bi$"Bquery", ]
         #   bi = bi[ bi$"Aref.start" <= bi$"Bref.start" & bi$"Aref.end" >= bi$"Bref.end", ]
            bi$"fragment" = paste(bi$"Bquery" ,bi$"Bref.start",bi$"Bref.end",sep=":")
            
            bis = split(bi  ,bi$"Bquery")
            #bbb = split(bbb,bbb$"query")
            
#              print(length(bis))
            if( length(bis)>0) { bis.stat = data.frame("query"= names(bis),"n.overlap.fragment"  =as.numeric(sapply(bis,function(x) length(unique(x$"fragment")))), "overlapp.align.bp"= as.numeric(sapply(bis,function(x) sum(x[ !duplicated(x$"fragment"),]$"overlap") )+1), stringsAsFactors=FALSE) }
          # if( length(bbb)>0) { bbb.stat = data.frame("query"= names(bbb),"n.all.fragment"      =as.numeric(sapply(bbb,nrow))                                               ,"all.align.bp"= as.numeric(sapply(bbb,function(x) sum(x$"ref.width"))), stringsAsFactors=FALSE) }
            
            # if no overlap available, create empty df
            if( length(bis)==0) { bis.stat = data.frame("query"= names(bbb),"n.overlap.fragment"  =0, "overlapp.align.bp"= 0, stringsAsFactors=FALSE) }
            
            }   else  {  
                if(length(bbb)>0 ) { bbb.stat = data.frame("query"=names(bbb) ,"n.overlap.fragment"  =0, "overlapp.align.bp"= 0, stringsAsFactors=FALSE) }  ;
                if(length(bbb)==0) { bbb.stat = bbb.empty }  ;
                }  
         
#         print(bbb.stat)
#         print(bis.stat)
#          
        if( nrow(bbb.stat) > 0 & nrow(bis.stat)> 0)
        {
            fin.stat = merge(bbb.stat,bis.stat,all.x=T, by="query")
            fin.stat$"perc.bp" = round(100 *  fin.stat$"overlapp.align.bp" /  fin.stat$"all.align.bp", 2) 
            fin.stat$"chr" = i 
            fin.stat$"sample" = ii
            }  else { fin.stat = data.frame("query"= 0, "n.all.fragment" = 0, "all.align.bp" = 0 , "n.overlap.fragment" = 0 , "overlapp.align.bp" = 0 , "perc.bp" = 0 , "chr" = 0 , "sample" = 0, stringsAsFactors=FALSE ) ; fin.stat = fin.stat[ fin.stat$"query">0, ]
          }

#         print(fin.stat)

 #      text( -max.len/60 ,  aa$"ya", adj=1, paste(aa$"query","vs Ref."), cex=0.5, col= aa$"apg.col")
        text( -max.len/60 ,  aa$"ya", adj=1, paste(aa$"query"          ), cex=1.0, col= aa$"apg.col")
        text( -max.len/60 ,  x      , adj=1, paste("Ref.",a$"ref.chr"),cex=1.5)
        text( -max.len/60 , -x      , adj=1, paste("Scaffold",a$"ref.chr"),cex=1.5)
        
        # CTG stat
        xx = max.len/40
        b.data = data.frame( "query"=names(as),stringsAsFactors=F)
        b.data$"len"  =   as.numeric(sapply(as, function(x) x$"q.len"[1]))   
        b.data$"LEN"  =   b.data$"len"
        ww = which(b.data$"len" >1e3 & b.data$"len" < 1e6) ; if(length(ww)>0) { b.data$"LEN"[ww] = paste( round(b.data$"len"[ww]/1e3) , "Kb") }
        ww = which(b.data$"len" >1e3 & b.data$"len" >=1e6) ; if(length(ww)>0) { b.data$"LEN"[ww] = paste( round(b.data$"len"[ww]/1e6) , "Mb") }
        b.data$"ya"  =   as.numeric(sapply(as, function(x) x$"ya"[1])) - 0.000  
        b.data$"xa1"  =  chr.max + 0.20*xx
        b.data$"xa2"  =  chr.max + 2.70*xx
        b.data$"xa3"  =  chr.max + 3.70*xx
        b.data$"xa4"  =  chr.max + 4.80*xx
        b.data$"xa5"  =  chr.max + 5.80*xx
        b.data$"xa6"  =  chr.max + 6.80*xx
        b.data$"xa7"  =  chr.max + 7.70*xx
        b.data$"xa8"  =  chr.max + 8.60*xx
        b.data$"xa9"  =  chr.max + 9.50*xx
        b.data$"col" = as.character(sapply(as, function(x) x$"col"[1]))
        b.data$"min.win" = as.numeric(sapply(as, function(x) min(x$"ref.width")))              
        b.data$"max.win" = as.numeric(sapply(as, function(x) max(x$"ref.width")))              
        b.data$"tot.cov" = as.numeric(sapply(as, function(x) sum(x$"ref.width")))              
        b.data$"perc.cov"= round(100*b.data$"tot.cov"/b.data$"len",0)              
        b.data$"perc.min"= round(100*b.data$"min.win"/b.data$"len",1)              
        b.data$"perc.max"= round(100*b.data$"max.win"/b.data$"len",1)              
        b.data$"perc.cov2"= paste0(b.data$"perc.cov","%")              
        b.data = merge(b.data, cnf, by="query",all.x=TRUE)
        b.data$"string"= paste0(b.data$"query"," %")  
        b.data$"MAX.WIN"  =   b.data$"max.win"
        ww = which(b.data$"max.win" >1e3 & b.data$"max.win" < 1e6) ; if(length(ww)>0) { b.data$"MAX.WIN"[ww] = paste( round(b.data$"max.win"[ww]/1e3) , "Kb") }
        ww = which(b.data$"max.win" >1e3 & b.data$"max.win" >=1e6) ; if(length(ww)>0) { b.data$"MAX.WIN"[ww] = paste( round(b.data$"max.win"[ww]/1e6) , "Mb") }
        b.data$"TOT.COV"  =   b.data$"tot.cov"
        ww = which(b.data$"tot.cov" >1e3 & b.data$"tot.cov" < 1e6) ; if(length(ww)>0) { b.data$"TOT.COV"[ww] = paste( round(b.data$"tot.cov"[ww]/1e3) , "Kb") }
        ww = which(b.data$"tot.cov" >1e3 & b.data$"tot.cov" >=1e6) ; if(length(ww)>0) { b.data$"TOT.COV"[ww] = paste( round(b.data$"tot.cov"[ww]/1e6) , "Mb") }
        b.data = merge(b.data, fin.stat, by="query",all.x=TRUE)
        B.DATA[[i]] = b.data
            
        text(  b.data$"xa1" ,  b.data$"ya", adj=0, b.data$"query"                              , cex=1.0, col= b.data$"col")
        text(  b.data$"xa2" ,  b.data$"ya", adj=1, b.data$"LEN"                                , cex=1.0, col= b.data$"col")
        text(  b.data$"xa3" ,  b.data$"ya", adj=1, b.data$"TOT.COV"                            , cex=1.0, col= b.data$"col")
        text(  b.data$"xa4" ,  b.data$"ya", adj=1, b.data$"MAX.WIN"                            , cex=1.0, col= b.data$"col")
        text(  b.data$"xa5" ,  b.data$"ya", adj=1, b.data$"perc.cov2"                          , cex=1.0, col= b.data$"col")
        text(  b.data$"xa6" ,  b.data$"ya", adj=1, substr(b.data$"grouping_confidence",1,4)    , cex=1.0, col= b.data$"col")
        text(  b.data$"xa7" ,  b.data$"ya", adj=1, substr(b.data$"location_confidence",1,4)    , cex=1.0, col= b.data$"col")
        text(  b.data$"xa8" ,  b.data$"ya", adj=1, substr(b.data$"orientation_confidence",1,4) , cex=1.0, col= b.data$"col")
        text(  b.data$"xa9" ,  b.data$"ya", adj=1, round(b.data$"perc.bp"/100,2)        ,       , cex=1.0, col= b.data$"col")
        text(  as.numeric(b.data[1, grep("xa", names(b.data))]) , rep( max(b.data$"ya") +x/2, sum( grepl("xa", names(b.data))) ), adj=1  , c("","CTG.len","cov.bp","max.win","%cov","grp","loc","orie","%ovlp") , cex=1.2 )
        text(  chr.max + 7.00*xx                                ,      max(b.data$"ya")  +x                                     , adj=0.5, "confidence scores" , cex=1.4 )
   
        # --> SCAFFOLD    
            polygon( c(               0,      chr.max  ,        chr.max,               0) ,c(    x-x/2,     x-x/2,     x+x/2,    x+x/2), lwd=2 , col=gray(0.80))      
n=nrow(g) ; polygon( c(g$"ref.start"[1], g$"ref.end"[n], g$"ref.end"[n],g$"ref.start"[1]) ,c(g$"za"[1], g$"za"[n], g$"zb"[n],g$"zb"[1]), lwd=2 , col=gray(0.90))
           segments( g$"ref.start" , g$"za", g$"ref.start"  , g$"zb", lwd=2 , col="black") 
    #   for (j in 1:nrow(g)) { polygon( c(g$"ref.start"[j], g$"ref.end"[j], g$"ref.end"[j],g$"ref.start"[j]) ,c(g$"za"[j], g$"za"[j], g$"zb"[j],g$"zb"[j]), lwd=2 , col=gray(0.90))     }  # very slow

        # --> SCAFFOLD - ctg names
        CTG.name = FALSE
        if( CTG.name==TRUE) { text(g$"ref.med", g$"zm", g$"QUERY",cex=0.8, col=g$"col",adj=0)  }
                               segments( g$"ref.ctg_start"   ,  g$"ZS"         , g$"ref.ctg_end"   , g$"ZS"          , lwd=2, col=g$"col"   ) ; text(g$"ref.ctg_start"   +chr.txt/9 ,  g$"ZS"         , g$"query"   , cex=1.0, col=g$"col",adj=1  ) ; text(g$"ref.ctg_end"  -chr.txt/30 ,  g$"ZS" , g$"Q.LEN" , cex=1.0, col=g$"col", adj=0) 
# for (j in 1:nrow(g)) { segments( g$"ref.ctg_start"[j],  g$"zs"[j]-j/100, g$"ref.ctg_end"[j], g$"zs"[j]-j/100 , lwd=2, col=g$"col"[j]) ; text(g$"ref.ctg_start"[j]+chr.txt/4 ,  g$"zs"[j]-j/100, g$"query"[j], cex=0.7, col=g$"col"[j]) }
      # for (j in 1:nrow(g)) { segments( g$"ref.ctg_start"[j],  g$"zs"[j]-j/100, g$"ref.ctg_end"[j], g$"zs"[j]-j/100 , lwd=2, col=g$"col"[j]) } # text(g$"ref.ctg_start"[j]+chr.txt/4 ,  g$"zs"[j]-j/100, g$"query"[j], cex=0.7, col=g$"col"[j]) }
      
      legend("topleft"   , "Alignment of contigs on the reference", bty="n", cex=3)
      legend("bottomleft", "Contig position on the Scaffold",bty="n", cex=3)
      legend("topright"   , c("FW contigs fragments", "REV contigs fragments"), fill=c("blue4","darkred"), cex=1.75)
#      legend("bottomright", c("FW contigs","FW contigs fragments","FW contigs fragments oriented","REV contigs","REV contigs fragments ","REV contigs fragments oriented") , fill=c("blue4","cornflowerblue","lightskyblue","darkred","firebrick3","pink"), cex=1.2)
      legend("bottomright", c("FW contigs","FW contigs fragments","REV contigs","REV contigs fragments ") , fill=c("blue4","cornflowerblue","darkred","firebrick3"), cex=1.75)
      legend("top", paste( "min plotted window =",MIN.WIND,"bp"), cex=1.0, bty="n")
    
    if( SINGLE_CHR_PLOT == TRUE )    dev.off()
    if( SINGLE_CHR_PLOT == TRUE )    system( paste0( " evince ",PDF) )
        
    
      # A = a[ a$"query" %in% "utg000010l",]
        A=a
        A = A[ A$"ctg.len" > MIN.WIND , ]
        A = A[ order(A$"ref.start") , ]
        A$"COL" = gsub("darkred","firebrick3",A$"col") 
        A$"COL" = gsub("blue4","cornflowerblue",A$"COL") 
        y1 = x 
        y2 = x/2 
        if (nrow(A)> 0)
        {
        if (INVERT_MINUS_CTG ==FALSE)
        {
            segments( A$"ref.start",  rep(y2 ,nrow(A)) , A$"scaf.pos", rep(-y2,nrow(A)) , lwd=1, col=gray(0.30))
            segments( A$"ref.end"  ,  rep(y2 ,nrow(A)) , A$"scaf.end", rep(-y2,nrow(A)) , lwd=1, col=gray(0.30))
        
            for (j in 1:nrow(A)) { polygon( c(A$"ref.start"[j], A$"ref.end"[j], A$"ref.end"[j],A$"ref.start"[j])  ,c( y2 , y2, y1, y1)   , lwd=1 , col=A$"COL"[j])     }
            for (j in 1:nrow(A)) { polygon( c(A$"scaf.pos" [j], A$"scaf.end"[j], A$"scaf.end"[j],A$"scaf.pos"[j]) ,c( -y1 , -y1, -y2, -y2), lwd=1 , col=A$"COL"[j])     }
            } else {
           
            # invert Strand of  MINUS contigs
            
            min.contigs = g[ g$"strand"=="-", ]$"query"
          
            B = A
           # B = B[ B$"query" == "utg000024l" , ]
            B$"Q.START"  = B$"q.start" 
            B$"Q.END"    = B$"q.end" 
            B$"SCAF.POS" = B$"scaf.pos"
            B$"SCAF.END" = B$"scaf.end" 
            w.min= which(B$"query" %in% min.contigs)
            
            if ( length(min.contigs)>0 )
            {  
                B$"Q.START" [w.min]  = B$"q.len"[w.min] - B$"q.end"[w.min] +1
                B$"Q.END"   [w.min]  = B$"q.len"[w.min] - B$"q.start"[w.min] +1
                B$"SCAF.POS"[w.min] = B$"scaf.start"[w.min] + B$"Q.START"[w.min] -1
                B$"SCAF.END"[w.min] = B$"scaf.start"[w.min] + B$"Q.END"[w.min] -1
#               B$"COL"[w.min]  = gsub("firebrick3","pink",B$"COL"[w.min] ) 
#                B$"COL"[w.min]  = gsub("cornflowerblue","lightskyblue",B$"COL"[w.min] ) 
                }
            segments( B$"ref.start",  rep(y2 ,nrow(B)) , B$"SCAF.POS", rep(-y2,nrow(B)) , lwd=1, col=gray(0.30))
            segments( B$"ref.end"  ,  rep(y2 ,nrow(B)) , B$"SCAF.END", rep(-y2,nrow(B)) , lwd=1, col=gray(0.30))
            for (j in 1:nrow(A)) { polygon( c(A$"ref.start"[j], A$"ref.end"[j], A$"ref.end"[j],A$"ref.start"[j])  ,c( y2 , y2, y1, y1)   , lwd=1 , col=A$"COL"[j])     }
            for (j in 1:nrow(B)) { polygon( c(B$"SCAF.POS" [j], B$"SCAF.END"[j], B$"SCAF.END"[j],B$"SCAF.POS"[j]) ,c( -y1 , -y1, -y2, -y2), lwd=1 , col=B$"COL"[j])     }
          }
          }
       }
     }
    CTG.DATA = do.call(rbind,B.DATA)
    for ( cc in grep("xa", names(CTG.DATA), value=T)) CTG.DATA[[cc]]=NULL
    CTG.DATA[[cc]]=NULL
    CTG.DATA[[cc]]=NULL
    write.table( CTG.DATA , paste0(out.dir,"/Chromosome_details.",ii,".",OUT,"ALL_CHR.overlapping_ctg.txt")  ,col.names=TRUE, row.names=F, quote=F, sep="\t")
    dev.off()
cat("\n")    
}

# system (  paste0( "evince ../../../RagTag_scaffold_check.",ii,".ALL_CHR.oriented.pdf"))



cat("\n")
cat("############################################################################\n")
cat("### 6) Single Chromosome reports including unplaced                     ###\n")
cat("############################################################################\n")
cat("\n")

# MIN.WIND = 

for( ii in names(CHR))
{

# ii =  "hifi200k";

# ii =  names(CHR)[2];
# ii =  names(CHR)[1];

cat(ii,"")

chr = data.frame(CHR[[ii]]) ; 
agp = setDT(AGP[[ii]]) ; 
cnf = setDT(CNF[[ii]]) ; 
asm = setDT(ASM[[ii]]) ; 
chr = chr[ rev(order(chr$"chr")) ,]
chr$"chr" = gsub("_RagTag","",chr$"chr")

CCmax.max = max(CCmax, chr$"len")

# chr = chr[ chr$"chr" != "Chr0", ]
bar=barplot(chr$"len", names.arg=chr$"chr", las=1, horiz=T, main=i, xlim=c(0, CCmax.max*1.10),col=gray(0.95) )
b2 = data.frame(chr.pos = bar[,1])
b2$"chr" = chr$"chr"
rrr = merge(ref, b2,by="chr", all.x=TRUE)

# separate GAPs annotation from AGP
agp$"chr" = gsub("_RagTag","",agp$"ref.chr")
gap = agp[ agp$"query" =="100",]
agp = agp[ agp$"query" !="100",]

# ASM2 select in ASM the same contigs as AGP 
# ASM3 select in ASM the same contigs as AGP 


agp2 = split( agp, agp$"chr")
asm2 = split( asm, asm$"ref.chr")
asm3 = split( asm, asm$"ref.chr")


CTG_unpl_filt =100000

sapply(asm2,nrow)
for (i in names(asm2)) { asm2[[i]] = asm2[[i]] [  asm2[[i]]$"query" %in% agp2[[i]]$"query",    ]  }
for (i in names(asm3)) { asm3[[i]] = asm3[[i]] [ !asm3[[i]]$"query" %in% agp2[[i]]$"query",    ]  }
# cat("\n")
# print(sapply(asm2,nrow))
# print(sapply(asm3,nrow))

for (i in names(asm3)) { asm3[[i]] = asm3[[i]] [ !asm3[[i]]$"ctg.len" >= CTG_unpl_filt/2,    ]  }
ASM3 = do.call(rbind,asm3)
ASM2 = do.call(rbind,asm2)

ASM2 =rbind(ASM2,ASM3)
asm2 = split( ASM2, ASM2$"ref.chr")

# print(sapply(asm3,nrow))
# print(sapply(asm2,nrow))
# 

#asm2 =asm3

INVERT_MINUS_CTG = TRUE
SINGLE_CHR_PLOT = FALSE

unplaced.ctg = CTG[[ii]]
unplaced.ctg = unplaced.ctg[ unplaced.ctg$"len">= CTG_unpl_filt,]
unplaced.ctg = unplaced.ctg[ rev(order(unplaced.ctg$"len")) ,]



 if (INVERT_MINUS_CTG ==FALSE) pdf( paste0(out.dir,"/GS-viewer_06.Chromosome_details_with_unplaced.",ii,".",OUT,"ALL_CHR.pdf"), width=plot.w, height=plot.h )
#if (INVERT_MINUS_CTG ==TRUE ) pdf( paste0(out.dir,"/GS-viewer_06.Chromosome_details_with_unplaced.",ii,".",OUT,"ALL_CHR.oriented.pdf"),width=plot.w, height=plot.h )
 if (INVERT_MINUS_CTG ==TRUE ) pdf( paste0(out.dir,"/GS-viewer_06.Chromosome_details_with_unplaced.",ii,".",OUT,"ALL_CHR.pdf"), width=plot.w, height=plot.h )

ctg_min_ass = min(agp$"q.len")


FINAL.AGP = list()
B.DATA   = list()
for (i in names(asm2))  { 
#   i = names(asm2)[1]
    cat(i,"")
    a = asm2[[i]]  
    g = agp2[[i]]
    a$"query" = as.character(a$"query")
    g$"query" = as.character(g$"query")
    
    a.unpl = intersect(unplaced.ctg$"chr", a$"query")
    if(length(a.unpl)==0) a.unpl=NA
    
   # print( a.unpl)
   # cat("\n")
#     print(a [ ! a$"query" %in%  unplaced.ctg$"chr" , ])
#     print(a [   a$"query" %in%  unplaced.ctg$"chr" , ])
#     
#     print( a.unpl)
#     print( a)
#    print( g)
#     print(nrow(g))
#     print(length(g))
#     print(dim(g))
#   print(is.null(g))
#     print(class(g))
#     print(unlist(g))
#     print(class(unlist(g)))
#     print(length(unlist(g)))
# 

    if (nrow(a)> 0 & length(g)> 0 & length(unlist(g))>0 )  # nrow(g) makes a error, when null is a empty list, so using length() is good for both data,frame and list
    {
        u = data.frame("ref.chr"= g$"ref.chr"[1],"ref.start"=NA,"ref.end"=NA,"part_number"=NA,"component_type"=NA,"query"=a.unpl,"q.start"=NA,"q.end"=NA,"strand"=NA,"q.len"=NA,"chr"=NA, stringsAsFactors=FALSE)
    
        # contigs on ref
        x = 0.02
        
        a$"ya" = x*2
        a$"yb" = a$"ya"
        
        a$"col"="blue4"
        if( any(a$"strand"=="-")) a[ a$"strand"=="-",]$"yb"= a[ a$"strand"=="-",]$"yb" *-1
        if( any(a$"strand"=="-")) a[ a$"strand"=="-",]$"col"="darkred"
        if( any(a$"strand"=="+")) a[ a$"strand"=="+",]$"col"="blue4"
        if( any(is.na(a$"strand"))) a[ is.na(a$"strand"),]$"col"="black"
        
        # sort by position of biggest fragment
        a$"ref.width" = a$"ref.end"- a$"ref.start"+1
        a$"q.left"    = a$"q.len"- a$"q.end"
        a$"strand.query" = paste(a$"query", a$"strand")
        a$"chr.ctg"      = paste0(a$"ref.chr",":", a$"query")
        a$"ctg.perc"     = round(100* a$"ctg.len"/ a$"q.len",2)
     
        as= split(a,a$"query")
        ass= lapply(as, function(x) x[ which.max(x$"ctg.len"), ])
        ass= sapply(ass, function(x) x$"ref.start")   
        as = as [ names(sort(ass)) ]
        
        # select only unplaced with a window > 10kb or cumulated >50% of contigs
      # as = as [ g$"query"]
        as = as [ c(g$"query",u$"query") ]
        as1 = as [ g$"query" ]
        as2 = as [ u$"query" ]
        as2.cum = sapply(as2, function(x) sum(x$"ctg.len") )
        as2.per = sapply(as2, function(x) sum(x$"ctg.perc") )
        as2.max = sapply(as2, function(x) max(x$"ctg.len") )

    #    w_max10k = which( as2.cum > 50)
        w_cum50p = which( as2.per > 60)
        w_max10K = which( as2.max > 20000)

        w_merge = union(names(w_cum50p),names(w_max10K))
        
        as = as [ c(g$"query",w_merge) ]
               
        # add cordinates for plot
        AS = list()
        for (j in 1:length(as)) { aa = as[[j]] ; aa=split(aa,aa$"strand.query") ; for (kk in names(aa)) AS[[kk]]=aa[[kk]]  }
        for (j in 1:length(AS)) { AS[[j]]$"ya"= AS[[j]]$"ya"+j/100  }
        a= do.call(rbind,AS)
        
         
        length(unique(a$query))
        length(unique(a$ya))
        length(unique(a$strand.query))
    
        n.a = length(as)
        col.a = rainbow(n.a)
        chr.max =a$"ref.len"[1]
        max.len  =max( chr.max, g$"ref.end" )
        MAX.LEN = max.len
        if( max(g$"ref.end") / chr.max < 1.20 ) { max.len = chr.max*1.25}
        chr.txt = -max.len/8
        n.ctg = length(unique(a$query))
        for (j in 1:length(as)) { as[[j]]$"ya"= as[[j]]$"ya"+j/100  }
        
        a= do.call(rbind,as)
          
        # scaffold
        g$"za"= -x/2
        g$"zb"= g$"za" -x
        g$"zm"=g$"zb"-g$"za"
     #  g$"zs"=g$"zb"-0.02
        g$"zs"=g$"zb" -x/2
        g$"ref.med"=(g$"ref.start"+g$"ref.end")/2
        g$"ref.ctg_start"=g$"ref.med"-g$"q.len"/2
        g$"ref.ctg_end"  =g$"ref.med"+g$"q.len"/2
        g$"ZS"=g$"zs" -(1:nrow(g))/100
        
        g$"col"="blue4"
        g$"QUERY"=g$"query"
        if( any(g$"strand"=="-")) g[ g$"strand"=="-",]$"col"="darkred"
        if( any(g$"strand"=="+")) g[ g$"strand"=="+",]$"col"="blue4"
        if( any(g$"strand"=="-")) g[ g$"strand"=="-",]$"QUERY"= paste0(g[ g$"strand"=="-",]$"QUERY", ".rev")

        # invs 
         as2= lapply( split(a$"strand",a$"query"), unique )
         qg= split(a,a$"strand.query")
       
        gg = g[, c(2,3,6)] 
        names(gg) = c("scaf.start","scaf.end","query")
        a = merge(a, gg, by="query",all.x=T)
        a$"scaf.pos"  =a$"scaf.start"+a$"q.start"-1
        a$"scaf.end"  =a$"scaf.start"+a$"q.end"-1
        
        gg3.col= c("query","strand","col")
        gg3 = g[,..gg3.col] 
        names(gg3) = c("query","agp.str","apg.col")
        a = merge(a,gg3, by="query", all.x=T)
        a = a[ order(a$"ref.chr",a$"query",a$"ref.start",a$"ref.end") , ]
        as= split(a,a$"query")
    
        # contig names
        col.aa = c("query","q.len","strand.query","ya","yb","col","agp.str","apg.col")
        aa = a[ ,..col.aa]
        aa = aa[ ! duplicated(aa),]
        all( g$"q.len" %in% aa$"q.len" )
            
        #CTG real start-ends
        a$"ref.st.ext" = NA
        a$"ref.e.ext"   = NA
        w.minu = which(a$"agp.str"=="-")
        w.plus = which(a$"agp.str"=="+")
        if( length(w.plus)>0 ) { a$"ref.st.ext"[w.plus]=a$"ref.start"[w.plus]-a$"q.start"[w.plus] ; a$"ref.e.ext"[w.plus]=a$"ref.end"[w.plus]+a$"q.left" [w.plus] }
        if( length(w.minu)>0 ) { a$"ref.st.ext"[w.minu]=a$"ref.start"[w.minu]-a$"q.left" [w.minu] ; a$"ref.e.ext"[w.minu]=a$"ref.end"[w.minu]+a$"q.start"[w.minu] }
        as= split(a,a$"query")
        
        # unplaces
     #  print("ok")
        
               
        # PLOT   
        b = do.call(rbind, as) ;  b=b[ b$"ref.width">MIN.WIND*1.0 ,] 
        Yax1 = max( b$"ya") + 2*x
        Yax2 = min( g$"ZS") - 2*x*1.5
        if (Yax1 <= 0.1) Yax1 = 0.1
        if (Yax2 >= -0.1) Yax2 = -0.1
        Yax2 = -Yax1
        ##   SINGLE_CHR_PLOT = TRUE
        if( SINGLE_CHR_PLOT == TRUE )  dir.create(paste0("../../../ragtag_scaffold.",ii, "/PDF/"))
        if( SINGLE_CHR_PLOT == TRUE )  PDF =  paste0("../../../ragtag_scaffold.",ii, "/PDF/Scaffold_vs_reference_",i,".pdf")
        if( SINGLE_CHR_PLOT == TRUE )  pdf( PDF, 30,16)
            # --> REF TRACK
            
            
        # create new X max
        
        xx = max.len/40
        b.data = data.frame( "query"=names(as),stringsAsFactors=F)
        b.data$"len"  =   as.numeric(sapply(as, function(x) x$"q.len"[1]))   
         b.data$"LEN"  =   b.data$"len"
         ww = which(b.data$"len" >1e3 & b.data$"len" < 1e6) ; if(length(ww)>0) { b.data$"LEN"[ww] = paste( round(b.data$"len"[ww]/1e3) , "Kb") }
         ww = which(b.data$"len" >1e3 & b.data$"len" >=1e6) ; if(length(ww)>0) { b.data$"LEN"[ww] = paste( round(b.data$"len"[ww]/1e6) , "Mb") }
        u.data = b.data[ b.data$"query" %in% unplaced.ctg$"chr", ]
       # u.data$"col" = "black"  
       
        if (nrow(u.data)>0 )
        {
         u.data$"xa10"  =  chr.max + 10.10*xx
         max.unplaced = max(u.data$"xa10") + max(u.data$"len")+2.0*xx
         max.len = max(max.len,max.unplaced )
        # text( u.data$"xa10"+ max(u.data$"len")+0.4*xx  , mean(u.data$"ya"), adj=0, "unplaced contigs"                        , cex=1.0, col= "black")
         }
         
         
         
            
        plot(c(0, chr.max ), c(x,x), type="l", lwd=5, xlim=c(-max.len/10,max.len), ylim=c(Yax2, Yax1), lend=1, cex.axis=2)
        abline(h=0, col="gray")
        
        g$"Q.LEN"  = as.character(g$"q.len")
        if( any(g$"q.len">=1000   ) )  g[ g$"q.len">=1000    ,]$"Q.LEN"= paste( as.character(round( g[ g$"q.len">=1000    ,]$"q.len"/1000   )),"Kb")
        if( any(g$"q.len">=1000000) )  g[ g$"q.len">=1000000 ,]$"Q.LEN"= paste( as.character(round( g[ g$"q.len">=1000000 ,]$"q.len"/1000000)),"Mb")
      #  print(head(data.frame(g)))
        
        # --> CONTIG ALIGN on REF
       #  MIN.WIND = 5000
        # for (j in 1:length(as)){ b= as[[j]];  b=b[ b$"ref.width">0000     ,] ; segments( b$"ref.start", b$"ya"-0.000, b$"ref.end", b$"ya"-0.000, lwd=1, col=b$"col") }
    #      for (j in 1:length(as)){ b= as[[j]];  b=b[ b$"ref.width">MIN.WIND*2.0 ,] ; segments( b$"ref.st.ext", b$"ya"-0.000, b$"ref.e.ext", b$"ya"-0.000, lwd=0.5, col="gray" ) ; points( b$"ref.st.ext", b$"ya", pch="|", cex=0.3,col="gray" ); points( b$"ref.e.ext", b$"ya", pch="|", cex=0.3,col="gray" ) }
    #      for (j in 1:length(bs)){ b= bs[[j]];  b=b[ b$"ref.width">BED.WIND*1.5 ,] ; segments( b$"ref.start" , b$"ya"-0.000, b$"ref.end"  , b$"ya"-0.000, lwd=5.0, col="gold") }
    #      for (j in 1:length(as)){ b= as[[j]];  b=b[ b$"ref.width">MIN.WIND*1.0 ,] ; segments( b$"ref.start" , b$"ya"-0.000, b$"ref.end"  , b$"ya"-0.000, lwd=2.0 , col=b$"col") }
         rs = as
         rs = lapply(rs, function(x) head(x[ rev( order(x$"ref.width")),]) )
     #   b = do.call(rbind, rs) ;  b=b[ b$"ref.width">MIN.WIND*2.0 ,] ; segments( b$"ref.st.ext", b$"ya"-0.000, b$"ref.e.ext", b$"ya"-0.000, lwd=0.5, col="gray" ) ; points( b$"ref.st.ext", b$"ya", pch="|", cex=0.3,col="gray" ); points( b$"ref.e.ext", b$"ya", pch="|", cex=0.3,col="gray" ) 
     #   b = do.call(rbind, bs) ;  b=b[ b$"ref.width">BED.WIND*1.5 ,] ; segments( b$"ref.start" , b$"ya"-0.000, b$"ref.end"  , b$"ya"-0.000, lwd=5.0, col="gold") 
         b = do.call(rbind, as) ;  b=b[ b$"ref.width">MIN.WIND*1.0 ,] ; segments( b$"ref.start" , b$"ya"-0.000, b$"ref.end"  , b$"ya"-0.000, lwd=2.0 , col=b$"col") 
         
         bbb.col=c("ref.chr","ref.start","ref.end","strand","query","q.len","ref.width")
         bbb = b[ , ..bbb.col]
         
        # check overlaps
        if( dir.exists(paste(out.dir,"/ctg_coordinates_on_ref",sep="")) == FALSE ) {  dir.create(paste(out.dir,"/ctg_coordinates_on_ref",sep=""))}
        bname = paste(out.dir,"/ctg_coordinates_on_ref/",ii,".",i,".ctg_coordinates_on_ref.txt",sep="")
        b.int = paste(out.dir,"/ctg_coordinates_on_ref/",ii,".",i,".ctg_coordinates_on_ref.intersect.txt",sep="")
        write.table( bbb , bname  ,col.names=FALSE, row.names=F, quote=F, sep="\t")
        system( paste0( " bedtools intersect -wao -a ",bname," -b ",bname, " > ",b.int) )
        system( paste0( "wc -l ",b.int ), inter=TRUE) -> bi.length
        bi.length = as.numeric(as.character(unlist(strsplit(bi.length," ")))) [1]

        bbb = split(bbb,bbb$"query")
#         cat("print bbb \n")
#         print(length(bbb))
        if( length(bbb)>0) { bbb.stat = data.frame("query"= names(bbb),"n.all.fragment"      =as.numeric(sapply(bbb,nrow))                                               ,"all.align.bp"= as.numeric(sapply(bbb,function(x) sum(x$"ref.width"))), stringsAsFactors=FALSE) }
        bbb.empty = data.frame("query"=0,"n.overlap.fragment"  =0, "overlapp.align.bp"= 0, stringsAsFactors=FALSE)
        bbb.empty = bbb.empty[ bbb.empty$"query">0, ]
        if( bi.length > 0 & is.na(bi.length)==FALSE)
        {        
            bi = fread(b.int, stringsAsFactors=FALSE)
            names(bi) = c(   paste0("A", bbb.col) , paste0("B", bbb.col), "overlap" )
            bi = bi[ bi$"Aquery" != bi$"Bquery", ]
         #   bi = bi[ bi$"Aref.start" <= bi$"Bref.start" & bi$"Aref.end" >= bi$"Bref.end", ]
            bi$"fragment" = paste(bi$"Bquery" ,bi$"Bref.start",bi$"Bref.end",sep=":")
            
            bis = split(bi  ,bi$"Bquery")
            #bbb = split(bbb,bbb$"query")
            
#              print(length(bis))
            if( length(bis)>0) { bis.stat = data.frame("query"= names(bis),"n.overlap.fragment"  =as.numeric(sapply(bis,function(x) length(unique(x$"fragment")))), "overlapp.align.bp"= as.numeric(sapply(bis,function(x) sum(x[ !duplicated(x$"fragment"),]$"overlap") )+1), stringsAsFactors=FALSE) }
          # if( length(bbb)>0) { bbb.stat = data.frame("query"= names(bbb),"n.all.fragment"      =as.numeric(sapply(bbb,nrow))                                               ,"all.align.bp"= as.numeric(sapply(bbb,function(x) sum(x$"ref.width"))), stringsAsFactors=FALSE) }
            
            # if no overlap available, create empty df
            if( length(bis)==0) { bis.stat = data.frame("query"= names(bbb),"n.overlap.fragment"  =0, "overlapp.align.bp"= 0, stringsAsFactors=FALSE) }
            
            }   else  {  
                if(length(bbb)>0 ) { bbb.stat = data.frame("query"=names(bbb) ,"n.overlap.fragment"  =0, "overlapp.align.bp"= 0, stringsAsFactors=FALSE) }  ;
                if(length(bbb)==0) { bbb.stat = bbb.empty }  ;
                }  
         
#         print(bbb.stat)
#         print(bis.stat)
#          
        if( nrow(bbb.stat) > 0 & nrow(bis.stat)> 0)
        {
            fin.stat = merge(bbb.stat,bis.stat,all.x=T, by="query")
            fin.stat$"perc.bp" = round(100 *  fin.stat$"overlapp.align.bp" /  fin.stat$"all.align.bp", 2) 
            fin.stat$"chr" = i 
            fin.stat$"sample" = ii
            }  else { fin.stat = data.frame("query"= 0, "n.all.fragment" = 0, "all.align.bp" = 0 , "n.overlap.fragment" = 0 , "overlapp.align.bp" = 0 , "perc.bp" = 0 , "chr" = 0 , "sample" = 0, stringsAsFactors=FALSE ) ; fin.stat = fin.stat[ fin.stat$"query">0, ]
          }

#         print(fin.stat)

 #      text( -max.len/60 ,  aa$"ya", adj=1, paste(aa$"query","vs Ref."), cex=0.5, col= aa$"apg.col")
        text( -max.len/60 ,  aa$"ya", adj=1, paste(aa$"query"          ), cex=1.0, col= aa$"apg.col")
        text( -max.len/60 ,  x      , adj=1, paste("Ref.",a$"ref.chr"),cex=1.5)
        text( -max.len/60 , -x      , adj=1, paste("Scaffold",a$"ref.chr"),cex=1.5)
        
        # CTG stat
        xx = max.len/40
        b.data = data.frame( "query"=names(as),stringsAsFactors=F)
        b.data$"len"  =   as.numeric(sapply(as, function(x) x$"q.len"[1]))   
        b.data$"LEN"  =   b.data$"len"
        ww = which(b.data$"len" >1e3 & b.data$"len" < 1e6) ; if(length(ww)>0) { b.data$"LEN"[ww] = paste( round(b.data$"len"[ww]/1e3) , "Kb") }
        ww = which(b.data$"len" >1e3 & b.data$"len" >=1e6) ; if(length(ww)>0) { b.data$"LEN"[ww] = paste( round(b.data$"len"[ww]/1e6) , "Mb") }
        b.data$"ya"  =   as.numeric(sapply(as, function(x) x$"ya"[1])) - 0.000  
        b.data$"xa1"  =  chr.max + 0.20*xx
        b.data$"xa2"  =  chr.max + 2.70*xx
        b.data$"xa3"  =  chr.max + 3.70*xx
        b.data$"xa4"  =  chr.max + 4.80*xx
        b.data$"xa5"  =  chr.max + 5.80*xx
        b.data$"xa6"  =  chr.max + 6.80*xx
        b.data$"xa7"  =  chr.max + 7.70*xx
        b.data$"xa8"  =  chr.max + 8.60*xx
        b.data$"xa9"  =  chr.max + 9.50*xx
        b.data$"col" = as.character(sapply(as, function(x) x$"col"[1]))
        b.data$"min.win" = as.numeric(sapply(as, function(x) min(x$"ref.width")))              
        b.data$"max.win" = as.numeric(sapply(as, function(x) max(x$"ref.width")))              
        b.data$"tot.cov" = as.numeric(sapply(as, function(x) sum(x$"ref.width")))              
        b.data$"perc.cov"= round(100*b.data$"tot.cov"/b.data$"len",0)              
        b.data$"perc.min"= round(100*b.data$"min.win"/b.data$"len",1)              
        b.data$"perc.max"= round(100*b.data$"max.win"/b.data$"len",1)              
        b.data$"perc.cov2"= paste0(b.data$"perc.cov","%")              
        b.data = merge(b.data, cnf, by="query",all.x=TRUE)
        b.data$"string"= paste0(b.data$"query"," %")  
        b.data$"MAX.WIN"  =   b.data$"max.win"
        ww = which(b.data$"max.win" >1e3 & b.data$"max.win" < 1e6) ; if(length(ww)>0) { b.data$"MAX.WIN"[ww] = paste( round(b.data$"max.win"[ww]/1e3) , "Kb") }
        ww = which(b.data$"max.win" >1e3 & b.data$"max.win" >=1e6) ; if(length(ww)>0) { b.data$"MAX.WIN"[ww] = paste( round(b.data$"max.win"[ww]/1e6) , "Mb") }
        b.data$"TOT.COV"  =   b.data$"tot.cov"
        ww = which(b.data$"tot.cov" >1e3 & b.data$"tot.cov" < 1e6) ; if(length(ww)>0) { b.data$"TOT.COV"[ww] = paste( round(b.data$"tot.cov"[ww]/1e3) , "Kb") }
        ww = which(b.data$"tot.cov" >1e3 & b.data$"tot.cov" >=1e6) ; if(length(ww)>0) { b.data$"TOT.COV"[ww] = paste( round(b.data$"tot.cov"[ww]/1e6) , "Mb") }
        b.data = merge(b.data, fin.stat, by="query",all.x=TRUE)

        B.DATA[[i]] = b.data
        
#       print(b.data)
        
        
        text(  b.data$"xa1" ,  b.data$"ya", adj=0, b.data$"query"                              , cex=1.0, col= b.data$"col")
        text(  b.data$"xa2" ,  b.data$"ya", adj=1, b.data$"LEN"                                , cex=1.0, col= b.data$"col")
        text(  b.data$"xa3" ,  b.data$"ya", adj=1, b.data$"TOT.COV"                            , cex=1.0, col= b.data$"col")
        text(  b.data$"xa4" ,  b.data$"ya", adj=1, b.data$"MAX.WIN"                            , cex=1.0, col= b.data$"col")
        text(  b.data$"xa5" ,  b.data$"ya", adj=1, b.data$"perc.cov2"                          , cex=1.0, col= b.data$"col")
        text(  b.data$"xa6" ,  b.data$"ya", adj=1, substr(b.data$"grouping_confidence",1,4)    , cex=1.0, col= b.data$"col")
        text(  b.data$"xa7" ,  b.data$"ya", adj=1, substr(b.data$"location_confidence",1,4)    , cex=1.0, col= b.data$"col")
        text(  b.data$"xa8" ,  b.data$"ya", adj=1, substr(b.data$"orientation_confidence",1,4) , cex=1.0, col= b.data$"col")
        text(  b.data$"xa9" ,  b.data$"ya", adj=1, round(b.data$"perc.bp"/100,2)        ,       , cex=1.0, col= b.data$"col")
        text(  as.numeric(b.data[1, grep("xa", names(b.data))]) , rep( max(b.data$"ya") +x/2, sum( grepl("xa", names(b.data))) ), adj=1  , c("","CTG.len","cov.bp","max.win","%cov","grp","loc","orie","%ovlp") , cex=1.2 )
        text(  chr.max + 7.00*xx                                ,      max(b.data$"ya")  +x                                     , adj=0.5, "confidence scores" , cex=1.4 )

       u.data = b.data[ b.data$"query" %in% unplaced.ctg$"chr", ]
   #    u.data$"col" = "black"  


       
        if (nrow(u.data)>0 )
        {
         u.data$"xa10"  =  chr.max + 10.10*xx
         segments( u.data$"xa10" , u.data$"ya", u.data$"xa10"+u.data$"len" , u.data$"ya", lwd=2 , col=u.data$"col" ) 
         segments( u.data$"xa10"+ max(u.data$"len")+0.2*xx , min(u.data$"ya"), u.data$"xa10"+ max(u.data$"len")+0.2*xx , max(u.data$"ya"), lwd=2 , col= "black" )  
         segments( -max.len/10 , min(u.data$"ya")-0.005, u.data$"xa10"+ max(u.data$"len")+1.2*xx ,min(u.data$"ya")-0.005, lwd=1 , col= "black")  
         text( u.data$"xa10"+ max(u.data$"len")+0.4*xx  , mean(u.data$"ya"), adj=0, "unplaced contigs"                        , cex=1.0, col= "black")
         }
         
         
        # --> SCAFFOLD    
            polygon( c(               0,      chr.max  ,        chr.max,               0) ,c(    x-x/2,     x-x/2,     x+x/2,    x+x/2), lwd=2 , col=gray(0.80))      
n=nrow(g) ; polygon( c(g$"ref.start"[1], g$"ref.end"[n], g$"ref.end"[n],g$"ref.start"[1]) ,c(g$"za"[1], g$"za"[n], g$"zb"[n],g$"zb"[1]), lwd=2 , col=gray(0.90))
           segments( g$"ref.start" , g$"za", g$"ref.start"  , g$"zb", lwd=2 , col="black") 
    #   for (j in 1:nrow(g)) { polygon( c(g$"ref.start"[j], g$"ref.end"[j], g$"ref.end"[j],g$"ref.start"[j]) ,c(g$"za"[j], g$"za"[j], g$"zb"[j],g$"zb"[j]), lwd=2 , col=gray(0.90))     }  # very slow

        # --> SCAFFOLD - ctg names
        CTG.name = FALSE
        if( CTG.name==TRUE) { text(g$"ref.med", g$"zm", g$"QUERY",cex=0.8, col=g$"col",adj=0)  }
                               segments( g$"ref.ctg_start"   ,  g$"ZS"         , g$"ref.ctg_end"   , g$"ZS"          , lwd=2, col=g$"col"   ) ; text(g$"ref.ctg_start"   +chr.txt/9 ,  g$"ZS"         , g$"query"   , cex=1.0, col=g$"col",adj=1  ) ; text(g$"ref.ctg_end"  -chr.txt/30 ,  g$"ZS" , g$"Q.LEN" , cex=1.0, col=g$"col", adj=0) 
# for (j in 1:nrow(g)) { segments( g$"ref.ctg_start"[j],  g$"zs"[j]-j/100, g$"ref.ctg_end"[j], g$"zs"[j]-j/100 , lwd=2, col=g$"col"[j]) ; text(g$"ref.ctg_start"[j]+chr.txt/4 ,  g$"zs"[j]-j/100, g$"query"[j], cex=0.7, col=g$"col"[j]) }
      # for (j in 1:nrow(g)) { segments( g$"ref.ctg_start"[j],  g$"zs"[j]-j/100, g$"ref.ctg_end"[j], g$"zs"[j]-j/100 , lwd=2, col=g$"col"[j]) } # text(g$"ref.ctg_start"[j]+chr.txt/4 ,  g$"zs"[j]-j/100, g$"query"[j], cex=0.7, col=g$"col"[j]) }
      
      legend("topleft"   , "Alignment of contigs on the reference", bty="n", cex=3)
      legend("bottomleft", "Contig position on the Scaffold",bty="n", cex=3)
      legend("topright"   , c("FW contigs fragments", "REV contigs fragments"), fill=c("blue4","darkred"), cex=1.75)
#      legend("bottomright", c("FW contigs","FW contigs fragments","FW contigs fragments oriented","REV contigs","REV contigs fragments ","REV contigs fragments oriented") , fill=c("blue4","cornflowerblue","lightskyblue","darkred","firebrick3","pink"), cex=1.2)
      legend("bottomright", c("FW contigs","FW contigs fragments","REV contigs","REV contigs fragments ") , fill=c("blue4","cornflowerblue","darkred","firebrick3"), cex=1.75)
      legend("top", paste( "min plotted window =",MIN.WIND,"bp"), cex=1.0, bty="n")
    
    if( SINGLE_CHR_PLOT == TRUE )    dev.off()
    if( SINGLE_CHR_PLOT == TRUE )    system( paste0( " evince ",PDF) )
        
    
      # A = a[ a$"query" %in% "utg000010l",]
        A=a
        A = A[ A$"ctg.len" > MIN.WIND , ]
        A = A[ order(A$"ref.start") , ]
        A$"COL" = gsub("darkred","firebrick3",A$"col") 
        A$"COL" = gsub("blue4","cornflowerblue",A$"COL") 
        y1 = x 
        y2 = x/2 
        if (nrow(A)> 0)
        {
        if (INVERT_MINUS_CTG ==FALSE)
        {
            segments( A$"ref.start",  rep(y2 ,nrow(A)) , A$"scaf.pos", rep(-y2,nrow(A)) , lwd=1, col=gray(0.30))
            segments( A$"ref.end"  ,  rep(y2 ,nrow(A)) , A$"scaf.end", rep(-y2,nrow(A)) , lwd=1, col=gray(0.30))
        
            for (j in 1:nrow(A)) { polygon( c(A$"ref.start"[j], A$"ref.end"[j], A$"ref.end"[j],A$"ref.start"[j])  ,c( y2 , y2, y1, y1)   , lwd=1 , col=A$"COL"[j])     }
            for (j in 1:nrow(A)) { polygon( c(A$"scaf.pos" [j], A$"scaf.end"[j], A$"scaf.end"[j],A$"scaf.pos"[j]) ,c( -y1 , -y1, -y2, -y2), lwd=1 , col=A$"COL"[j])     }
            } else {
           
            # invert Strand of  MINUS contigs
            
            min.contigs = g[ g$"strand"=="-", ]$"query"
          
            B = A
           # B = B[ B$"query" == "utg000024l" , ]
            B$"Q.START"  = B$"q.start" 
            B$"Q.END"    = B$"q.end" 
            B$"SCAF.POS" = B$"scaf.pos"
            B$"SCAF.END" = B$"scaf.end" 
            w.min= which(B$"query" %in% min.contigs)
            
            if ( length(min.contigs)>0 )
            {  
                B$"Q.START" [w.min]  = B$"q.len"[w.min] - B$"q.end"[w.min] +1
                B$"Q.END"   [w.min]  = B$"q.len"[w.min] - B$"q.start"[w.min] +1
                B$"SCAF.POS"[w.min] = B$"scaf.start"[w.min] + B$"Q.START"[w.min] -1
                B$"SCAF.END"[w.min] = B$"scaf.start"[w.min] + B$"Q.END"[w.min] -1
#               B$"COL"[w.min]  = gsub("firebrick3","pink",B$"COL"[w.min] ) 
#                B$"COL"[w.min]  = gsub("cornflowerblue","lightskyblue",B$"COL"[w.min] ) 
                }
            segments( B$"ref.start",  rep(y2 ,nrow(B)) , B$"SCAF.POS", rep(-y2,nrow(B)) , lwd=1, col=gray(0.30))
            segments( B$"ref.end"  ,  rep(y2 ,nrow(B)) , B$"SCAF.END", rep(-y2,nrow(B)) , lwd=1, col=gray(0.30))
            for (j in 1:nrow(A)) { polygon( c(A$"ref.start"[j], A$"ref.end"[j], A$"ref.end"[j],A$"ref.start"[j])  ,c( y2 , y2, y1, y1)   , lwd=1 , col=A$"COL"[j])     }
            for (j in 1:nrow(B)) { polygon( c(B$"SCAF.POS" [j], B$"SCAF.END"[j], B$"SCAF.END"[j],B$"SCAF.POS"[j]) ,c( -y1 , -y1, -y2, -y2), lwd=1 , col=B$"COL"[j])     }
          }
          }
       }
     }
    CTG.DATA = do.call(rbind,B.DATA)
    for ( cc in grep("xa", names(CTG.DATA), value=T)) CTG.DATA[[cc]]=NULL
    CTG.DATA[[cc]]=NULL
    CTG.DATA[[cc]]=NULL
    write.table( CTG.DATA , paste0(out.dir,"/Chromosome_details_with_unplaced.",ii,".",OUT,"ALL_CHR.overlapping_ctg.txt")  ,col.names=TRUE, row.names=F, quote=F, sep="\t")
    dev.off()
cat("\n")    
}

# system (  paste0( "evince ../../../RagTag_scaffold_check.",ii,".ALL_CHR.oriented.pdf"))




cat("\n\n\n")
cat("############################################################################\n")
cat("### 7) UNPLACED CONTIGS                                                  ###\n")
cat("############################################################################\n")


### UNPLACED CONTIGS

#PLOT_TYPES = c("single", "multiple","separate")

for (plot_type in PLOT_TYPES)
{
if ( plot_type == "single"  ) { cat("\n ==> Plot type:   ",plot_type,"") ; pdf( paste0(out.dir,"/GS-viewer_07.Unplaced_contigs.single.",OUT,"pdf"  ) ,width=plot.w, height=plot.h) ; par(mfcol=c(1    ,1    ), mar=c(2,8,4,1) , oma=c(1,4,1,1) ); LG.size=1 }
if ( plot_type == "multiple") { cat("\n ==> Plot type: "  ,plot_type,"") ; pdf( paste0(out.dir,"/GS-viewer_07.Unplaced_contigs.multiple.",OUT,"pdf") ,width=plot.w, height=plot.h) ; par(mfcol=c(n.row,n.col), mar=c(2,9,4,1) , oma=c(1,4,1,1) ); LG.size=0.5 }

for( ii in names(CHR))
{

# ii =  names(CHR)[2];
if ( plot_type %in% c("single","multiple")) {  cat(ii,"")  }

# recreate bar --> barplot to get REF X.axis (before it was used for SCAFFOLD instead of reference and it may miss 1 or more chr )
bar=barplot(ref$"ref.len", names.arg=ref$"chr", plot=F)
b2 = data.frame(chr.pos = bar[,1])
b2$"chr" = ref$"chr"
rrr = merge(ref,b2, by="chr")



chr = data.frame(CHR[[ii]]) ; 
agp = setDT(AGP[[ii]]) ; 
cnf = setDT(CNF[[ii]]) ; 
asm = setDT(ASM[[ii]]) ; 
chr = chr[ rev(order(chr$"chr")) ,]
chr$"chr" = gsub("_RagTag","",chr$"chr")

CCmax.max = max(CCmax, chr$"len")


# separate GAPs annotation from AGP
agp$"chr" = gsub("_RagTag","",agp$"ref.chr")
gap = agp[ agp$"query" =="100",]
agp = agp[ agp$"query" !="100",]

# ASM2 select in ASM the same      contigs as AGP 
# ASM3 select in ASM the different contigs as AGP 

agp2 = split( agp, agp$"chr")
asm2 = split( asm, asm$"ref.chr")
asm3 = split( asm, asm$"ref.chr")

# sapply(asm2,nrow)
for (i in names(asm2)) { asm2[[i]] = asm2[[i]] [  asm2[[i]]$"query" %in% agp2[[i]]$"query",    ]  }
for (i in names(asm3)) { asm3[[i]] = asm3[[i]] [ !asm3[[i]]$"query" %in% agp2[[i]]$"query",    ]  }
# sapply(asm2,nrow)


CTG_unpl_filt = 500000 #0.5e+06
unplaced.ctg = CTG[[ii]]
unplaced.ctg = unplaced.ctg[ unplaced.ctg$"len">= CTG_unpl_filt,]
unplaced.ctg = unplaced.ctg[ rev(order(unplaced.ctg$"len")) ,]

# top 10 only
unplaced.ctg = head(unplaced.ctg,20)

n.ctg = nrow(unplaced.ctg)
n.nnn = n.ctg - length(COL.LIST)
if (n.nnn <=0 ) n.nnn=1

pos = 1.5 
if( n.ctg>6)  pos = 1
if( n.ctg>9)  pos = 1.1
if( n.ctg>18) pos = 1.3
if( n.ctg>24) pos = 4

unplaced.ctg$"pos" = -( 1:nrow(unplaced.ctg)/pos) -0.5
unplaced.ctg$"start" = 0
unplaced.ctg$"col" = c(COL.LIST, rev(rainbow(n.nnn)))[ 1:n.ctg]
unplaced.ctg$"name" = paste0(unplaced.ctg$"chr"," (", unplaced.ctg$"len", "bp)")
unplaced.ctg$"query" = unplaced.ctg$"chr"


asm3 = lapply(asm3, function(x) { x$"ctg.perc"=round( x$"ctg.len"/ x$"q.len",2); return(x) })
ASM3 = do.call(rbind,asm3)
ASM4 = ASM3[ ASM3$"query" %in% unplaced.ctg$"chr", ]
ASM4 = ASM4[ ASM4$"ctg.len">=MIN.WIND, ]
#ASM4 = ASM4[ ASM4$"ctg.perc">=0.20, ]

ASM4 = merge(ASM4,unplaced.ctg, by="query", all.x=TRUE)
ASM4$"pos" = NULL

rr2= rrr
names(rr2) = c("ref.chr","ref.len","pos")
rr2$"ref.len" = NULL

ASM4 = merge(ASM4,rr2, by="ref.chr", all.x=TRUE)

# 
# print(chr) # fai file of scaffold  -> 13 ctg
# print(bar) # barplot hight created on ref lengths --> 15 ctgs
# print(rrr) # ref fai + barplot coordinates 
# print(rr2)
# print(ref)
# print(unplaced.ctg)

unpl_ctg = data.frame( "query"=unplaced.ctg$"query","ctg.pos"=unplaced.ctg$"pos", stringsAsFactors=FALSE )
ASM4 = merge(ASM4,unpl_ctg, by="query", all.x=TRUE)

if (nrow(unplaced.ctg)>0)
{
    if ( plot_type == "separate")  { cat("\n ==> Plot type: ",plot_type,"") ; cat(ii,"") ;  pdf( paste0(out.dir,"/Unplaced_contigs.",ii,".",OUT,"pdf"),width=plot.w, height=plot.h ); par(mfcol=c(1    ,1    ), mar=c(2,5,4,1), oma=c(1,3,1,1) ) }
#   pdf( paste0(out.dir,"/RagTag_scaffold__Unplaced_contigs.",ii,".",OUT,"pdf"),width=plot.w, height=plot.h )
#   par( mar=c(2,8,4,1))
    lll =  0.7+(0:(nrow(ref)-0)*1.2)
 #   ll2 = min(low_conf.ctg$"pos")
    ll2 = -15
    bar=barplot(ref$"ref.len", names.arg=ref$"chr", las=1, horiz=T, main="top 20 unplaced contigs alignment on Reference",col=gray(0.95), cex.main=2*LG.size,cex.axis=2,  xlim=c(0, max(ref$"ref.len")*1.15), cex.main=2*LG.size, cex.names=2*LG.size, ylim=c(ll2,max(lll)))
    abline( h=-0.3, col="black", lwd=5)
    text( max(ref$"ref.len")/2, max(lll),  ii, cex=3*LG.size )
    
    segments( unplaced.ctg$"start" , unplaced.ctg$"pos", unplaced.ctg$"len"  , unplaced.ctg$"pos", lwd=10*LG.size , col= unplaced.ctg$"col",lend=1) 
    text( unplaced.ctg$"len"+CCmax.max/100 , unplaced.ctg$"pos", unplaced.ctg$"name" , cex=1., adj=0, col=unplaced.ctg$"col")   
    segments( ASM4$"ref.start" , ASM4$"pos"    , ASM4$"ref.end"  , ASM4$"pos"    , lwd=16 , col= ASM4$"col",lend=1) 
    segments( ASM4$"q.start"   , ASM4$"ctg.pos", ASM4$"q.end"    , ASM4$"ctg.pos", lwd=20 , col= ASM4$"col",lend=1) 
    
    if ( plot_type == "separate")  dev.off()  
    }
    }
    
if ( plot_type %in% c("single","multiple")) {   dev.off() }
 }
 
cat("\n")

cat("\n\n\n")
cat("############################################################################\n")
cat("### 8) MISPLACED CONTIGS                                                  ###\n")
cat("############################################################################\n")


#PLOT_TYPES = c("single", "multiple","separate")

for (plot_type in PLOT_TYPES)
{
if ( plot_type == "single"  ) { cat("\n ==> Plot type:   ",plot_type,"") ; pdf( paste0(out.dir,"/GS-viewer_08.Misplaced_contigs.single.",OUT,"pdf"  ) ,width=plot.w, height=plot.h) ; par(mfcol=c(1    ,1    ), mar=c(2,8,4,1) , oma=c(1,4,1,1) ); LG.size=1  }
if ( plot_type == "multiple") { cat("\n ==> Plot type: "  ,plot_type,"") ; pdf( paste0(out.dir,"/GS-viewer_08.Misplaced_contigs.multiple.",OUT,"pdf") ,width=plot.w, height=plot.h) ; par(mfcol=c(n.row,n.col), mar=c(2,9,4,1) , oma=c(1,4,1,1) ); LG.size=0.5  }

for( ii in names(CHR))
{


# ii =  names(CHR)[2];
if ( plot_type %in% c("single","multiple")) {  cat(ii,"")  }

# recreate bar --> barplot to get REF X.axis (before it was used for SCAFFOLD instead of reference and it may miss 1 or more chr )
bar=barplot(ref$"ref.len", names.arg=ref$"chr", plot=F)
b2 = data.frame(chr.pos = bar[,1])
b2$"chr" = ref$"chr"
rrr = merge(ref,b2, by="chr")



chr = data.frame(CHR[[ii]]) ; 
agp = setDT(AGP[[ii]]) ; 
cnf = setDT(CNF[[ii]]) ; 
asm = setDT(ASM[[ii]]) ; 
chr = chr[ rev(order(chr$"chr")) ,]
chr$"chr" = gsub("_RagTag","",chr$"chr")

CCmax.max = max(CCmax, chr$"len")


# separate GAPs annotation from AGP
agp$"chr" = gsub("_RagTag","",agp$"ref.chr")
gap = agp[ agp$"query" =="100",]
agp = agp[ agp$"query" !="100",]

agp2 = split( agp, agp$"chr")
asm2 = split( asm, asm$"ref.chr")
asm3 = split( asm, asm$"ref.chr")

# ASM2 select in ASM the same      contigs as AGP 
# ASM3 select in ASM the different contigs as AGP 

#sapply(asm2,nrow)
for (i in names(asm2)) { asm2[[i]] = asm2[[i]] [  asm2[[i]]$"query" %in% agp2[[i]]$"query",    ]  }
for (i in names(asm3)) { asm3[[i]] = asm3[[i]] [ !asm3[[i]]$"query" %in% agp2[[i]]$"query",    ]  }
#sapply(asm2,nrow)


CTG_filt = 100000 #0.5e+05

# remove unplaced contigs
placed.ctg = agp[, c("query","q.len","ref.chr")]
names(placed.ctg) = c("chr","len","ref.chr")
placed.ctg = placed.ctg[ placed.ctg$"len">= CTG_filt,]
placed.ctg = placed.ctg[ rev(order(placed.ctg$"len")) ,]
placed.ctg = placed.ctg[ placed.ctg$"chr" != placed.ctg$"ref.chr",]



# top 20 only
# placed.ctg = head(  placed.ctg,20)
# calculate in a different way

# n.ctg = nrow(placed.ctg)
# n.nnn = n.ctg - length(COL.LIST)
# if (n.nnn <=0 ) n.nnn=1
# 
# pos = 1.5 
# if( n.ctg>6)  pos = 1
# if( n.ctg>9)  pos = 1.1
# if( n.ctg>18) pos = 1.3
# if( n.ctg>24) pos = 4


MIN_WIND = 1000

asm3 = lapply(asm3, function(x) { x$"ctg.perc"=round( x$"ctg.len"/ x$"q.len",2); return(x) })
ASM3 = do.call(rbind,asm3)
ASM4 = ASM3[ ASM3$"query" %in% placed.ctg$"chr", ]
ASM4 = ASM4[ ASM4$"ctg.len">= MIN_WIND, ]


# calculate in a different way
# calculate total bo covered on other chromomsomes
asm4 = split( ASM4, ASM4$"query")
asm5 = lapply(asm4, function(x) sum(x$"ctg.len"))
asm5 = setDT(data.frame("query"=names(asm5), "tot.cov"=as.numeric(asm5),stringsAsFactors=FALSE))
asm5 = asm5[ asm5$"tot.cov">= CTG_filt/2 , ]
asm5 = asm5[ rev(order(asm5$"tot.cov")) ,]
asm5 = head(  asm5,20)

# filter 
ASM4 = ASM4[ ASM4$"query" %in% asm5$"query", ]
ASM4 = ASM4[ ASM4$"ctg.len">= CTG_filt/2, ]



placed.ctg$"query" = placed.ctg$"chr"
placed.ctg = placed.ctg[ placed.ctg$"query" %in% asm5$"query", ]
# placed.ctg = head(  placed.ctg, 20)



placed.ctg = placed.ctg[ placed.ctg$"query" %in% ASM4$"query", ]



n.ctg = nrow(placed.ctg)
n.nnn = n.ctg - length(COL.LIST)
if (n.nnn <=0 ) n.nnn=1
pos = 1.5 
if( n.ctg>6)  pos = 1
if( n.ctg>9)  pos = 1.1
if( n.ctg>18) pos = 1.3
if( n.ctg>24) pos = 4


placed.ctg$"pos" = -( 1:nrow(placed.ctg)/pos) -0.5
placed.ctg$"start" = 0
placed.ctg$"col" = c(COL.LIST, rev(rainbow(n.nnn)))[ 1:n.ctg]
placed.ctg$"name" = paste0(placed.ctg$"chr"," (",placed.ctg$"ref.chr"," - ", placed.ctg$"len", "bp)")

placed.ctg$"ref.chr" = NULL

#ASM4 = ASM4[ ASM4$"ctg.perc">=0.20, ]



 
ASM4 = merge(ASM4,placed.ctg, by="query", all.x=TRUE)
ASM4$"pos" = NULL


rr2= rrr
names(rr2) = c("ref.chr","ref.len","pos")
rr2$"ref.len" = NULL

ASM4 = merge(ASM4,rr2, by="ref.chr", all.x=TRUE)



# 
# print(chr) # fai file of scaffold  -> 13 ctg
# print(bar) # barplot hight created on ref lengths --> 15 ctgs
# print(rrr) # ref fai + barplot coordinates 
# print(rr2)
# print(ref)
# 

pl_ctg = data.frame( "query"=placed.ctg$"query","ctg.pos"=placed.ctg$"pos", stringsAsFactors=FALSE )

ASM4 = merge(ASM4,pl_ctg, by="query", all.x=TRUE)


if (nrow(placed.ctg)>0)
{
    if ( plot_type == "separate")  { cat("\n ==> Plot type: ",plot_type,"") ; cat(ii,"") ;  pdf( paste0(out.dir,"/Misplaced_contigs.",ii,".",OUT,"pdf"),width=plot.w, height=plot.h ); par(mfcol=c(1    ,1    ), mar=c(2,5,4,1), oma=c(1,3,1,1) ) }
#   pdf( paste0(out.dir,"/RagTag_scaffold__Unplaced_contigs.",ii,".",OUT,"pdf"),width=plot.w, height=plot.h )
#   par( mar=c(2,8,4,1))
    lll =  0.7+(0:(nrow(ref)-0)*1.2)
 #   ll2 = min(low_conf.ctg$"pos")
    ll2 = -15
    bar=barplot(ref$"ref.len", names.arg=ref$"chr", las=1, horiz=T, main="(Top20) Scaffolded contigs that have alignment also on other chromosomes of the Reference [>50kb fragments only]",col=gray(0.95), cex.main=2*LG.size, cex.axis=2,  xlim=c(0, max(ref$"ref.len")*1.15), cex.names=2*LG.size, ylim=c(ll2,max(lll)))
    abline( h=-0.3, col="black", lwd=5)
    text( max(ref$"ref.len")/2, max(lll),  ii, cex=3*LG.size )
    
    segments( placed.ctg$"start" , placed.ctg$"pos", placed.ctg$"len"  , placed.ctg$"pos", lwd=10*LG.size , col= placed.ctg$"col",lend=1) 
    text( placed.ctg$"len"+CCmax.max/100 , placed.ctg$"pos", placed.ctg$"name" , cex=1.6, adj=0, col=placed.ctg$"col")   
    segments( ASM4$"ref.start" , ASM4$"pos"    , ASM4$"ref.end"  , ASM4$"pos"    , lwd=16 , col= ASM4$"col",lend=1) 
    segments( ASM4$"q.start"   , ASM4$"ctg.pos", ASM4$"q.end"    , ASM4$"ctg.pos", lwd=20 , col= ASM4$"col",lend=1) 
#    print(ASM4)
#    cat("\n printed ASM4 \n")

    if ( plot_type == "separate")  dev.off()  
    }
    }
    
if ( plot_type %in% c("single","multiple")) {   dev.off() }
 }
 
cat("\n")




cat("\n\n\n")
cat("############################################################################\n")
cat("### 9) Low confidence Contigs                                            ###\n")
cat("############################################################################\n")

cat("\n  --> Confidence Threshold",thr.conf,"\n")

#PLOT_TYPES = c("single", "multiple","separate")

for (plot_type in PLOT_TYPES)
{
if ( plot_type == "single"  ) { cat("\n ==> Plot type:   ",plot_type,"") ; pdf( paste0(out.dir,"/GS-viewer_09.Low_confidence_contigs.single.",OUT,"pdf"  ) ,width=plot.w, height=plot.h) ; par(mfcol=c(1    ,1    ), mar=c(2,8,4,1) , oma=c(0,5,0,0) ) }
if ( plot_type == "multiple") { cat("\n ==> Plot type: "  ,plot_type,"") ; pdf( paste0(out.dir,"/GS-viewer_09.Low_confidence_contigs.multiple.",OUT,"pdf") ,width=plot.w, height=plot.h) ; par(mfcol=c(n.row,n.col), mar=c(2,9,4,1) , oma=c(0,5,0,0) ) }

for( ii in names(CHR))
{
# i =  names(CHR)[3];
# long_i = 
# if ( plot_type == "separate")  { cat("\n ==> Plot type: ",plot_type,"") ; pdf( paste0(out.dir,"/RagTag_scaffold__Low_confidence_contigs.",ii,".",OUT,"pdf") ,width=plot.w, height=plot.h) ; par(mfcol=c(1    ,1    ), mar=c(2,5,4,1), oma=c(0,3,0,0) ) }
# pdf( paste0(out.dir,"/RagTag_scaffold__Low_confidence_contigs.",ii,".",OUT,"pdf"),width=plot.w, height=plot.h )

# for( ii in names(CHR)) {

# ii =  names(CHR)[2];
if ( plot_type %in% c("single","multiple")) {  cat(ii,"")  }

chr = data.frame(CHR[[ii]]) ; 
agp = setDT(AGP[[ii]]) ; 
cnf = setDT(CNF[[ii]]) ; 
asm = setDT(ASM[[ii]]) ; 
chr = chr[ rev(order(chr$"chr")) ,]
chr$"chr" = gsub("_RagTag","",chr$"chr")

CCmax.max = max(CCmax, chr$"len")

# separate GAPs annotation from AGP
agp$"chr" = gsub("_RagTag","",agp$"ref.chr")
gap = agp[ agp$"query" =="100",]
agp = agp[ agp$"query" !="100",]

# create barplot and get coordinates of each chromosome in the plot 
# lll =  0.7+(0:(nrow(chr)-0)*1.2)
# bar=barplot(chr$"len", names.arg=chr$"chr", las=1, horiz=T, main=i, cex.main=2.5 ,xlim=c(0, CCmax*1.20), ylim=c(0,max(lll)*0.98 ), cex.names=2, plot=F)
# b2 = data.frame(chr.pos = bar[,1])
# b2$"chr" = chr$"chr"
# rrr = merge(ref, b2,by="chr", all.x=TRUE) # add reference chr length

# AGP file
agp=  agp[ grep("_RagTag",agp$"ref.chr") ,] # rm GAP lines
agp$"ref.chr" = as.factor( agp$"ref.chr")
gap=  agp[ agp$"query"=="100" ,] # rm GAP lines
agp=  agp[ agp$"query"!="100" ,] # rm GAP lines
agp$"chr" = gsub("_RagTag","",agp$"ref.chr")
agp$"ref.mid" = (agp$"ref.start" + agp$"ref.end")/2  # mid point of the contig on the reference 
agp = merge(agp, cnf, by="query", all.x=TRUE)        # add confidence columns
agp = merge(agp, b2, by="chr", all.x=TRUE)           # add coordinated on the plot
agp = agp[ order(agp$"ref.chr"), ]

low_conf.ctg = agp [ agp$"grouping_confidence" <= 0.5      | agp$"location_confidence" <= 0.5,]
low_conf.ctg = agp [ agp$"grouping_confidence" <= thr.conf | agp$"location_confidence" <= thr.conf,]

if (nrow(low_conf.ctg)>0)
{
    n.ctg = nrow(low_conf.ctg)
    n.nnn = n.ctg - length(COL.LIST)
    if (n.nnn <=0 ) n.nnn=1
    pos = 1.5 
    if( n.ctg>6)  pos = 3
    if( n.ctg>12) pos = 4.5
    if( n.ctg>18) pos = 4.5
    if( n.ctg>24) pos = 6
    
    low_conf.ctg$"pos" = -( 1:nrow(low_conf.ctg)/pos) -0.5
    low_conf.ctg$"start" = 0
    low_conf.ctg$"col" = c(COL.LIST[ 1:length(CHR)], rev(rainbow(n.nnn)))[ 1:n.ctg]
    low_conf.ctg$"name"=""
    w.grp = which(low_conf.ctg$"grouping_confidence"<= thr.conf) ;   if ( length(w.grp)>0)  low_conf.ctg$"name"[w.grp] = paste(low_conf.ctg$"name"[w.grp],"grouping", round(low_conf.ctg$"grouping_confidence"[w.grp],2))
    w.loc = which(low_conf.ctg$"location_confidence"<= thr.conf) ;   if ( length(w.loc)>0)  low_conf.ctg$"name"[w.loc] = paste(low_conf.ctg$"name"[w.loc],"location", round(low_conf.ctg$"location_confidence"[w.loc],2))
    low_conf.ctg$"name" = gsub("groupinglocation","grouping/location",low_conf.ctg$"name")
    
    low_conf.ctg$"name" = paste0(low_conf.ctg$"query"," (", low_conf.ctg$"q.len", "bp on ",low_conf.ctg$"chr",") ",low_conf.ctg$"name")
    low_conf.ctg$"len" = low_conf.ctg$"q.len"
    
   
    col_umpl = names(unplaced.ctg)
    low_conf.ctg = low_conf.ctg[, ..col_umpl ]
    
    # asm2 has not bueen fioltered for conting on same chr or not
    asm2 = lapply(asm2, function(x) { x$"ctg.perc"=round( x$"ctg.len"/ x$"q.len",2); return(x) })
    ASM2 = do.call(rbind,asm2)
    ASM4 = ASM2[ ASM2$"query" %in% low_conf.ctg$"query", ]
    ASM4 = ASM4[ ASM4$"ctg.len">=MIN.WIND, ]
    #ASM4 = ASM4[ ASM4$"ctg.perc">=0.30, ]
    
    ASM4 = merge(ASM4,low_conf.ctg, by="query", all.x=TRUE)
    ASM4$"pos" = NULL
    
    rr2= rrr
    names(rr2) = c("ref.chr","ref.len","pos")
    rr2$"ref.len" = NULL
    
    ASM4 = merge(ASM4,rr2, by="ref.chr", all.x=TRUE)
    
    
   # pdf( paste0(out.dir,"/RagTag_scaffold__Low_confidence_contigs.",ii,".",OUT,"pdf"),width=plot.w, height=plot.h )   
   # par( mar=c(2,8,4,1))
    if ( plot_type == "separate")  { cat("\n ==> Plot type: ",plot_type,"") ; cat(ii,"") ; pdf( paste0(out.dir,"/Low_confidence_contigs.",ii,".",OUT,"pdf") ,width=plot.w, height=plot.h) ; par(mfcol=c(1    ,1    ), mar=c(2,5,4,1), oma=c(0,3,0,0) ) }

    lll =  0.7+(0:(nrow(ref)-0)*1.2)
 #   ll2 = min(low_conf.ctg$"pos")
    ll2 = -5    
    bar=barplot(ref$"ref.len", names.arg=ref$"chr", las=1, horiz=T, main="Low confidence contigs alignment on Scaffolds",col=gray(0.95), cex.main=2, cex.axis=2,xlim=c(0,max(ref$"ref.len")*1.15), cex.names=2, ylim=c(ll2,max(lll)))
    abline( h=-0.3, col="black", lwd=1)
    text( max(ref$"ref.len")/2, max(lll),  ii, cex=3*LG.size )
    
    segments( low_conf.ctg$"start" , low_conf.ctg$"pos", low_conf.ctg$"len"  , low_conf.ctg$"pos", lwd=12 , col= low_conf.ctg$"col",lend=1) 
    text( low_conf.ctg$"len"+CCmax.max/100 , low_conf.ctg$"pos", low_conf.ctg$"name" , cex=1.2, adj=0, col=low_conf.ctg$"col")   
    segments( ASM4$"ref.start" , ASM4$"pos", ASM4$"ref.end"  , ASM4$"pos", lwd=12 , col= ASM4$"col",lend=1 )
    

    if ( plot_type == "separate")  dev.off()  
    }
    }
    
if ( plot_type %in% c("single","multiple")) {   dev.off() }
 }
 
cat("\n\n")

colors2$"PERC"=  colors2$"perc"/max(colors2$"perc")
colors2$"GAPS"= (colors2$"gaps"/max(colors2$"gaps"))^-1
colors2$"perc.gaps" = colors2$"PERC" * colors2$"GAPS"

# 
# print(colors)
# print(colors2)
# print(colors3)
# print( colors2[ which.max(colors2$"perc.gaps"),]$"assembly")

# if ( bionano == TRUE) 
# {
# 
# 
# system ( paste0( "mkdir ", out.dir,"/Bionano_fasta"))
# 
# cat("\n\n\n")
# cat("############################################################################\n")
# cat("### 10a) BIONANO SUPERSCAFFOLD                                             ###\n")
# cat("############################################################################\n")
# 
# 
# #PLOT_TYPES = c("single", "multiple","separate")
# 
# for (plot_type in PLOT_TYPES)
# {
# if ( plot_type == "single"  ) { cat("\n ==> Plot type:   ",plot_type,"") ; pdf( paste0(out.dir,"/GS-viewer_10a.Bionano_superscaffold.single.",OUT,"pdf"  ) ,width=plot.w, height=plot.h) ; par(mfcol=c(1    ,1    ), mar=c(2,8,4,1) , oma=c(1,5,1,1) ) }
# if ( plot_type == "multiple") { cat("\n ==> Plot type: "  ,plot_type,"") ; pdf( paste0(out.dir,"/GS-viewer_10a.Bionano_superscaffold.multiple.",OUT,"pdf") ,width=plot.w, height=plot.h) ; par(mfcol=c(n.row,n.col), mar=c(2,9,4,1) , oma=c(1,5,1,1) ) }
# 
# for( ii in names(CHR))
# {
# 
# # create barplot and get coordinates of each chromosome in the plot 
# lll =  0.7+(0:(nrow(chr)-0)*1.2)
# bar=barplot(chr$"len", names.arg=chr$"chr", las=1, horiz=T, main=i, cex.main=2.5 ,xlim=c(0, CCmax*1.20), ylim=c(0,max(lll)*0.98 ), cex.names=2, plot=F)
# b2 = data.frame(chr.pos = bar[,1])
# b2$"chr" = chr$"chr"
# rrr = merge(ref, b2,by="chr", all.x=TRUE) # add reference chr length
# 
# 
# # ii =  names(CHR)[2];
# if ( plot_type %in% c("single","multiple")) {  cat(ii,"")  }
# 
# # ii =  names(CHR)[1];
# chr = data.frame(CHR[[ii]]) ; 
# ctg = data.frame(CTG[[ii]]) ; 
# agp = setDT(AGP[[ii]]) ; 
# cnf = setDT(CNF[[ii]]) ; 
# asm = setDT(ASM[[ii]]) ; 
# bio = setDT(BIO[[ii]]) ; 
# 
# chr = chr[ rev(order(chr$"chr")) ,]
# chr$"chr" = gsub("_RagTag","",chr$"chr")
# 
# CCmax.max = max(CCmax, chr$"len")
# 
# # separate GAPs annotation from AGP
# agp$"chr" = gsub("_RagTag","",agp$"ref.chr")
# gap = agp[ agp$"query" =="100",]
# agp = agp[ agp$"query" !="100",]
# agp$"ctg" = agp$"query"
# agp$"chr.start" = agp$"ref.start"
# agp$"chr.end" = agp$"ref.end"
# 
# agp2.names = c("chr","chr.start","chr.end","ctg")
# agp2 = agp [, ..agp2.names]
# 
# 
# # add data to bionano
# names(ctg) = c("ctg","length")
# bio = bio[ bio$"q.start" !="scaffold",]
# bio$"sup.start" =  as.numeric(bio$"sup.start")
# bio$"sup.end"   =  as.numeric(bio$"sup.end")
# 
# bio$"ctg" = do.call(rbind,strsplit( as.character(bio$"query"),"_subseq_"))[,1]
# bio$"ctg.frag" = as.numeric(bio$"q.end") - as.numeric(bio$"q.start") +1 
# bio = merge(bio, ctg, by="ctg", all.x=TRUE)
# bio = merge(bio, agp2, by="ctg", all.x=TRUE)
# bio$"ctg.perc" = 100*bio$"ctg.frag"/bio$"length"
# 
# bio = bio [ order(bio$"sup.start") , ]
# 
# # SPLIT DATA
# # a)  1 ctgs  -->  1 superscaffold
# # b) >1 ctgs  --> >1 superscaffold
# 
# bio_single = bio[ bio$"ctg" == gsub("_obj","",bio$"superscaffold" ),]
# bio_multi  = bio[ bio$"ctg" != gsub("_obj","",bio$"superscaffold" ),]
# 
# # b) >1 ctgs  --> >1 superscaffold
# # b)   --> split by CHR
# # b)   --> split by SCAFFOLD
# 
# bio.all.sup = split(bio, bio$"superscaffold")
# bios_chr = split(bio_multi, bio_multi$"chr")
# bios_sup = split(bio_multi, bio_multi$"superscaffold")
# agp__chr = split(agp      ,       agp$"chr")
# gap__chr = split(gap      ,       gap$"chr")
# 
# bios_sup.same___chr = which( sapply(bios_sup, function(x) { all( nrow(x)> 1 & x$"chr" %in% ref$"chr" & length(unique(x$"chr"))==1)   }   ))
# bios_sup.single_chr = which( sapply(bios_sup, function(x) { all( nrow(x)==1 & x$"chr" %in% ref$"chr")   }   ))
# bios_sup.single_ctg = which( sapply(bios_sup, function(x) { all( nrow(x)==1 & x$"chr" %in% ctg$"ctg")   }   ))
# 
# # b)   --> split by SCAFFOLD
# # b)         --> 1 SUPERSCAFFOLD vs  1 CTG
# # b)         --> 1 SUPERSCAFFOLD vs >1 CTG
# bios_sup.singl = bios_sup[  c( bios_sup.single_chr,bios_sup.single_ctg,bios_sup.same___chr) ]
# bios_sup.multi = bios_sup[ -c( bios_sup.single_chr,bios_sup.single_ctg,bios_sup.same___chr) ]
# 
# 
# # b)   --> split by CHR/ctg
# # b)         --> 1 CHR/ctg vs  1 SUPERSCAFFOLD
# # b)         --> 1 CHR/ctg vs >1 SUPERSCAFFOLD
# bios_chr.single_sca = which( sapply(bios_chr, function(x) { length(unique(x$"chr"))==1 & length(unique(x$"superscaffold"))==1  }   ))
# 
# bios_chr.singl = bios_chr[  c( bios_chr.single_sca) ]
# bios_chr.multi = bios_chr[ -c( bios_chr.single_sca) ]
# 
# bios_chr.multi.names = names(bios_chr.multi)
# agp__chr [ bios_chr.multi.names ]
# 
# 
# print(bios_sup.multi)
# print(bios_chr.multi)
# print(agp__chr [ bios_chr.multi.names ])
# 
# 
# # SINGLE DATA
# 
# bios_sup.SINGL = do.call(rbind,bios_sup.singl)
# bios_chr.SINGL = do.call(rbind,bios_chr.singl)
# 
# length(unique(bios_sup.SINGL$"chr"))
# length(unique(bios_sup.SINGL$"ctg"))
# length(unique(bios_sup.SINGL$"superscaffold"))
# length(unique(bios_chr.SINGL$"chr"))
# length(unique(bios_chr.SINGL$"ctg"))
# length(unique(bios_chr.SINGL$"superscaffold"))
# 
# 
# ##  filter REF.chr only
# 
# bios_chr.singl.ref = bios_chr.singl[  names(bios_chr.singl) %in% ref$"chr" ]
# bios_chr.multi.ref = bios_chr.multi[  names(bios_chr.multi) %in% ref$"chr" ]
# 
# 
# # Add scaffold completeners to CHR 
# 
# for ( jj in names(bios_chr.singl.ref) ) { a=bios_chr.singl.ref [[ jj ]] ; for( sca in 1:nrow(a) ) { SCA=a$"superscaffold"[sca] ; bios_chr.singl.ref [[ jj ]]$"superscaffold.len"[sca] = max( bios_sup [[ SCA ]] $"sup.end") } }
# for ( jj in names(bios_chr.multi.ref) ) { a=bios_chr.multi.ref [[ jj ]] ; for( sca in 1:nrow(a) ) { SCA=a$"superscaffold"[sca] ; bios_chr.multi.ref [[ jj ]]$"superscaffold.len"[sca] = max( bios_sup [[ SCA ]] $"sup.end") } }
# 
# 
# if ( plot_type == "separate")  { cat("\n ==> Plot type: ",plot_type,"") ; cat(ii,"") ;  pdf( paste0(out.dir,"/Bionano_vs_ragtag.",ii,".",OUT,"pdf"),width=plot.w, height=plot.h ); par(mfcol=c(1    ,1    ), mar=c(2,6,4,1), oma=c(1,3,1,1) ) }
# #   pdf( paste0(out.dir,"/Bionano_vs_ragtag.",ii,".",OUT,"pdf"),width=plot.w, height=plot.h )
# #   par( mar=c(2,9,4,1))
# 
# 
# col.ctg1 = c("skyblue","skyblue4")
# col.ctg2 = c("lightpink","lightpink4")
# 
# lll =  0.7+(0:(nrow(ref)-0)*1.2)
# ll2 = -14
# 
# max.plot = max(ref$"ref.len")*1.20
# plot.shift = max.plot/500
# 
# 
# bar=barplot(ref$"ref.len", names.arg=ref$"chr", las=1, horiz=T, main="Ragtag Chromosomes vs Bionano Superscaffolds",col=gray(0.95), cex.main=3.0,cex.axis=2,  xlim=c(0, max(ref$"ref.len")*1.15), cex.names=2, ylim=c(ll2,max(lll)*1.025))
# b2 = data.frame(chr.pos = bar[,1])
# b2$"chr" = ref$"chr"
# rrr = merge(ref, b2,by="chr", all.x=TRUE) # add reference chr length
# 
# # abline( h=-0.3, col="black", lwd=3)
# text( max(ref$"ref.len")*1.15/2, max(lll)*1.05,  ii, cex=2 )
# 
# ploted.scaffolds = c( )
# perc.plot.scaffolds = c( )
# 
# #  1 Superscaffold == 1 Contig
# for ( jj in names(bios_chr.singl.ref) ) 1==1
# for ( jj in names(bios_chr.singl.ref) )
# {
#   print(jj)
#   a = bios_chr.singl.ref [[ jj ]]
#   a = data.frame(a)
#   a$"sup.len" = a$"sup.end" - a$"sup.start"+1
#   a$"superscaffold.perc" = NA
#   for ( sca in  unique(  a$"superscaffold" )) {  b=a[ a$"superscaffold" == sca, ]; a[ a$"superscaffold" == sca, ]$"superscaffold.perc" = round( 100*sum(b$"sup.len")/b$"superscaffold.len"[1], 0) }
#   a$"SupSca.stat" = paste0(a$"superscaffold", "  (",a$"superscaffold.perc","%)")
#   
# 
#   a.pos = rrr[ rrr$"chr" == jj, ]$"chr.pos"
#   a.len = rrr[ rrr$"chr" == jj, ]$"ref.len"
# 
#   #  print(a.pos)
#   #  print(unique(a$ctg))
# 
#    gg  = gap__chr[[jj]]
#    if( is.null(gg)==FALSE) ggg = merge(gg, b2,by="chr", all.x=TRUE) # add reference chr length
# 
#   if (nrow(a)==1 ) 
#   {     polygon( c(0,a.len,a.len,0), c(a.pos-0.50,a.pos-0.50,a.pos+0.50,a.pos+0.50), lwd=3 , col="lightblue") #, border=FALSE) 
#         text( a.len + plot.shift, a.pos+0.20 , a$"SupSca.stat", cex=1, adj=0, col="black")   
#         text( a.len + plot.shift, a.pos-0.20 , a$"ctg", cex=1, adj=0, col="gray")  
#       if( is.null(gg)==FALSE) points( ggg$"ref.start" , ggg$"chr.pos" , pch=23, cex=1.4, adj=0, col="black", bg="yellow", lwd=2)   
#       }
# 
#   if (nrow(a)>1 ) 
#   { 
#       a$"query.coord"=                         do.call(rbind,strsplit( as.character(a$"query")      ,"_subseq_"))[,2]
#       a$"query.start"= as.numeric(as.character(do.call(rbind,strsplit( as.character(a$"query.coord"),":"))       [,1]))
#       a$"query.end"  = as.numeric(as.character(do.call(rbind,strsplit( as.character(a$"query.coord"),":"))       [,2]))
#       a = a[ order(a$"query.start") , ]
#       
#       a.qs = a$"query.start"[-1]      -1 
#       a.qe = a$"query.end"  [-nrow(a)] 
#       a.nn = length(a.qs)
# 
#       if( all(a.qs == a.qe) & length(unique(a$"ctg") ==1) )  
#       {
#         polygon( c(0,a.len,a.len,0), c(a.pos-0.50,a.pos-0.50,a.pos+0.50,a.pos+0.50), lwd=3 , col="lightblue") #, border=FALSE) 
#         text( a.len + plot.shift, a.pos+0.20 , a$"SupSca.stat"[1], cex=1, adj=0, col="black")   
#         text( a.len + plot.shift, a.pos-0.20 , a$"ctg"[1], cex=1, adj=0, col="gray")   
#        # for (A.QS in a.qs) segments(A.QS, a.pos-0.3, A.QS  , a.pos+0.3, lwd=3 , col="blue",lend=1) 
#         for (A.QS in a.qs)    text(A.QS, a.pos, "x", cex=2.5, col="blue",lend=1) 
#       if( is.null(gg)==FALSE) points( ggg$"ref.start" , ggg$"chr.pos" , pch=23, cex=1.4, adj=0, col="black", bg="yellow", lwd=2)   
#         }
# 
# }
# ploted.scaffolds = c(ploted.scaffolds, unique(a$"superscaffold"))
# perc.plot.scaffolds = c(perc.plot.scaffolds, unique(a$"SupSca.stat"))
# }
# 
# remaining_chr = setdiff( ref$"chr" , names(bios_chr.singl.ref) )
#  
# for ( jj in remaining_chr ) 1==1
# for ( jj in names(bios_chr.multi.ref) )
# {
#     # jj = names(bios_chr.multi.ref)[1]
#     print(jj)
#     a = bios_chr.multi.ref [[ jj ]]
#     a = data.frame(a)
#     a.pos = rrr[ rrr$"chr" == jj, ]$"chr.pos"
#     a.len = rrr[ rrr$"chr" == jj, ]$"ref.len"
#     
#     dd = chr
#     
#     gg =gap__chr[[jj]]
#     if( is.null(gg)==FALSE) ggg = merge(gg, b2,by="chr", all.x=TRUE) # add reference chr length
# 
#   a$"sup.len" = a$"sup.end" - a$"sup.start"+1
#   a$"superscaffold.perc" = NA
#   for ( sca in  unique(  a$"superscaffold" )) {  b=a[ a$"superscaffold" == sca, ]; a[ a$"superscaffold" == sca, ]$"superscaffold.perc" = round( 100*sum(b$"sup.len")/b$"superscaffold.len"[1], 1) }
#  
# 
#     a$"SupScaff"=gsub("Super-Scaffold_","", a$"superscaffold")
#     a$"SupScaff"= sapply( strsplit( as.character(a$"SupScaff")   ,"_") , function(x) { y=x[1]; if ( any(x=="subseq" )) { y= paste( head(x,1),"subseq", tail(x,1), sep="_") } ; return(y) }  )
#     a$"query.coord"=                         do.call(rbind,strsplit( as.character(a$"query")      ,"_subseq_"))[,2]
#     a$"query.start"= as.numeric(as.character(do.call(rbind,strsplit( as.character(a$"query.coord"),":"))       [,1]))
#     a$"query.end"  = as.numeric(as.character(do.call(rbind,strsplit( as.character(a$"query.coord"),":"))       [,2]))
#     a = a[ order(a$"chr.start") , ]
#     
#     a.qs = a$"query.start"[-1]      -1 
#     a.qe = a$"query.end"  [-nrow(a)] 
#     a.nn = length(a.qs)
#     print(a)
#     print(a.nn)
#     print(a.qs)
#    
#     n.ctg = length(unique(a$"ctg"))
#     n.sup = length(unique(a$"SupScaff"))
#     a$"SupSca.stat"  = paste0(a$"SupScaff", "  (",a$"superscaffold.perc","%)")
#     a$"SupSca.stat2" = paste0(a$"superscaffold", "  (",a$"superscaffold.perc","%)")
#     a$"chr.pieces.len" = a$"chr.end" - a$"chr.start"  +1
#     a$"sup.pieces.len" = a$"sup.end" - a$"sup.start"  +1
#     a$"pieces.perc" = abs(log10(a$"chr.pieces.len" / a$"sup.pieces.len"))
#    
#     
#     SUP_lab = paste( "Super-Scaffolds :", paste(unique(a$"SupSca.stat"), collapse=" , ") )
#     #CTG_lab = paste( "Contigs:", paste(unique(a$"ctg"), collapse=" , ") )
#     CTG_lab =                     paste(unique(a$"ctg"), collapse=" , ") 
#     
#     a$"col" = col.ctg1[1]
#     
#     as = split(a, a$"superscaffold")
#     
#     even.pos = seq(2,nrow(a),by=2)
#     
#     as[ even.pos ] = lapply( as[ even.pos ] , function(x) {x$"col" = col.ctg1[2] ; return(x) } )
# 
#   
#  if (n.ctg == 1) 
#   {  
#      #   polygon( c(0,a.len,a.len,0), c(a.pos-0.50,a.pos-0.50,a.pos+0.50,a.pos+0.50), lwd=3 , col="lightyellow") #, border=FALSE) 
#        # for (A.QS in a.qs) segments(A.QS, a.pos-0.5, A.QS  , a.pos+0.5, lwd=3 , col="black",lend=1) 
#         for (ai in names(as)) { b=as[[ai]]; polygon( c(b$"query.start",b$"query.end",b$"query.end",b$"query.start"), c(a.pos-0.50,a.pos-0.50,a.pos+0.50,a.pos+0.50), lwd=3 , col=b$"col",border = NA) }
#         text( a.len + plot.shift, a.pos+0.20 ,SUP_lab, cex=1, adj=0, col="blue")   
#      #  text( a.len + plot.shift, a.pos+0.20 ,"Super-Scaffolds :", cex=1, adj=0, col="black")   
#         text( a.len + plot.shift, a.pos-0.20 , CTG_lab, cex=1, adj=0, col="gray")  
#         for (A.QS in a.qs)     text(A.QS, a.pos, "x", cex=2.5, col="black",lend=1) 
#       #  for (A.QS in a.qs) segments(A.QS, a.pos-0.4, A.QS  , a.pos+0.4, lwd=3 , col="black",lend=1) 
#        if( is.null(gg)==FALSE) points( ggg$"ref.start" , ggg$"chr.pos" , pch=23, cex=1.4, adj=0, col="black", bg="yellow", lwd=2)   
#   }
# polygon( c(0,a.len,a.len,0), c(a.pos-0.50,a.pos-0.50,a.pos+0.50,a.pos+0.50), lwd=3 , col=NA)
#   
#   CTG_lab = paste( "Contigs:", paste(unique(a$"ctg"), collapse=" , ") )
#   if ( n.ctg >1 ) 
#   { 
#   # overlapping?
#   ovrlp.test = sum(a$"chr.pieces.len" / max(a$"chr.end"))
#   if( ovrlp.test  >95 &  ovrlp.test <105 )
#   {
#   #no overalp -->ok
#         polygon( c(0,a.len,a.len,0), c(a.pos-0.50,a.pos-0.50,a.pos+0.50,a.pos+0.50), lwd=3 , col="lightpink") #, border=FALSE) 
#         text( a.len + plot.shift, a.pos+0.20 ,SUP_lab, cex=1, adj=0, col="darkred")   
#         text( a.len + plot.shift, a.pos+0.20 ,"Super-Scaffolds :", cex=1, adj=0, col="black")   
#         text( a.len + plot.shift, a.pos-0.20 , CTG_lab, cex=1, adj=0, col="gray")  
#      #   for (A.QS in a.qs) segments(A.QS, a.pos-0.5, A.QS  , a.pos+0.5, lwd=3 , col="black",lend=1) 
#         for (A.QS in a.qs)     text(A.QS, a.pos, "+", cex=2.4, col="black",lend=1) 
#         if( is.null(gg)==FALSE) points( ggg$"ref.start" , ggg$"chr.pos" , pch=23, cex=1.4, adj=0, col="black", bg="yellow", lwd=2)   
#         } else {
#             # overlap of different bionano --> need to separate in 2.
#             # identify main
#             a$"chr.coord"     = paste(a$"chr.start" , a$"chr.end" )
#             a$"chr.coord.grp" = NA
#             kkk=0
#             for ( grp in unique(a$"chr.coord")) { kkk=kkk+1 ;  a[ a$"chr.coord" == grp,]$"chr.coord.grp" = kkk }
#             
#             as = split(a, a$"chr.coord.grp")
#             
#             # any single scaffold?
#             best = as [[  which(sapply(as,nrow)==1) ]]$"superscaffold"
#             if( length(best)>0) best = unique(grep("Super-Scaffold",a$"superscaffold", value=TRUE))
#             if( length(best)>1) best = a[ which( grepl("Super-Scaffold",a$"superscaffold") & a$"superscaffold.perc"==max( a$"superscaffold.perc")), ]$"superscaffold"
#             if( length(best)>0)
#             {        
#                 a.best = a[ a$"superscaffold" == best,]
#                 missing.best = setdiff( names(as) , a.best$"chr.coord.grp" )
#                 
#                 a.miss = as[ missing.best ]       
#             
#                 for ( kkk in names(a.miss))
#                 {
#                     am = a.miss[[kkk]]
#                     if( any( grepl("Super-Scaffold",am$"superscaffold") )) am = am[ which( grepl("Super-Scaffold",am$"superscaffold")), ]
#                     if( nrow(am)>1)
#                     {
#                         am = am [ am$"superscaffold.perc" == max(am$"superscaffold.perc"), ]
#                         am = am [ am$"pieces.perc"        == min(am$"pieces.perc"), ]
#                         a.miss[[kkk]] = head(am,1)
#                     }
#                     }
#                 a.new = rbind(a.best,do.call(rbind,a.miss))
#                 a.new = a.new[ order(a.new$"chr.start") , ]
#                 a=a.new
#                 a.qs = a$"chr.start"[-1]      -1 
#                 a.qe = a$"chr.end"  [-nrow(a)] 
#                 SUP_lab = paste( "Super-Scaffolds :", paste(unique(a$"SupSca.stat"), collapse=" , ") )
#                 CTG_lab = paste( "Contigs:", paste(unique(a$"ctg"), collapse=" , ") )
# 
#                 a$"col" = col.ctg2[1]
#                 a$"pch.break.col" = NA
#                 as = split(a, a$"superscaffold")                
#                 even.pos = seq(2,nrow(a),by=2)               
#                 as[ even.pos ] = lapply( as[ even.pos ] , function(x) {x$"col" = col.ctg2[2] ; return(x) } )
#                 
#                 as = do.call(rbind,as)
#                 as = as[ order(as$"chr.start") , ]
#                 if( length(unique(as$superscaffold)) >1 ) for ( ai in 2:nrow(as) ) { aj=ai-1 ; if( as$"superscaffold"[ai]==as$"superscaffold"[aj]) {as$"pch.break.col"[ai]="blue" } else  as$"pch.break.col"[ai]="black"  }
#                 as1 = split(as, as$"ctg")                
#                 as2 = split(as, as$"superscaffold")                
# 
#                for (ai in names(as1)) { b=as1[[ai]]; polygon( c(b$"chr.start",b$"chr.end",b$"chr.end",b$"chr.start"), c(a.pos-0.50,a.pos-0.50,a.pos+0.50,a.pos+0.50), lwd=3 , col=b$"col", border=NA) }
#            #   for (A.QS in a.qs)     text(A.QS, a.pos, "x", cex=2.5, col="black",lend=1) 
#            #   for (A.QS in a.qs) segments(A.QS, a.pos-0.4, A.QS  , a.pos+0.4, lwd=3 , col="black",lend=1) 
#                text(as$"chr.start"[-1], rep(a.pos, nrow(as)-1), "x", cex=2.5, col=as$"pch.break.col"[-1],lend=1) 
#                text( a.len + plot.shift, a.pos+0.20 ,SUP_lab, cex=1, adj=0, col="darkred")   
#            #   text( a.len + plot.shift, a.pos+0.20 ,"Super-Scaffolds :", cex=1, adj=0, col="black")   
#                text( a.len + plot.shift, a.pos-0.20 , CTG_lab, cex=1, adj=0, col="gray")                    
#                if( is.null(gg)==FALSE) points( ggg$"ref.start" , ggg$"chr.pos" , pch=23, cex=1.5, adj=0, col="black", bg="yellow", lwd=2)  
#                # break points same scaffold
# polygon( c(0,a.len,a.len,0), c(a.pos-0.50,a.pos-0.50,a.pos+0.50,a.pos+0.50), lwd=3 , col=NA)
#                
#             }
#   
#         }
# 
# }
# ploted.scaffolds = c(ploted.scaffolds, unique(a$"superscaffold"))
# perc.plot.scaffolds = c(perc.plot.scaffolds, unique(a$"SupSca.stat2"))
# 
# }
# 
# ###  extra scaffold 
# 
# other.scaffolds = grep( "Super-Scaffold", setdiff( names(bios_sup), ploted.scaffolds), value=TRUE)
# incomplete.scaffolds = grep( "Super-Scaffold", perc.plot.scaffolds, value=TRUE)
# incomplete.scaffolds = data.frame( do.call(rbind,strsplit( gsub('%)',"",incomplete.scaffolds),"  [(]")))
# incomplete.scaffolds$"X1" =            as.character(incomplete.scaffolds$"X1")
# incomplete.scaffolds$"X2" = as.numeric(as.character(incomplete.scaffolds$"X2"))
# incomplete.scaffolds = incomplete.scaffolds[ incomplete.scaffolds$"X2" < 98, ]$"X1"
# 
# 
# OTH.SCAF = bios_sup [ other.scaffolds ]
# INC.SCAF = bios_sup [ incomplete.scaffolds ]
# 
# inc.scaf = data.frame( scaffold="broken Scaffold", name=names(INC.SCAF), stringsAsFactors=FALSE)
# oth.scaf = data.frame( scaffold="other  Scaffold" , name=names(OTH.SCAF), stringsAsFactors=FALSE)
# 
# new.scaf = rbind(inc.scaf,oth.scaf)
# NEW.SCAF = c(INC.SCAF,OTH.SCAF)
# 
# 
# ##  AAAAa remove the ctgs from the OTH.SCAF and INC.SCAF
# 
# ctg.scaffolfd = setdiff( unique(bio$superscaffold),  c(ploted.scaffolds,other.scaffolds,incomplete.scaffolds))
# CTG.SCAF = bio.all.sup [ ctg.scaffolfd ]
# ctg.scaf = data.frame( scaffold="unplaced contig", name=names(CTG.SCAF), "scaff.len"=sapply(CTG.SCAF, function(x) max(x$"sup.end")), stringsAsFactors=FALSE)
# ctg.scaf = ctg.scaf [ rev(order(ctg.scaf$"scaff.len")) , ]
# ctg.scaf$"ctg" = as.character(do.call(rbind,strsplit(ctg.scaf$"name","_"))[,1])
# ctg.scaf$"bionano" = "no"
# rownames(ctg.scaf) = NULL
# 
# # mark also the _subseq: CTGS
# if( nrow(new.scaf)>0) for( jj in new.scaf$"name") {   jj.ctg = NEW.SCAF [[jj]]$"ctg" ; if( any(ctg.scaf$"ctg" %in% jj.ctg) ) {ctg.scaf [ ctg.scaf$"ctg" %in% jj.ctg,]$"bionano" = jj} }
# 
# # add info to the subseq
# 
# ctg.temp = do.call(rbind, CTG.SCAF[ ctg.scaf[ctg.scaf$"bionano" != "no", ]$"name"] )
# ctg.temp = data.frame(ctg.temp)
# ctg.temp = ctg.temp[ , c("superscaffold","ctg.frag") ]
# ctg.temp = merge(ctg.scaf,ctg.temp, by.y="superscaffold", by.x="name", all.x=TRUE)
# ctg.temp = merge(ctg.temp,ctg, by="ctg", all.x=TRUE)
# ctg.temp$"ctg.perc" = 100*ctg.temp$"ctg.frag"/ctg.temp$"length"
# 
# ctg.scaf= ctg.temp
# # keep only unplaced by ragtag
# ctg.scaf = ctg.scaf [ ctg.scaf$"ctg" %in% ctg$"ctg", ]
# #ctg.scaf = merge(ctg.scaf,ctg, by="ctg", all.x=TRUE)
# ctg.scaf = ctg.scaf [ rev(order(ctg.scaf$"length")) , ]
# 
# 
# n.ctg_to_add = 10-nrow(new.scaf)
# 
# 
# new.scaf$"length"=sapply(NEW.SCAF, function(x) max(x$"sup.end"))
# new.scaf$"pos" = - lll [ 2: (nrow(new.scaf)+1) ] + 0.2
# ctg.scaf$"pos" = NA
# 
# for ( jj in setdiff(names(ctg.scaf), names(new.scaf)) ) new.scaf[[jj]]=NA
# new.scaf = rbind(new.scaf, head(ctg.scaf, n.ctg_to_add ) )
# NEW.SCAF = c(NEW.SCAF, CTG.SCAF [ head(ctg.scaf$"name", n.ctg_to_add)] )
# 
# 
# 
# new.scaf$"pos" = - lll [ 2: (nrow(new.scaf)+1) ] + 0.1
# w.unpl = which(new.scaf$"scaffold"=="unplaced contig")
# if (length(w.unpl)>0) new.scaf$"pos"[w.unpl] =  new.scaf$"pos"[w.unpl] -1.8
# 
# for ( jj in names(NEW.SCAF))
# { 
#   # jj= names(NEW.SCAF)[1]
#     a = NEW.SCAF[[jj]] 
#     b = new.scaf[ new.scaf$"name" == jj, ]
#     
#     a.pos  = new.scaf[ new.scaf$"name" == jj, ]$"pos"
#     a.type = new.scaf[ new.scaf$"name" == jj, ]$"scaffold"
#     
#     if( a.type == "broken Scaffold")  { a.col ="lightpink"  ;  a.col2 ="darkred" ;  lwd.polyg = 3}
#     if( a.type == "other  Scaffold")  { a.col ="gold"       ;  a.col2 ="brown"   ;  lwd.polyg = 1}
#     if( a.type == "unplaced contig")  { a.col ="cornsilk"   ;  a.col2 ="black"   ;  lwd.polyg = 1 }
#     
#     a.len = b$"length"
#     a.qs = a$"sup.start"[-1]      -1 
#     
#     SUP_lab = paste( "Super-Scaffolds :", paste(unique(a$"superscaffold"), collapse=" , ") )
#     CTG_lab = paste(          "Contigs:", paste(unique(a$"ctg"), collapse=" , ") )
#     
#     
#     polygon( c(0,a.len,a.len,0), c(a.pos-0.50,a.pos-0.50,a.pos+0.50,a.pos+0.50), lwd=1 , col= a.col) #, border=FALSE) 
#     for (A.QS in a.qs) segments(A.QS, a.pos-0.4, A.QS  , a.pos+0.4, lwd=lwd.polyg , col="black",lend=1) 
# 
#     if( a.type == "broken Scaffold")   text( a.len + plot.shift, a.pos+0.20 ,SUP_lab, cex=1, adj=0, col= a.col2 )   
#     if( a.type == "other  Scaffold")   text( a.len + plot.shift, a.pos+0.20 ,SUP_lab, cex=1, adj=0, col= a.col2 )   
#     if( a.type == "unplaced contig")   text( a.len + plot.shift, a.pos+0.20 ,CTG_lab, cex=1, adj=0, col= a.col2 )   
#     if( a.type == "broken Scaffold")   text( a.len + plot.shift, a.pos-0.20 , CTG_lab, cex=1, adj=0, col="gray")  
#     if( a.type == "other  Scaffold")   text( a.len + plot.shift, a.pos-0.20 , CTG_lab, cex=1, adj=0, col="gray")  
#              
#     if( is.na(b$"bionano")==FALSE ) if( b$"bionano" != "no")   
#     {          
#         sss = bio[ bio$"ctg" == b$"ctg" , ]
#         sss = sss[ sss$"superscaffold" != b$"name" , ]
#         ssss = data.frame(sss)
#         
#         sss$"subseq" = do.call(rbind,strsplit(sss$"query" ,"_subseq_"))[,2]
#         sss$"subseq.start" = as.numeric(as.character(do.call(rbind,strsplit(sss$"subseq" ,":"))[,1]))
#         sss$"subseq.end"   = as.numeric(as.character(do.call(rbind,strsplit(sss$"subseq" ,":"))[,2]))
#         b.start = sss$"subseq.start" #  invert !coordinates smaller box 
#         b.end = sss$"subseq.end"      
# 
#         polygon( c(b.start,b.end,b.end,b.start), c(a.pos-0.50,a.pos-0.50,a.pos+0.50,a.pos+0.50), lwd=3 , col= "lightpink", border=NA) #, border=FALSE) 
#         polygon( c(0,a.len,a.len,0), c(a.pos-0.50,a.pos-0.50,a.pos+0.50,a.pos+0.50), lwd=1 , col= NA) #, border=FALSE) 
#         CTG_lab =  paste(CTG_lab, "   --> pink region is on", b$bionano )
#         if( a.type == "unplaced contig")   text( a.len + plot.shift, a.pos+0.20 ,CTG_lab, cex=1, adj=0, col="darkred"  )   
#         # fix 
#         # add 
#     }
# 
#     
#     }
# 
# }
# 
# legend("right", c("1 Chr =  1 Superscaffold", "1 Chr =  >1 Superscaffolds","1 Chr =   incomplete Superscaffold","Bionano break points","Gaps"), col=c("black","black","black","blue","black"),pt.bg=c("lightblue",col.ctg1[1],"lightpink","blue","yellow"), pch=c(22,22,22,4,23), cex=1.3, bg="white")
# legend("right", c("1 Chr =  1 Superscaffold", "1 Chr =  >1 Superscaffolds","1 Chr =  incomplete Superscaffold(s)","*same color --> same Superscaffold","","Bionano break points (same Superscaffold)","Bionano break points (diff. Superscaffolds)","Gaps in RagTag" ), col=c("black","black","black","white","white","blue","black","black"),pt.bg=c("lightblue",col.ctg1[2],"lightpink","white","white","blue","black","yellow"), pch=c(22,22,22,22,4,4,4,23), cex=2, bg="white",  pt.cex=c(3.5,3.5,3.5,2.5,2.5,1.8,1.8,2 ), pt.lwd=c(1,1,1,1,1,4,4,1))
# 
# 
# 
# 
# legend("bottomright", c("Complete Superscaffolds with partial match with Ragtag Chr","Superscaffolds absent in Ragtag Chr","unplaced CTGs on reference (RagTag) "), pch=c(22,22,22), col=c("black","black","black"), pt.bg=c("lightpink","gold","cornsilk"), cex=2, bg="white",pt.cex=3.5)
# 
# unpl_ctg_text.pos = new.scaf$"pos"[w.unpl][1] +1
# supersca_text.pos = new.scaf$"pos"[1] +1
# 
# text( 0 , max(bar)+1        , "RagTag chromosomes + Bionano superscaffold alignment", cex=1.8, adj=0 )  
# text( 0 , supersca_text.pos , "Bionano Superscaffold that have partial or not match with RagTag chromosomes", cex=1.8, adj=0 )  
# text( 0 , unpl_ctg_text.pos , "RagTag unplaced contigs on reference ", cex=1.8, adj=0 )  
# 
# 
# 
# extra.scaf = rbind(inc.scaf,oth.scaf)
# assembly.file     = LOG[ LOG$"short_name" ==i,]$"assembly_file"
# 
# cat("Creating fasta files with Bionano Superscaffold that have partial or not match with RagTag chromosomes: \n\n\n")    
# 
# for ( jj in extra.scaf$"name")
# {
#     # jj = extra.scaf$"name"[1]
#     cat(jj)
#     a = NEW.SCAF[[jj]] 
#     a.ctgs = a$"ctg"
#     jj.ctg.list  = paste0(out.dir,"/Bionano_fasta/",jj,".ctg.list"  ) 
#     jj.ctg.fasta = paste0(out.dir,"/Bionano_fasta/",jj,".ctg.fasta"  ) 
#     jj.ctg.fai   = paste0(out.dir,"/Bionano_fasta/",jj,".ctg.fasta.fai"  ) 
#     jj.ctg.nnnnn = paste0(out.dir,"/Bionano_fasta/",jj,".NNN.txt"  ) 
#     jj.ctg.command = paste0(out.dir,"/Bionano_fasta/",jj,".command_merge.sh"  )  ;# add "#!/bin/bash\n"
#     jj.sca.fasta = paste0(out.dir,"/Bionano_fasta/",jj,".scaffold.fasta"  ) 
#     
#     
#     a.nnnn = rep(  paste(  rep("N",100), collapse=""), length(a.ctgs)-1 )
#     command = paste0( " cat <( echo '>",jj,"') <(paste <(grep -v '^>' ",jj.ctg.fasta,") <(cat ",jj.ctg.nnnnn,") -d '' | sed  ':a;N;$!ba;s/WW//1g' ) > ", jj.sca.fasta ) # cannot put \n because lionuex will generate new line, add WW and sed later
#     command = c("#!/bin/bash\n", paste( "cd",  getwd() ), command)
#     command
#     write.table(        a.ctgs   , file = jj.ctg.list     , col.names=FALSE, row.names=F, quote=F, sep="\t")
#     write.table(        a.nnnn   , file = jj.ctg.nnnnn    , col.names=FALSE, row.names=F, quote=F, sep="\t")
#     write.table(        command  , file = jj.ctg.command  , col.names=FALSE, row.names=F, quote=F, sep="\t")
#     cat("\n")    
# 
#       system( paste( "seqtk subseq", assembly.file,   jj.ctg.list  ,">", jj.ctg.fasta , "; samtools faidx ", jj.ctg.fasta , ";  echo ' ---> Assembly:' ",i, ";  echo ' ---> contig list:'",jj.ctg.fai, ";  cut -f1 ",jj.ctg.fai ) )   
# #     system( paste( "sed -i 's/WW/xxxxxxxxxxxxn/g' ", jj.ctg.command  ) )   
# #     system( paste( "sed    's/xxxxxxxxxxxx/\/g' ", jj.ctg.command  ) )   ###EDIT HERE
# #     system( paste( "bash ", jj.ctg.command  ) )   
#     
#    # https://www.biostars.org/p/294920/
#    # cat multifasta.fa | sed -e '1!{/^>.*/d;}' | sed  ':a;N;$!ba;s/\n//2g' > output.fa
#    
#    #  CTG
#    
#    library(seqinr)
#    jj.ctg.seq = read.fasta(jj.ctg.fasta)
#    jj.ctg.seq = lapply(jj.ctg.seq, paste, collapse="")
#    cat(" ---> Contig lengths: \n")    
#    print( sapply(jj.ctg.seq, nchar) )
#    polyN = paste( rep("N",100), collapse="")
#    jj.ctg.seq [-1 ]= lapply( jj.ctg.seq [-1 ], function(x) {  x=paste(polyN,x, sep="") }  )
#                 lapply(jj.ctg.seq,  function(x)  substr(x,1,150))
#    cat(" ---> Adding 100N: \n")    
#    jj.scaff.seq = sapply(jj.ctg.seq, as.character)
#    cat(" ---> Check 100N at the beginning of Contigs excluding the first one: \n")    
#      print( substr(jj.scaff.seq,1,150) )
#    jj.scaff.seq = paste(jj.scaff.seq, collapse="")
#      print( substr(jj.scaff.seq,1,150) )
#    jj.scaff.seq = c( paste0(">",jj),jj.scaff.seq)  
#    cat("\n ---> Check New scaffold name: \n")    
#      print( substr(jj.scaff.seq,1,150) )
# 
#    write.table( jj.scaff.seq  , file = jj.sca.fasta  , col.names=FALSE, row.names=F, quote=F, sep="\t")
#    }
#    
#    
# 
# # 
# #   seqtk subseq /ibex/scratch/projects/c2067/celiim/Ajwa_dp/Ajwa_hifiasm.100Kb.fasta <( echo "ptg000180l")  > "ptg000180l".fasta  ; samtools faidx "ptg000180l".fasta 
# # 
# #    
# #     cat <( echo '>Super-Scaffold_100071') 
# #         <(paste 
# #              <(grep -v '^>' test/Bionano_fasta/Super-Scaffold_100071.ctg.fasta) 
# #              <(cat test/Bionano_fasta/Super-Scaffold_100071.NNN.txt) -d '' | sed  ':a;N;$!ba;s/\n//1g' ) 
# #              > test/Bionano_fasta/Super-Scaffold_100071.scaffold.fasta
# #              
# #   cat <( echo '>Super-Scaffold_100071') <(paste <(grep -v '^>' multifasta.fa                                     ) <(cat test.n.fa                                       ) -d '' | sed  ':a;N;$!ba;s/\n//1g' ) > tes.merged.fasta        
# #   cat <( echo '>Super-Scaffold_100071') <(paste <(grep -v '^>' Super-Scaffold_100071.ctg.fasta                   ) <(cat Super-Scaffold_100071.NNN.txt                   ) -d '' | sed  ':a;N;$!ba;s/\n//1g' ) > Super-Scaffold_100071.scaffold.fasta          
# #   cat <( echo '>Super-Scaffold_100071') <(paste <(grep -v '^>' test/Bionano_fasta/Super-Scaffold_100071.ctg.fasta) <(cat test/Bionano_fasta/Super-Scaffold_100071.NNN.txt) -d '' | sed  ':a;N;$!ba;s/\n//1g' ) > test/Bionano_fasta/Super-Scaffold_100071.scaffold.fasta
# #   cat <( echo '>Super-Scaffold_100071') <(paste <(grep -v '^>' test/Bionano_fasta/Super-Scaffold_100071.ctg.fasta) <(cat test/Bionano_fasta/Super-Scaffold_100071.NNN.txt) -d '' | sed  ':a;N;$!ba;s/\n//1g' ) > test/Bionano_fasta/Super-Scaffold_100071.scaffold.fasta                 
# #     cat("\n")
#  
#  
# 
#      
# if ( plot_type %in% c("single","multiple")) {   dev.off() }
#  }
#  
# cat("\n")
# 
# } else {
#     cat("\n\n\n")
#     cat("############################################################################\n")
#     cat("### 10a) no BIONANO AGP provided                                         ###\n")
#     cat("############################################################################\n")
#     
#     }





if ( ctg.gap == TRUE) 
{
cat("\n\n\n")
cat("############################################################################\n")
cat("### 11a) Select contigs for Gap Filling of the best scaffold              ###\n")
cat("###         --> still a work in progress                                 ###\n")
cat("############################################################################\n")

# Select Scaffold by base.name
# TO BE CORRECTED IN ABSENCE OF MAIN
if ( (main %in% SAMPLE$"short_name")==FALSE &  sum(grepl(main,SAMPLE$"short_name"))==1  ) main = grep(main,SAMPLE$"short_name", value=TRUE)

# if scaffold not present, take the one with less gaps
if ( length(main)==0 ) main = colors2[ which.min(colors2$"gaps"),]$"assembly"

# if base.name not present, take the one with less gaps
if ( length(main)>0 & (main %in% SAMPLE$"short_name")==FALSE & sum(grepl(main,SAMPLE$"short_name"))!=1 ) main = colors2[ which.max(colors2$"perc.gaps"),]$"assembly"


if ( dir.exists(paste0(out.dir,"/contig_selection"))==FALSE ) dir.create(paste0(out.dir,"/contig_selection"))

names.ass = c( main, paste(main,"Ctgs"), "gap",  setdiff(colors2$"assembly", main) )
NAMES.ass =1:length(names.ass)
names(NAMES.ass) = names.ass
OTHER.ASS=       setdiff(colors2$"assembly", main) 
other.ass=paste( setdiff(colors2$"assembly", main) , collapse="/")
other.as2=paste( setdiff(colors2$"assembly", main) , collapse="+")
other.as3=paste( "{", paste(OTHER.ASS,collapse=","),  "}" , sep="")


#   CTG.min = 000
pdf( paste0(out.dir,"/GS-viewer_11.",toupper(main),"_GAP_filling.",OUT,"pdf"),width=plot.w, height=plot.h )
par( mar=c(2,8,4,1))
# chr = chr[ chr$"chr" != "Chr0", ]

# Select best Scaffold

bar=barplot(ref$"ref.len", names.arg=ref$"chr", las=1, horiz=T, border="black", col="white", main= paste("Gaps in",toupper(main),"filled with", other.ass,"contigs"), cex.main=2.5, xlim=c(0,max(ref$"ref.len")*1.20), cex.names=2,cex.axis=2)
b2 = data.frame(chr.pos = bar[,1])
b2$"chr" = ref$"chr"
rrr = merge(ref,b2, by="chr")

text( max(ref$"ref.len")/2, max(b2$"chr.pos")*1.06,  paste( "min plotted window =",CTG.min,"bp"), cex=0.5 )

sh = max(ref$"ref.len")*0.01/2

n.ass = length( A.list[[1]] ) +2

shifts = rev(seq(-0.5,0.5,length.out=n.ass+1))

for( p in shifts ) segments(rep(0, nrow(rrr)),rrr$"chr.pos"+p, rrr$"ref.len", rrr$"chr.pos"+p, col="gray", cex=2)

for (i in names(A.list)) 1==1

names.ass = c( main, paste(main,"Ctgs"), "gap",  setdiff(colors2$"assembly", main) )
NAMES.ass =1:length(names.ass)
names(NAMES.ass) = names.ass

CTG.fill=list()

overalp.ctg = 0

for (i in names(A.list)) { 
#  i=  names(A.list)[1]
#  i=  names(A.list)[5]
    A.list [[i]]  -> AA
    B.list [[i]]  -> BB
    g.list [[i]]  -> gg
        
    bb = rrr[ rrr$"chr" == i, ]
    bb$"chr.pos"
    
    CEX.txt = 1
    if(length(CHR)>5) CEX.txt = 0.6

    x = 0.02
    cat("\n",i,"")
    
    # k = 1 -- > MAIN
    kkk = which( names(NAMES.ass) == main)
    kk  = NAMES.ass[[kkk]]
    k   = which(names(AA)==main)
    
    A=AA[[k]] ;
    ass=names.ass[kkk]
    as2=names(AA)[k]
   # cat(ass,"+", as2,"")
    cat("main =", as2,"")
    COL =  unique(col.data [ col.data$"assembly" == ass ,]$"col" )
    s1=shifts[kk];
    s2=shifts[kk+1];
    y1 = bb$"chr.pos" + s1
    y2 = bb$"chr.pos" + s2
    y3 = (y1+y2)/2
   # cat(k,"")
    if( nrow(A)>0)
    {
         
        A.coord = list()
        for (j in 1:nrow(A)) A.coord[[j]]  = c(A$"ref.start"[j]: A$"ref.end"[j])           
        A.coord = unique(do.call(c,A.coord))
        A.perc =round(100* length(A.coord)/A$"ref.len"[1])
        A.coo.list [[ass]] [[i]] = A.coord
                
        A$"chr.ctg"= paste(A$"ref.chr",A$"query",sep=".")
        A = A[ A$"ctg.len" > CTG.min , ]
        A = A[ order( A$"chr.ctg",A$"ref.start") , ]
    
        for (j in 1:nrow(A)) { polygon( c(A$"ref.start"[j], A$"ref.end"[j], A$"ref.end"[j],A$"ref.start"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5 , col=COL , border=NA)    }
        text( A$"ref.len"[1]+sh  , y3 ,  paste0(A.perc,"% of ref. covered"), cex=CEX.txt ,col=COL, adj=0)
        
        
        col.int = c("chr.ctg","ref.start","ref.end","query","q.len","ctg.len","scaf.perc","ctg.ref")
    
        ac = A[,..col.int]
        ac = ac[ order(ac$"chr.ctg",ac$"ref.start",ac$"ref.end"), ]
    
        INPUT.m  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.bed")
        OUT.MRG  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.merged,bed")
        OUT.M.I  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.merged.intersect.bed")
        write.table( ac , INPUT.m  ,col.names=FALSE, row.names=F, quote=F, sep="\t")
        BED.WIND = 50000
        system( paste0( " bedtools merge -i ",INPUT.m," -d ",BED.WIND, "  > ",OUT.MRG) )
        system( paste0( " bedtools intersect -wao -a ",OUT.MRG," -b ",INPUT.m, " > ",OUT.M.I) )
        
        bi= fread(OUT.M.I, stringsAsFactors=FALSE)
        bm= fread(OUT.MRG, stringsAsFactors=FALSE)
        
        names(bi) = c("CHR","START","END", col.int, "overlap")
        names(bm) = c("CHR","START","END")
        
        bm = split(bm, bm$"CHR")
        b3 = lapply(bm, function(x) { y=x[1,] ; y$"START"=min(x$"START") ;  y$"END"=max(x$"END") ; return(y)  })
        b3 = do.call(rbind,b3)
        b3 = b3[ order(b3$"START"),]  
        n3 =nrow(b3)
        b3$"gap.start"= c(1,b3$"END"[-n3])
        b3$"gap.end"   =c(b3$"START"-1)
        b3$"gap.chr" =i
        if( max(b3$"END") < A$"ref.len"[1])  b3 = rbind(b3, data.frame("CHR"=paste0(i,".end"),"START"=NA, "END"=NA, "gap.start"=max(b3$"END")+1, "gap.end"=A$"ref.len"[1]-1, "gap.chr" =i, stringsAsFactors=TRUE ) )
        b3$"gap.len" =  b3$"gap.end" -b3$"gap.start" +1
        b3$"ctg.overlap"="no"
        
        # overlapping ctg on main assembly
        w = which(b3$"gap.len"< 0)
        if ( length(w)>0) { b3$"ctg.overlap"[w]="yes" ;   x=b3$"gap.start"[w] ; b3$"gap.start"[w]=b3$"gap.end"[w]+1 ;  b3$"gap.end"[w]=x-1  }
       # if ( length(w)>0) { b3$"ctg.overlap"[w]="yes" ;  b3 = b3[ b3$"ctg.overlap" == "no", ] }
        b3$"gap.len" =  b3$"gap.end" -b3$"gap.start" +1
        b3$"gap.mid" = (b3$"gap.end" +b3$"gap.start")/2
        
        # k = 2 -- > MAIN CONTIGS
        k = which( names.ass == paste(main,"Ctgs"))
         
        ass=names.ass[k]
      #  cat(ass,"")
        s1=shifts[k];
        s2=shifts[k+1];
        y1 = bb$"chr.pos" + s1
        y2 = bb$"chr.pos" + s2
        y3 = (y1+y2)/2
       
        for (j in 1:nrow(b3)) { polygon( c(b3$"START"[j], b3$"END"[j], b3$"END"[j],b3$"START"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5   ,  col=adjustcolor(COL,alpha.f=0.5) , border="black")    }
        text( A$"ref.len"[1]+sh  , y3 ,  paste(main,"contigs"), cex=CEX.txt ,col="black", adj=0)
       
        # k = 3 -- > GAPS
        k = which( names.ass == "gap")
         
        ass=names.ass[k]
        cat( paste0( " -->  CTGs vs ",ass,"s :  ") ) 
        s1=shifts[k];
        s2=shifts[k+1];
        y1 = bb$"chr.pos" + s1
        y2 = bb$"chr.pos" + s2
        y3 = (y1+y2)/2
       
        for (j in 1:nrow(b3)) { polygon( c(b3$"gap.start"[j], b3$"gap.end"[j], b3$"gap.end"[j],b3$"gap.start"[j]) ,c( y2 , y2, y1, y1) ,   col="black" , border="black")    }
        ovlp.col ="yellow"
        ovlp.col ="seagreen1"
        if ( length(w)>0 ) { points(b3$"gap.mid"[w],  rep(y3,length(w)), type="p", pch=21, cex=1.5 ,col="black", bg=ovlp.col, lwd=2.5 )  }  
        if ( length(w)==0) { text(  A$"ref.len"[1]+sh  , y3 , paste( nrow(b3),"gaps"), cex=CEX.txt ,col="black", adj=0) }
        if ( length(w)> 0) { text(  A$"ref.len"[1]+sh  , y3 , paste( nrow(b3),"gaps", length(w), " on overlapping CTGa"), cex=CEX.txt ,col="black", adj=0) }
        if ( length(w)> 0) { overalp.ctg = overalp.ctg + length(w) }

        b3.ov = b3
        b3 = b3[ b3$"ctg.overlap" == "no",]
       
        #  other k --> other assemblies
      #  for ( kkk in  setdiff(colors2$"assembly", main) ) print(kkk)
        for ( kkk in  setdiff(colors2$"assembly", main) ) 
        {
            # kkk = setdiff(colors2$"assembly", main)[1] ;
            kk = NAMES.ass[[kkk]]
            k = which(names(AA)==kkk)
            A=AA[[k]] ;
            ass=names.ass[kk]
            as2=names(AA)[k]
        #    cat(ass,"+",as2,"")
            cat(ass,"  ")
            COL =  unique(col.data [ col.data$"assembly" == ass ,]$"col" )
            s1=shifts[kk];
            s2=shifts[kk+1];
            y1 = bb$"chr.pos" + s1
            y2 = bb$"chr.pos" + s2
            y3 = (y1+y2)/2
            
            if( nrow(A)>0)
            { 
                # plot assembly coverage ib gray
                A$"chr.ctg"= paste(A$"ref.chr",A$"query",sep=".")
                A = A[ A$"ctg.len" > CTG.min , ]
                A = A[ order(A$"chr.ctg",A$"ref.start") , ]
                for (j in 1:nrow(A)) { polygon( c(A$"ref.start"[j], A$"ref.end"[j], A$"ref.end"[j],A$"ref.start"[j]) ,c( y2 , y2, y1, y1) , lwd=1 , col=adjustcolor("gray70",alpha.f=1) , border=NA)    }
                        
                col.int = c("chr.ctg","ref.start","ref.end","query","q.len","ctg.len","scaf.perc","ctg.ref","ref.chr")
                col.b4  = c("gap.chr","gap.start","gap.end")
            
                ac =  A[,..col.int]
                b4 = b3[,..col.b4]  # gaps coordinates
                b4 = b4[which(b4$"gap.end" > b4$"gap.start" ), ]  # gaps coordinates
                
               
                # merge ass CTGs 
                GAP.fil  = paste0(out.dir,"/contig_selection/",main, ".gap_coordinates.bed")
                INPUT.m  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.bed")
                OUT.MRG  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.merged,bed")
                OUT.M.I  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.merged.intersect.bed")
                OUT.M_2  = paste0(out.dir,"/contig_selection/",ass, ".contig_max_coordinates.bed")
                OUT.GAP  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.merged.intersect.vs_gaps.bed")
                write.table( ac , INPUT.m  ,col.names=FALSE, row.names=F, quote=F, sep="\t")
                write.table( b4 , GAP.fil  ,col.names=FALSE, row.names=F, quote=F, sep="\t")
                BED.WIND = 50000
                system( paste0( " bedtools merge -i ",INPUT.m," -d ",BED.WIND, "  > ",OUT.MRG) )
                system( paste0( " bedtools intersect -wao -a ",OUT.MRG," -b ",INPUT.m, " > ",OUT.M.I) )
              # system( paste0( " paste <(cut -f12 ",OUT.M.I,") <(cut -f2-30 ",OUT.M.I, ") > ",OUT.M_2) ) #does not work         
                
                bi = fread(OUT.M.I, stringsAsFactors=FALSE)
                bm = fread(OUT.MRG, stringsAsFactors=FALSE)            
                names(bm) = c("CHR","START","END")
                names(bi) = c("CHR","START","END", col.int, "overlap")
                bi$"CHR" = bi$"ref.chr" 
                bm = split(bm, bm$"CHR")
                b33 = lapply(bm, function(x) { y=x[1,] ; y$"START"=min(x$"START") ;  y$"END"=max(x$"END") ; return(y)  })
                b33 = do.call(rbind,b33)
                b33 = b33[ order(b33$"START"),]  
                b33$"query" = as.character( do.call(rbind,strsplit(b33$"CHR",".", fixed=T))[,2]) 
                b33$"CHR"   = as.character( do.call(rbind,strsplit(b33$"CHR",".", fixed=T))[,1])
    
                write.table( b33 , OUT.M_2  ,col.names=FALSE, row.names=F, quote=F, sep="\t")
    
                # whole.CTG vs GAPs
                system( paste0( " bedtools intersect -wao -a ",OUT.M_2," -b ",GAP.fil, " > ",OUT.GAP) )
                
                bg = fread(OUT.GAP, stringsAsFactors=FALSE)
                if(nrow(bg) >1 )
                {
                    names(bg) = c("CHR","START","END", "query","gap.chr", "gap.start","gap.end" ,"overlap")
                    bg = bg[ bg$"gap.start">0, ]
                    bg$"gap.len" =  bg$"gap.end" -bg$"gap.start" +1
                   
                    bii = bi
                   
                    # CTG filling GAPS ( 100% or less)
                    
                    bg4 = bg
                    bg4$"perc.gap" = round(100*(bg4$"overlap"+1)/bg4$"gap.len",2)
                    bg4$"col.adj" = "-"
                    if ( any(bg4$"perc.gap"< 100)) { bg4[bg4$"perc.gap"< 100,]$"col.adj" = adjustcolor( COL, alpha.f=0.25) }
                    if ( any(bg4$"perc.gap"==100)) { bg4[bg4$"perc.gap"==100,]$"col.adj" = adjustcolor( COL, alpha.f=0.75) }
                   
                    # if 1 ctg covers 2 gaps, get the max percentage for the plot
                    col.bg5 = c("query","perc.gap","col.adj")
                    bg5 = bg4[,..col.bg5]
                    bg5 = split(bg5,bg5$"query")
                    bg5 = lapply(bg5, function(x) x[which.max(x$"perc.gap"),])
                    bg5 = do.call(rbind, bg5)
                                    
                    b44 = b33[   b33$"query" %in% bg4$"query",]
                    b32 = b33[ ! b33$"query" %in% bg4$"query",]
                    
                    if( nrow(b32)>0) for (j in 1:nrow(b32)) { polygon( c(b32$"START"[j], b32$"END"[j], b32$"END"[j],b32$"START"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5 ,  col=adjustcolor( "white", alpha.f=0.5) , border="black")    }
                    if( nrow(bg4)>0)
                    {
                        b44 = merge(b44, bg5, by="query", all.x=TRUE)     
                  #     for (j in 1:nrow(b33)) { polygon( c(b33$"START"[j], b33$"END"[j], b33$"END"[j],b33$"START"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5 ,  col="white" , border="black")    }         
                  #     for (j in 1:nrow(b32)) { polygon( c(b32$"START"[j], b32$"END"[j], b32$"END"[j],b32$"START"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5 ,  col=adjustcolor( "white", alpha.f=0.5) , border="black")    }
                        for (j in 1:nrow(b44)) { polygon( c(b44$"START"[j], b44$"END"[j], b44$"END"[j],b44$"START"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5 ,  col=b44$"col.adj"[j] , border="black")    }
                       # for (j in 1:nrow(bg4)) { polygon( c(bg4$"START"[j], bg4$"END"[j], bg4$"END"[j],bg4$"START"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5 ,  col=adjustcolor( COL, alpha.f=0.5) , border="black")    }
                        gap.stat = paste0("-> ", ass,": ", sum(bg4$perc.gap ==100), " closed, ", sum(bg4$perc.gap < 100), " partial.")
                        text( A$"ref.len"[1]+sh  , y3 , gap.stat , cex=CEX.txt ,col="black", adj=0)
                        
                        CTG.fill[[ass]][[i]] = bg4
                        } else { gap.stat = paste0("-> ", ass,": 0 closed, 0 partial.") ; text( rrr[ rrr$"chr" ==i, ]$"ref.len" +sh  , y3 , gap.stat, cex=CEX.txt ,col="black", adj=0) }
                    }
               } else { gap.stat = paste0("-> ", ass,": 0 closed, 0 partial.") ; text( rrr[ rrr$"chr" ==i, ]$"ref.len" +sh  , y3 , gap.stat, cex=CEX.txt ,col="black", adj=0)  }
          }
        
        
        }
        Y1 = bb$"chr.pos" - 0.5 -0.02
        Y2 = bb$"chr.pos" + 0.5 +0.02
        CL = bb$"ref.len"         
        polygon( c(0 , CL,CL,0) ,c( Y2 , Y2, Y1, Y1) , border="black", lwd=4 , col=NA )   

}


CTG.FILL = lapply(CTG.fill, function(x) do.call(rbind,x))
CTG.FILL = lapply(CTG.FILL, function(x) { x$"gap" = paste(x$"gap.chr" ,x$"gap.start" ,x$"gap.end" ,sep=".") ; return(x) } )
for( i in names(CTG.FILL)) CTG.FILL[[i]]$"scaffold"=i
#CTG.all[["all"]] = do.call(rbind,CTG.FILL)

#CTG.FILL = lapply(CTG.FILL, function(x) {split(x, x$"gap") } )



# merging ctgs
# - prioritiy to CTGS spanning over more gaps
# - prioritiy to assembly with more gaps filled

if ( length(CTG.FILL) >0)
{
    CTG.MRG = do.call(rbind,CTG.FILL)
    CTG.MRG$"gap" = paste(CTG.MRG$"gap.chr" ,CTG.MRG$"gap.start" ,CTG.MRG$"gap.end" ,sep=".")
    
    sort.ass = rev(sort(sapply(CTG.FILL,nrow)))
    sort.ass = data.frame( scaffold= names(sort.ass), "gap_per_ass"=as.numeric(sort.ass),stringsAsFactors=FALSE )
    
    gap.ctg = rev(sort(table(CTG.MRG$"query")))
    gap.ctg = data.frame( "query"= names(gap.ctg), "gaps_per_ctg"=as.numeric(gap.ctg),stringsAsFactors=FALSE )
    
    CTG.all = CTG.FILL
    CTG.all[["all_merged"]] = CTG.MRG
    CTG.all = lapply(CTG.all, function(x) { x$"gap" = paste(x$"gap.chr" ,x$"gap.start" ,x$"gap.end" ,sep=".") ; return(x) } )
    CTG.all = lapply(CTG.all, function(x) { x = merge( x, sort.ass ,by="scaffold", all.x=T) ; return(x) } )
    CTG.all = lapply(CTG.all, function(x) { x = merge( x, gap.ctg  ,by="query", all.x=T) ; return(x) } )
    CTG.all = lapply(CTG.all, function(x) { x = split( x , x$"gap"); return(x) } )
    CTG.all = lapply(CTG.all, lapply,   function(x) { x [ rev(order(x$"perc.gap",x$"gaps_per_ctg",x$"gap_per_ass")),] } )
    CTG.all = lapply(CTG.all, lapply,   head,1  )
    CTG.all = lapply(CTG.all, function(x) do.call(rbind,x)  )
    
    n_fill = sapply(CTG.all, function(x) length(unique(x[x$"perc.gap"==100,]$"gap")))
    n_part = sapply(CTG.all, function(x) length(unique(x[x$"perc.gap" <100,]$"gap")))
    n_tota = sapply(CTG.all, function(x) length(unique(x[x$"perc.gap"<=100,]$"gap")))
    
    n_fill = data.frame( "assembly"= names(n_fill), "filled"=as.numeric(n_fill),stringsAsFactors=FALSE )
    n_part = data.frame( "assembly"= names(n_part), "partially.filled"=as.numeric(n_part),stringsAsFactors=FALSE )
    n_tota = data.frame( "assembly"= names(n_tota), "total"      =as.numeric(n_tota),stringsAsFactors=FALSE )
    
    n_fill = merge(n_fill,n_part,by="assembly",all.x=TRUE, all.y=TRUE)
    n_fill = merge(n_fill,n_tota,by="assembly",all.x=TRUE, all.y=TRUE)
    n_fill = merge(n_fill,col.data,by="assembly",all.x=TRUE)
    
    n_fill$"nchar" =  trunc( (-nchar( n_fill$"assembly") + max(nchar( n_fill$"assembly"))) *1.80, 1)
    wnc = which(n_fill$"nchar" >=4); if(length(wnc>0)) n_fill$"nchar"[wnc] = n_fill$"nchar"[wnc] +1
    n_fill$"assembly" = paste0( n_fill$"assembly",": ")
    
    for ( nr in 1:nrow(n_fill)) n_fill$"assembly"[nr] = paste0(n_fill$"assembly"[nr], paste( rep(" ", n_fill$"nchar"[nr]), collapse=""))
    n_fill$"comment" = paste0( n_fill$"assembly", n_fill$"filled", " filled + ", n_fill$"partially.filled"," incomplete")
    n_fill$"assembly" = gsub(" |:","", n_fill$"assembly")
    n_fill[ n_fill$"assembly" =="all_merged",]$"col" ="black"
    
    rownames(n_fill) = n_fill$"assembly"
    n_fill = n_fill[ c(setdiff(colors2$"assembly", main), "all_merged"),]
}    

colors4 = colors3
colors4 = colors4[ ! colors4$"assembly" %in% c("Reference","Gaps"," ") ,]
colors4$"nchar" =  trunc( (-nchar( colors4$"assembly") + max(nchar( colors4$"assembly"))) *1.80, 1)
wnc = which(colors4$"nchar" >=4); if(length(wnc>0)) colors4$"nchar"[wnc] = colors4$"nchar"[wnc] +1
rownames(colors4) = colors4$"assembly"
colors4 = colors4[ order(colors4$"gaps"),]
colors4 = colors4[ c(main,setdiff(colors2$"assembly", main)),]

for ( nr in 1:nrow(colors4)) colors4$"assembly"[nr] = paste0(colors4$"assembly"[nr], paste( rep(" ", colors4$"nchar"[nr]), collapse=""))
colors4$"n.gap" = paste( colors4$"assembly","  ", colors4$"gaps", "gaps")
colors4$"assembly" = gsub(" ","", colors4$"assembly")


colors4  = colors4[, c("assembly","col","gaps","n.gap")]
colors22 = colors2[, c("assembly","ta3")]
colors4$"order"= 1:nrow(colors4)
colors4 = merge(colors4,colors22,by="assembly",all.x=TRUE)
colors4$"n.gap" = paste0( colors4$"n.gap",", ", colors4$"ta3", " Ref. covered")
colors4 = colors4[ order(colors4$"order"),]
colors4$"ta3" = NULL

colors4$"bg" = colors4$"col"
#colors3$"col" = "black"
colors4$"pch" = 22
#colors4 = rbind(colors4, c("all_merged","na","na","na CTGs",ovlp.col,23))
colors4 = rbind(         c("","white","na","SCAFFOLDS:",0,"white",NA,NA),colors4)
colors4 = rbind(colors4, c("","white","na"," ",nrow(colors4),"white",NA)        )


#colors4=merge(colors4,n_fill, by="assembly", all.x=TRUE)
colors4[ colors4$"assembly" ==main,]$"n.gap" = paste(colors4[ colors4$"assembly" ==main,]$"n.gap", "--> main!")

if ( length(CTG.FILL) >0)
{
    n_fill$"bg" = n_fill$"col"
    n_fill$"n.gap" = n_fill$"comment"
    n_fill$"gaps" = NA
    n_fill$"order" = ( nrow(colors4)+1): (nrow(colors4)+nrow(n_fill))
    n_fill$"pch" = 22
    
    colors5 = n_fill[ , names(colors4)]
    colors4 = rbind(colors4, c("GAP FILLINGS:","white","na","GAP filling:","white",nrow(colors4)+1,NA,NA))  
    colors6 = rbind(colors4,colors5)

    #colors6 = rbind(colors6, c("","white","na"," ","white",NA,nrow(colors6)+1))
    CTG.overlap_n = paste0( "Overlaping CTGs in the same assembly (",overalp.ctg, ")")
    
    colors6 = rbind(colors6, c(CTG.overlap_n                       ,"black","na",  CTG.overlap_n,nrow(colors6)+1,ovlp.col,21))
    colors6 = rbind(colors6, c("CTG Dark color = Gap filled 100%"  ,"white","na","CTG Dark color = Gap filled 100%",nrow(colors6)+1,"white",22))
    colors6 = rbind(colors6, c("CTG Light color = Gap part. filled","white","na","CTG Light color = Gap part. filled",nrow(colors6)+1,"white",22))
    
    colors6
    } else  {     colors6 = rbind(colors4, c("Overlaping CTGs in the same assembly","black","na","Overlaping CTGs in the same assembly",nrow(colors4)+1,ovlp.col,21)) }

legend("right", colors6$"n.gap", cex=CEX.leg*0.8, pt.bg=colors6$"bg", col=colors6$"col" , pch=as.numeric(colors6$"pch") , bg=adjustcolor("white",alpha.f=1), pt.lwd =3 )


 dev.off()       

 cat("\n")

    
# STATISTICS

cat("\n\n\n")
cat("############################################################################\n")
cat("### 11b) Create Fasta file with selected Contigs for filling GAPS     ###\n")
cat("############################################################################\n")

if ( ctg.fa == TRUE & length(CTG.FILL) >0) 
{ 
cat("\n")
cat(" ==> Saving contigs lists:")
cat("\n")
for( i in names(CTG.all)) write.table( CTG.all[[i]]        , file = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.",i,".",OUT,"tsv"  )  , col.names=TRUE, row.names=F, quote=F, sep="\t")
for( i in names(CTG.all)) write.table( CTG.all[[i]]$"query", file = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.",i,".",OUT,"list"  ) , col.names=FALSE, row.names=F, quote=F, sep="\t")

# crate separate files for the merged ctg list
CTG.merg_split = split(CTG.all[["all_merged"]],CTG.all[["all_merged"]]$"scaffold")
for( i in names(CTG.merg_split)) write.table(         CTG.merg_split[[i]]        , file = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.",i,".","for_merged_ctg.",OUT,"tsv"  )  , col.names=TRUE, row.names=F, quote=F, sep="\t")
for( i in names(CTG.merg_split)) write.table( unique(CTG.merg_split[[i]]$"query"), file = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.",i,".","for_merged_ctg.",OUT,"list"  ) , col.names=FALSE, row.names=F, quote=F, sep="\t")

cat("\n")
cat(" ==> Contigs lists:")
cat("\n")
system( paste( "ls -lth" , paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.*tsv"  ) ) )
cat("\n")


cat(" ==> Contigs .fasta and .fai:\n")
cat("\n")

# SAVE CTG LISTS and create CTG.fasta
for (i in OTHER.ASS) 1
for (i in OTHER.ASS) 
{
    assembly.file     = LOG[ LOG$"short_name" ==i,]$"assembly_file"
    ctgs_for_gaps     = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.",i,".",OUT,"list"  ) 
    ctgs_for_gaps.fa  = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.",i,".",OUT,"fasta"  ) 
    ctgs_for_gaps.fai = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.",i,".",OUT,"fasta.fai"  ) 
    
    system( paste( "seqtk subseq", assembly.file,   ctgs_for_gaps  ,">", ctgs_for_gaps.fa , "; samtools faidx ", ctgs_for_gaps.fa , ";  echo ",i, ";  echo ",ctgs_for_gaps.fai, ";  cut -f1 ",ctgs_for_gaps.fai ) )   
    cat("\n")
    
    }

#merge
cat(" ==> Contigs .fasta and .fai MERGED :\n")
cat("\n")

ii = "all_merged"

for (i in names(CTG.merg_split)) 
{
    assembly.file     = LOG[ LOG$"short_name" ==i,]$"assembly_file"
    ctgs_for_gaps     = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.",i,".","for_merged_ctg.",OUT,"list"  ) 
    ctgs_for_gaps.fa  = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.",i,".","for_merged_ctg.",OUT,"fasta"  ) 
    ctgs_for_gaps.fai = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.",i,".","for_merged_ctg.",OUT,"fasta.fai"  ) 
    
    system( paste( "seqtk subseq", assembly.file,   ctgs_for_gaps  ,">", ctgs_for_gaps.fa , "; samtools faidx ", ctgs_for_gaps.fa , ";  echo ",i, ";  echo ",ctgs_for_gaps.fai, ";  cut -f1 ",ctgs_for_gaps.fai ) )   
    cat("\n")
    
    }
    

all.ctgs_for_gaps.fa      = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.","*","for_merged_ctg.",OUT,"fasta"  ) 
all.ctgs_for_gaps.fai     = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.","*","for_merged_ctg.",OUT,"fasta.fai"  ) 
all.ctgs_for_gaps.out.fa  = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.","MERGED_CTG_LIST.",OUT,"fasta"  ) 
all.ctgs_for_gaps.out.fai = paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.","MERGED_CTG_LIST.",OUT,"fasta.fai"  ) 

system( paste( "ls -lth"       , all.ctgs_for_gaps.fa) )   
system( paste( "cat"           , all.ctgs_for_gaps.fa, ">", all.ctgs_for_gaps.out.fa ) )   
system( paste( "samtools faidx", all.ctgs_for_gaps.out.fa ) )   
system( paste( "cut -f1"       , all.ctgs_for_gaps.out.fai ) )   


# check
ctgs.fai = system( paste( "ls " , paste0(out.dir,"/Ragtag_scaffolds_gaps_in_",main,"_filled_by.",other.as3,".",OUT,"list"  ) ) , inter=TRUE)   
CTGS.FAI = list()
for( i in ctgs.fai ) {  j=gsub( paste0(out.dir,"|Ragtag_scaffolds_gaps_in_|",main,"|_filled_by.|",OUT,"|list")   ,"",i) ; j=gsub( ".|/","",j,fixed=TRUE); if(j=="") {j=SS} ; CTGS.FAI[[j]]=system( paste( "cat ", i ), inter=T)}
CTGS.FAI[["all_merge"]] = as.character(fread(all.ctgs_for_gaps.out.fai )$"V1")


#length
cat("\n")
cat(" ==> check number of contigs in fasta.fai file\n")
print( sapply(CTGS.FAI,length) == sapply(CTG.all,nrow))
cat("\n")

cat(" ==> check all of contigs names in fasta.fai file\n")
for (i in names(CTGS.FAI)) cat("    ", i, all(sort(CTGS.FAI[[i]])==sort(CTG.all[[i]]$"query") ) , "\n")
cat("\n")

} else {  cat("\n ==> NO Overalpping Contigs! \n")  }



cat("############################################################################\n")
cat("### 12)  Ref. Coverage + Scaffold Length/GAPs + GAP Filling in 1 Figure  ###\n")
cat("############################################################################\n")

cat("\n  --> Scaffold length","\n")

CHR2 = ref
CHR2$"chr" = paste0(CHR2$"chr","_RagTag")
for (i in names(CHR)) { x=CHR[[i]] ; names(x)=c("chr",i) ; CHR2 = merge(CHR2,x,by="chr", all.x=TRUE) }
CHR2 = as.data.frame(CHR2)
CHR2 [ is.na(CHR2) ] = 0
CHR2$"chr" = gsub("_RagTag","",CHR2$"chr")
rownames(CHR2) = CHR2$"chr" 
CHR2$"ref.len" = NULL
CHR2$"chr" = NULL
CHR2 = t(as.matrix(CHR2))
CHR2.colnames= colnames(CHR2)
CHR2.rownames= rownames(CHR2)

if(nrow(CHR2)>1) CHR2 = CHR2[  rev(rownames(CHR2)) , ]
if(ncol(CHR2)>1) CHR2 = CHR2[ , rev(colnames(CHR2)) ]

if( is.null(dim(CHR2)) )  { CHR2= as.matrix(t(CHR2)) ;  dimnames(CHR2)[[1]] = CHR2.rownames } 


# print(CHR2)
# print(max(CHR2))
# print(max(c(ref$"ref.len",CHR2)))

#   CTG.min = 000
# chr = chr[ chr$"chr" != "Chr0", ]
bar=barplot(ref$"ref.len", names.arg=ref$"chr", las=1, horiz=T, main="Scaffold length vs Reference length",col=gray(0.95), cex.main=3, xlim=c(0, max(c(ref$"ref.len",CHR2,CCmax)*1.20 )), cex.names=CEX.names,cex.axis=CEX.names)
b2 = data.frame(chr.pos = bar[,1])
b2$"chr" = ref$"chr"
rrr = merge(ref,b2, by="chr")

sh = max(ref$"ref.len")*0.01

n.ass = length( A.list[[1]] )

shifts = rev(seq(-0.5,0.5,length.out=n.ass+1))

#for( p in shifts ) segments(rep(0, nrow(rrr)),rrr$"chr.pos"+p, rrr$"ref.len", rrr$"chr.pos"+p, col="gray", cex=2)


GAP2 = list()
for (i in names(A.list)) 1==1
for (i in names(A.list)) { 
#  i=  names(A.list)[4]
    AA = CHR2[,i]
    if( nrow(CHR2)==1 ) { names(AA) = rownames(CHR2)}
    AA = AA[ rev(names(AA))]
    bb = rrr[ rrr$"chr" == i, ]
    bb$"chr.pos"
    
       
    CEX.txt = 1.1
    if(length(CHR)>7) CEX.txt = 0.6

    x = 0.02
    cat("\n",i,"")
     
    for (k in length(AA):1) 
    {  
        A=AA[[k]] ;
        ass=names(AA)[k]
        cat(ass,"")
        COL =  unique(col.data [ col.data$"assembly" == ass ,]$"col" )
        s1=shifts[k];
        s2=shifts[k+1];
        y1 = bb$"chr.pos" + s1
        y2 = bb$"chr.pos" + s2
        y3 = bb$"chr.pos" + (s1+s2)/2
        Y1 = bb$"chr.pos" - 0.5 -0.01
        Y2 = bb$"chr.pos" + 0.5 +0.01
        CL = bb$"ref.len" 
    #    cat(y3,"")

        polygon( c(0 , A, A,0) ,c( y2 , y2, y1, y1) , lwd=2 , col=COL , border="black", lwd=0.2)   
        
        gap = GAPS[[ ass]]
        gap$"y3" = y3
        gap = gap[ gap$"chr" == i, ]
        GAP2[[i]][[ ass]] = gap
        
        if ( nrow(gap)> 0)  points(gap$"ref.end",  gap$"y3", type="p", pch=18, cex=1 ,col="black")
        text( max(c(CL,A))+sh  , y3 ,  paste0(nrow(gap)," gaps"), cex=CEX.txt ,col=COL,adj=0)
        
       # segments(rep(0, nrow(rrr)),rrr$"chr.pos"+p, rrr$"ref.len", rrr$"chr.pos"+p, col="gray", cex=2)

        }          
        polygon( c(0 , CL,CL,0) ,c( Y2 , Y2, Y1, Y1) , lwd=1 , border="blue", lwd=2, lty=1 , col=adjustcolor("gray",alpha.f=0.5) )   
       
}
cat("\n")
cat("\n  --> Adding GAPs","\n")
for (i in names(A.list)) { AA = CHR2[,i] ; if( nrow(CHR2)==1 ) { names(AA) = rownames(CHR2)}; AA = AA[ rev(names(AA))] ; cat("\n",i,"") ; for (k in length(AA):1) { ass=names(AA)[k] ; cat(ass,"") ; gap = GAP2[[i]][[ ass]] ; if ( nrow(gap)> 0) points(gap$"ref.end",  gap$"y3", type="p", pch=23, cex=2 ,col="black", bg="yellow", lwd=2 )}}

GAP.tot  = sapply(GAPS,nrow)
GAP.tot  = data.frame("assembly"= names(GAP.tot), "gaps"=as.numeric(GAP.tot), stringsAsFactors=TRUE)

colors2 = colors
colors2 = merge(colors2, GAP.tot, by="assembly", all.x=TRUE)
colors2$"n.gap" = paste0(colors2$"assembly","   ",colors2$"gap"," gaps")

colors3 = colors2[, c("assembly","col","gaps","n.gap")]

colors3$"bg" = colors3$"col"
#colors3$"col" = "black"
colors3$"pch" = 22

rownames(colors3) = colors3$"assembly"
#cat("\n"); print(data.frame(colors3))

colors3 = colors3[   (SAMPLE$"short_name") , ]
#cat("\n");print(data.frame(colors3))

#colors3 = colors3[  rev(SAMPLE$"short_name") , ]

colors3 = rbind(colors3, c(" ","white"," "," ","white",22))
colors3 = rbind(colors3, c("Reference","blue","na","Reference","gray",22))
colors3 = rbind(colors3, c("Gaps","black","na","Gaps","yellow",23))

rownames(colors3) = colors3$"assembly"
#colors3 = colors3[  names(A.list[[1]]) , ]

#cat("\n");print(data.frame(colors3))


colors3$"tag" = paste0(colors3$"assembly","        ",colors3$"gaps"," gaps")
colors3$"ta2" = paste0(colors3$"assembly","                 ")
colors3$"ta3" = paste0(colors3$"gaps"," gaps")
colors3[ !colors3$"assembly" %in% GAP.tot$"assembly", ]$"ta3"=""

w_1gap = which(colors3$"ta3" ==1)
if(length(w_1gap)>0) colors3$"ta3"[w_1gap] = gsub("gaps","gap",colors3$"ta3"[w_1gap])


legend("right", colors3$"ta2", cex=CEX.leg*0.8, pt.bg=colors3$"bg", col=colors3$"col" , pch=as.numeric(colors3$"pch") , bg=adjustcolor("white",alpha.f=0.5), pt.lwd =3 )
legend("right", colors3$"ta3", cex=CEX.leg*0.8, pt.lwd =3 , bty="n"  )

cat("\n")

cat("\n  -->  Reference Coverage","\n")

#   CTG.min = 000
pdf( paste0(out.dir,"/GS-viewer_12.Coverage_on_Reference_plus_Scaffold_length_gaps.",OUT,"min.",CTG.min,"bp_plus_gap_filling.pdf"),width=plot.w, height=plot.h )
#pdf( paste0(out.dir,"/GS-viewer_11.Coverage_on_Reference_plus_Scaffold_length_gaps.",OUT,","_plus_gap_filling.pdf"),width=plot.w, height=plot.h )
par( mfcol=c(1,3) ,  mar=c(2,7,4,1))
# chr = chr[ chr$"chr" != "Chr0", ]
CEX.names = 2
if ( nchar(max(ref$"chr"))>6 ) CEX.names = 1.4

bar=barplot(ref$"ref.len", names.arg=ref$"chr", las=1, horiz=T, main="Scaffold coverage on Reference Genome",col=gray(0.95), cex.main=3, xlim=c(0,max(ref$"ref.len")*1.20), cex.names=CEX.names,cex.axis=CEX.names)
b2 = data.frame(chr.pos = bar[,1])
b2$"chr" = ref$"chr"
rrr = merge(ref,b2, by="chr")

text( max(ref$"ref.len")/2, max(b2$"chr.pos")*1.06,  paste( "min plotted window =",CTG.min,"bp"), cex=0.5 )

sh = max(ref$"ref.len")*0.01

n.ass = length( A.list[[1]] )

shifts = rev(seq(-0.5,0.5,length.out=n.ass+1))

for( p in shifts ) segments(rep(0, nrow(rrr)),rrr$"chr.pos"+p, rrr$"ref.len", rrr$"chr.pos"+p, col="gray", cex=2)

for (i in names(A.list)) 1==1

A.coo.list=list()

for (i in names(A.list)) { 
#  i=  names(A.list)[1]
    A.list [[i]]  -> AA
    B.list [[i]]  -> BB
    g.list [[i]]  -> gg
        
    bb = b2[ b2$"chr" == i, ]
    bb$"chr.pos"
    
    CEX.txt = 1
    if(length(CHR)>7) CEX.txt = 0.6

    x = 0.02
    cat("\n",i,"")
     
    for (k in length(AA):1) 
    {  
        A=AA[[k]] ;
        ass=names(AA)[k]
        cat(ass,"")
        COL =  unique(col.data [ col.data$"assembly" == ass ,]$"col" )
        s1=shifts[k];
        s2=shifts[k+1];
        y1 = bb$"chr.pos" + s1
        y2 = bb$"chr.pos" + s2
        y3 = (y1+y2)/2
       # cat(k,"")
        if( nrow(A)>0)
        {
         
        A.coord = list()
        for (j in 1:nrow(A)) A.coord[[j]]  = c(A$"ref.start"[j]: A$"ref.end"[j])           
        A.coord = unique(do.call(c,A.coord))
        A.perc =round(100* length(A.coord)/A$"ref.len"[1])
        A.coo.list [[ass]] [[i]] = A.coord
                
        A = A[ A$"ctg.len" > CTG.min , ]

        for (j in 1:nrow(A)) { polygon( c(A$"ref.start"[j], A$"ref.end"[j], A$"ref.end"[j],A$"ref.start"[j]) ,c( y2 , y2, y1, y1) , lwd=2 , col=COL , border=NA)    }
        text( A$"ref.len"[1]+sh*2  , y3 ,  paste0(A.perc,"%"), cex=CEX.txt ,col=COL)
        }   else { A.coo.list [[ass]] [[i]] = data.frame()}       
        }
}


#COL.LIST=c("firebrick3","black","dodgerblue4","darkorange3","darkturquoise","blueviolet","forestgreen","lightcoral")
col.list = COL.LIST[ 1:length(CHR)]

colors= data.frame(assembly=names(A.list [[1]] ), col=col.list[1:length(A.list [[1]] )], row.names=names(A.list [[1]] ), stringsAsFactors=FALSE)

colors = col.data

chr.sum = sapply(A.coo.list,sapply,length)
chr.sum = as.data.frame(chr.sum)
gen.sum = sapply(chr.sum,sum)
gen.sum = data.frame( assembly=names(gen.sum), len=gen.sum, row.names=names(gen.sum), stringsAsFactors=FALSE)
colors  = merge( colors, gen.sum, by="assembly")
colors$"perc" = round( 100* colors$"len" / sum(ref$"ref.len") ,2)
colors$"tag" = paste0(colors$"assembly","        ",colors$"perc","%")
colors$"ta2" = paste0(colors$"assembly","                 ")
colors$"ta3" = paste0(colors$"perc","%")

colors$"bg" = colors$"col"
colors$"pch" = 22

rownames(colors) = colors$"assembly"
colors = colors[  names(A.list[[1]]) , ]
CEX.leg=2.5
legend("right", gsub(".ch0","",colors$"ta2"), cex=CEX.leg*0.8, pt.bg=colors$"bg", col=colors$"col",  pch=as.numeric(colors$"pch")  )
#legend("bottomright", gsub(".ch0","",colors$"tag"), fill=colors$col  )
legend("right", gsub(".ch0","",colors$"ta3"),cex=CEX.leg*0.8, bty="n" )


cat("\n")

cat("\n  --> Filling GAPs","\n")

# Select best Scaffold

bar=barplot(ref$"ref.len", names.arg=ref$"chr", las=1, horiz=T, border="black", col="white", main= paste("Gaps in",toupper(main),"filled with", other.ass,"contigs"), cex.main=2.5, xlim=c(0,max(ref$"ref.len")*1.20), cex.names=CEX.names,cex.axis=CEX.names)
b2 = data.frame(chr.pos = bar[,1])
b2$"chr" = ref$"chr"
rrr = merge(ref,b2, by="chr")

text( max(ref$"ref.len")/2, max(b2$"chr.pos")*1.06,  paste( "min plotted window =",CTG.min,"bp"), cex=0.5 )

sh = max(ref$"ref.len")*0.01/2

n.ass = length( A.list[[1]] ) +2

shifts = rev(seq(-0.5,0.5,length.out=n.ass+1))

for( p in shifts ) segments(rep(0, nrow(rrr)),rrr$"chr.pos"+p, rrr$"ref.len", rrr$"chr.pos"+p, col="gray", cex=2)

for (i in names(A.list)) 1==1

names.ass = c( main, paste(main,"Ctgs"), "gap",  setdiff(colors2$"assembly", main) )
NAMES.ass =1:length(names.ass)
names(NAMES.ass) = names.ass

CTG.fill=list()
overalp.ctg = 0

for (i in names(A.list)) { 
#  i=  names(A.list)[1]
#  i=  names(A.list)[5]
    A.list [[i]]  -> AA
    B.list [[i]]  -> BB
    g.list [[i]]  -> gg
        
    bb = rrr[ rrr$"chr" == i, ]
    bb$"chr.pos"
    
    CEX.txt = 1
    if(length(CHR)>5) CEX.txt = 0.6

    x = 0.02
    cat("\n",i,"")
    
    # k = 1 -- > MAIN
    kkk = which( names(NAMES.ass) == main)
    kk  = NAMES.ass[[kkk]]
    k   = which(names(AA)==main)
    
    A=AA[[k]] ;
    ass=names.ass[kkk]
    as2=names(AA)[k]
    cat("main =",as2,"")
    COL =  unique(col.data [ col.data$"assembly" == ass ,]$"col" )
    s1=shifts[kk];
    s2=shifts[kk+1];
    y1 = bb$"chr.pos" + s1
    y2 = bb$"chr.pos" + s2
    y3 = (y1+y2)/2
   # cat(k,"")
    if( nrow(A)>0)
    {
         
        A.coord = list()
        for (j in 1:nrow(A)) A.coord[[j]]  = c(A$"ref.start"[j]: A$"ref.end"[j])           
        A.coord = unique(do.call(c,A.coord))
        A.perc =round(100* length(A.coord)/A$"ref.len"[1])
        A.coo.list [[ass]] [[i]] = A.coord
                
        A$"chr.ctg"= paste(A$"ref.chr",A$"query",sep=".")
        A = A[ A$"ctg.len" > CTG.min , ]
        A = A[ order( A$"chr.ctg",A$"ref.start") , ]
    
        for (j in 1:nrow(A)) { polygon( c(A$"ref.start"[j], A$"ref.end"[j], A$"ref.end"[j],A$"ref.start"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5 , col=COL , border=NA)    }
        text( A$"ref.len"[1]+sh  , y3 ,  paste0(A.perc,"% of ref. covered"), cex=CEX.txt ,col=COL, adj=0)
        
        
        col.int = c("chr.ctg","ref.start","ref.end","query","q.len","ctg.len","scaf.perc","ctg.ref")
    
        ac = A[,..col.int]
    
        INPUT.m  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.bed")
        OUT.MRG  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.merged,bed")
        OUT.M.I  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.merged.intersect.bed")
        write.table( ac , INPUT.m  ,col.names=FALSE, row.names=F, quote=F, sep="\t")
        BED.WIND = 50000
        system( paste0( " bedtools merge -i ",INPUT.m," -d ",BED.WIND, "  > ",OUT.MRG) )
        system( paste0( " bedtools intersect -wao -a ",OUT.MRG," -b ",INPUT.m, " > ",OUT.M.I) )
        
        bi= fread(OUT.M.I, stringsAsFactors=FALSE)
        bm= fread(OUT.MRG, stringsAsFactors=FALSE)
        
        names(bi) = c("CHR","START","END", col.int, "overlap")
        names(bm) = c("CHR","START","END")
        
        bm = split(bm, bm$"CHR")
        b3 = lapply(bm, function(x) { y=x[1,] ; y$"START"=min(x$"START") ;  y$"END"=max(x$"END") ; return(y)  })
        b3 = do.call(rbind,b3)
        b3 = b3[ order(b3$"START"),]  
        n3 =nrow(b3)
        b3$"gap.start"= c(1,b3$"END"[-n3])
        b3$"gap.end"   =c(b3$"START"-1)
        b3$"gap.chr" =i
        if( max(b3$"END") < A$"ref.len"[1])  b3 = rbind(b3, data.frame("CHR"=paste0(i,".end"),"START"=NA, "END"=NA, "gap.start"=max(b3$"END")+1, "gap.end"=A$"ref.len"[1]-1, "gap.chr" =i, stringsAsFactors=TRUE ) )
        b3$"gap.len" =  b3$"gap.end" -b3$"gap.start" +1
        b3$"ctg.overlap"="no"
        
        # overlapping ctg on main assembly
        w = which(b3$"gap.len"< 0)
        if ( length(w)>0) { b3$"ctg.overlap"[w]="yes" ;   x=b3$"gap.start"[w] ; b3$"gap.start"[w]=b3$"gap.end"[w]+1 ;  b3$"gap.end"[w]=x-1  }
       # if ( length(w)>0) { b3$"ctg.overlap"[w]="yes" ;  b3 = b3[ b3$"ctg.overlap" == "no", ] }
        b3$"gap.len" =  b3$"gap.end" -b3$"gap.start" +1
        b3$"gap.mid" = (b3$"gap.end" +b3$"gap.start")/2
        
        # k = 2 -- > MAIN CONTIGS
        k = which( names.ass == paste(main,"Ctgs"))
         
        ass=names.ass[k]
      #  cat(ass,"")
        s1=shifts[k];
        s2=shifts[k+1];
        y1 = bb$"chr.pos" + s1
        y2 = bb$"chr.pos" + s2
        y3 = (y1+y2)/2
       
        for (j in 1:nrow(b3)) { polygon( c(b3$"START"[j], b3$"END"[j], b3$"END"[j],b3$"START"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5   ,  col=adjustcolor(COL,alpha.f=0.5) , border="black")    }
        text( A$"ref.len"[1]+sh  , y3 ,  paste(main,"contigs"), cex=CEX.txt ,col="black", adj=0)
       
        # k = 3 -- > GAPS
        k = which( names.ass == "gap")
         
        ass=names.ass[k]
        cat( paste0( " -->  CTGs vs ",ass,"s :  ") ) 
        s1=shifts[k];
        s2=shifts[k+1];
        y1 = bb$"chr.pos" + s1
        y2 = bb$"chr.pos" + s2
        y3 = (y1+y2)/2
       
        for (j in 1:nrow(b3)) { polygon( c(b3$"gap.start"[j], b3$"gap.end"[j], b3$"gap.end"[j],b3$"gap.start"[j]) ,c( y2 , y2, y1, y1) ,   col="black" , border="black")    }
        ovlp.col ="yellow"
        ovlp.col ="seagreen1"
        if ( length(w)>0 ) { points(b3$"gap.mid"[w],  rep(y3,length(w)), type="p", pch=21, cex=1.5 ,col="black", bg=ovlp.col, lwd=2.5 )  }  
        if ( length(w)==0) { text(  A$"ref.len"[1]+sh  , y3 , paste( nrow(b3),"gaps"), cex=CEX.txt ,col="black", adj=0) }
        if ( length(w)> 0) { text(  A$"ref.len"[1]+sh  , y3 , paste( nrow(b3),"gaps", length(w), " on overlapping CTGa"), cex=CEX.txt ,col="black", adj=0) }
        if ( length(w)> 0) { overalp.ctg = overalp.ctg + length(w) }
                
        b3.ov = b3
        b3 = b3[ b3$"ctg.overlap" == "no",]
       
        #  other k --> other assemblies
      #  for ( kkk in  setdiff(colors2$"assembly", main) ) print(kkk)
        for ( kkk in  setdiff(colors2$"assembly", main) ) 
        {
            # kkk = setdiff(colors2$"assembly", main)[1] ;
            kk = NAMES.ass[[kkk]]
            k = which(names(AA)==kkk)
            A=AA[[k]] ;
            ass=names.ass[kk]
            as2=names(AA)[k]
         # cat(ass,as2,"")
            cat(ass," ")
            COL =  unique(col.data [ col.data$"assembly" == ass ,]$"col" )
            s1=shifts[kk];
            s2=shifts[kk+1];
            y1 = bb$"chr.pos" + s1
            y2 = bb$"chr.pos" + s2
            y3 = (y1+y2)/2
            
            if( nrow(A)>0)
            { 
                # plot assembly coverage ib gray
                A$"chr.ctg"= paste(A$"ref.chr",A$"query",sep=".")
                A = A[ A$"ctg.len" > CTG.min , ]
                A = A[ order(A$"chr.ctg",A$"ref.start") , ]
                for (j in 1:nrow(A)) { polygon( c(A$"ref.start"[j], A$"ref.end"[j], A$"ref.end"[j],A$"ref.start"[j]) ,c( y2 , y2, y1, y1) , lwd=1 , col=adjustcolor("gray70",alpha.f=1) , border=NA)    }
                        
                col.int = c("chr.ctg","ref.start","ref.end","query","q.len","ctg.len","scaf.perc","ctg.ref","ref.chr")
                col.b4  = c("gap.chr","gap.start","gap.end")
            
                ac =  A[,..col.int]
                b4 = b3[,..col.b4]  # gaps coordinates
                b4 = b4[which(b4$"gap.end" > b4$"gap.start" ), ]  # gaps coordinates
               
                # merge ass CTGs 
                GAP.fil  = paste0(out.dir,"/contig_selection/",main, ".gap_coordinates.bed")
                INPUT.m  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.bed")
                OUT.MRG  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.merged,bed")
                OUT.M.I  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.merged.intersect.bed")
                OUT.M_2  = paste0(out.dir,"/contig_selection/",ass, ".contig_max_coordinates.bed")
                OUT.GAP  = paste0(out.dir,"/contig_selection/",ass, ".contig_coordinates.merged.intersect.vs_gaps.bed")
                write.table( ac , INPUT.m  ,col.names=FALSE, row.names=F, quote=F, sep="\t")
                write.table( b4 , GAP.fil  ,col.names=FALSE, row.names=F, quote=F, sep="\t")
                BED.WIND = 50000
                system( paste0( " bedtools merge -i ",INPUT.m," -d ",BED.WIND, "  > ",OUT.MRG) )
                system( paste0( " bedtools intersect -wao -a ",OUT.MRG," -b ",INPUT.m, " > ",OUT.M.I) )
              # system( paste0( " paste <(cut -f12 ",OUT.M.I,") <(cut -f2-30 ",OUT.M.I, ") > ",OUT.M_2) ) #does not work         
                
                bi = fread(OUT.M.I, stringsAsFactors=FALSE)
                bm = fread(OUT.MRG, stringsAsFactors=FALSE)            
                names(bm) = c("CHR","START","END")
                names(bi) = c("CHR","START","END", col.int, "overlap")
                bi$"CHR" = bi$"ref.chr" 
                bm = split(bm, bm$"CHR")
                b33 = lapply(bm, function(x) { y=x[1,] ; y$"START"=min(x$"START") ;  y$"END"=max(x$"END") ; return(y)  })
                b33 = do.call(rbind,b33)
                b33 = b33[ order(b33$"START"),]  
                b33$"query" = as.character( do.call(rbind,strsplit(b33$"CHR",".", fixed=T))[,2]) 
                b33$"CHR"   = as.character( do.call(rbind,strsplit(b33$"CHR",".", fixed=T))[,1])
    
                write.table( b33 , OUT.M_2  ,col.names=FALSE, row.names=F, quote=F, sep="\t")
    
                # whole.CTG vs GAPs
                system( paste0( " bedtools intersect -wao -a ",OUT.M_2," -b ",GAP.fil, " > ",OUT.GAP) )
                bg = fread(OUT.GAP, stringsAsFactors=FALSE)
                if(nrow(bg) >1 )
                {
                names(bg) = c("CHR","START","END", "query","gap.chr", "gap.start","gap.end" ,"overlap")
                    bg = bg[ bg$"gap.start">0, ]
                    bg$"gap.len" =  bg$"gap.end" -bg$"gap.start" +1
                   
                    bii = bi
                   
                    # CTG filling GAPS ( 100% or less)
                    
                    bg4 = bg
                    bg4$"perc.gap" = round(100*(bg4$"overlap"+1)/bg4$"gap.len",2)
                    bg4$"col.adj" = "-"
                    if ( any(bg4$"perc.gap"< 100)) { bg4[bg4$"perc.gap"< 100,]$"col.adj" = adjustcolor( COL, alpha.f=0.25) }
                    if ( any(bg4$"perc.gap"==100)) { bg4[bg4$"perc.gap"==100,]$"col.adj" = adjustcolor( COL, alpha.f=0.75) }
                   
                    # if 1 ctg covers 2 gaps, get the max percentage for the plot
                    col.bg5 = c("query","perc.gap","col.adj")
                    bg5 = bg4[,..col.bg5]
                    bg5 = split(bg5,bg5$"query")
                    bg5 = lapply(bg5, function(x) x[which.max(x$"perc.gap"),])
                    bg5 = do.call(rbind, bg5)
                                    
                    b44 = b33[   b33$"query" %in% bg4$"query",]
                    b32 = b33[ ! b33$"query" %in% bg4$"query",]
                    
                    if( nrow(b32)>0) for (j in 1:nrow(b32)) { polygon( c(b32$"START"[j], b32$"END"[j], b32$"END"[j],b32$"START"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5 ,  col=adjustcolor( "white", alpha.f=0.5) , border="black")    }
                    if( nrow(bg4)>0)
                    {
                        b44 = merge(b44, bg5, by="query", all.x=TRUE)     
                  #     for (j in 1:nrow(b33)) { polygon( c(b33$"START"[j], b33$"END"[j], b33$"END"[j],b33$"START"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5 ,  col="white" , border="black")    }         
                  #     for (j in 1:nrow(b32)) { polygon( c(b32$"START"[j], b32$"END"[j], b32$"END"[j],b32$"START"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5 ,  col=adjustcolor( "white", alpha.f=0.5) , border="black")    }
                        for (j in 1:nrow(b44)) { polygon( c(b44$"START"[j], b44$"END"[j], b44$"END"[j],b44$"START"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5 ,  col=b44$"col.adj"[j] , border="black")    }
                       # for (j in 1:nrow(bg4)) { polygon( c(bg4$"START"[j], bg4$"END"[j], bg4$"END"[j],bg4$"START"[j]) ,c( y2 , y2, y1, y1) , lwd=0.5 ,  col=adjustcolor( COL, alpha.f=0.5) , border="black")    }
                        gap.stat = paste0("-> ", ass,": ", sum(bg4$perc.gap ==100), " closed, ", sum(bg4$perc.gap < 100), " partial.")
                        text( A$"ref.len"[1]+sh  , y3 , gap.stat , cex=CEX.txt ,col="black", adj=0)
                        
                        CTG.fill[[ass]][[i]] = bg4
                        } else { gap.stat = paste0("-> ", ass,": 0 closed, 0 partial.") ; text( rrr[ rrr$"chr" ==i, ]$"ref.len" +sh  , y3 , gap.stat, cex=CEX.txt ,col="black", adj=0) }
                   }
               } else { gap.stat = paste0("-> ", ass,": 0 closed, 0 partial.") ; text( rrr[ rrr$"chr" ==i, ]$"ref.len" +sh  , y3 , gap.stat, cex=CEX.txt ,col="black", adj=0)  }
          }
        
        
        }
        Y1 = bb$"chr.pos" - 0.5 -0.02
        Y2 = bb$"chr.pos" + 0.5 +0.02
        CL = bb$"ref.len"         
        polygon( c(0 , CL,CL,0) ,c( Y2 , Y2, Y1, Y1) , border="black", lwd=4 , col=NA )   

}
cat("\n")

colors7 = colors6
colors7$"n.gap" = gsub("covered| −","",colors7$"n.gap")

legend("right", colors7$"n.gap", cex=CEX.leg*0.45, pt.bg=colors7$"bg", col=colors7$"col" , pch=as.numeric(colors7$"pch") , bg=adjustcolor("white",alpha.f=1), pt.lwd =3 )


 dev.off()       
 

}


system(" rm Rplots.pdf" )
# Zippung
system(paste0(" zip -rm  ",out.dir,"/ctg_coordinates_on_ref.zip  " , out.dir,"/ctg_coordinates_on_ref > /dev/null 2>&1 " ) )

END.time =Sys.time() 
TOT.time =difftime(END.time, START.time, units="mins")


cat( paste( "==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n\n" ) )
cat( paste( "==>> Output folder: ", out.dir, " \n" ) )
cat( paste( "==>>     full path:  ", getwd(),"/",out.dir, " \n", sep=""))
cat("\n")
cat("\n")









