#!/bin/R

# R
options(width=340)
options(scipen=9999)

args <- commandArgs(TRUE)



   if("--help" %in% args) 
    {
        cat('\n')
        cat(
        " Indentify and visualize Telomeres, Tandem repeats, rDNA and transposable elements in a genome \n",
        " \n",
        "USAGE: Rscript Telomere_finder.R --fasta=<query_scaffold.fasta> --lib=<query_TE_lib.fasta> [ options ] \n",
        " \n",
        "Before the scritp it will run Tandem Repeat Finder and Repeat Masker \n",
        " --> e.g. : \n",
        "   # trf <query_scaffold.fasta>  2 3 5 80 10 30 2000 -l 6 -dat -h  : \n",
        "   # RepeatMasker -qq -no_is --pa 64 -lib <query_TE_lib.fasta> <query_scaffold.fasta> -gff -xsmall   # Do not use -nolow \n",
        " \n",
        " TRF parameters :\n",
       "mandatory arguments:\n",
        "  --fasta=            |  scaffold of assembly .fasta or .fa file ] ","\n",
        "  --lib=              |  TE library, possibily including rDNA  ","\n",
        " \n",
        "external rDNA repeat profile: \n",
        "  --rDNA.fa =         |  rDNA library, if not included in TE library] ","\n",
        "  --rDNA.out=         |  repeat masker out with rDNA library        ] ","\n",
        " \n",
        " other parameters:\n",
        "  --sort.size=        |  sort by size                                                [ default = FALSE  ]  ","\n",
        "  --tel=              |  identify a specififc telomeric patter                       [ default = absent ]  ","\n",
        "  --diploid=          |  group chromosomes by name ( e.g Chr1.A + Chr1.B  )          [ default = absent ]  ","\n",
        "  --polyploid=        |  group chromosomes by name ( e.g Chr1.A + Chr1.B  )          [ default = absent ]  ","\n",
        "  --ref.name=         |  replace GenBank IDs with chromosome names from reference.fna[ default = absent ]  ","\n",
        "  --gene.gff=         |  gene prediction in gff3 format                              [ default = absent ]  ","\n",
        "  --add.dots=         |  .txt file --> Chr, pos,      ,name, color, pch, size        [ default = absent ]  ","\n",
        "  --add.segm=         |  .bed file --> Chr, start, end,name, color                   [ default = absent ]  ","\n",
        "  --chr.filt=         |  Integer, minimal Chr size plotted                           [ default = 1Mb / 1000Kb / 1000000 ; numeric and Kb,Mb character accepted  ]  ","\n",
        "  --sample.name=      |  get samples name from main folder                           [ default = TRUE   ]  ","\n",
        "  --plot.h=           |  Hight in inches of PDF output files                         [ default 18       ]  ","\n",
        "  --plot.w=           |  Width in inches of PDF output files                         [ default 35       ]  ","\n",
        "  --name.space=       |  space for chromosome names in inches                        [ default 6        ]  ","\n",
        "  --leg.size=         |  increase (>0) or decrease (<0) legend size by               [ default 1        ]  ","\n",
        "  --name.del=         |  strings to be deleted from all names separated by ,         [ default absent   ] ", "\n",
        "  --run=              |  plot, all (overwrite intermediate files)                    [ default absent   ] ", "\n",
        " \n",
        " TRF parameters :\n",
        "  --match=            |  matching weight                     [ default = 2    ]  ","\n",
        "  --mismatch=         |  mismatching penalty                 [ default = 3    ]  ","\n",
        "  --delta=            |  indel penalty                       [ default = 5    ]  ","\n",
        "  --pm=               |  match probability (whole number)    [ default = 80   ]  ","\n",
        "  --pi=               |  indel probability (whole number)    [ default = 10   ]  ","\n",
        "  --minscore=         |  minimum alignment score to report   [ default = 30   ]  ","\n",
        "  --maxperiod=        |  maximum period size to report       [ default = 2000 ]  ","\n",
        "  --l=                |  matching weight                     [ default = 6    ]  ","\n",
        "  --consensus.size=   |  consensus.size in bp                [ default = 7    ]  ","\n",
        "  --match=            |  matching weight                     [ default = 2    ]  ","\n",
     " \n",
        "Modules required: \n",
        " - R \n",
        " - R packages: data.table  \n",
      #  " - R packages: data.table , seqinr \n",
        " - trf \n",
        " - samtools \n",
        " - bedtools \n",
        " - Repeat Masker \n",
        " - RagTag (if you need to create AGP files to plot gaps) \n",
         "\n"
        )

        cat('\n')
        q(save="no")
    }



cat('\n')
cat('\n')
cat("############################################################################\n")
cat("############################################################################\n")
cat("###                          Telomere_finder.R                          ###\n")
cat("############################################################################\n")
cat("############################################################################\n")

START.time=Sys.time()
system ( "rm Rplots.pdf")



cat('\n')
cat("############################################################################\n")
cat('-------------------------------OPTIONS--------------------------------------\n')
cat("############################################################################\n")
if (length(args)>0 ) print(as.data.frame(args))
cat('\n')


# COMMANDLINE ARGS

    if (length(args)==0) ARG = data.frame("X1"="--option", "X2"="argument")
    if (length(args)>0 ) ARG = data.frame(do.call(rbind,strsplit(args,"=")))
    #print(ARG)
    
cat('\n')
cat("############################################################################\n")
cat('-------------------------------OPTIONS check--------------------------------\n')
cat("############################################################################\n")
cat('\n')
if (length(args)==0) cat("--> all default \n" )
if (length(args)>0 ) print(ARG, justify = 3 ) 
    
	fasta           = as.character(ARG[ ARG$X1=="--fasta"      ,]$"X2")
	TElib           = as.character(ARG[ ARG$X1=="--lib"        ,]$"X2")
	rDNA.fa         = as.character(ARG[ ARG$X1=="--rDNA.fa"    ,]$"X2")
	rDNA.out        = as.character(ARG[ ARG$X1=="--rDNA.out"   ,]$"X2")
	Match           = as.character(ARG[ ARG$X1=="--match"      ,]$"X2")
	Mismatch        = as.character(ARG[ ARG$X1=="--mismatch"   ,]$"X2")
	Delta           = as.character(ARG[ ARG$X1=="--delta"      ,]$"X2")
	PM              = as.character(ARG[ ARG$X1=="--pm"         ,]$"X2")
	PI              = as.character(ARG[ ARG$X1=="--pi"         ,]$"X2")
	Minscore        = as.character(ARG[ ARG$X1=="--minscore"   ,]$"X2")
	MaxPeriod       = as.character(ARG[ ARG$X1=="--maxperiod"  ,]$"X2")
	L               = as.character(ARG[ ARG$X1=="--l"         ,]$"X2")
	consensus.size  = as.character(ARG[ ARG$X1=="--consensus.size"  ,]$"X2")
	chr.filt        = as.character(ARG[ ARG$X1=="--chr.filt"   ,]$"X2")
	plot.h          = as.character(ARG[ ARG$X1=="--plot.h"     ,]$"X2")
	plot.w          = as.character(ARG[ ARG$X1=="--plot.w"     ,]$"X2")
	sample.name     = as.character(ARG[ ARG$X1=="--sample.name",]$"X2")
	GeneGFF         = as.character(ARG[ ARG$X1=="--gene.gff"   ,]$"X2")
	polyploid       = as.character(ARG[ ARG$X1=="--polyploid"  ,]$"X2")
	diploid         = as.character(ARG[ ARG$X1=="--diploid"    ,]$"X2")
	add.dots        = as.character(ARG[ ARG$X1=="--add.dots"   ,]$"X2")
	add.segm        = as.character(ARG[ ARG$X1=="--add.segm"   ,]$"X2")
	name.space      = as.character(ARG[ ARG$X1=="--name.space" ,]$"X2")
	leg.size        = as.character(ARG[ ARG$X1=="--leg.size"   ,]$"X2")
	name.del        = as.character(ARG[ ARG$X1=="--name.del"   ,]$"X2")
	run             = as.character(ARG[ ARG$X1=="--run"        ,]$"X2")
	ref.name        = as.character(ARG[ ARG$X1=="--ref.name"   ,]$"X2")
	tel.seq         = as.character(ARG[ ARG$X1=="--tel.seq"    ,]$"X2")
	sort.size       = as.character(ARG[ ARG$X1=="--sort.size"  ,]$"X2")
	
	plot.h     = as.numeric(plot.h)
	plot.w     = as.numeric(plot.w)
	name.space = as.numeric(name.space)
	leg.size   = as.numeric(leg.size)

#  DEFAULT VALUES

    if ( length(Match         )==0  ) Match = 2  
    if ( length(Mismatch      )==0  ) Mismatch = 3  
    if ( length(Delta         )==0  ) Delta = 5
    if ( length(PM            )==0  ) PM = 80
    if ( length(PI            )==0  ) PI = 10 
    if ( length(Minscore      )==0  ) Minscore = 30  
    if ( length(MaxPeriod     )==0  ) MaxPeriod = 2000  
    if ( length(L             )==0  ) L = 6  
    if ( length(consensus.size)==0  ) consensus.size = 7  
    if ( length(chr.filt      )==0  ) chr.filt = 1 * 1e6
    if ( length(TElib)==0           ) TElib=NA  
    if ( length(polyploid)==0       ) polyploid=FALSE  else { polyploid=TRUE }
    if ( length(diploid  )==0       ) diploid  =FALSE  else { diploid  =TRUE }
    if ( length(rDNA.fa)==0         ) rDNA.fa=NA  
    if ( length(rDNA.out)==0        ) rDNA.out=NA  
    if ( length(plot.h)==0          ) plot.h      = 18
    if ( length(plot.w)==0          ) plot.w      = 35
    if ( length(sample.name)==0     ) sample.name = TRUE
    if ( length(GeneGFF)==0         ) GeneGFF     = NA
    if ( length(add.dots)==0        ) add.dots    = NA
    if ( length(add.segm)==0        ) add.segm    = NA
    if ( length(name.space)==0      ) name.space  = 6
    if ( length(leg.size)==0        ) leg.size  = 1
    if ( length(name.del )==0       ) name.del    = FALSE
    if ( length(run )==0            ) run    = "plot"
    if ( length(ref.name )==0       ) ref.name    = FALSE
    if ( ref.name =="--ref.name"    ) ref.name = fasta
    if ( is.na(GeneGFF)==FALSE & GeneGFF=="yes" ) GeneGFF = system("ls *gff3", intern=TRUE) # FIND gene GFF   
    if ( length(tel.seq)==0         ) tel.seq     = NA
    if ( length(sort.size)==0       ) sort.size   = FALSE
    
## if file exists, take its parameters    
    trf.out  =  paste(fasta,".", Match,".", Mismatch,".",Delta,".",PM,".",PI,".",Minscore,".",MaxPeriod,".dat",sep="")
    if ( file.exists(trf.out)== FALSE  & length( system("ls *dat", intern=TRUE))==1   ) # FIND *.dat file
    {
         trf.out = system("ls *dat", intern=TRUE)
         trf.split = strsplit(trf.out,".", fixed=T)[[1]]
         Match     = trf.split [4]  
         Mismatch  = trf.split [5]   
         Delta     = trf.split [6] 
         PM        = trf.split [7] 
         PI        = trf.split [8]  
         Minscore  = trf.split [9]  
         MaxPeriod = trf.split [10]  
        }


# TEST PARAMETERS for interactive session
if( 2 ==1 )
{
    setwd("/ibex/scratch/projects/c2067/celiim/Ajwa_dp/ragtag_scaffold.hifiasm.100Kb.Ajw3")
    
    TElib="/ibex/scratch/projects/c2067/celiim/Ajwa_dp/EDTA/Ajwa_hifiasm.fasta.mod.EDTA.TElib_No_Helitron.plus_rDNA.fa"
#   TElib="../EDTA/Coffea_arabica_hifiasm.all.fasta.mod.EDTA.TElib_No_Helitron.plus_rDNA.fa"

# sbatch -t 20:00:00 --mem=50G  -N 1 -c 64  -J "TELfinder" --wrap 'Rscript $TELfinder  --fasta=Acacia_flavaa_4_subgenomes.fasta --lib=../EDTA/Acacia_flavaa_hifiasm.l2.primary.all.fasta.mod.EDTA.TElib_No_Helitron.plus_rDNA.fa  '

    fasta="ragtag.scaffold.fasta" ; TElib="../../Acacia_flavaa_hifiasm.l2.primary.all.fasta.mod.EDTA.TElib_No_Helitron.plus_rDNA.fa" 
#   TElib="../EDTA/Acacia_flavaa_hifiasm.l2.primary.all.fasta.mod.EDTA.TElib_No_Helitron.plus_rDNA.fa" ; fasta="Acacia_flavaa_4_subgenomes.fasta"
    fasta="4species.fasta" ; TElib="4species_TElib_noHEL_LINEs.edit.sativa_rRNA.fasta" 

    GeneGFF="Barhee_BC4_genes_on_Ajwa_220x_18seq_nogap.gff3"
    GeneGFF = grep("genes", grep("gff3", dir(), value=TRUE), value=TRUE)
    GeneGFF = NA
    Match = 2  
    Mismatch = 3  
    polyploid = NA
    diploid   = NA
    Delta = 5
    PM = 80
    PI = 10 
    Minscore = 30  
    MaxPeriod = 2000  
    L = 6  
    consensus.size = 7  
    plot.h         = 18
    plot.w         = 35
    chr.filt   =  1 * 1e5
    sample.name   =  TRUE
    add.segm    = NA
    add.dots    = NA
    name.space    = 6
    leg.size = 1
    name.del       = FALSE
    run       = "all"
    ref.name= FALSE ;# ref.name = "../GCA_030873655.1_IGA_Cara_2.4_genomic.fna"
    rDNA.fa  = NA
    rDNA.out = NA
    diploid  = FALSE
    polyploid=FALSE
    tel.seq     = NA
    sort.size   = FALSE
 }


cat('\n')
cat("############################################################################\n")
cat('------------------------ Summary of Parameters -----------------------------\n')
cat("############################################################################\n")
cat('\n')
    cat("   directory      = ", getwd(),"\n")
    cat("   fasta          = ", fasta,"\n")
    cat("   TE.library     = ", TElib,"\n")
    cat("   ref.name       = ", ref.name,"\n")
    cat("   rDNA.fa        = ", rDNA.fa,"\n")
    cat("   rDNA.out       = ", rDNA.out,"\n")
    cat("   GeneGFF        = ", GeneGFF,"\n")
    cat("   add.segm       = ", add.segm,"\n")
    cat("   add.dots       = ", add.dots,"\n")
    cat("   tel.seq        = ", tel.seq,"\n")
    cat("   diploid        = ", diploid,"\n")
    cat("   polyploid      = ", polyploid,"\n")
    cat("   chr.filt       = ", chr.filt,"\n")
    cat("   sort.size      = ", sort.size,"\n")
    cat("   sample.name    = ", sample.name,"\n")
    cat("   plot.h         = ", plot.h,"\n")
    cat("   plot.w         = ", plot.w,"\n")
    cat("   name.del       = ", name.del,"\n")
    cat("   name.space     = ", name.space,"\n")
    cat("   run            = ", run,"\n")
    cat("   TRF PARAMETERS: \n")
    cat("   - Match         = ", Match,"\n")
    cat("   - Mismatch      = ", Mismatch,"\n")
    cat("   - Delta         = ", Delta,"\n")
    cat("   - PM            = ", PM,"\n")
    cat("   - PI            = ", PI,"\n")
    cat("   - Minscore      = ", Minscore,"\n")
    cat("   - MaxPeriod     = ", MaxPeriod,"\n")
    cat("   - L             = ", L,"\n")
    cat("   - consensus.size= ", consensus.size,"\n")
cat('\n')
cat('\n')


if( 2 == 1 )
{ 

# TELfinder="/ibex/scratch/projects/c2067/celiim/tools/Telomere_finder.v1.R"


# for other
# sbatch -t 25:00:00 --mem=50G  -N 1 -c 32  -J "TRF"          --wrap ' module load trf          ; trf ragtag.scaffold.fasta 2 3 5 80 10 30 2000 -l 6 -dat -h ' ;
# sbatch -t 25:00:00 --mem=50G  -N 1 -c 32  -J "RepeatMasker" --wrap ' module load repeatmasker ; RepeatMasker -qq -no_is --pa 32 ragtag.scaffold.fasta -lib /ibex/scratch/projects/c2067/celiim/Ajwa_dp/EDTA/Ajwa_hifiasm.fasta.mod.EDTA.TElib_No_Helitron.plus_rDNA.fa -gff '

# Tamarind
# sbatch -t 10:00:00 --mem=50G  -N 1 -c 32  -J "TRF"          --wrap ' module load trf          ; trf ragtag.scaffold.fasta 2 3 5 80 10 30 2000 -l 6 -dat -h '
# sbatch -t 20:00:00 --mem=50G  -N 1 -c 32  -J "RepeatMasker" --wrap ' module load repeatmasker ; RepeatMasker -qq -no_is --pa 32 ragtag.scaffold.fasta -lib ../Tamarind_EDTA_rDNA_lib.fasta -gff '

# Zizi rDNA nly
#   sbatch -t 10:00:00 --mem=50G  -N 1 -c 32  -J "TRF"          --wrap ' module load trf          ; trf ragtag.scaffold.fasta 2 3 5 80 10 30 2000 -l 6 -dat -h '
#   sbatch -t 10:00:00 --mem=50G  -N 1 -c 32  -J "RepeatMasker" --wrap ' module load repeatmasker ; RepeatMasker -qq -no_is --pa 32 ragtag.scaffold.fasta -lib /ibex/scratch/projects/c2067/celiim/Ziziphus_spina-christi/EDTA/Ziziphus_spina-christi_hifiasm.fasta.mod.EDTA.TElib_No_Helitron_plus_rDNA.fa -gff '
#   rDNA only
#   sbatch -t 20:00:00 --mem=50G  -N 1 -c 32  -J "RepeatMasker" --wrap ' module load repeatmasker ; RepeatMasker -qq -no_is --pa 32 ragtag.scaffold.fasta -lib /ibex/scratch/projects/c2067/celiim/Ziziphus_spina-christi/EDTA/Ziziphus_spina-chirsti_rDNA.fasta -gff -dir rDNA_rm '


# sbatch -t 10:00:00 --mem=100G  -N 1 -c 32  -J "TRF"          --wrap ' module load trf          ; trf *.softmasked.fasta 2 3 5 80 10 30 2000 -l 6 -dat -h '
# sbatch -t 10:00:00 --mem=100G  -N 1 -c 32  -J "RepeatMasker" --wrap ' module load repeatmasker ; RepeatMasker -qq -no_is --pa 32 ragtag.scaffold.fasta -lib ../EDTA/Ajwa_hifiasm.fasta.mod.EDTA.TElib_No_Helitron.plus_rDNA.fa -gff '



# sbatch -t 03:00:00 --mem=50G  -N 1 -c 32  -J "RepeatMasker" --wrap ' module load repeatmasker ; RepeatMasker -qq -no_is --pa 32 omicsbox_transcripts.fasta -lib /ibex/scratch/projects/c2067/celiim/Ziziphus_spina-christi/EDTA/Ziziphus_spina-christi_hifiasm.fasta.mod.EDTA.TElib_No_Helitron_plus_rDNA.fa -gff '


# sbatch -t 15:00:00 --mem=50G  -N 1 -c 32  -J "TElfinder.hap1" --wrap ' Rscript $TELfinder --fasta=ragtag.scaffold.fasta --lib=../Tamarind_EDTA_rDNA_lib.fasta --chr.filt=70Kb'

}



# R
suppressWarnings(library(data.table))

options(width=400)
options(scipen=9999) 

base.name=gsub(".fasta$|.fa$|.fna$","",fasta)

# RM
raw.data =  paste(fasta,".", Match,".", Mismatch,".",Delta,".",PM,".",PI,".",Minscore,".",MaxPeriod,".dat",sep="")
rep.data =  paste(fasta,".out",sep="")
ref.data =  paste(fasta,".out.gff",sep="")
rem.data =  paste(fasta,".out.merged.bed",sep="")
rei.data =  paste(fasta,".out.intersect.bed",sep="")

rem.dataTE =  paste(fasta,".out.TE.merged.bed",sep="")
rei.dataTE =  paste(fasta,".out.TE.intersect.bed",sep="")
rem.dataSI =  paste(fasta,".out.SIMPLE.merged.bed",sep="")
rei.dataSI =  paste(fasta,".out.SIMPLE.intersect.bed",sep="")
rem.dataRD =  paste(fasta,".out.rDNA.merged.bed",sep="")
rei.dataRD =  paste(fasta,".out.rDNA.intersect.bed",sep="")
rem.dataGE =  paste(fasta,".out.genes.merged.bed",sep="")
rei.dataGE =  paste(fasta,".out.genes.intersect.bed",sep="")

# TDR
trf.out  =  paste(fasta,".", Match,".", Mismatch,".",Delta,".",PM,".",PI,".",Minscore,".",MaxPeriod,".dat",sep="")
pro.data =  paste(raw.data,".with_chr.txt",sep="")
mer.data =  paste(raw.data,".merged.txt",sep="")
int.data =  paste(raw.data,".intersect.txt",sep="")

#AGP
agp.data = gsub(".fasta$|.fa$|.fna$",".agp",fasta)
ass.data = paste(base.name,".split_assembly.fasta",sep="")
 
rag=NA

# CREATE AGP if not present using Ragtag

cat ("\n\n")
if ( file.exists(agp.data                   )== FALSE ) 
{ cat (" --> no AGP file with fasta file base name \n ") ;
  system("module load ragtag; which ragtag.py",intern=T) -> rag ; 
  RAGTAG_SPLITASM = paste( "ragtag.py splitasm ",fasta," -o ",agp.data," > ",ass.data )
  if (is.na(rag)==TRUE ) cat("--> you need AGP files for locating GAPS, e.g. create it with: ",RAGTAG_SPLITASM,"  \n") ;
  if (is.na(rag)==FALSE) cat("--> Creating AGP files with: ",RAGTAG_SPLITASM," \n") ; system ( paste(RAGTAG_SPLITASM, "\n"))  ;  system (paste("rm ",ass.data, "\n") ) ; 
  if ( file.exists(agp.data) ) { cat (" --> AGP file present: ", agp.data,"\n") }
} else { cat (" --> AGP file present: ", agp.data,"\n") }

# check presence of other files
#if ( file.exists(paste(fasta,".fai",sep="") )== FALSE ) { cat (" --> generating FAIDX           \n\n");  system( paste("samtools faidx", fasta ) ) } else { cat (" --> FAIDX file present \n") }
cat (" --> generating FAIDX           \n\n");  system( paste("samtools faidx", fasta ) ) 
if ( file.exists(rep.data                   )== FALSE ) { cat (" --> running RepeatMasker       \n\n");  system( paste("RepeatMasker -qq -no_is --pa 32 -lib",TElib,fasta, "-gff -xsmall" ) )    } else { cat (" --> REPEAT MASKER output present \n") }
if ( file.exists(trf.out                    )== FALSE ) { cat (" --> running TRF                \n\n");  system( paste("trf", fasta, Match , Mismatch ,Delta ,PM,PI,Minscore,MaxPeriod,"-l",L,"-dat -h" ) ) } else { cat (" --> TRF output present \n") }
if ( file.exists(pro.data                   )== FALSE ) { cat (" --> TRF output to be processed \n\n");                                                                                                                    } else { cat (" --> TRF output processed present \n") }
cat ("\n\n")

# Chr filtering 
chr.filt.char = chr.filt
if ( grepl("Kb|Mb|Gb",chr.filt.char) ==TRUE ) 
{   chr.filt = as.numeric(gsub("bp","",(gsub("Kb","000",gsub("Mb","000000",gsub("Gb","000000000",chr.filt))))))
  } else {   
    chr.filt.char = gsub("0000000","0Mb",chr.filt)
    chr.filt.char = gsub("000000" , "Mb",chr.filt.char)
    chr.filt.char = gsub("00000" ,"00Kb",chr.filt.char)
    chr.filt.char = gsub("0000"   ,"0Kb",chr.filt.char)
    chr.filt.char = gsub("000"     ,"Kb",chr.filt.char)
    chr.filt = as.numeric(chr.filt)
    }
    

  
# import FAI, AGP, repeatmasker
chr = fread( paste(fasta,".fai",sep=""),fill=TRUE,header=F)[,1:2]
#rep = fread( paste(rep.data   ,sep=""),fill=TRUE,header=F)
rep = fread(cmd=paste("grep 'C\\|' ", rep.data, " | sed 's/ [*]//g'"))


agp = setDT(read.table( agp.data ))
names(agp)=c("ref.chr","ref.start","ref.end","part_number","component_type","query","q.start","q.end","strand")
names(chr)=c("chr","chr.len")
chr = data.frame(chr)
input.bed = gsub(".fasta$|.fa$|.fna$",".bed",fasta)
chh =  data.frame("chr"=chr$"chr", "start"=1, end=chr$"chr.len")
write.table( chh, file =input.bed , col.names=F, row.names=F, quote=F, sep="\t")

# sorting chromoomes
chr$"chr" = as.character(chr$"chr")
chr = chr[ order(chr$"chr" ) , ]

# Sorting chromosomes by name, then contings by size, if present
if( length(grep("CHR|Chr|LG|CM0|C|RagTag",chr$"chr"))>0)
{
    ctg = chr[ -grep("CHR|Chr|LG|CM0|C|RagTag",chr$"chr") , ] # Identify Chromosomes or Linkage Groups or 
    chh = chr[  grep("CHR|Chr|LG|CM0|C|RagTag",chr$"chr") , ]
    ctg = ctg[ rev(order(ctg$"chr.len" )) , ]
    
    # Filt by CHR/contig size
    chr = rbind(chh,ctg)
    }

CHR = chr
chr = CHR[ which(CHR$"chr.len" > chr.filt) , ]
# CHR --> full     CHR/contig dataset
# chr --> filtered CHR/contig dataset
CHR.bbb = CHR
CHR.bbb$"chr" = gsub("_RagTag","",CHR.bbb$"chr" )



agp$"ref.chr" = gsub("_RagTag","",agp$"ref.chr"  )
gap = agp[ agp$"strand" =="align_genus",]
agp = agp[ agp$"strand" !="align_genus",]


GAP = gap
AGP = agp

#cat ("\n\n")

cat ("\n --> Chromosomes/Contig names and length \n\n")
print(chr)
# print(ctg)
# print(chh)
cat ("\n --> AGP Contigs \n\n")
print(agp)
cat ("\n --> GAP positions \n\n")
print(gap)
cat ("\n\n")

# import gene gff if present
if ( is.na(GeneGFF)==FALSE ) { gene = fread( GeneGFF,fill=TRUE,header=F)[,1:2]  ; if ( file.exists(GeneGFF))  cat (" --> GENE GFF3 output present \n") else cat (" --> GENE GFF3 output missing \n") }
if ( is.na(GeneGFF)==TRUE  ) {  cat (" --> no GENE GFF3 output loaded \n")  }
cat ("\n")

if ( ref.name != FALSE )
{
cat("\n\n")
cat("   ==> Chromosome complete names:\n\n")

# input names
fas.name = gsub("_RagTag","",chr$"chr")

# select only input fasta names in the reference
ref.names = system( paste(' grep ">" ',ref.name), inter=TRUE) 
ref.names =  grep(   paste(fas.name,collapse="|") , ref.names , value = TRUE )


ref.df= data.frame( name=ref.names, stringsAsFactors=FALSE     )

print(ref.df)

 # remove common words
ref.words = strsplit(ref.names, " ")
ref.all.words = sort(unique(unlist(ref.words)))

 #count words occurencies for each ">" line
 REF.ALL.WORDS = list()
 for (i in ref.all.words) REF.ALL.WORDS[[i]] = sapply(ref.words, function(x) { sum(x==i) })
 
 REF.ALL.WORDS.sin = REF.ALL.WORDS [ which(sapply(REF.ALL.WORDS, function(x) ! all(x==1))) ] 
 REF.ALL.WORDS.dup = REF.ALL.WORDS [ which(sapply(REF.ALL.WORDS, function(x)   all(x>=1))) ] 
 ref.dup = names(REF.ALL.WORDS.dup)
 ref.dup = grep("Chr|chr|LG|CTG",ref.dup, invert=T, value = TRUE)
 ref.dup = paste(ref.dup, collapse= "|")

ref.df$"simp_name" = ref.df$"name"
ref.df$"simp_name" = gsub(ref.dup, "", ref.df$"simp_name" )
ref.df$"simp_name" = gsub("       |     |   |  |,|;", "  ", ref.df$"simp_name" )

 ref.df$"simp_name" = gsub("whole|shotgun|sequence|,|;", "  ", ref.df$"simp_name" )
 ref.df$"simp_name" = gsub("\\s{2,}", "  ", ref.df$"simp_name")

ref.edit = suppressWarnings( data.frame(do.call(rbind, strsplit(ref.df$"simp_name","  " ))))
for ( ii in names(ref.edit)) if ( all( as.character(ref.edit[,ii])==" ") | all(as.character(ref.edit[,ii])=="")) ref.edit[ii] = NULL


CHR_PREFIX = paste(fas.name, collapse="|")

w1.col = names( which(sapply(ref.edit,function(x) length(grep(CHR_PREFIX,x))) == nrow(ref.df)))
w2.col = names( which(sapply(ref.edit,function(x) length(grep("Chr|chr|LG",x))) > 0))

name.ref.edit = names(ref.edit)
name.ref.edit[ name.ref.edit == w1.col] = "chr"
name.ref.edit[ name.ref.edit == w2.col] = "chr.rename"
names(ref.edit) = name.ref.edit
ref.edit$"chr"        = gsub(">","", ref.edit$"chr")
ref.edit$"chr.rename" = gsub("Chromosome |Chromosome|chromosome |chromosome","Chr", ref.edit$"chr.rename")
ref.edit$"chr.rename" = gsub("Chr: ","Chr", ref.edit$"chr.rename")
ref.edit$"chr.rename" = gsub(": ",":", ref.edit$"chr.rename")
ref.edit$"chr.rename" = gsub("contig:","", ref.edit$"chr.rename")

w_chl = grep("plastid|chloroplast",ref.edit$"chr.rename",ignore.case = TRUE)
w_mit = grep("mitochondrion"      ,ref.edit$"chr.rename",ignore.case = TRUE)

grep("plastid|chloroplast",ref.edit$"chr.rename",ignore.case = TRUE)
grep("mitochondrion"      ,ref.edit$"chr.rename",ignore.case = TRUE)
 
if (length(w_chl) == 1 ) ref.edit$"chr.rename"[w_chl]="ChrC"
if (length(w_mit) == 1 ) ref.edit$"chr.rename"[w_mit]="ChrM"

ref.edit$"nchar" = nchar(ref.edit$"chr.rename")
chr_lines = grep("Chr[1-9]", ref.edit$"chr.rename")
min_char = min(ref.edit$"nchar"[chr_lines] )
max_char = max(ref.edit$"nchar"[chr_lines ])
if( min_char == max_char-1  ) ref.edit[ ref.edit$"nchar" == min_char & grepl("Chr[1-9]", ref.edit$"chr.rename"), ]$"chr.rename" = gsub("Chr","Chr0",ref.edit[ ref.edit$"nchar" == min_char & grepl("Chr[1-9]", ref.edit$"chr.rename"), ]$"chr.rename")

ref.edit$"nchar" = NULL
ref.edit.full = ref.edit 
ref.edit = ref.edit [ , 1:2]

names(ref.edit) = c("chr","chr.rename")
ref.edit$"chr" = gsub(">","",ref.edit$"chr")

# Contigs named "ptg" may create confusion with unplaced contigs from ragtag, so keep the NCBI name or add ">"
w_ptg = grep("ptg|h1tg|h2tg|tig",ref.edit$"chr.rename"  )
 
if ( length(w_ptg)>0   ) {  ref.edit$"chr.rename"[w_ptg] = ref.edit$"chr"[w_ptg]                    ; cat("\n\n   --> unplaced contig names not replaced:         \n"); print(ref.edit[w_ptg,])   }
#if ( length(w_ptg)>0   ) {  ref.edit$"chr.rename"[w_ptg] = paste0(">",ref.edit$"chr.rename"[w_ptg]) ; cat("   --> unplaced contig names replaced with extra >:\n"); print(ref.edit[w_ptg,])   }


cat("\n\n")
  
cat("   --> New Chromosome Names:\n")    
print(head(ref.edit,50))
cat("\n\n")
}

#
TElib.fai = paste0(TElib,".fai")

if ( file.exists(TElib.fai)== FALSE ) { cat (" --> generating FAIDX  TElibrary         \n");  system( paste("samtools faidx", TElib ) ) } else { cat (" --> FAIDX TElibrary file present \n") }

cat ("\n --> TE library \n\n")

# print(TElib)
# print(TElib.fai)


# IMPORT TE and rDNA LIBRARY if present

seq.len = fread( TElib.fai)[,1:2]
names(seq.len) = c("repeat.class","rep.len")

if ( is.na(rDNA.fa)==FALSE ) 
{
    rDNA.fai = paste0(rDNA.fa,".fai")
    if ( file.exists(rDNA.fai)== FALSE ) { cat (" --> generating FAIDX  rDNA library        \n\n");  system( paste("samtools faidx", rDNA.fa ) ) } else { cat (" --> FAIDX rDNA library file present \n\n") }
    rDNA.len = fread( rDNA.fai)[,1:2]
    names(rDNA.len) = c("repeat.class","rep.len")
    
    seq.len = seq.len[ ! seq.len$"repeat.class" %in% rDNA.len$"repeat.class", ]
    seq.len = rbind(seq.len, rDNA.len)
   cat (" --> rDNA library added to TE library      \n\n")

}


rDNA.presence = grep("rRNA|rDNA|ribos|Ribos", seq.len$"repeat.class")
print(seq.len)

rDNA.len=seq.len[rDNA.presence, ]


if ( length(rDNA.presence)>0){ cat("\n   --> rDNA repeats detected \n\n"); print( rDNA.len )} else  cat("\n   --> no rDNA repeats detected \n\n")


# Processing TRF and Repealt Masker output, create 100Kb windows profile
cat('\n\n')
cat("############################################################################\n")
cat('--------------- TRF and Repeat Masker outputs processing -------------------\n')
cat("############################################################################\n")
cat('\n')


if ( file.exists(pro.data)==TRUE & file.exists(rei.data)==TRUE & run=="plot" ) {  cat("\n   --> TRF processed files present:",pro.data," \n\n");
    } else {
   cat (" --> Loading TRF data      \n\n")

        trf = fread( paste(fasta,".", Match,".", Mismatch,".",Delta,".",PM,".",PI,".",Minscore,".",MaxPeriod,".dat",sep=""),fill=TRUE,header=F)
        
        names(trf)=c("start","end","period","copy.n","consensus.size","perc.match","prc.indels","score","A","C","G","T","entropy","consensus.seq","tot.seq")
        trf[ 1:3, 1:14]
        # line with chr
        # line with data
        
        trf = trf [   trf$"start" != "", ]
        trf = trf [ - grep("Tandem|Gary|Benson|Program|Boston|Version|Parameters", trf$"start"), ]
        trf[ 1:3, 1:14]
              
        new.names= c("chr", names(trf) )
        trf$"chr" = "-"
        
        trf = trf [,..new.names]
        w.chrom = which(trf$"start" == "Sequence:")
        trf[w.chrom,]
        
   cat (" --> Adding Chromosome/Contig columns to files     \n\n")
        www = c(w.chrom, nrow(trf))
        cat(" --> N° of chr/contigs: ",length(www)-1,"\n")
        cat ("     ")
        for ( i in 2:length(www))  { cat( i-1,""); a=www[i-1] ;  b=www[i] ; chrom = trf$"end"[a] ;  trf$"chr"[a:b] = chrom }
   cat ("\n\n")
   cat (" --> Merging Data     \n\n")

        trf[ 1:3, 1:14]
        trf = trf [,..new.names]
        
        # check 
#         print("check 1 \n" )
        
        ts = split(trf,trf$"chr")
        
        aa = sapply(ts, function(x) as.character(x$"end"[1]  ))
        bb = sapply(ts, function(x) unique(x$"chr" ))
        
        all(aa ==bb )
        trf = trf[- w.chrom,] 
        trf = trf[ ! duplicated(trf),] 
#           print("check 2 \n" )
             
        trf$"end"   = as.numeric(as.character( trf$"end" ))
        trf$"start" = as.numeric(as.character( trf$"start" ))
        trf$"period" = as.numeric(as.character( trf$"period" ))
        trf$"copy.n" = as.numeric(as.character( trf$"copy.n" ))
        trf$"consensus.size" = as.numeric(as.character( trf$"consensus.size" ))
        trf$"perc.match" = as.numeric(as.character( trf$"perc.match" ))
        trf$"prc.indels" = as.numeric(as.character( trf$"prc.indels" ))
        trf$"score" = as.numeric(as.character( trf$"score" ))
        trf$"A" = as.numeric(as.character( trf$"A" ))
        trf$"C" = as.numeric(as.character( trf$"C" ))
        trf$"G" = as.numeric(as.character( trf$"G" ))
        trf$"T" = as.numeric(as.character( trf$"T" ))
        trf$"entropy" = as.numeric(as.character( trf$"entropy" ))
        
        trf$"len" = trf$"end" -trf$"start"  +1
        trf$"rep.tag" =paste0("rep.",1:nrow(trf))
        
   cat (" --> Saving TRF processed output    \n\n")
        trf = trf[ order(trf$"chr", trf$"start") , ]
        write.table( trf, file = pro.data , col.names=TRUE, row.names=F, quote=F, sep="\t")    
#         print(head(trf))
#         print("check 3 \n" )
#         
        #base.name=gsub(".fasta|.fa|.fna","",fasta)

# raw.data =  paste(fasta,".", Match,".", Mismatch,".",Delta,".",PM,".",PI,".",Minscore,".",MaxPeriod,".dat",sep="")
# ref.data =  paste(fasta,".out.gff",sep="")
# rep.data =  paste(fasta,".out",sep="")
# pro.data =  paste(raw.data,".with_chr.txt",sep="")
# mer.data =  paste(raw.data,".merged.txt",sep="")
# int.data =  paste(raw.data,".intersect.txt",sep="")
#    
#         
   cat (" --> Computing Repeat density per 100kb window    \n\n")
        system( paste0( " grep -v entropy ",pro.data," >", pro.data ,".temp.txt" ) )
        system( paste0( " sort -k1,1 -k4,4n  ",ref.data," >", ref.data ,".temp.gff" ) )        
        system( paste0( " grep    'Motif:TE_\\|rDNA\\|ribosom' ", ref.data ,".temp.gff" , " > ",ref.data ,".TE.temp.gff" ) )        
        system( paste0( " grep -v 'Motif:TE_\\|rDNA\\|ribosom' ", ref.data ,".temp.gff" , " > ",ref.data ,".SIMPLE.temp.gff" ) )        
        system( paste0( " grep    'rDNA\\|ribosom'             ", ref.data ,".temp.gff" , " > ",ref.data ,".rDNA.temp.gff" ) )        
        
        system( paste0( " bedtools merge -i ",pro.data,".temp.txt > "       ,mer.data) )
        system( paste0( " bedtools merge -i ",ref.data,".temp.gff > "       ,rem.data) )
        system( paste0( " bedtools merge -i ",ref.data,".TE.temp.gff > "    ,rem.dataTE) )
        system( paste0( " bedtools merge -i ",ref.data,".SIMPLE.temp.gff > ",rem.dataSI) )
        system( paste0( " bedtools merge -i ",ref.data,".rDNA.temp.gff > "  ,rem.dataRD) )
            
        system( paste0( " bedtools makewindows -b ",base.name,".bed  -w 100000 -s 100000 > ",base.name,".100Kb_windows.bed"))
        system( paste0( " bedtools intersect -wao -a ",base.name,".100Kb_windows.bed -b ",mer.data, " > ",int.data) )
        print ( paste0( " bedtools intersect -wao -a ",base.name,".100Kb_windows.bed -b ",mer.data, " > ",int.data) )
#         system( paste0( " wc -l  ",pro.data) )
#         system( paste0( " head  ",pro.data) )
#         system( paste0( " wc -l  ",base.name,".100Kb_windows.bed")) 
#         system( paste0( " head  ",base.name,".100Kb_windows.bed")) 
#         system( paste0( " wc -l  ",int.data) )
#         system( paste0( " head  ",int.data) )
#         system( paste0( " wc -l  ",mer.data) )
#         system( paste0( " head  ",mer.data) )
        system( paste0( " bedtools intersect -wao -a ",base.name,".100Kb_windows.bed -b ",rem.data, " > ",rei.data) )
        system( paste0( " bedtools intersect -wao -a ",base.name,".100Kb_windows.bed -b ",rem.dataSI, " > ",rei.dataSI) )
        system( paste0( " bedtools intersect -wao -a ",base.name,".100Kb_windows.bed -b ",rem.dataRD, " > ",rei.dataRD) )
        system( paste0( " bedtools intersect -wao -a ",base.name,".100Kb_windows.bed -b ",rem.dataTE, " > ",rei.dataTE) )
}

# Create Genes 100Kb windows profile 


if ( is.na(GeneGFF)==FALSE ) 

      {
        # to be fixed
#     cat("test gene bed \n")
#          paste0( " bedtools merge -i <( cut -f1-5 ",GeneGFF, " | grep Liftoff | grep gene |  cut -f1,4,5) > ",rem.dataGE)
#          paste0( " bedtools intersect -wao -a <( sed 's/_RagTag//g' ",base.name,".100Kb_windows.bed )  -b  ",rem.dataGE, " > ",rei.dataGE)
#         system( paste0( " bedtools merge -i <( cut -f1-5 ",GeneGFF, " | grep Liftoff | grep gene |  cut -f1,4,5) > ",rem.dataGE) )
#         system( paste0( " bedtools intersect -wao -a <( sed 's/_RagTag//g' ",base.name,".100Kb_windows.bed )  -b  ",rem.dataGE, " > ",rei.dataGE) )
         cat("Generating gene density .bed file \n")
        rem.dataGE_input = gsub("merged","edit",rem.dataGE)
        system( paste0( " bedtools makewindows -b ",base.name,".bed  -w 100000 -s 100000 > ",base.name,".100Kb_windows.bed"))
        system( paste0( " cut -f1-5 ",GeneGFF, " | grep -v '#' | grep gene |  cut -f1,4,5 > ",rem.dataGE_input) )                 
        system( paste0( " sed 's/_RagTag//g' ",base.name,".100Kb_windows.bed > ",base.name,".100Kb_windows.edit.bed "))      
        system( paste0( " bedtools intersect -wao -a ",base.name,".100Kb_windows.edit.bed  -b  ",rem.dataGE_input, " > ",rei.dataGE) )

        }




chr$"chr" = gsub("_RagTag","",chr$"chr"  )
chh$"chr" = gsub("_RagTag","",chh$"chr"  )
CHR$"chr" = gsub("_RagTag","",CHR$"chr"  )
agp$"chr" = gsub("_RagTag","",agp$"ref.chr"  )
AGP$"chr" = gsub("_RagTag","",AGP$"ref.chr"  )
gap$"chr" = gsub("_RagTag","",gap$"ref.chr"  )
GAP$"chr" = gsub("_RagTag","",GAP$"ref.chr"  )





# CHR renaming 
if ( ref.name != FALSE) 
{

cat("\n")
cat("############################################################################\n")
cat(" ==> replacing GenBank indentifier with Chromosomename from the fasta header  \n")
cat("############################################################################\n")
cat("\n")

    chr = merge(chr,ref.edit, by="chr",                   all.x=T) ; w_ren = which(! is.na(chr$"chr.rename")) ;   chr$"chr"[w_ren]     =chr$"chr.rename"[w_ren];  chr$"chr.rename"=NULL 
    chh = merge(chh,ref.edit, by="chr",                   all.x=T) ; w_ren = which(! is.na(chh$"chr.rename")) ;   chh$"chr"[w_ren]     =chh$"chr.rename"[w_ren];  chh$"chr.rename"=NULL 
    CHR = merge(CHR,ref.edit, by="chr",                   all.x=T) ; w_ren = which(! is.na(CHR$"chr.rename")) ;   CHR$"chr"[w_ren]     =CHR$"chr.rename"[w_ren];  CHR$"chr.rename"=NULL 
    agp = merge(agp,ref.edit, by.x="ref.chr", by.y="chr", all.x=T) ; w_ren = which(! is.na(agp$"chr.rename")) ;   agp$"ref.chr"[w_ren] =agp$"chr.rename"[w_ren];  agp$"chr.rename"=NULL 
    AGP = merge(AGP,ref.edit, by.x="ref.chr", by.y="chr", all.x=T) ; w_ren = which(! is.na(AGP$"chr.rename")) ;   AGP$"ref.chr"[w_ren] =AGP$"chr.rename"[w_ren];  AGP$"chr.rename"=NULL 
    gap = merge(gap,ref.edit, by.x="ref.chr", by.y="chr", all.x=T) ; w_ren = which(! is.na(gap$"chr.rename")) ;   gap$"ref.chr"[w_ren] =gap$"chr.rename"[w_ren];  gap$"chr.rename"=NULL 
    GAP = merge(GAP,ref.edit, by.x="ref.chr", by.y="chr", all.x=T) ; w_ren = which(! is.na(GAP$"chr.rename")) ;   GAP$"ref.chr"[w_ren] =GAP$"chr.rename"[w_ren];  GAP$"chr.rename"=NULL 
      }






cat('\n')
cat("############################################################################\n")
cat('----------------------- Chromosomes/Contig recap ---------------------------\n')
cat("############################################################################\n")
cat('\n')
    
# Chr disposition: haploid or polyploid 
# create matrix with CHR/ctg names
      
mat = chr$"chr.len"
names(mat) = chr$chr
names(mat) = gsub("_RagTag","",names(mat) )
mat = mat [ rev(names(mat)) ]

print(1)

## sort bu size
if(sort.size==TRUE) mat = sort(mat)


# Sorting first CHR, then plastidian, mitochondrion and unplaced contigs if present
#print(mat)

# Chloroplast, Mitoichondrion and small contigs at the end
if ( any(grepl("contig|Contig|ChrC|Chl|ChrM|Mit|hloroplast|itochondrion",names(mat)) )==TRUE )
{  
    main_chr= grep("contig|Contig|ChrC|Chl|ChrM|Mit|hloroplast|itochondrion",names(mat), invert=TRUE) 
    smal_chr= grep("contig|Contig|ChrC|Chl|ChrM|Mit|hloroplast|itochondrion",names(mat)) 
    
    sort_chr = c(smal_chr,main_chr)
    mat = mat [ sort_chr ]

}

space    = 0.2
mat.spac = 0.2

# Chr disposition: 
#   haploid --> all CHR/ctg separate
#   polyploid --> CHR paired, plastif, mitochondrion, ctg separate

print(2)


if ( diploid==TRUE | polyploid==TRUE   )
{  
   
   chr_names= data.frame(chr=names(mat), stringsAsFactors=FALSE)
   
print(3)
    # mt and chr present?
    w_org = grep("Mit|Chl|Plast|ChrC|ChrM",chr_names$"chr" )
    if(length(w_org)>0) { org_names = data.frame(chr= chr_names[w_org, ], stringsAsFactors=FALSE) ;   chr_names = data.frame(chr= chr_names[-w_org, ], stringsAsFactors=FALSE) }

print(4)
    
    # Homologous chr separated by "." "-" e.g Chr01.1 or not (eg Chr01a)
    #  
    if( all(grepl( "[.]|-|_",chr_names$chr)) ==TRUE )
    { 
      name_split= strsplit(  gsub("_",".",chr_names$"chr"),".", fixed=TRUE)
      x_names =  data.frame( do.call(rbind , strsplit(  gsub("_",".",chr_names$"chr"),".", fixed=TRUE)), stringsAsFactors=FALSE)
      names(x_names) = gsub("X","x", names(x_names))
      chr_names = cbind(chr_names,x_names)
      ww_ctg = grep("chr|Chr|CHR", chr_names$"chr", invert=T)
      ctg_names =  data.frame(chr= chr_names[ww_ctg,] , stringsAsFactors=FALSE)
print(5)
   
      } else {
print(6)
         
          ww_chr = grep("chr|Chr|CHR", chr_names$"chr")
          ww_ctg = grep("chr|Chr|CHR", chr_names$"chr", invert=T)
          ctg_names =  data.frame(chr= chr_names[ww_ctg,] , stringsAsFactors=FALSE)
          chr_names =  data.frame(chr= chr_names[ww_chr,] , stringsAsFactors=FALSE)
          Nchar = nchar(chr_names$"chr")                                # Nchar Chr names
          base  = unique( substr(chr_names$"chr",1,3))                  # First 3 character of Chr names
          sub  = unique( substr(chr_names$"chr",Nchar,Nchar))           # uniq 3-char names
#          num  = as.numeric( substr(chr_names$"chr",Nchar-2,Nchar-1))   # Chr numeric part
#          chrom  =           substr(chr_names$"chr",Nchar-5,Nchar-1)
          num  = as.numeric( substr(chr_names$"chr",4,5) )  # Chr numeric part
          chrom  =           substr(chr_names$"chr",1,5)
print(7)
          
          chr_names$"x1" = chrom
          chr_names$"x2" = sub 
          
       #   ctg_names$"x1" = substr(ctg_names$"chr",1,3)
       #  ctg_names$"x2" = "ctg" 
                  
          
      }
print(8)
    # N of fields
    n_fields = ncol(chr_names)
    if(n_fields>2 ) { chr_names$"x2" =  apply(chr_names[,-1],1,paste,collapse=".") ; chr_names=chr_names[,1:2]}       
    
    # add Mit and Chloroplast
    if(length(w_org)>0)   { org_names$"x1" = org_names$"chr"; org_names$"x2" = 1; }
    if(length(w_org)>0)   { chr_names = rbind(org_names,chr_names); }

    # identiry chr-pair base name
    if( any(grepl("Chr|chr|tig",chr_names$"x1" ) ) == TRUE )  chr_names$"pair" =  chr_names$"x1"
    if( any(grepl("Chr|chr|tig",chr_names$"x2" ) ) == TRUE )  chr_names$"pair" =  chr_names$"x2"
print(9)
print(chr_names)
print(9.1)

    chr_pairs = as.data.frame(table(chr_names$"pair"))
    chr_names$"x1" =   chr_names$"x2"  = NULL
    chr_paired = as.character(chr_pairs[ which(chr_pairs$"Freq">=2) ,]$"Var1")
    chr_unpair = as.character(chr_pairs[ which(chr_pairs$"Freq"==1) ,]$"Var1")
    chr_names$"type" =NA
print(9.2)
    if(length(chr_paired)>0 ) chr_names[  chr_names$"pair" %in% chr_paired, ]$"type" = "polyploid"
    if(length(chr_unpair)>0 ) chr_names[  chr_names$"pair" %in% chr_unpair, ]$"type" = "haploid"

    chr_names = rbind(chr_names[  chr_names$"type"  =="haploid",  ] , chr_names[ chr_names$"type"  =="polyploid", ] )

print(10)
    # add ploidy > 4?
    names(chr_pairs)=c("pair","n")
    chr_names = merge(chr_names,chr_pairs, by="pair", all.x=T)
    # add space
    chr_names$"barplot.space"=1
    for( ii in 2:nrow(chr_names)) if( chr_names$"pair"[ii]==chr_names$"pair"[ii-1] ) { chr_names$"barplot.space"[ii]=0 }
    
    rownames(chr_names) = chr_names$"chr"
 print(11)
   
    # add ctgs
    if (nrow(ctg_names) > 0)
    {
        ctg_names$"pair" = ctg_names$"chr" 
        ctg_names$"type" = "haploid" 
        ctg_names$"n" = 1
        ctg_names$"barplot.space" = 1
        rownames(ctg_names) = ctg_names$"chr"
        ctg_names = ctg_names[ , names(chr_names)]
        chr_names = rbind(chr_names,ctg_names)
        }
print(12)
    
    mat.spac = chr_names[  names(mat) , ]$"barplot.space"
    if (any(chr_names$"n" ==2)) chr_names[  chr_names$"n" ==2 , ]$"type" = "diploid"    
    rownames(chr_names) = NULL
    cat('\n ---> Chromosome grouping \n\n')
    print(chr_names)
       
    } else  print(as.data.frame(mat))
    
    
print(13)
    


#print(mat)
par(mfcol=c(1    ,1    ), mar=c(2,7,4,1) , oma=c(1,4,1,1) )
bar = barplot(mat, horiz=T, beside=T, las=1, col="white",xlim=c(-max(mat)*0.07,max(mat)*1.15), main="Telomeres and Tandem Repeat profile",space=mat.spac )
bbb = data.frame( chr=names(mat), "chr.len"=mat, pos = bar[,1])



cat('\n\n')
cat("############################################################################\n")
cat('------------------------ Output file suffix----------------------------\n')
cat("############################################################################\n")



chr.filt.char = gsub("0000000","0Mb",chr.filt)
chr.filt.char = gsub("000000" , "Mb",chr.filt)
chr.filt.char = gsub("00000" ,"00Kb",chr.filt.char)
chr.filt.char = gsub("0000"   ,"0Kb",chr.filt.char)
chr.filt.char = gsub("000"     ,"Kb",chr.filt.char)

if (polyploid==TRUE) chr.filt.char = paste0(chr.filt.char,".polyploid")
if (diploid  ==TRUE) chr.filt.char = paste0(chr.filt.char,".diploid")

# add sample names and haplotype
DIR = getwd()
DIR_split = unlist(strsplit(DIR,"/"))
DIR_sampl = DIR_split[ length(DIR_split)-1  ]
DIR_other = DIR_split[ length(DIR_split)    ] ; 
DIR_other= grep("hap|HAP|Hap", unlist(strsplit(DIR_other,".", fixed=TRUE)), value=TRUE)
if (length(DIR_other)==1) DIR_sampl = paste(DIR_sampl,DIR_other,sep=".")

    # FILE SUFFIX 
MAIN.NAME="" 
if ( sample.name== TRUE              ) chr.filt.char = paste(DIR_sampl  ,chr.filt.char,sep=".")
if ( is.character(sample.name)==TRUE ) chr.filt.char = paste(sample.name,chr.filt.char,sep=".")
if ( is.character(sample.name)==TRUE ) MAIN.NAME = sample.name

cat(' \n --->   FILE SUFFIX =',chr.filt.char,' \n')
cat(' \n --->   MAIN NAME   =',MAIN.NAME,' \n')

RRR =""
if ( ref.name != FALSE) RRR=".RefName"

DIR_out=paste0("TELfinder.",chr.filt.char,RRR,"/")
cat(' --->   Output Folder =',DIR_out,' \n\n')
dir.create(DIR_out)


cat('\n\n')
cat("############################################################################\n")
cat('------------------------ Telomere identification----------------------------\n')
cat("############################################################################\n")
cat('\n')

# Load TRF
trf = fread( pro.data ,fill=TRUE,header=T)
trf = trf[ ! duplicated(trf),] 
trf$"chr" =  gsub("_RagTag","", trf$"chr" )

#print(trf)
# rename
if ( ref.name != FALSE)  {  trf = merge(trf,ref.edit, by="chr", all.x=T) ; w_ren = which(! is.na(trf$"chr.rename")) ;   trf$"chr"[w_ren] =trf$"chr.rename"[w_ren];  trf$"chr.rename"=NULL }
#print(trf)

trf$"chr" =  gsub("_RagTag","", trf$"chr" )


# Load TRF data by 100Kb windows

int = fread( int.data,fill=TRUE,header=F)
names(int)=c("chr","start","end","chrom","s","e","n")
int$"chr"   = gsub("_RagTag","",int$"chr"  )
int$"chrom" = gsub("_RagTag","",int$"chrom"  )
int$"wind" =paste0(int$"chr",":",int$"start","-",int$"end")
cols.fil2 = setdiff( names(trf) , c("tot.seq","consensus.seq"))
cols.filt = setdiff( names(trf) , c("tot.seq"))
#trf[,..cols.filt]

# rename
if ( ref.name != FALSE)  {  int = merge(int,ref.edit, by="chr", all.x=T) ; w_ren = which(! is.na(int$"chr.rename")) ;   int$"chr"[w_ren] =int$"chr.rename"[w_ren];  int$"chr.rename"=NULL }


cat('\n --->  TRF data per 100Kb windows: \n')

trf = merge(trf, CHR, by="chr")

ts = split(trf,trf$"chr")
is = split(int,int$"wind")

is2 = lapply(is, function(x) { y=data.frame(chr=x$"chr"[1],start=x$"start"[1],end=x$"end"[1],n=sum(x$"n"))    ;  setDT(y); return(y) } )
is2 = do.call(rbind,is2)
is2 = is2[ order(is2$"chr", is2$"start") , ]

#last window is not 100Kb, so normalize
is2$"perc.rep" = is2$"n" / (is2$"end" -is2$"start"+1)

#print(is2)

is3 = split(is2,is2$"chr")
is4 = lapply(is3, function(x) { y=data.frame( chr=x$"chr"[1], len=max(x$"end"), rep=sum(x$"n"), perc_rep=sum(x$"n")/max(x$"end") )    ;  setDT(y); return(y) } )
is4 = do.call(rbind,is4)
#head(is4,20)

is5 = lapply(is3, function(x) { head(x[ rev(order(x$"perc.rep")), ])  } )
is5 = lapply(is5, function(x) { y=merge(x, bbb, by="chr"); return(y)} )
is5 = is5 [ which(sapply(is5, nrow)>0) ]

# 
# summary(trf$"len")
#     sum(trf$"len" > 100000)
#     sum(trf$"len" >  50000)
#     sum(trf$"len" >  10000)
#     sum(trf$"len" >   5000)
#     sum(trf$"len" >   1000)
#     sum(trf$"len" >    500)
#     sum(trf$"len" >    100)
    
ncol.trf  = ncol(trf)
ncol.trf2 = ncol(trf) -2

cols.fil2 = setdiff( names(trf) , c("tot.seq","consensus.seq"))
cols.filt = setdiff( names(trf) , c("tot.seq"))
#trf[,..cols.filt]


# Filter by >30 copies and consensu.size = 7
ts2 = lapply(ts, function(x) { x=data.frame(x); y = x[ x$"copy.n"> 30 , ] ;  setDT(y); return(y) } )
ts3 = lapply(ts, function(x) { x=data.frame(x); y = x[ which(x$"copy.n"> 30 & x$"consensus.size"== consensus.size) , ] ; setDT(y); return(y) } )

first = data.frame( do.call(rbind, lapply(ts3 , function(x) { n=1       ; x=x[n, ..cols.filt] ; return(x)} ) ))
last  = data.frame( do.call(rbind, lapply(ts3 , function(x) { n=nrow(x) ; x=x[n, ..cols.filt] ; return(x)} ) ))

rownames(first) = NULL
rownames(last) = NULL

if (any(is.na(first$"chr")) ==TRUE)  first = first[ - which(is.na(first$"chr")), ]
if (any(is.na( last$"chr")) ==TRUE)  last= last [ - which(is.na(last$"chr")), ]

# 
# print(first)
# print(last)

# watch also after line 1
first.top10 = lapply(ts2 , function(x) head(x[, 1:ncol.trf2], 10)  ) 
last.last10 = lapply(ts2 , function(x) tail(x[, 1:ncol.trf2], 10)  ) 

first.top10 = do.call(rbind, first.top10)
last.last10 = do.call(rbind, last.last10)

# loger telo?
first.telo.top10 = lapply(ts3 , function(x) head(x[, 1:ncol.trf2], 10)  ) 
last.telo.last10 = lapply(ts3 , function(x) tail(x[, 1:ncol.trf2], 10)  ) 


# telomere only
telo = lapply(ts, function(x) { x=data.frame(x); ; y = x[ x$"copy.n"> 30 & x$"period"== consensus.size , ] ; y=rbind(head(y,1), tail(y,1)); y=y[ ! duplicated(y), ] ;setDT(y); return(y) } )
telo = do.call(rbind, telo)
#telo[, ..cols.filt]


# rm telomeres

# other = lapply(ts, function(x) { y = x[ !( x$"copy.n"> 30 & x$"period"== consensus.size) , ] ; nn = nrow(y);  y=y[  11:(nn-10), ] ; return(y) } )
# ot_k1 = lapply(other , function(x) head(x[ rev(order(x$"len")),], 50)  ) 
# other = do.call(rbind, other)
# #other[, ..cols.filt]

# a = ts5 ; aa = sapply(a , nrow) ; aa [ aa>0 ] ;   a2 = a[ which( sapply(a , nrow) >0) ] ; a3 = data.frame(do.call(rbind,a2)) ; a3[, cols.filt ]
# a = ts6 ; aa = sapply(a , nrow) ; aa [ aa>0 ] ;   a2 = a[ which( sapply(a , nrow) >0) ] ; a3 = data.frame(do.call(rbind,a2)) ; a3[, cols.filt ]


# consensus of CHRs

 last$"consensus.seq.rc" = toupper( gsub('A', 't', gsub('C', 'g', gsub('G', 'c',  gsub('T', 'a',  sapply( strsplit(last$"consensus.seq","") , function(x) paste(rev(x), collapse="")) )))))


# if not Chr, sort by contig length and take the first 20
w_chr = grep("Chr|chr",first$"chr")
if (length(w_chr)==0 ) 
{
    first = first[ rev(order(first$"chr.len")), ]
    last  =  last[ rev(order( last$"chr.len")), ]
    w_chr = 1:20
    }



consesus.stat.first =  table(    as.character( first$"consensus.seq"[ w_chr] ) )
consesus.stat.last  =  table(    as.character(  last$"consensus.seq.rc"[w_chr ] ) )
consesus.stat       =  table( c( as.character( first$"consensus.seq"[w_chr ]) , as.character( last$"consensus.seq.rc"[w_chr ]) ) ) 

consesus.stat.top       = names( which.max(consesus.stat))[1]
# consesus.stat.first.top = names( which.max(consesus.stat.first))[1]
# consesus.stat.last.top  = names( which.max(consesus.stat.last))[1]


# is it AAACCCT?

consesus.stat.top.2x = paste(consesus.stat.top,consesus.stat.top,sep="")
consesus.stat.top.2x = gsub("AAACCCT", "aaaccct", consesus.stat.top.2x)
consesus.stat.top.2x = gsub("A|G|C|T", "", consesus.stat.top.2x )



if (consesus.stat.top.2x == "aaaccct") consesus.stat.top = "AAACCCT"

cat('\n --->  Most frequent Telomere sequence:\n')
cat('\n            ',consesus.stat.top,"\n\n")

# tel.seq provided
if ( is.na(tel.seq)==FALSE ) {  consesus.stat.top = toupper(tel.seq) ; cat('\n --->  user defined Telomere sequence: --> ', toupper(consesus.stat.top), ' \n\n') } ; 



consesus.stat.first.top = consesus.stat.top
consesus.stat.last.top =  toupper( gsub('A', 't', gsub('C', 'g', gsub('G', 'c',  gsub('T', 'a',  sapply( strsplit(consesus.stat.top,"") , function(x) paste(rev(x), collapse="")) )))))

cat("\n --->  5' Telomere sequence: --> ", toupper(consesus.stat.first.top), "")
cat("\n --->  3' Telomere sequence: --> ", toupper(consesus.stat.last.top), " [ reverse complement ]\n")

first$"consensus.2x" = paste(first$"consensus.seq",first$"consensus.seq",sep="")
last$"consensus.2x"  = paste( last$"consensus.seq", last$"consensus.seq",sep="")

first$"consensus.2x" = gsub(consesus.stat.first.top, tolower(consesus.stat.first.top), first$"consensus.2x" )
 last$"consensus.2x" = gsub(consesus.stat.last.top, tolower(consesus.stat.last.top), last$"consensus.2x" )
 
first$"consensus.2x" = gsub("A|G|C|T", "_", first$"consensus.2x" )
 last$"consensus.2x" = gsub("A|G|C|T", "_", last$"consensus.2x" )


consesus.stat.first.top.2x = paste0(tolower(consesus.stat.first.top),tolower(consesus.stat.first.top))
consesus.stat.last.top.2x  = paste0(tolower(consesus.stat.last.top) ,tolower(consesus.stat.last.top))

consesus.stat.first.top.1x = paste0(tolower(consesus.stat.first.top), gsub("a|g|c|t", "_", tolower(consesus.stat.first.top)))
consesus.stat.last.top.1x  = paste0(tolower(consesus.stat.last.top) , gsub("a|g|c|t", "_", tolower(consesus.stat.last.top)))

 
first$"consensus.2x" = gsub(consesus.stat.first.top.2x ,  consesus.stat.first.top.1x , first$"consensus.2x" )
 last$"consensus.2x" = gsub(consesus.stat.last.top.2x  ,  consesus.stat.last.top.1x ,  last$"consensus.2x" )
 
first$"CONSENSUS" = gsub("_","",first$"consensus.2x")
 last$"CONSENSUS" = gsub("_","", last$"consensus.2x")

if ( any( first$"CONSENSUS"=="" )) first[ which(first$"CONSENSUS"=="" ), ]$"CONSENSUS" = NA
if ( any(  last$"CONSENSUS"=="" ))  last[ which( last$"CONSENSUS"=="" ), ]$"CONSENSUS" = NA


cons.first = tolower(consesus.stat.first.top)
cons.lastt = tolower(consesus.stat.last.top)

# add telo start
FIRST = data.frame(first)
LASTT = data.frame(last)
# 
# 
# w_na = which(is.na(FIRST$"CONSENSUS")) ; if(length(w_na)>0) {  FIRST$"CONSENSUS"[w_na ]="--"}
# w_na = which(is.na(LASTT$"CONSENSUS")) ; if(length(w_na)>0) {  LASTT$"CONSENSUS"[w_na ]="--"}
# 

cat('\n --->  FIRST tandem repeat of each Chromosomes \n\n')
print(FIRST)
cat('\n --->  LAST  tandem repeat of each Chromosomes \n\n')
print(LASTT)
cat('\n')


cat('\n')
cat("############################################################################\n")
cat('------------------------ Repeat profile ------------------------------------\n')
cat("############################################################################\n")
cat('\n')

# rep = fread( paste(rep.data   ,sep=""),fill=TRUE,header=F)
if( ncol(rep)==15) names(rep)=c("SWscore","perc.div","perc.del","perc.ins","chr","start","end","left","strand","repeat","class","rep.start","rep.end","rep.left","ID")
if( ncol(rep)==16) names(rep)=c("SWscore","perc.div","perc.del","perc.ins","chr","start","end","left","strand","repeat","class","rep.start","rep.end","rep.left","ID","x")
#rep = rep[ 4:nrow(rep) , ]
#rep = rep[ ! duplicated(rep),] 

# rename
rep$"chr" = gsub("_RagTag","",rep$"chr"  )
if ( ref.name != FALSE)   {  rep = merge(rep,ref.edit, by="chr", all.x=T) ; w_ren = which(! is.na(rep$"chr.rename")) ;   rep$"chr"[w_ren] =rep$"chr.rename"[w_ren];  rep$"chr.rename"=NULL }


if ( is.na(rDNA.out)==FALSE)
{
    cat (" --> Adding rDNA.out to TE repeate masker .out    \n\n")
    rDNA = fread( paste(rDNA.out  ,sep=""),fill=TRUE,header=F)
    if( ncol(rDNA)==15) rDNA$"x"=NA
    names(rDNA)=c("SWscore","perc.div","perc.del","perc.ins","chr","start","end","left","strand","repeat","class","rep.start","rep.end","rep.left","ID","x")
    rDNA = rDNA[ 4:nrow(rDNA) , ]
    rDNA = rDNA[ grep("rDNA|ribos|Ribos", rDNA$"class"), ]
    rDNA$"chr" = gsub("_RagTag","",rDNA$"chr"  )
     if ( ref.name != FALSE)     rDNA = merge(rDNA,ref.edit, by="chr", all.x=T) ; w_ren = which(! is.na(rDNA$"chr.rename")) ;   rDNA$"chr"[w_ren] =rDNA$"chr.rename"[w_ren];  rDNA$"chr.rename"=NULL 
   # table(rep$"repeat" )
   # table(rDNA$"repeat" )
   # setdiff( unique(rDNA$"repeat"), unique(rep$"repeat") )
    rep = rep[ ! rep$"repeat" %in% rDNA$"repeat", ]
   # setdiff( unique(rDNA$"repeat"), unique(rep$"repeat") )
    rep = rbind(rep, rDNA)
    cat (" --> rDNA.out added to TE repeate masker .out    \n\n")

    }


cat('\n --->  Repeat Classes by total bp \n\n')
      as.data.frame( rev(sort(  table(rep$"class") )))
cat('\n --->  Low complexity and  Simple_repeat counts \n\n')
head( as.data.frame( rev(sort( table(rep[ rep$"class" %in% c("Low_complexity","Simple_repeat") ,]$"repeat" )))))
cat('\n')

rep$"end"      = as.numeric(as.character( rep$"end" ))
rep$"start"    = as.numeric(as.character( rep$"start" ))
rep$"rep.start"= as.numeric(as.character( rep$"rep.start" ))
rep$"rep.end"  = as.numeric(as.character( rep$"rep.end" ))
rep$"perc.ins" = as.numeric(as.character( rep$"perc.ins" ))
rep$"perc.del" = as.numeric(as.character( rep$"perc.del" ))
rep$"perc.div" = as.numeric(as.character( rep$"perc.div" ))

rep$"left"     = as.numeric(as.character( gsub("[(]|[)]","",rep$"left" )))
rep$"rep.left" = as.numeric(as.character( gsub("[(]|[)]","",rep$"rep.left" )))

rep$"len"     = rep$"end"- rep$"start" +1
rep$"repeat.class" =  paste(rep$"repeat", rep$"class",sep="#")
#rep$"rep.len" = rep$"rep.end" +rep$"rep.left"
rep = merge(rep, seq.len, by="repeat.class", all.x=TRUE)
rep$"perc.rep" = round(100*rep$"len" / rep$"rep.len",2)

rep$"REP.START" = rep$"rep.start" 
rep[ is.na( rep$"REP.START") , ]$"REP.START" = rep[ is.na( rep$"REP.START")  , ]$"rep.left"

 #  print(rep)
    
    
sim = rep[ rep$"class" == "Simple_repeat", ]
sim = sim [ order(sim$"chr",sim$"start")]
sim$"consensus.seq" = gsub("[)]|[(]|n","",sim$"repeat")
sim$"consensus.size" = nchar(sim$"consensus.seq")
sim$"copy.n" = round(    (sim$"rep.end"- sim$"rep.start"+1 ) / sim$"consensus.size" , 1)
sim$"repeat.class" = sim$"class" = sim$"rep.len" = sim$"rep.left" = sim$"perc.rep" = sim$"x" = sim$"REP.START"  = NULL

sim = sim [ sim$"copy.n"> 30 ,]
sim = sim [ sim$"len"   > 200 ,]
sim = sim [ sim$"consensus.size"   > 4 ,]



sim$"consensus.2x" = paste(sim$"consensus.seq",sim$"consensus.seq",sep="")
sim$"consensus.2x" = gsub(consesus.stat.first.top, tolower(consesus.stat.first.top), sim$"consensus.2x" )
sim$"consensus.2x" = gsub(consesus.stat.last.top , tolower(consesus.stat.last.top), sim$"consensus.2x" )
sim$"consensus.2x" = gsub("A|G|C|T", "_", sim$"consensus.2x" )
sim$"CONSENSUS" = toupper(gsub("_","",sim$"consensus.2x"))


w_ccc = which(sim$"consensus.seq" ==consesus.stat.first.top |sim$"consensus.seq" ==consesus.stat.last.top )
if ( length(w_ccc)>0) sim$"CONSENSUS"[w_ccc] = sim$"consensus.seq"[w_ccc]



si2 = sim[ sim$"CONSENSUS" %in% c(consesus.stat.first.top,consesus.stat.last.top ), ]

ss= split(si2, si2$"chr")

first_ss = lapply(ss , function(x) { y=x[ x$"start" < 200000 , ] ; return(y)} ) 
last__ss = lapply(ss , function(x) { y=x[ x$"left"  < 200000 , ] ; return(y)} )

first_ss = lapply(first_ss , function(x) { y=x[ which.min(x$"start" ) , ] ; return(y)} ) 
last__ss = lapply(last__ss , function(x) { y=x[ which.min(x$"left")  , ] ; return(y)} )


first_ss = do.call(rbind, first_ss)
last__ss = do.call(rbind, last__ss)

## Correct TRF telomere calls with repeat masker if possibile
cat('\n --->  Correct TRD with repeat masker if possibile\n\n')

FIRST$"tool" = "TRF"
LASTT$"tool" = "TRF"


first_ss$"tool" = "RepeatMasker"
last__ss$"tool" = "RepeatMasker"

first_ss$"CONSENSUS" =tolower( first_ss$"CONSENSUS" )
last__ss$"CONSENSUS" =tolower( last__ss$"CONSENSUS" )


F_ss = data.frame(chr=first_ss$"chr",left=first_ss$"left", "Repeat"=first_ss$"repeat",  "RepMask.consensus"=first_ss$"CONSENSUS","rm.start"=first_ss$"start", "rm.end"=first_ss$"end",  "RepMask.copy.n"=first_ss$"copy.n",stringsAsFactors=FALSE  )
L_ss = data.frame(chr=last__ss$"chr",left=last__ss$"left", "Repeat"=last__ss$"repeat",  "RepMask.consensus"=last__ss$"CONSENSUS","rm.start"=last__ss$"start", "rm.end"=last__ss$"end","RepMask.copy.n"=last__ss$"copy.n",stringsAsFactors=FALSE  )


FIRST = merge(FIRST,F_ss,by="chr", all.x=TRUE)
LASTT = merge(LASTT,L_ss,by="chr", all.x=TRUE)

FIRST$"Telomere_added" ="--"
LASTT$"Telomere_added" ="--"

 w_F_na = which(is.na(FIRST$"CONSENSUS") & FIRST$"RepMask.consensus"==tolower(consesus.stat.first.top))
 w_L_na = which(is.na(LASTT$"CONSENSUS") & LASTT$"RepMask.consensus"==tolower(consesus.stat.last.top))

if( length(w_F_na)> 0) FIRST$"Telomere_added"[w_F_na]="Repeat_masker"
if( length(w_L_na)> 0) LASTT$"Telomere_added"[w_L_na]="Repeat_masker"

if( length(w_F_na)> 0) FIRST$"CONSENSUS"[w_F_na]=FIRST$"RepMask.consensus"[w_F_na]
if( length(w_L_na)> 0) LASTT$"CONSENSUS"[w_L_na]=LASTT$"RepMask.consensus"[w_L_na]


FIRST$"start_shift" = FIRST$"start" - FIRST$"rm.start"  
LASTT$"start_shift" = LASTT$"start" - LASTT$"rm.start"  
FIRST$"end_shift" = FIRST$"end" - FIRST$"rm.end"  
LASTT$"end_shift" = LASTT$"end" - LASTT$"rm.end"  

check_cols = c("chr", "start","end","CONSENSUS","copy.n","Repeat","RepMask.consensus", "rm.start","rm.end","left","RepMask.copy.n","Telomere_added","start_shift","end_shift")
# FIRST[, check_cols]
# LASTT[, check_cols]

if( length(w_F_na)>0 ) { cat("\n ---> ",length(w_F_na)," Telomere(s) added at 5' from Repeat Masker \n\n") ; print(FIRST[, check_cols]) }
if( length(w_L_na)>0 ) { cat("\n ---> ",length(w_L_na)," Telomere(s) added at 3' from Repeat Masker \n\n") ; print(LASTT[, check_cols]) }

cat('\n ---> PLOT size estimation \n\n')



# PLOT coordinates

par(mfcol=c(1    ,1    ), mar=c(2,7,4,1) , oma=c(1,4,1,1) )
bar = barplot(mat, horiz=T, beside=T, las=1, col="white",xlim=c(-max(mat)*0.07,max(mat)*1.15), main=paste(MAIN.NAME,"Telomeres and Tandem Repeat profile"),space=mat.spac )
bbb = data.frame( chr=names(mat), "chr.len"=mat, pos = bar[,1])

is2 = merge(is2, bbb, by="chr")
agp = merge(AGP, bbb, by.y="chr",by.x="ref.chr")
gap = merge(GAP, bbb, by.y="chr",by.x="ref.chr")

MAX = max(mat)



FIRST$"chr" =  gsub("_RagTag","", FIRST$"chr" )
LASTT$"chr" =  gsub("_RagTag","", LASTT$"chr" )


FIRST = merge(FIRST, bbb, by=c("chr","chr.len"))
LASTT = merge(LASTT, bbb, by=c("chr","chr.len"))


## Differential Telomeres
FIRST$"TEL_seq"=FIRST$"CONSENSUS"
LASTT$"TEL_seq"=LASTT$"CONSENSUS"

w_first_diff = which( is.na(FIRST$"CONSENSUS"))
w_lastt_diff = which( is.na(LASTT$"CONSENSUS"))

for ( i in w_first_diff) { if (FIRST$"copy.n"[i]>300 & FIRST$"start"[i]<1000 )  {FIRST$"TEL_seq"[i]= FIRST$"consensus.seq"[i]   }}
for ( i in w_lastt_diff) { if (LASTT$"copy.n"[i]>300 & LASTT$"start"[i]<1000 )  {LASTT$"TEL_seq"[i]= LASTT$"consensus.seq"[i]   }}



if (nrow(FIRST)>0) FIRST$"chr.beg"=1
FIRST$"pos.tel"=NA
LASTT$"pos.tel"=NA

FIRST$"col.tel"=NA
LASTT$"col.tel"=NA


if ( any( grepl(cons.first,FIRST$"CONSENSUS" ))) FIRST[ grep(cons.first,FIRST$"CONSENSUS" ), ]$"pos.tel" = 1
if ( any( grepl(cons.lastt,LASTT$"CONSENSUS" ))) LASTT[ grep(cons.lastt,LASTT$"CONSENSUS" ), ]$"pos.tel" = LASTT[ grep(cons.lastt,LASTT$"CONSENSUS" ), ]$"chr.len"

w_blue = which(FIRST$"CONSENSUS" == FIRST$"TEL_seq")
w_gray = which(is.na(FIRST$"CONSENSUS") & FIRST$"consensus.seq" == FIRST$"TEL_seq")

if( length(w_blue)>0) FIRST$"col.tel"[w_blue]="blue"
if( length(w_gray)>0) FIRST$"col.tel"[w_gray]="cornflowerblue"
if( length(w_gray)>0) FIRST$"pos.tel"[w_gray] = 1


w_blue = which(LASTT$"CONSENSUS" == LASTT$"TEL_seq")
w_gray = which(is.na(LASTT$"CONSENSUS") & LASTT$"consensus.seq" == LASTT$"TEL_seq")

if( length(w_blue)>0) LASTT$"col.tel"[w_blue]="blue"
if( length(w_gray)>0) LASTT$"col.tel"[w_gray]="cornflowerblue"
if( length(w_gray)>0) LASTT$"pos.tel"[w_gray] = 1



if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.5, FIRST$"pos.tel",FIRST$"pos"+0.5, col=FIRST$"col", lwd=10 , lend=1)
if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.5, LASTT$"pos.tel",LASTT$"pos"+0.5, col=LASTT$"col", lwd=10 , lend=1)

if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", tolower(FIRST$"TEL_seq"), cex=1, adj=1, col=FIRST$"col")  
if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", tolower(LASTT$"TEL_seq"), cex=1, adj=0, col=LASTT$"col")  



write.table( first[ !is.na(first$"CONSENSUS"), ] , file = paste0(DIR_out,"/Chr_contigs_with_telomere_5prime.txt")  , col.names=TRUE, row.names=F, quote=F, sep="\t")    
write.table(  last[ !is.na( last$"CONSENSUS"), ] , file = paste0(DIR_out,"/Chr_contigs_with_telomere_3prime.txt")  , col.names=TRUE, row.names=F, quote=F, sep="\t")    

cat("\n\n")
END.time =Sys.time() 
TOT.time =difftime(END.time, START.time, units="mins")
cat( paste( "==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n\n" ) )




cat('\n --->  Internal Telomere sequences \n\n')



## Internal Telomere sequences 


telo.reptag = c( FIRST[ !is.na(FIRST$"CONSENSUS"), ]$"rep.tag", LASTT[ !is.na(LASTT$"CONSENSUS"), ]$"rep.tag"   )

miss.first = FIRST[ which(is.na(FIRST$"CONSENSUS")),]$"chr"
miss.lastt = LASTT[ which(is.na(LASTT$"CONSENSUS")),]$"chr"
miss.chr = sort(unique(c(miss.first,miss.lastt)))


trf.filt = data.frame(trf)
trf.filt = trf.filt[  which(trf.filt$"copy.n"> 30 & trf.filt$"consensus.size"== consensus.size)  ,]
trf.filt$"consensus.2x" = paste(trf.filt$"consensus.seq",trf.filt$"consensus.seq",sep="")
trf.filt = trf.filt[  grep( toupper(paste0(cons.first,"|",cons.lastt)), trf.filt$"consensus.2x")  ,]
trf.filt$"chr" = gsub("_RagTag","",trf.filt$"chr" )

#trf.filt[,cols.filt]

# --> internal telomeres in telomere-devoid Chr only
#trf.filt = trf.filt[ trf.filt$"chr" %in% miss.chr , ]

# --> remove real telomere annotations
trf.filt = trf.filt[ !( trf.filt$"rep.tag" %in% telo.reptag) , ]
#trf.filt[,cols.filt]


trf.filt = merge(trf.filt, bbb, by=c("chr","chr.len"))
trf.filt$"mid" = round(trf.filt$"start" + trf.filt$"len"/2)
trf.filt$"dist_to_end" = trf.filt$"chr.len" - trf.filt$"end"+1

# rm very close to telomeres
if(nrow(trf.filt)>0) trf.filt = trf.filt[ trf.filt$"start"       > 100000, ]
if(nrow(trf.filt)>0) trf.filt = trf.filt[ trf.filt$"dist_to_end" > 100000, ]
if(nrow(trf.filt)>0) trf.filt$"lab" = paste0(trf.filt$"consensus.seq"," ", trf.filt$"len" ,"bp")

write.table( trf.filt , file =  paste0(DIR_out,"/Chr_contigs_with_internal_telomeres.txt" ) , col.names=TRUE, row.names=F, quote=F, sep="\t")    



if(nrow(trf.filt)>0)  print(trf.filt[,cols.filt])

if(nrow(trf.filt)>0)  segments(trf.filt$"mid", trf.filt$"pos"-0.5, trf.filt$"mid",trf.filt$"pos"+0.1, col="red", lwd=5 , lend=1)
# text(    trf.filt$"mid" +MAX/100, trf.filt$"pos"+0.25, trf.filt$"consensus.seq", cex=0.8, adj=0, col="red")  
# text(    trf.filt$"mid" +MAX/100, trf.filt$"pos"-0.25, paste0(trf.filt$"len"," bp"), cex=0.8, adj=0, col="red")  
if(nrow(trf.filt)>0)  trf.filt.text = split(trf.filt ,trf.filt$"pos")

if(nrow(trf.filt)>0)  trf.filt.text = data.frame( POS = as.numeric(names(trf.filt.text)) , chr.len = sapply(trf.filt.text,function(x) x$"chr.len"[1]) ,LAB = sapply(trf.filt.text,function(x) paste(x$"lab", collapse=", "))   )
if(nrow(trf.filt)>0)  trf.filt.text$"LAB" = paste0(" internal ", trf.filt.text$"LAB" )
#text( trf.filt.text$"chr.len" +MAX/100, trf.filt.text$"POS"-0.25, trf.filt.text$"LAB", cex=0.8, adj=0, col="red")  


is2$"POS" = is2$"pos" -0.5 + (is2$"n"/100000)
is2$"mid" = round( (is2$"start" + is2$"end")/2 ,0)

segments(is2$"start", is2$"pos"-0.5, is2$"start",is2$"POS", col="gray", lwd=1 , lend=1)


if ( nrow(gap)> 0)  segments(gap$"ref.end", gap$"pos"-0.5, gap$"ref.end",gap$"pos"+0.5, col="black", lwd=1)
if ( nrow(gap)> 0)  points(gap$"ref.end",  gap$"pos", type="p", pch=23, cex=1 ,col="black", bg="yellow", lwd=2)


#SAVE 
# trf.filt
# is2
# FIRST
# LASTT

cat('\n')
END.time =Sys.time() 
TOT.time =difftime(END.time, START.time, units="mins")
cat( paste( "==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n\n" ) )

# 
# cat('\n')
# cat("############################################################################\n")
# cat('------------------------ Repeat profile ------------------------------------\n')
# cat("############################################################################\n")
# cat('\n')
# 
# rep = fread( paste(rep.data   ,sep=""),fill=TRUE,header=F)
# if( ncol(rep)==15) rep$"x"=NA
# names(rep)=c("SWscore","perc.div","perc.del","perc.ins","chr","start","end","left","strand","repeat","class","rep.start","rep.end","rep.left","ID","x")
# rep = rep[ 4:nrow(rep) , ]
# rep = rep[ ! duplicated(rep),] 
# 
# # rename
# rep$"chr" = gsub("_RagTag","",rep$"chr"  )
# if ( ref.name != FALSE)   {  rep = merge(rep,ref.edit, by="chr", all.x=T) ; w_ren = which(! is.na(rep$"chr.rename")) ;   rep$"chr"[w_ren] =rep$"chr.rename"[w_ren];  rep$"chr.rename"=NULL }
# 
# 
# if ( is.na(rDNA.out)==FALSE)
# {
#     cat (" --> Adding rDNA.out to TE repeate masker .out    \n\n")
#     rDNA = fread( paste(rDNA.out  ,sep=""),fill=TRUE,header=F)
#     if( ncol(rDNA)==15) rDNA$"x"=NA
#     names(rDNA)=c("SWscore","perc.div","perc.del","perc.ins","chr","start","end","left","strand","repeat","class","rep.start","rep.end","rep.left","ID","x")
#     rDNA = rDNA[ 4:nrow(rDNA) , ]
#     rDNA = rDNA[ grep("rDNA|ribos|Ribos", rDNA$"class"), ]
#     rDNA$"chr" = gsub("_RagTag","",rDNA$"chr"  )
#      if ( ref.name != FALSE)     rDNA = merge(rDNA,ref.edit, by="chr", all.x=T) ; w_ren = which(! is.na(rDNA$"chr.rename")) ;   rDNA$"chr"[w_ren] =rDNA$"chr.rename"[w_ren];  rDNA$"chr.rename"=NULL 
#    # table(rep$"repeat" )
#    # table(rDNA$"repeat" )
#    # setdiff( unique(rDNA$"repeat"), unique(rep$"repeat") )
#     rep = rep[ ! rep$"repeat" %in% rDNA$"repeat", ]
#    # setdiff( unique(rDNA$"repeat"), unique(rep$"repeat") )
#     rep = rbind(rep, rDNA)
#     cat (" --> rDNA.out added to TE repeate masker .out    \n\n")
# 
#     }
# 
# 
# cat('\n --->  Repeat Classes by total bp \n\n')
#       as.data.frame( rev(sort(  table(rep$"class") )))
# cat('\n --->  Low complexity and  Simple_repeat counts \n\n')
# head( as.data.frame( rev(sort( table(rep[ rep$"class" %in% c("Low_complexity","Simple_repeat") ,]$"repeat" )))))
# cat('\n')
# 
# rep$"end"      = as.numeric(as.character( rep$"end" ))
# rep$"start"    = as.numeric(as.character( rep$"start" ))
# rep$"rep.start"= as.numeric(as.character( rep$"rep.start" ))
# rep$"rep.end"  = as.numeric(as.character( rep$"rep.end" ))
# rep$"perc.ins" = as.numeric(as.character( rep$"perc.ins" ))
# rep$"perc.del" = as.numeric(as.character( rep$"perc.del" ))
# rep$"perc.div" = as.numeric(as.character( rep$"perc.div" ))
# 
# rep$"left"     = as.numeric(as.character( gsub("[(]|[)]","",rep$"left" )))
# rep$"rep.left" = as.numeric(as.character( gsub("[(]|[)]","",rep$"rep.left" )))
# 
# rep$"len"     = rep$"end"- rep$"start" +1
# rep$"repeat.class" =  paste(rep$"repeat", rep$"class",sep="#")
# #rep$"rep.len" = rep$"rep.end" +rep$"rep.left"
# rep = merge(rep, seq.len, by="repeat.class", all.x=TRUE)
# rep$"perc.rep" = round(100*rep$"len" / rep$"rep.len",2)
# 
# rep$"REP.START" = rep$"rep.start" 
# rep[ is.na( rep$"REP.START") , ]$"REP.START" = rep[ is.na( rep$"REP.START")  , ]$"rep.left"
# 
#     print(rep)


r80 = rep[   rep$"perc.rep" >= 80, ]
r80 = r80[ ! r80$"class" %in% c("Low_complexity","Simple_repeat"), ]
r80 = r80[   grep("rDNA|ribo", r80$"class" , invert=T), ]
r80 = r80[   r80$"len" > 500, ]
r80 = r80[ ! which(is.na(r80$"rep.start")), ]

  
# rename
if ( ref.name != FALSE)  {  CHR.bbb = merge(CHR.bbb,ref.edit, by="chr", all.x=T) ; w_ren = which(! is.na(CHR.bbb$"chr.rename")) ;   CHR.bbb$"chr"[w_ren] =CHR.bbb$"chr.rename"[w_ren];  CHR.bbb$"chr.rename"=NULL }

  
if( length(rDNA.presence)>0 )
{
     
    rib = rep[   grep("rDNA|ribo", rep$"class" ), ]
    rib = rib[   rib$"perc.rep" >= 20, ]
    rib = rib[   rib$"len" >= min(rDNA.len$"rep.len")*0.60, ]
#     print(rib)
#     print(CHR.bbb)
#     print(bbb)

    rib$"chr" = gsub("_RagTag","",rib$"chr"  )
    rib = merge(rib, CHR.bbb, by="chr",all.x=TRUE )          # all.x=TRUE  --> important for statistics of contigs < chr.filt
    rib = merge(rib, bbb, by=c("chr","chr.len"),all.x=TRUE ) # all.x=TRUE  --> important for statistics of contigs < chr.filt
    
    ris = split(rib,rib$"repeat")
    rxx = lapply( ris, function(x) rev(sort(  table(x$"chr") )) )
    rxx = lapply( rxx, function(x) data.frame( "Var1" = names(x),"Freq" = as.numeric(x)   )) 
    
 #   print(ris)

    
    for( i in names(rxx)) names(rxx[[i]])=c("chr",i)
    
    RXX = data.frame( chr=CHR$"chr", stringsAsFactors=FALSE)
    RXX$"chr" = gsub("_RagTag","",RXX$"chr"  )
    for( i in names(rxx)) RXX=merge(RXX, rxx[[i]], by="chr", all.x=T)
    
    RXX [ is.na(RXX) ] =0
    
    rownames(RXX) = RXX$"chr"
    RXX$"chr" = NULL
        
    cat('\n --->  Ribosomial RNA counts per Chr/contig \n\n')
    print(RXX)
    cat('\n --->  Chr/contigs with Ribosomial RNA hits \n\n')
    print(RXX[RXX$"rDNA_28S" >0 | RXX$"rDNA_5S" >0 ,])
    
    cat('\n --->  Ribosomial bp per Chr/Contig \n\n')
    ric  = split(rib,rib$"chr")
    ric = do.call(rbind, lapply( ric, function(x) data.frame( chr=x$"chr"[1], "rDNA.bp"=sum(x$"len"), chr.len=x$"chr.len"[1],stringsAsFactors=TRUE) ))
    ric$"perc.rDNA"=round(100*ric$"rDNA.bp"/ric$"chr.len",1)
   
    print(ric[,-1])
    END.time =Sys.time() 
    TOT.time =difftime(END.time, START.time, units="mins")
    cat( paste( "\n==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n\n" ) )
    }
    
    
    

     
cat('\n --->  TE class counts (>80% length of query) \n\n')
as.data.frame( rev(sort(  table(r80$"class") )))
cat('\n')

re2 = rep[rep$"len" > 5000, ]
re2 = re2[re2$"class" %in% c("Low_complexity","Simple_repeat"), ]


rrr = fread( rei.data,fill=TRUE,header=F)
names(rrr)=c("chr","start","end","chrom","s","e","n")
rrr$"chr" = gsub("_RagTag","",rrr$"chr"  )
if ( ref.name != FALSE)  {  rrr = merge(rrr,ref.edit, by="chr", all.x=T) ; w_ren = which(! is.na(rrr$"chr.rename")) ;   rrr$"chr"[w_ren] =rrr$"chr.rename"[w_ren];  rrr$"chr.rename"=NULL }
rrr$"chrom" = rrr$"chr"
rrr$"wind" =paste0(rrr$"chr",":",rrr$"start","-",rrr$"end")

RRR = list()
RRR[["all"]]    = fread( paste(rei.data   ,sep=""),fill=TRUE,header=F)
RRR[["TEs"]]    = fread( paste(rei.dataTE   ,sep=""),fill=TRUE,header=F)
RRR[["simple"]] = fread( paste(rei.dataSI   ,sep=""),fill=TRUE,header=F)
RRR[["rDNA"]]   = fread( paste(rei.dataRD   ,sep=""),fill=TRUE,header=F)


# if simple repeat or rDNA are asbsent, remove them
for ( ii in names(RRR)) if (length(RRR[[ii]])==0) {RRR[[ii]] = NULL  }


RRR = lapply(RRR, function(x) { names(x)= c("chr","start","end","chrom","s","e","n") ; return(x)})
RRR = lapply(RRR, function(x) { x$"chr" = gsub("_RagTag","",x$"chr"  ) ; return(x)})

#chr rename
if ( ref.name != FALSE) RRR = lapply(RRR, function(x) { x=merge(x,ref.edit, by="chr"    , all.x=T) ; w_rep = which( ! is.na(x$"chr.rename")) ; x$"chr"     [w_rep] =x$"chr.rename"[w_rep] ; x$"chr.rename"=NULL ; return(x) } )
if ( ref.name != FALSE) RRR = lapply(RRR, function(x) { x$"chrom" = x$"chr" ; return(x) } )



RRR = lapply(RRR, function(x) { x$"wind" = paste0(x$"chr",":",x$"start","-",x$"end"); return(x)})
RRs = lapply(RRR, function(x) { x= split(x,x$"wind"); return(x)})
RR2 = lapply(RRs, lapply, function(x) { y=data.frame(chr=x$"chr"[1],start=x$"start"[1],end=x$"end"[1],n=sum(x$"n")) ;  setDT(y); return(y) } )
RR3 = lapply(RR2, function(x) { y=do.call(rbind,x); y = y[ order(y$"chr", y$"start") , ] ;  return(y)   } )
RR3 = lapply(RR3, function(x) { x$"perc.rep" = x$"n" / (x$"end" -x$"start"+1) ; return(x)})

Rs2 = lapply(RR3, function(x) { y=merge(x, bbb, by="chr"); y$"POS" = y$"pos" -0.5 + (y$"n"/100000) ; return(y)})


rs = split(rrr,rrr$"wind")

rs2 = lapply(rs, function(x) { y=data.frame(chr=x$"chr"[1],start=x$"start"[1],end=x$"end"[1],n=sum(x$"n"))    ;  setDT(y); return(y) } )
rs2 = do.call(rbind,rs2)
rs2 = rs2[ order(rs2$"chr", rs2$"start") , ]

#last window is not 100Kb, so normalize
rs2$"perc.rep" = rs2$"n" / (rs2$"end" -rs2$"start"+1)


rs3 = split(rs2,rs2$"chr")
rs4 = lapply(rs3, function(x) { y=data.frame( chr=x$"chr"[1], len=max(x$"end"), rep=sum(x$"n"), perc_rep=sum(x$"n")/max(x$"end") )    ;  setDT(y); return(y) } )
rs4 = do.call(rbind,rs4)
#head(rs4)



rs2 = merge(rs2, bbb, by="chr")
rs2$"POS" = rs2$"pos" -0.5 + (rs2$"n"/100000)
rs2$"mid" = round( (is2$"start" + rs2$"end")/2 ,0)


r80$"chr" = gsub("_RagTag","",r80$"chr"  )
r80 = merge(r80, bbb, by="chr")

###  rm _LTR _INT
r80$"repeat" = gsub("_LTR|_INT","",r80$"repeat"  )


TTT = list()
TTT[["TEs"]]    = fread( paste(rei.dataTE   ,sep=""),fill=TRUE,header=F)

names(TTT[["TEs"]] )=c("chr","start","end","chrom","s","e","n")
TTT[["TEs"]]$"chr" = gsub("_RagTag","",TTT[["TEs"]]$"chr"  )
if ( ref.name != FALSE) TTT = lapply(TTT, function(x) { x=merge(x,ref.edit, by="chr"    , all.x=T) ; w_rep = which( ! is.na(x$"chr.rename")) ; x$"chr"     [w_rep] =x$"chr.rename"[w_rep] ; x$"chr.rename"=NULL ; return(x) } )
TTT[["TEs"]]$"chrom" = TTT[["TEs"]]$"chr"
#rrr$"wind" =paste0(rrr$"chr",":",rrr$"start","-",rrr$"end")


te_fam_col=data.frame( rev(sort(  table(rep$"class") )))
names(te_fam_col) = c("class","n")
te_fam_col$"col" = rainbow(nrow(te_fam_col))
#r80 = merge(r80, te_fam_col, by="class")

END.time =Sys.time() 
TOT.time =difftime(END.time, START.time, units="mins")
cat( paste( "==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n\n" ) )



# top TEs 80%
cat('\n --->  top TEs by 80-80 rule  \n\n')

top.te = setDT(as.data.frame( rev(sort(  table(r80$"repeat") ))))
x = setDT ( data.frame( "repeat" =r80$"repeat", "class"=r80$"class", stringsAsFactors=FALSE))
names(top.te) = c("repeat","freq")
names(x)      = c("repeat","class")
x = x[ ! duplicated(x), ]
top.te = merge(top.te,x, vy="repeat")
top.te = top.te[ rev(order(top.te$"freq")), ]

top10.te = head(top.te,10)
top10.te$"col" = rainbow(nrow(top10.te))

top10.te

r80_top10 = r80[ r80$"repeat" %in%  top.te$"repeat", ]
r80_top10 = merge(r80_top10,top10.te, by=c("repeat","class"))
r80_top10$"chr_te" = paste(r80_top10$"chr" ,r80_top10$"repeat",r80_top10$"class",sep="::")

bed_cols = c("chr","start","end")
new_cols = c(bed_cols ,  setdiff( names(r80_top10), bed_cols) )

r80_top10 = r80_top10[ , ..new_cols]

# top TEs all
cat('\n --->  top TEs total annotations \n\n')

rep$"TE" = gsub("_LTR|_INT","",rep$"repeat")
te.tot_bp = sapply( split(rep$"len",rep$"TE") , sum )
te.tot_bp = rev(sort(  te.tot_bp ))
top.TE    = setDT ( data.frame( "TE" =names(te.tot_bp), "freq"=te.tot_bp, stringsAsFactors=FALSE))


col.rep=c("repeat","class","TE")
r2 = rep[,..col.rep]
r2 = r2[ ! duplicated(r2), ]
top.TE = merge(top.TE, r2,by="TE")
top.TE = top.TE[ rev(order(top.TE$"freq")), ]
top.TE = top.TE[ grep( "TE" , top.TE$"TE"  ), ]

top10.TE = head(top.TE,10)
top10.TE$"col" = rainbow(nrow(top10.TE))

#top10.TE

cat('\n')

R80_TOP10 = rep[ rep$"repeat" %in%  top10.TE$"repeat", ]
R80_TOP10$"chr" = gsub("_RagTag","",R80_TOP10$"chr"  )

R80_TOP10 = merge(R80_TOP10,top10.TE,by=c("TE","repeat","class"))
R80_TOP10$"chr_te" = paste(R80_TOP10$"chr" ,R80_TOP10$"TE",R80_TOP10$"class",sep="::")

bed_cols = c("chr","start","end")
new_cols = c(bed_cols ,  setdiff( names(R80_TOP10), bed_cols) )

R80_TOP10 = R80_TOP10[ , ..new_cols]


top10_bed = "Top10_TE_repeatfile.bed"
top10_mer = "Top10_TE_repeatfile.bed"
top10_bed = "Top10_TE_repeatfile.bed"
 
r80s = split(r80_top10, r80_top10$"repeat")
r80s = split(R80_TOP10, R80_TOP10$"repeat")
#print(r80s[1])

# we will intersect r80 with original files, that requires NCBI chromosomes names if they were converted before
if ( ref.name != FALSE)
{
ref.re_edit = ref.edit
r80ss = r80s

ref.re_edit$"chr" ->  old_chr
ref.re_edit$"chr" = ref.re_edit$"chr.rename" 
ref.re_edit$"chr.rename" = old_chr


cat("check ref-edit \n")
# print(ref.edit)
# print(ref.re_edit)
}

# print(ref.name)
# print(r80s[1])

#print(r80ss[[1]])

if ( ref.name != FALSE) r80ss = lapply(r80ss, function(x) { x=merge(x,ref.re_edit, by="chr"    , all.x=T) ; w_rep = which( ! is.na(x$"chr.rename")) ; x$"chr"     [w_rep] =x$"chr.rename"[w_rep] ; x$"chr.rename"=NULL ; return(x) } )

#print(r80ss[[1]])
  


 if ( file.exists( paste0(DIR_out,"/Top10_TE_REPEATFILE.ALL.intersect_100Kb_wind.bed")) == FALSE | run=="all"  )
 {
       
    system( paste0( " sed 's|_RagTag||g' ",base.name,".100Kb_windows.bed  > ",base.name,".100Kb_windows.v2.bed  " ))

    for (i in names(r80s)) {
           # print(i)
            j  =paste0(DIR_out,"/Top10_TE_repeatfile.",i,".bed")   
            jjj=paste0(DIR_out,"/Top10_TE_repeatfile.",i,".intersect_100Kb_wind.bed")   
            write.table( r80s[[i]], file = j , col.names=FALSE, row.names=F, quote=F, sep="\t")    
            system( paste0( " bedtools intersect -wao -a ",base.name,".100Kb_windows.v2.bed -b ",j, " > ",jjj) )
            }
    
    system( paste0( " cat ",DIR_out,"/Top10_TE_repeatfile*intersect_100Kb_wind.bed > ",DIR_out,"/Top10_TE_REPEATFILE.ALL.intersect_100Kb_wind.bed" ),ignore.stdout = FALSE )
    system( paste0( " cat ",DIR_out,"/Top10_TE_repeatfile*intersect_100Kb_wind.bed > ",DIR_out,"/Top10_TE_REPEATFILE.ALL.intersect_100Kb_wind.bed" ),ignore.stdout = FALSE)
 }       

cat("check TOP TEs\n")
# chr is the A file (windows)
# CHR is the B file (TE files))

TOP = list()
for (i in names(r80s)) { print(i); jjj=paste0(DIR_out,"Top10_TE_repeatfile.",i,".intersect_100Kb_wind.bed") ; TOP[[i]]=fread(jjj,fill=TRUE,header=F) }

#cat("check TOP rename \n")

top_names = names(R80_TOP10)
top_names[ top_names=="chr"] = "CHR"


TOP = lapply(TOP, function(x) { names(x)=c("chr","START","END",top_names,"overlap") ; return(x) }) 
#print(TOP[[1]])
if ( ref.name != FALSE) TOP = lapply(TOP, function(x) { x=merge(x,ref.edit, by="chr"    , all.x=T) ; w_rep = which( ! is.na(x$"chr.rename")) ; x$"chr"     [w_rep] =x$"chr.rename"[w_rep] ; x$"chr.rename"=NULL ; return(x) } )
if ( ref.name != FALSE) TOP = lapply(TOP, function(x) { x$"CHR"  = x$"chr" ; return(x) } )
#print(TOP[[1]])



TOP = lapply(TOP, function(x) { names(x)=c("chr","START","END",top_names,"overlap") ; return(x) }) 
# print(TOP[[1]])
TOP = lapply(TOP, function(x) { x$"CHR"  = x$"chr"; return(x) }) 
TOP = lapply(TOP, function(x) { x$"wind" =paste0(x$"CHR",":",x$"START","-",x$"END"); return(x) }) 
TOP = lapply(TOP, function(x) { xs = split(x,x$"wind"); return(xs) }) 
TOP = lapply(TOP, lapply,function(x) { y=data.frame( "CHR" =x$"CHR"[1],  "START" =x$"START"[1], "END" =x$"END"[1], "overlap" =sum(as.numeric(x$"overlap")), stringsAsFactors=FALSE); return(y) }) 
TOP = lapply(TOP, function(x) { do.call(rbind,x) }) 
TOP = lapply(TOP, function(x) { setDT(x) }) 
TOP = lapply(TOP, function(x) { x$"wind" =paste0(x$"CHR",":",x$"START","-",x$"END"); return(x) }) 
for (i in names(TOP) )  names(TOP[[i]]) = c("CHR","START","END",i,"wind")  
TOP2 = lapply(TOP, function(x) { y = data.frame(x[,4]); rownames(y)=x$"wind"; return(y) }) 
# print(lapply(TOP,head))
# print(lapply(TOP2,head))
TOP3 = setDT(do.call(cbind,TOP2))
#print(TOP3)
END.time =Sys.time()
TOT.time =difftime(END.time, START.time, units="mins")
cat( paste( "==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n\n" ) )

TOP4 = cbind(TOP[[1]][,-4], TOP3)
#print(TOP4)

TOP5 = merge(TOP4, bbb,  by.x="CHR", by.y="chr")
TOP5$"POS" = TOP5$"pos"-0.5
#print(TOP5)

MAX.TE = max(TOP3)
#tops = split(TOP, TOP$"repeat")

# TE superfamilies
cat('\n --->  TE superfamilies  \n\n')

if( file.exists( paste0(DIR_out,"/TE_merged_annotations_by_class.edit.txt'")) == FALSE |  run=="all" )
{
cat(  '      ->  Merging annotations by class \n\n')
    re3 = rep
#re3$"chr" = gsub("_RagTag","",re3$"chr")
    re3$"chr.class" = paste(re3$"chr",re3$"class",sep="___")
    re3_col = c("chr.class","start","end","strand","repeat")
    re3 = re3[,..re3_col]
    re3 = re3[ order(re3$"chr.class", re3$"start") ,]
#cat("check 6 \n")
    write.table( re3, file = paste0(DIR_out,"/TE_annotations_by_class.txt") , col.names=FALSE, row.names=F, quote=F, sep="\t")    
    system( paste0(" bedtools merge -i ",DIR_out,"/TE_annotations_by_class.txt | sed 's/___/\t/g' > ",DIR_out,"TE_merged_annotations_by_class.merged.txt " ))
    system( paste0(" awk '{ t=$2 ; $2=$3; $3=$4; $4=t; print }' ",DIR_out,"/TE_merged_annotations_by_class.merged.txt | sed 's/ /\t/g'> ",DIR_out,"/TE_merged_annotations_by_class.merged.edit.txt  " ))
 #  system( paste0(" bedtools intersect -wao -a ",base.name,".100Kb_windows.bed -b ",DIR_out,"/TE_merged_annotations_by_class.merged.edit.txt | cut -f1-8 > ",DIR_out,"/TE_merged_annotations_by_class.intersect.txt" ))
    system( paste0(" bedtools intersect -wao -a <( cat ",base.name,".100Kb_windows.bed | sed 's/_RagTag//g' ) -b ",DIR_out,"/TE_merged_annotations_by_class.merged.edit.txt | cut -f1-8 > ",DIR_out,"/TE_merged_annotations_by_class.intersect.txt " ))
     
#cat("check 7 \n")
    re4 = fread(paste0(DIR_out,"/TE_merged_annotations_by_class.intersect.txt"),fill=TRUE)
    names(re4)=c("chr","start","end","chrom","s","e","class","n")
# print(head(re4))
    re4$"wind" =paste0(re4$"chr",":",re4$"start","-",re4$"end")
    re4$"n" =as.numeric(re4$"n" )
    re4[ which(is.na(re4$"n")) ,]
# cat("check 8 \n")
    
    re4 = split(re4,re4$"wind")
    re4 = lapply(re4, function(x) split(x,x$"class") )
    
    re5 = lapply(re4, lapply, function(x) { y=data.frame(chr=x$"chr"[1],start=x$"start"[1],end=x$"end"[1],class=x$"class"[1],wind=x$"wind"[1],n=sum(x$"n"))    ;  setDT(y); return(y) } )
    re6 = lapply(re5,function(x) do.call(rbind,x) )
    re6 = do.call(rbind,re6)
    re6$"chr" = gsub("_RagTag","",re6$"chr"  )
#  cat("check 9 \n")
#  print(head(re4))
#  print(head(re5))
#  print(head(re6))
    
    re6 = re6[ order(re6$"chr", re6$"start") , ]
    re6$"perc.rep" = re6$"n" / (re6$"end" -re6$"start"+1)
    
    re6 = merge(re6, bbb, by="chr")
    re6$"POS" = re6$"pos" -0.5 + (re6$"n"/100000)
    re6$"mid" = round( (re6$"start" + re6$"end")/2 ,0)
    
#cat("check 10 \n")
    write.table( re6, file = paste0(DIR_out,"/TE_merged_annotations_by_class.edit.txt") , col.names=TRUE, row.names=F, quote=F, sep="\t")   
} else { cat(  '      -> Merged annotations by class present! \n\n')   }


re6 = fread( paste0(DIR_out,"/TE_merged_annotations_by_class.edit.txt"),fill=TRUE, header=TRUE)
#print(re6)
re7 = split(re6,re6$"class")
re7 = re7 [ names(re7) != "" ]
cat('\n ')
        
    
END.time =Sys.time() 
TOT.time =difftime(END.time, START.time, units="mins")
cat( paste( "\n==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n" ) )

if ( is.na(GeneGFF)==FALSE ) 
{
    
    cat('\n --->  Gene Density \n\n')
    
    genes = fread( rei.dataGE,fill=TRUE,header=F)
    names(genes)=c("chr","start","end","chrom","s","e","n")
    genes$"wind" =paste0(genes$"chr",":",genes$"start","-",genes$"end")
    gs = split(genes,genes$"wind")
    gs2 = lapply(gs, function(x) { y=data.frame(chr=x$"chr"[1],start=x$"start"[1],end=x$"end"[1],n=sum(x$"n"))    ;  setDT(y); return(y) } )
    gs2 = do.call(rbind,gs2)
    gs2 = gs2[ order(gs2$"chr", gs2$"start") , ]
    gs2$"perc.rep" = gs2$"n" / (gs2$"end" -gs2$"start"+1)
    gs2 = merge(gs2, bbb, by="chr")
    gs2$"POS" = gs2$"pos" -0.5 + (gs2$"n"/100000)
    gs2$"mid" = round( (gs2$"start" + gs2$"end")/2 ,0)
    } else {  cat('\n --->  no Gene Density data \n\n') }


if ( is.na(add.dots)==FALSE ) 
{
    
    cat('\n --->  adding DOTS from FILE ',add.dots,' \n\n')
    
    dots = fread( add.dots,fill=TRUE,header=F)
  #  print(dots)
    dots.names=c("chr","start","name","color","pch", "size","legend.name")
    if( length(dots.names) == ncol(dots) ) { names(dots)=dots.names ;
    } else {
        dots.names_1 = dots.names[ 1:ncol(dots)]
        dots.names_2 = setdiff(dots.names, dots.names_1) 
        names(dots)=dots.names_1 ; 
        for ( ddd in dots.names_2 ) dots[[ ddd ]]= NA;  
        if ( all(is.na(dots$"color")) == TRUE ) dots$"color" = "red"
        if ( all(is.na(dots$"pch"  )) == TRUE ) dots$"pch"   =  19
        if ( all(is.na(dots$"size" )) == TRUE ) dots$"size"  = 2
        if ( all(is.na(dots$"legend.name" )) == TRUE ) dots$"size"  = "custom.region"
        }
 #   print(dots)
    
    dots.leg = unique(dots[, -(1:3) ])
    dots$"chr" = gsub("_RagTag","", dots$"chr")
    dots = merge(dots, bbb, by="chr")
 #   print(bbb)
    print(dots)
    cat('\n --->  adding to Legend \n\n')
    print(dots.leg)
  
       }


if ( is.na(add.segm)==FALSE ) 
{
    
    cat('\n --->  adding SEGMENTS from FILE',add.segm,' \n\n')
    
    segm = fread( add.dots,fill=TRUE,header=F)
    segm.names=c("chr","start","end","name","color","pch", "size","legend.name")
    if( length(segm.names) == ncol(segm) ) { names(dots)=segm.names ;
    } else {
        segm.names_1 = segm.names[ 1:ncol(dots)]
        segm.names_2 = setdiff(segm.names, segm.names_1) 
        names(dots)=segm.names_1 ; 
        for ( ddd in segm.names_2 ) segm[[ ddd ]]= NA;  
        if ( all(is.na(segm$"color")) == TRUE ) segm$"color" = "red"
        if ( all(is.na(segm$"pch"  )) == TRUE ) segm$"pch"   =  19
        if ( all(is.na(segm$"size" )) == TRUE ) segm$"size"  = 2
        }
    }


END.time =Sys.time() 
TOT.time =difftime(END.time, START.time, units="mins")
cat( paste( "==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n\n" ) )

#print(rDNA.presence)

if( length(rDNA.presence)>0 )
{
    
  # print(RXX)
    rDNA_contigs = RXX[ RXX$"rDNA_28S" >1 & grepl(".tg|tig", rownames(RXX)),]
    print(rDNA_contigs)
    rDNA_contigs_length = sum(chr[ chr$"chr" %in% rownames(rDNA_contigs), ]$"chr.len")
    rDNA_chr = ris[[ "rDNA_28S" ]][  grep(".tg|tig",  ris[[ "rDNA_28S" ]]$"chr", invert=T) ,]
  # cat("test\n")
  # print(class(rDNA_chr))
  # print((rDNA_chr))
   
if( length(rDNA_chr)>0 )
   {
    rDNA_chr = split(rDNA_chr,rDNA_chr$"chr")
    rDNA_chr = lapply(rDNA_chr, function(x) data.frame(chr=x$"chr"[1],start=min(x$"start"), end=max(x$"end"), n.loci=nrow(x), stringsAsFactors=FALSE ))
    rDNA_chr[["contigs"]] = data.frame(chr="contigs",start=1, end=rDNA_contigs_length, n.loci=sum(rDNA_contigs$"rDNA_28S"), stringsAsFactors=FALSE )
    rDNA_chr = do.call(rbind,rDNA_chr)
    rDNA_chr = rbind( data.frame(chr="rDNA",start=0, end=0, n.loci=0, stringsAsFactors=FALSE ) ,rDNA_chr)
    rDNA_chr$"length" = rDNA_chr$"end"-rDNA_chr$"start"+1
    rDNA_chr$"LEN" =  round( rDNA_chr$"length"/1e6,1)
    rDNA_chr$"CHR" =  rDNA_chr$"chr"
    rDNA_chr$"SIZE"=  rDNA_chr$"LEN"
    rDNA_chr$"LOCI"=  rDNA_chr$"n.loci"
    y=max(nchar( rDNA_chr$"chr"))    ; for (i in 1:nrow(rDNA_chr)) { z = y-nchar(rDNA_chr$"chr"[i])   ; rDNA_chr$"CHR"[i]  = paste0( paste0(rep(" ",z),collapse=""), rDNA_chr$"chr"[i],collapse="")}
    y=max(nchar( rDNA_chr$"LEN"))    ; for (i in 1:nrow(rDNA_chr)) { z = y-nchar(rDNA_chr$"LEN"[i])   ; rDNA_chr$"SIZE"[i] = paste0( paste0(rep(" ",z),collapse=""), rDNA_chr$"LEN"[i],collapse="")}
    y=max(nchar( rDNA_chr$"n.loci")) ; for (i in 1:nrow(rDNA_chr)) { z = y-nchar(rDNA_chr$"n.loci"[i]); rDNA_chr$"LOCI"[i] = paste0( paste0(rep(" ",z),collapse=""), rDNA_chr$"n.loci"[i],collapse="")}
    
    rDNA_chr$"legend" = paste0( rDNA_chr$"CHR", " ", rDNA_chr$"SIZE"," Mb ", rDNA_chr$"LOCI"," loci")
    rDNA_chr$"pch"=19 
    rDNA_chr$"col"="orange"
    rDNA_chr$"pt.bg"="black"
    rDNA_chr$"legend"[1]="rDNA"
    rDNA_chr$"pch"[1]=NA
    rDNA_chr$"pch"[1]=NA
    
    rDNA_CHR = rDNA_chr[ grep(".tg|tig", rDNA_chr$"chr", invert=TRUE), ]

    
    RXX$"chr" = rownames(RXX)
    RXX = RXX[,c(ncol(RXX),1: (ncol(RXX)-1))]
    
write.table( rDNA_chr[-1,], file = paste0(DIR_out,"/rDNA_chr.stats") , col.names=TRUE, row.names=F, quote=F, sep="\t")    
write.table( RXX          , file = paste0(DIR_out,"/rDNA_loci_per_chromosome.txt"  ), col.names=TRUE, row.names=F, quote=F, sep="\t")    
write.table( ric          , file = paste0(DIR_out,"/rDNA_percentage_per_chromosome.txt" ), col.names=TRUE, row.names=F, quote=F, sep="\t")    

    }
    }


# relead files:
# rs2
# gs2
# is2
# gap
# ris
# FIRST LAST


MGP = 1 - round(3/name.space,1)
if (name.space >10 ) MGP=0

PT.SIZE = 2.2 * leg.size
LG.CEX  = 2   * leg.size

# simplify names with




if( name.del != FALSE  ) 
{ 
    cat(' \n --->   Chromosome renaming: \n')      
    cat(' \n        --> deleting strings: ',gsub(",","  ", name.del),'  \n\n')      
    name.del=gsub(",","|", name.del)
    names.mat = gsub(name.del,"", names(mat) )  ;
    xx = data.frame( old = names(mat) , new =names.mat )
    print(xx)
    names(mat)= names.mat
    }




################################################################
####  PLOT 0)  Contig names only
################################################################
cat(' \n --->   PLOT 0) Contig position \n\n')

pdf( paste0(DIR_out,"/TELfinder_00.Contig_position.",chr.filt.char,".pdf" )  ,width=plot.w, height=plot.h) ; 
par(mfcol=c(1    ,1    ), mar=c(2,6,4,1) ,  oma=c(1,3+round(name.space/5),1,1) ,mgp = c(3,MGP, 0) )

bg.col = rgb(0  , 0  , 0, alpha=0.15)
bar = barplot(mat, horiz=T, beside=T, las=1, col=bg.col, xlim=c(-max(mat)*0.05,max(mat)*1.15), main=paste(MAIN.NAME,"Contigs") , cex.names=1.5, cex.main=2.5,space=mat.spac )

if ( nrow(gap)> 0)  segments(gap$"ref.end", gap$"pos"-0.5, gap$"ref.end",gap$"pos"+0.5, col="black", lwd=2)
agp$"ref.mid"  =  (agp$"ref.start" + agp$"ref.end")/2

# if scaffold
if ( all(grepl("Chr|Ragtag|LG|[.]",agp$"ref.chr"))==TRUE  ) agp$"ctg_name" =  agp$"query"
# if assembly
if ( all(grepl("seq",agp$"query"))==TRUE &  nrow(gap) == 0) agp$"ctg_name" =  agp$"ref.chr"

 text( agp$"ref.mid"  , agp$"pos", agp$"ctg_name", cex=1.5, col="black")  
 
ctg.df = split( agp, agp$"pos") 
#print(ctg.df)
ctg.df = lapply(ctg.df , function(x) data.frame( "pos" = unique(x$"pos") , "ref.end" = max(x$"ref.end"), "ctg_name" = paste(" --> sorted contigs:", paste(x$"ctg_name", collapse= " "))   ))
#print(ctg.df)
ctg.df = do.call(rbind,ctg.df)

 text( ctg.df$"ref.end"+max(mat)*0.02  , ctg.df$"pos", ctg.df$"ctg_name", cex=1.2, col="black", adj=0)  
# 
# print(agp)
# print(ctg.df)
cat('\n')
dev.off()


################################################################
####  PLOT 1.b) Telomeres, tandem and all repeats
################################################################
cat(' \n --->   PLOT 1.a) Telomeres, tandem and all repeats   \n\n')

pdf( paste0(DIR_out,"TELfinder_01a.Telomeres_rDNA_repeat_profile.",chr.filt.char,".pdf" )  ,width=plot.w, height=plot.h) ; 
par(mfcol=c(1    ,1    ), mar=c(2,6,4,1) ,  oma=c(1,3+round(name.space/5),1,1) ,mgp = c(3,MGP,0) )

bg.col = rgb(0  , 0  , 0, alpha=0.15)
bar = barplot(mat, horiz=T, beside=T, las=1, col=bg.col, xlim=c(-max(mat)*0.05,max(mat)*1.15), main=paste(MAIN.NAME,"Telomeres and Tandem Repeat profile"), cex.names=1.5, cex.main=2.5,space=mat.spac )


GREEN = rgb(0  , 0.5, 0,   alpha=0.50)
GRAY  = rgb(0  , 0.0, 0,   alpha=0.50)
RED   = rgb(0  , 1  , 1,   alpha=0.50)
RED   = rgb(1  , 1  , 0,   alpha=0.60)

                             segments(rs2$"end"  , rs2$"pos"-0.5, rs2$"end"  ,rs2$"POS", col= GREEN  , lwd=4 , lend=1, border=NULL)
if ( is.na(GeneGFF)==FALSE ) segments(gs2$"start", gs2$"pos"-0.5, gs2$"start",gs2$"POS", col= RED   , lwd=4 , lend=1, border=NULL)
                             segments(is2$"start", is2$"pos"-0.5, is2$"start",is2$"POS", col= GRAY   , lwd=4 , lend=1, border=NULL)



# rDNA
if( length(rDNA.presence)>0 ) a = "rDNA_28S" ; points(ris[[ a ]]$"start", ris[[ a ]]$"pos"-0.4, col="orange", cex=0.8, pch=19, lend=1)
if( length(rDNA.presence)>0 ) a = "rDNA_5S"  ; points(ris[[ a ]]$"start", ris[[ a ]]$"pos"-0.4, col="magenta", cex=0.8, pch=19,lend=1)

# redraw a box
par(new = T,lwd = 2)
barplot(mat, horiz=T, beside=T, las=1, col=rgb(0, 0, 0, alpha=0) ,xlim=c(-max(mat)*0.05,max(mat)*1.15),cex.names=1.5 ,space=mat.spac )


if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.5, FIRST$"pos.tel",FIRST$"pos"+0.5, col=FIRST$"col", lwd=10 , lend=1)
if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.5, LASTT$"pos.tel",LASTT$"pos"+0.5, col=LASTT$"col", lwd=10 , lend=1)

if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", tolower(FIRST$"TEL_seq"), cex=1, adj=1, col=FIRST$"col")  
if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", tolower(LASTT$"TEL_seq"), cex=1, adj=0, col=LASTT$"col")  

# Gaps

if ( nrow(gap)> 0)  segments(gap$"ref.end", gap$"pos"-0.5, gap$"ref.end",gap$"pos"+0.5, col="black", lwd=2)
if ( nrow(gap)> 0)  points(gap$"ref.end",  gap$"pos", type="p", pch=23, cex= 2  ,col="black", bg="yellow", lwd=3)

if ( is.na(add.dots)==FALSE )  points(dots$"start",  dots$"pos", type="p", pch= dots$"pch", cex= dots$"size" ,col=dots$"color")


# legeng
LEGEND = data.frame( "name"=c("Telomere presence + seq","non-canonical Telomeres","Repeat %","Tandem Repeat %","GAP","rDNA 28S","rDNA 5S") , pch=c(22,22,22,22,23,22,22) , col=c("black","black","black","black","black","black","black"),  pt.bg=c("blue","cornflowerblue",GREEN,GRAY,"yellow","orange","magenta"), stringsAsFactors=FALSE)

if ( is.na(GeneGFF)==FALSE  ) LEGEND = rbind(         data.frame( "name"="Gene %"               , pch=22            , col="black",  pt.bg=RED             , stringsAsFactors=FALSE) , LEGEND )
if ( is.na(add.dots)==FALSE ) LEGEND = rbind( LEGEND, data.frame( "name"=dots.leg$"legend.name" , pch=dots.leg$"pch", col="black",  pt.bg=dots.leg$"color", stringsAsFactors=FALSE)  )

legend("bottomright", LEGEND$"name", pch=LEGEND$"pch",  col=LEGEND$col , pt.bg=LEGEND$"pt.bg", cex= LG.CEX, bg="white",pt.cex=PT.SIZE)
#legend("bottom",    rDNA_chr$"legend", pch=rDNA_chr$"pch",  col=rDNA_chr$col , pt.bg=rDNA_chr$"pt.bg", cex= LG.CEX, bg="white",pt.cex=PT.SIZE)


dev.off()



################################################################
####  PLOT 1.a) Telomeres, tandem and all repeats
################################################################
cat(' \n --->   PLOT 1.b bis) Telomeres, tandem and all repeats + telomere length and internal telomere sequences \n\n')

pdf( paste0(DIR_out,"TELfinder_01b.Telomeres_rDNA_repeat_profile.",chr.filt.char,".with_data.pdf" )  ,width=plot.w, height=plot.h) ; 
par(mfcol=c(1    ,1    ), mar=c(2,6,4,1) ,  oma=c(1,3+round(name.space/5),1,1) ,mgp = c(3,MGP, 0) )

bg.col = rgb(0  , 0  , 0, alpha=0.15)
bar = barplot(mat, horiz=T, beside=T, las=1, col=bg.col, xlim=c(-max(mat)*0.05,max(mat)*1.15), main=paste(MAIN.NAME,"Telomeres and Tandem Repeat profile"), cex.names=1.5, cex.main=2.5,space=mat.spac )


GREEN = rgb(0  , 0.5, 0,   alpha=0.50)
GRAY  = rgb(0  , 0.0, 0,   alpha=0.50)
RED   = rgb(0  , 1  , 1,   alpha=0.50)
RED   = rgb(1  , 1  , 0,   alpha=0.60)

                             segments(rs2$"end"  , rs2$"pos"-0.5, rs2$"end"  ,rs2$"POS", col= GREEN  , lwd=4 , lend=1, border=NULL)
if ( is.na(GeneGFF)==FALSE ) segments(gs2$"start", gs2$"pos"-0.5, gs2$"start",gs2$"POS", col= RED   , lwd=4 , lend=1, border=NULL)
                             segments(is2$"start", is2$"pos"-0.5, is2$"start",is2$"POS", col= GRAY   , lwd=4 , lend=1, border=NULL)

# 
# xxx=Rs2[["all"]]    ; segments(xxx$"end"  , xxx$"pos"-0.5, xxx$"end"  ,xxx$"POS", col= "pink"  , lwd=4 , lend=1, border=NULL)
# xxx=Rs2[["simple"]] ; segments(xxx$"end"  , xxx$"pos"-0.5, xxx$"end"  ,xxx$"POS", col= "blue"  , lwd=4 , lend=1, border=NULL)
# xxx=Rs2[["rDNA"]]   ; segments(xxx$"end"  , xxx$"pos"-0.5, xxx$"end"  ,xxx$"POS", col= "red"  , lwd=4 , lend=1, border=NULL)
# xxx=Rs2[["TEs"]]    ; segments(xxx$"end"  , xxx$"pos"-0.5, xxx$"end"  ,xxx$"POS", col= "black"  , lwd=4 , lend=1, border=NULL)
# 

# rDNA
if( length(rDNA.presence)>0 ) a = "rDNA_28S" ; points(ris[[ a ]]$"start", ris[[ a ]]$"pos"-0.4, col="orange", cex=0.8, pch=19, lend=1)
if( length(rDNA.presence)>0 ) a = "rDNA_5S"  ; points(ris[[ a ]]$"start", ris[[ a ]]$"pos"-0.4, col="magenta", cex=0.8, pch=19,lend=1)

# redraw a box
par(new = T,lwd = 2)
barplot(mat, horiz=T, beside=T, las=1, col=rgb(0, 0, 0, alpha=0) ,xlim=c(-max(mat)*0.05,max(mat)*1.15),cex.names=1.5 ,space=mat.spac )



# Telomers
# if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.52, FIRST$"pos.tel",FIRST$"pos"+0.52, col="blue", lwd=9 , lend=1)
# if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.52, LASTT$"pos.tel",LASTT$"pos"+0.52, col="blue", lwd=9 , lend=1)
# 
# if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", FIRST$"CONSENSUS", cex=1.2, adj=1, col="blue")  
# if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", LASTT$"CONSENSUS", cex=1.2, adj=0, col="blue")  

if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.5, FIRST$"pos.tel",FIRST$"pos"+0.5, col=FIRST$"col", lwd=10 , lend=1)
if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.5, LASTT$"pos.tel",LASTT$"pos"+0.5, col=LASTT$"col", lwd=10 , lend=1)

if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", tolower(FIRST$"TEL_seq"), cex=1, adj=1, col=FIRST$"col")  
if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", tolower(LASTT$"TEL_seq"), cex=1, adj=0, col=LASTT$"col")  

if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos"-0.5, FIRST$"len", cex=0.8, adj=1, col="blue")  
if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos"-0.5, LASTT$"len", cex=0.8, adj=0, col="blue")  

# gaps 

if ( nrow(gap)> 0)  segments(gap$"ref.end", gap$"pos"-0.5, gap$"ref.end",gap$"pos"+0.5, col="black", lwd=2)
if ( nrow(gap)> 0)  points(gap$"ref.end",  gap$"pos", type="p", pch=23, cex= 2  ,col="black", bg="yellow", lwd=3)

# internal tel
if(nrow(trf.filt)>0)  segments(trf.filt$"mid", trf.filt$"pos"-0.5, trf.filt$"mid",trf.filt$"pos"+0.0, col="blue", lwd=5 , lend=1)
if(nrow(trf.filt)>0)  text( trf.filt.text$"chr.len" +MAX/20, trf.filt.text$"POS"-0.25, trf.filt.text$"LAB", cex=1, adj=0, col="blue")  


# legeng
LEGEND = data.frame( "name"=c("Telomere presence + seq","Telomere-like internal seq","non-canonical Telomeres","Repeat %","Tandem Repeat %","GAP","rDNA 28S","rDNA 5S") , pch=c(22,22,22,22,22,23,22,22) , col=c("black","black","black","black","black","black","black","black"),  pt.bg=c("blue","blue","cornflowerblue",GREEN,GRAY,"yellow","orange","magenta"), stringsAsFactors=FALSE)

if ( is.na(GeneGFF)==FALSE ) LEGEND = rbind( data.frame( "name"="Gene %" , pch=22, col="black",  pt.bg=RED, stringsAsFactors=FALSE) , LEGEND )

# print(LEGEND)
# print(rDNA_chr)

if( length(rDNA_chr)>0 ) 
{
legend("bottom",    rDNA_chr$"legend", pch=rDNA_chr$"pch",  col=rDNA_chr$col , pt.bg=rDNA_chr$"pt.bg", cex= LG.CEX, bg="white",pt.cex=PT.SIZE)
}
legend("bottomright", LEGEND$"name", pch=LEGEND$"pch",  col=LEGEND$col , pt.bg=LEGEND$"pt.bg", cex= LG.CEX, bg="white",pt.cex=PT.SIZE)


dev.off()



LEGEND2 = LEGEND[ grep("rDNA|Telo",LEGEND$"name", invert=T) ,]
LEGEND0 = data.frame( "name"="" , pch=22 , col=NA,  pt.bg=NA, stringsAsFactors=FALSE)

################################################################
####  PLOT 2a) Top10 TEs overlaps
################################################################
cat(' \n --->   PLOT 2a) Top10 TEs overlaps  \n\n')


pdf(  paste0(DIR_out,"TELfinder_02a.Telomeres_top10_TE_profile.",chr.filt.char,".pdf" )  ,width=plot.w, height=plot.h) ; 
par(mfcol=c(1    ,1    ), mar=c(2,6,4,1) ,  oma=c(1,3+round(name.space/5),1,1) ,mgp = c(3,MGP, 0) )


bg.col = rgb(0  , 0  , 0, alpha=0.10)
bar = barplot(mat, horiz=T, beside=T, las=1, col=bg.col, xlim=c(-max(mat)*0.05,max(mat)*1.05), main=paste(MAIN.NAME,"Telomeres, Top 10 TEs profile"), cex.names=1.5, cex.main=2.5,space=mat.spac )


GREEN = rgb(0  , 0.5, 0,   alpha=0.20)
GRAY  = rgb(0  , 0.0, 0,   alpha=0.20)
RED   = rgb(0  , 1  , 1,   alpha=0.20)
RED   = rgb(1  , 1  , 0,   alpha=0.20)

                             segments(rs2$"end"  , rs2$"pos"-0.5, rs2$"end"  ,rs2$"POS", col= GREEN  , lwd=4 , lend=1, border=NULL)
if ( is.na(GeneGFF)==FALSE ) segments(gs2$"start", gs2$"pos"-0.5, gs2$"start",gs2$"POS", col= RED   , lwd=4 , lend=1, border=NULL)
                             segments(is2$"start", is2$"pos"-0.5, is2$"start",is2$"POS", col= GRAY   , lwd=4 , lend=1, border=NULL)


# 
# if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.5, FIRST$"pos.tel",FIRST$"pos"+0.5, col="blue", lwd=10 , lend=1)
# if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.5, LASTT$"pos.tel",LASTT$"pos"+0.5, col="blue", lwd=10 , lend=1)
# 
# if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", FIRST$"CONSENSUS", cex=1.5, adj=1, col="blue")  
# if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", LASTT$"CONSENSUS", cex=1.5, adj=0, col="blue")  

if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.5, FIRST$"pos.tel",FIRST$"pos"+0.5, col=FIRST$"col", lwd=10 , lend=1)
if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.5, LASTT$"pos.tel",LASTT$"pos"+0.5, col=LASTT$"col", lwd=10 , lend=1)

if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", tolower(FIRST$"TEL_seq"), cex=1, adj=1, col=FIRST$"col")  
if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", tolower(LASTT$"TEL_seq"), cex=1, adj=0, col=LASTT$"col")  





par(new = T,lwd = 2)
barplot(mat, horiz=T, beside=T, las=1, col=rgb(0, 0, 0, alpha=0) ,xlim=c(-max(mat)*0.05,max(mat)*1.05),cex.names=1.5 ,space=mat.spac )


if ( nrow(gap)> 0)  segments(gap$"ref.end", gap$"pos"-0.5, gap$"ref.end",gap$"pos"+0.5, col="black", lwd=2)
if ( nrow(gap)> 0)  points(gap$"ref.end",  gap$"pos", type="p", pch=23, cex= 2  ,col="black", bg="yellow", lwd=3)

k=-8000

TOP.COL = top10.te
TOP.COL = top10.TE
TOP.COL$"repeat2" = paste( TOP.COL$"repeat", TOP.COL$"class")
TOP.COL$"name" = paste( TOP.COL$"repeat", TOP.COL$"class")
TOP.COL$"pch" = 22
TOP.COL$"pt.bg" = TOP.COL$"col"
TOP.COL$"col" = "black"
for(i in names(TOP3)) { print(i); k=k+8000;  segments(TOP5$"START"+k  , TOP5$"POS", TOP5$"START"+k  ,TOP5$"POS"+ TOP5[[i]]/MAX.TE  , col= TOP.COL[ TOP.COL$"repeat" == i , ]$"pt.bg" , border=NULL)   }

TOP.COL2 = data.frame(TOP.COL)

LEGEND3 = TOP.COL2 [, names(LEGEND2)]
LEGEND33 = data.frame( "name"="*TE content zoomed, not % of 100Kb" , pch=22 , col=NA,  pt.bg=NA, stringsAsFactors=FALSE)
LEGEND4 = rbind(LEGEND2,LEGEND0,LEGEND3,LEGEND33)

legend("bottomright", LEGEND4$"name", pch=LEGEND4$"pch", col=LEGEND4$"col",  pt.bg=LEGEND4$"pt.bg", cex= LG.CEX  , bg="white",pt.cex=PT.SIZE)


dev.off()

################################################################
####  PLOT 2b) Top10 TEs one by one
################################################################
cat(' \n --->   PLOT 2b) Top10 TEs one by one  \n\n')


pdf( paste0(DIR_out,"TELfinder_02b.Telomeres_top10_TE_profile.single.",chr.filt.char,".pdf" )    ,width=plot.w, height=plot.h) ; 
par(mfcol=c(1    ,1    ), mar=c(2,6,4,1) ,  oma=c(1,3+round(name.space/5),1,1) ,mgp = c(3,MGP, 0) )

for(i in names(TOP3)) { print(i); 

bg.col = rgb(0  , 0  , 0, alpha=0.10)
bar = barplot(mat, horiz=T, beside=T, las=1, col=bg.col, xlim=c(-max(mat)*0.05,max(mat)*1.05), main=paste(MAIN.NAME,"Telomeres, Top 10 TEs profile"), cex.names=1.5, cex.main=2.5,space=mat.spac )

GREEN = rgb(0  , 0.5, 0,   alpha=0.15)
GRAY  = rgb(0  , 0.0, 0,   alpha=0.15)
RED   = rgb(0  , 1  , 1,   alpha=0.15)
RED   = rgb(1  , 1  , 0,   alpha=0.15)

                             segments(rs2$"end"  , rs2$"pos"-0.5, rs2$"end"  ,rs2$"POS", col= GREEN  , lwd=4 , lend=1, border=NULL)
if ( is.na(GeneGFF)==FALSE ) segments(gs2$"start", gs2$"pos"-0.5, gs2$"start",gs2$"POS", col= RED   , lwd=4 , lend=1, border=NULL)
                             segments(is2$"start", is2$"pos"-0.5, is2$"start",is2$"POS", col= GRAY   , lwd=4 , lend=1, border=NULL)


# if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.5, FIRST$"pos.tel",FIRST$"pos"+0.5, col="blue", lwd=10 , lend=1)
# if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.5, LASTT$"pos.tel",LASTT$"pos"+0.5, col="blue", lwd=10 , lend=1)
# 
# if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", FIRST$"CONSENSUS", cex=1.5, adj=1, col="blue")  
# if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", LASTT$"CONSENSUS", cex=1.5, adj=0, col="blue")  

if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.5, FIRST$"pos.tel",FIRST$"pos"+0.5, col=FIRST$"col", lwd=10 , lend=1)
if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.5, LASTT$"pos.tel",LASTT$"pos"+0.5, col=LASTT$"col", lwd=10 , lend=1)

if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", tolower(FIRST$"TEL_seq"), cex=1, adj=1, col=FIRST$"col")  
if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", tolower(LASTT$"TEL_seq"), cex=1, adj=0, col=LASTT$"col")  



if ( nrow(gap)> 0)  segments(gap$"ref.end", gap$"pos"-0.5, gap$"ref.end",gap$"pos"+0.5, col="black", lwd=2)
if ( nrow(gap)> 0)  points(gap$"ref.end",  gap$"pos", type="p", pch=23, cex= 2  ,col="black", bg="yellow", lwd=3)

k=0
segments(TOP5$"START"+k  , TOP5$"POS", TOP5$"START"+k  ,TOP5$"POS"+ TOP5[[i]]/MAX.TE  , col= TOP.COL[ TOP.COL$"repeat" == i , ]$"pt.bg" , border=NULL)   


par(new = T,lwd = 2)
barplot(mat, horiz=T, beside=T, las=1, col=rgb(0, 0, 0, alpha=0) ,xlim=c(-max(mat)*0.05,max(mat)*1.05),cex.names=1.5 ,space=mat.spac )


legend("bottomright", LEGEND4$"name", pch=LEGEND4$"pch", col=LEGEND4$"col",  pt.bg=LEGEND4$"pt.bg", cex= LG.CEX  , bg="white",pt.cex=PT.SIZE)
}


dev.off()

################################################################
####  PLOT 2c) Top10 TEs one by one  -- version 2, TE in black
################################################################

cat(' \n --->   PLOT 2c) Top10 TEs one by one  -- version 2, TE in black  \n\n')

pdf( paste0(DIR_out,"TELfinder_02c.Telomeres_top10_TE_profile.single_black.",chr.filt.char,".pdf" )    ,width=plot.w, height=plot.h) ; 
par(mfcol=c(1    ,1    ), mar=c(2,6,4,1) ,  oma=c(1,3+round(name.space/5),1,1) ,mgp = c(3,MGP, 0) )

for(i in names(TOP3)) { print(i); 

bg.col = rgb(0  , 0  , 0, alpha=0.10)
bar = barplot(mat, horiz=T, beside=T, las=1, col=bg.col, xlim=c(-max(mat)*0.05,max(mat)*1.05), main=paste(MAIN.NAME,"Telomeres, Top 10 TEs profile"), cex.names=1.5, cex.main=2.5,space=mat.spac )

GREEN = rgb(0  , 0.5, 0,   alpha=0.15)
GRAY  = rgb(0  , 0.0, 0,   alpha=0.15)
RED   = rgb(0  , 1  , 1,   alpha=0.15)
RED   = rgb(1  , 1  , 0,   alpha=0.15)

                             segments(rs2$"end"  , rs2$"pos"-0.5, rs2$"end"  ,rs2$"POS", col= GREEN  , lwd=4 , lend=1, border=NULL)
if ( is.na(GeneGFF)==FALSE ) segments(gs2$"start", gs2$"pos"-0.5, gs2$"start",gs2$"POS", col= RED   , lwd=4 , lend=1, border=NULL)
                             segments(is2$"start", is2$"pos"-0.5, is2$"start",is2$"POS", col= GRAY   , lwd=4 , lend=1, border=NULL)


# if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.5, FIRST$"pos.tel",FIRST$"pos"+0.5, col="blue", lwd=10 , lend=1)
# if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.5, LASTT$"pos.tel",LASTT$"pos"+0.5, col="blue", lwd=10 , lend=1)
# 
# if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", FIRST$"CONSENSUS", cex=1.5, adj=1, col="blue")  
# if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", LASTT$"CONSENSUS", cex=1.5, adj=0, col="blue")  
# 
if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.5, FIRST$"pos.tel",FIRST$"pos"+0.5, col=FIRST$"col", lwd=10 , lend=1)
if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.5, LASTT$"pos.tel",LASTT$"pos"+0.5, col=LASTT$"col", lwd=10 , lend=1)

if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", tolower(FIRST$"TEL_seq"), cex=1, adj=1, col=FIRST$"col")  
if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", tolower(LASTT$"TEL_seq"), cex=1, adj=0, col=LASTT$"col")  



if ( nrow(gap)> 0)  segments(gap$"ref.end", gap$"pos"-0.5, gap$"ref.end",gap$"pos"+0.5, col="black", lwd=2)
if ( nrow(gap)> 0)  points(gap$"ref.end",  gap$"pos", type="p", pch=23, cex= 2  ,col="black", bg="yellow", lwd=3)

k=0
segments(TOP5$"START"+k  , TOP5$"POS", TOP5$"START"+k  ,TOP5$"POS"+ TOP5[[i]]/MAX.TE  , col= "black", border=NULL)   


par(new = T,lwd = 2)
barplot(mat, horiz=T, beside=T, las=1, col=rgb(0, 0, 0, alpha=0) ,xlim=c(-max(mat)*0.05,max(mat)*1.05),cex.names=1.5 ,space=mat.spac )



LEGEND3bis = LEGEND3[ grep(i,LEGEND3$"name"),]
LEGEND3bis$"pt.bg" = "black"
LEGEND4bis = rbind(LEGEND2,LEGEND0,LEGEND3bis)

legend("bottomright", LEGEND4bis$"name", pch=LEGEND4bis$"pch", col=LEGEND4bis$"col",  pt.bg=LEGEND4bis$"pt.bg", cex= LG.CEX  , bg="white",pt.cex=PT.SIZE)
}


dev.off()




################################################################
####  PLOT 3) TE classes separate
################################################################
cat(' \n --->   PLOT 3) TE classes separate  \n\n')


# TE classes
for (i in names(re7)) { xx=re7[[i]]; print(i) ;  segments(xx$"start", xx$"pos"-0.5, xx$"start",xx$"POS", col= "black"   , lwd=4 , lend=1, border=NULL)} 


pdf( paste0(DIR_out,"TELfinder_03.Telomeres_TEclasses.single.",chr.filt.char,".pdf" )    ,width=plot.w, height=plot.h) ; 
par(mfcol=c(1    ,1    ), mar=c(2,6,4,1) ,  oma=c(1,3+round(name.space/5),1,1) ,mgp = c(3,MGP, 0) )

for(i in names(re7)) { print(i); 

bg.col = rgb(0  , 0  , 0, alpha=0.10)
bar = barplot(mat, horiz=T, beside=T, las=1, col=bg.col, xlim=c(-max(mat)*0.05,max(mat)*1.05), main=paste(MAIN.NAME,"Telomeres, Top 10 TEs profile"), cex.names=1.5, cex.main=2.5,space=mat.spac )

GREEN = rgb(0  , 0.5, 0,   alpha=0.30)
GRAY  = rgb(0  , 0.0, 0,   alpha=0.30)
RED   = rgb(0  , 1  , 1,   alpha=0.30)
RED   = rgb(1  , 1  , 0,   alpha=0.30)

                             segments(rs2$"end"  , rs2$"pos"-0.5, rs2$"end"  ,rs2$"POS", col= GREEN  , lwd=4 , lend=1, border=NULL)
if ( is.na(GeneGFF)==FALSE ) segments(gs2$"start", gs2$"pos"-0.5, gs2$"start",gs2$"POS", col= RED   , lwd=4 , lend=1, border=NULL)
                             segments(is2$"start", is2$"pos"-0.5, is2$"start",is2$"POS", col= GRAY   , lwd=4 , lend=1, border=NULL)

# 
# if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.5, FIRST$"pos.tel",FIRST$"pos"+0.5, col="blue", lwd=10 , lend=1)
# if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.5, LASTT$"pos.tel",LASTT$"pos"+0.5, col="blue", lwd=10 , lend=1)
# 
# if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", FIRST$"CONSENSUS", cex=1.5, adj=1, col="blue")  
# if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", LASTT$"CONSENSUS", cex=1.5, adj=0, col="blue")  
# 
if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.5, FIRST$"pos.tel",FIRST$"pos"+0.5, col=FIRST$"col", lwd=10 , lend=1)
if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.5, LASTT$"pos.tel",LASTT$"pos"+0.5, col=LASTT$"col", lwd=10 , lend=1)

if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", tolower(FIRST$"TEL_seq"), cex=1, adj=1, col=FIRST$"col")  
if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", tolower(LASTT$"TEL_seq"), cex=1, adj=0, col=LASTT$"col")  



if ( nrow(gap)> 0)  segments(gap$"ref.end", gap$"pos"-0.5, gap$"ref.end",gap$"pos"+0.5, col="black", lwd=2)
if ( nrow(gap)> 0)  points(gap$"ref.end",  gap$"pos", type="p", pch=23, cex= 2  ,col="black", bg="yellow", lwd=3)

 xx=re7[[i]]
 segments(xx$"start", xx$"pos"-0.5, xx$"start",xx$"POS", col= "black"   , lwd=4 , lend=1, border=NULL)

par(new = T,lwd = 2)
barplot(mat, horiz=T, beside=T, las=1, col=rgb(0, 0, 0, alpha=0) ,xlim=c(-max(mat)*0.05,max(mat)*1.05),cex.names=1.5 ,space=mat.spac )



LEGEND_i = data.frame( "name"=i , pch=22 , col="black",  pt.bg="black", stringsAsFactors=FALSE)
LEGEND_5 = rbind(LEGEND2,LEGEND0,LEGEND_i)


legend("bottomright", LEGEND_5$"name", pch=LEGEND_5$"pch", col=LEGEND_5$"col",  pt.bg=LEGEND_5$"pt.bg", cex= LG.CEX  , bg="white",pt.cex=PT.SIZE)
}


dev.off()



################################################################
####  PLOT 4) LTR 
################################################################

cat(' \n --->   PLOT 4) LTR superfamilies separate  \n\n')

# TE classes

re8 = re7[ grep("LTR|Simple|Low",names(re7)) ]
re8 = re7[ grep("LTR",names(re7)) ]


pdf( paste0(DIR_out,"TELfinder_04.Telomeres_LTRclasses.single.",chr.filt.char,".pdf" )    ,width=plot.w, height=plot.h) ; 
par(mfcol=c(1    ,1    ), mar=c(2,6,4,1) ,  oma=c(1,3+round(name.space/5),1,1) ,mgp = c(3,MGP, 0) )

bg.col = rgb(0  , 0  , 0, alpha=0.10)
bar = barplot(mat, horiz=T, beside=T, las=1, col=bg.col, xlim=c(-max(mat)*0.05,max(mat)*1.05), main=paste(MAIN.NAME,"Telomeres, Top 10 TEs profile"), cex.names=1.5, cex.main=2.5,space=mat.spac )

GREEN = rgb(0  , 0.5, 0,   alpha=0.30)
GRAY  = rgb(0  , 0.0, 0,   alpha=0.30)
RED   = rgb(0  , 1  , 1,   alpha=0.30)
RED   = rgb(1  , 1  , 0,   alpha=0.30)

                             segments(rs2$"end"  , rs2$"pos"-0.5, rs2$"end"  ,rs2$"POS", col= GREEN  , lwd=4 , lend=1, border=NULL)
if ( is.na(GeneGFF)==FALSE ) segments(gs2$"start", gs2$"pos"-0.5, gs2$"start",gs2$"POS", col= RED   , lwd=4 , lend=1, border=NULL)
                             segments(is2$"start", is2$"pos"-0.5, is2$"start",is2$"POS", col= GRAY   , lwd=4 , lend=1, border=NULL)

# 
# if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.5, FIRST$"pos.tel",FIRST$"pos"+0.5, col="blue", lwd=10 , lend=1)
# if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.5, LASTT$"pos.tel",LASTT$"pos"+0.5, col="blue", lwd=10 , lend=1)
# 
# if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", FIRST$"CONSENSUS", cex=1.5, adj=1, col="blue")  
# if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", LASTT$"CONSENSUS", cex=1.5, adj=0, col="blue")  

if (nrow(FIRST)>0) segments(FIRST$"pos.tel", FIRST$"pos"-0.5, FIRST$"pos.tel",FIRST$"pos"+0.5, col=FIRST$"col", lwd=10 , lend=1)
if (nrow(LASTT)>0) segments(LASTT$"pos.tel", LASTT$"pos"-0.5, LASTT$"pos.tel",LASTT$"pos"+0.5, col=LASTT$"col", lwd=10 , lend=1)

if (nrow(FIRST)>0) text( -FIRST$"pos.tel" -MAX/100 , FIRST$"pos", tolower(FIRST$"TEL_seq"), cex=1, adj=1, col=FIRST$"col")  
if (nrow(LASTT)>0) text(  LASTT$"pos.tel" +MAX/100 , LASTT$"pos", tolower(LASTT$"TEL_seq"), cex=1, adj=0, col=LASTT$"col")  


if ( nrow(gap)> 0)  segments(gap$"ref.end", gap$"pos"-0.5, gap$"ref.end",gap$"pos"+0.5, col="black", lwd=2)
if ( nrow(gap)> 0)  points(gap$"ref.end",  gap$"pos", type="p", pch=23, cex= 2  ,col="black", bg="yellow", lwd=3)

LEGEND_ii = data.frame( "name"=names(re8) , pch=22 , col="black",  pt.bg=rainbow(length(re8)), stringsAsFactors=FALSE)

if(nrow(LEGEND_ii)==3) LEGEND_ii$"pt.bg" = c("red","blue","black")

for(i in names(re8)) { print(i);  xx=re7[[i]] ;  segments(xx$"start", xx$"pos"-0.5, xx$"start",xx$"POS", col= LEGEND_ii[ LEGEND_ii$"name" == i , ]$"pt.bg"  , lwd=4 , lend=1, border=NULL) }

par(new = T,lwd = 2)
barplot(mat, horiz=T, beside=T, las=1, col=rgb(0, 0, 0, alpha=0) ,xlim=c(-max(mat)*0.05,max(mat)*1.05),cex.names=1.5 ,space=mat.spac )



LEGEND_6 = rbind(LEGEND2,LEGEND0,LEGEND_ii)


legend("bottomright", LEGEND_6$"name", pch=LEGEND_6$"pch", col=LEGEND_6$"col",  pt.bg=LEGEND_6$"pt.bg", cex= LG.CEX  , bg="white",pt.cex=PT.SIZE)



dev.off()

system ( "rm Rplots.pdf")

END.time =Sys.time() 
TOT.time =difftime(END.time, START.time, units="mins")


cat( paste( "==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n\n" ) )
cat( paste( "==>> Output folder: ", DIR_out, " \n" ) )
cat( paste( "==>>     full path:  ", getwd(),"/",DIR_out, "/TELfinder.*",chr.filt.char,"*pdf \n", sep=""))
cat("\n")
cat("\n")





