
#!/bin/R

# R

args <- commandArgs(TRUE)



   if("--help" %in% args) 
    {
        cat('\n')
        cat(
        " Plot contig coverage from HIFIASM gfa file \n",
        " \n",
        "USAGE: Rscript Coverage_plot.R --gfa=my_hifiasm_assembly.fasta.bp.p_ctg.noseq.gfa  --agp=my_scaffold.agp \n",
        " \n",
        "Arguments:\n",
        "  --gfa=              | *noseq.gfa from HIFIASM","\n",
        "  --agp=              | Scaffold AGP file                     [ optional, default absent ] if absent contgis will be plot one by one","\n",
        "  --ref=              | Reference FASTA file                  [ optional, default absent ] used for importing real Chr names","\n",
        "  --chr.filt=         | default = 1Mb                         [ Integer                  ] minimal Chr size plotted [ default = 1Mb / 1000Kb / 1000000 ; numeric-only and Kb,Mb character accepted  ]  ","\n",
        "  --out=              | no default                            [ string                   ] output file suffix             ","\n",
        "  --wind=             | default 100Kb/100000)                 [ string                   ] window size [ default = 1000Kb / 1000000 ; numeric-only, Mb,Kb,bp character accepted  ] ", "\n",
        "  --plot.h=           | Height in inches of PDF output files  [ default 18               ]  ","\n",
        "  --plot.w=           | Width in inches of PDF output files   [ default 35               ]  ","\n",
        "Modules required: \n",
        " - R \n",
        " - R packages: data.table  \n",
        " - samtools \n",
        " - bedtools \n",
         "\n"
        )

        cat('\n')
        q(save="no")
    }



cat('\n')
cat('\n')
cat("############################################################################\n")
cat("############################################################################\n")
cat("###                            Coverage_plot.r                            ###\n")
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
    
	gfa.file        = as.character(ARG[ ARG$X1=="--gfa"      ,]$"X2")
	agp.file        = as.character(ARG[ ARG$X1=="--agp"      ,]$"X2")
	ref.file        = as.character(ARG[ ARG$X1=="--ref"      ,]$"X2")
	chr.filt        = as.character(ARG[ ARG$X1=="--chr.filt" ,]$"X2")
	wind            = as.character(ARG[ ARG$X1=="--wind"     ,]$"X2")
	plot.h          = as.character(ARG[ ARG$X1=="--plot.h"   ,]$"X2")
	plot.w          = as.character(ARG[ ARG$X1=="--plot.w"   ,]$"X2")
	out             = as.character(ARG[ ARG$X1=="--out"      ,]$"X2")
		
	plot.h = as.numeric(plot.h)
	plot.w = as.numeric(plot.w)


#  DEFAULT VALUES

    if ( length(agp.file)==0  ) agp.file = NA  
    if ( length(ref.file)==0  ) ref.file = NA  
    if ( length(chr.filt)==0  ) chr.filt = 1 * 1e6
    if ( length(plot.h)==0    ) plot.h      = 18
    if ( length(plot.w)==0    ) plot.w      = 35
    if ( length(wind)==0      ) wind      = 100000
    if ( length(out)==0       ) out      = NA
    
 
# TEST PARAMETERS for interactive session
if( 2 ==1 )
{
    gfa.file="../HIFIASM/Ziziphus_spina-christi_hifiasm.bp.p_ctg.noseq.gfa"
    agp.file="ragtag.scaffold.agp"
    ref.file="../Z_jujuba_Genome_ref.fasta"
    chr.filt="15Mb"
    wind = 10000000
    plot.h = 18
    plot.w = 35    
    out    = NA
  }


cat('\n')
cat("############################################################################\n")
cat('------------------------ Summary of Parameters -----------------------------\n')
cat("############################################################################\n")
cat('\n')
    cat("   gfa.file   = ", gfa.file,"\n")
    cat("   agp.file   = ", agp.file,"\n")
    cat("   ref.file   = ", ref.file,"\n")
    cat("   chr.filt   = ", chr.filt,"\n")
    cat("   wind       = ", wind  ,"\n")
    cat("   plot.h     = ", plot.h,"\n")
    cat("   plot.w     = ", plot.w,"\n")
    cat("   out        = ", out,"\n")

cat('\n')


OUT = paste0(".",out,".")
if ( is.na(out) ==TRUE ) OUT ="."

library(data.table)
options(width=300)
options(scipen=999)

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

wind.char = wind
if ( grepl("Kb|Mb|Gb|bp",wind.char) ==TRUE ) 
{   wind = as.numeric(gsub("bp","",(gsub("Kb","000",gsub("Mb","000000",gsub("Gb","000000000",wind))))))
  } else {
    
    wind.char = gsub("0000000","0Mb",wind)
    wind.char = gsub("000000" , "Mb",wind.char)
    wind.char = gsub("00000" ,"00Kb",wind.char)
    wind.char = gsub("0000"   ,"0Kb",wind.char)
    wind.char = gsub("000"     ,"Kb",wind.char)
    wind = as.numeric(wind)
    }


# 
# http://gfa-spec.github.io/GFA-spec/GFA1.html
# 
# # Terminology
# Segment: a continuous sequence or subsequence.
# Link: an overlap between two segments. Each link is from the end of one segment to the beginning of another segment. The link stores the orientation of each segment and the amount of basepairs overlapping.
# Jump: (since v1.2) a connection between two oriented segments. Similar to link, but does not imply a direct adjacency between the segments, instead providing an estimated distance between the segments. Main use case is to specify segment relations across assembly gaps.
# Containment: an overlap between two segments where one is contained in the other.
# Path: an ordered list of oriented segments, where each consecutive pair of oriented segments is supported by a link or a jump record.
# Walk: (since v1.1) an ordered list of oriented segments, intended for pangenome use cases. Each consecutive pair of oriented segments must correspond to a 0-overlap link record.



# Line structure
# Each line in GFA has tab-delimited fields and the first field defines the type of line. 
# The type of the line defines the following required fields. The required fields are followed by optional fields.
# 
# Type	Description
# #	Comment
# H	Header
# S	Segment
# L	Link
# J	Jump (since v1.2)
# C	Containment
# P	Path
# W	Walk (since v1.1)

# Optional fields
# All optional fields follow the TAG:TYPE:VALUE format where TAG is a two-character string that matches /[A-Za-z][A-Za-z0-9]/. Each TAG can only appear once in one line. A TAG containing lowercase letters are reserved for end users. A TYPE is a single case-sensitive letter which defines the format of VALUE.
# 
# Type	Regexp	Description
# A	[!-~]	Printable character
# i	[-+]?[0-9]+	Signed integer
# f	[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?	Single-precision floating number
# Z	[ !-~]+	Printable string, including space
# J	[ !-~]+	JSON, excluding new-line and tab characters
# H	[0-9A-F]+	Byte array in hex format
# B	[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+	Array of integers or floats


#################################################
#####   HIFIASM
#################################################

# Output file formats
# Hifiasm broadly follows the specification for GFA 1.0. There are several fields that are specifically used by hifiasm. For S segment line:
# 
# rd:i:: read coverage. 
# It is calculated by the reads coming from the same contig/unitig.
# 
# Hifiasm outputs A lines including the information of reads which are used to construct contig/unitig. 
# Each A line is plain-text, tab-separated, and the columns appear in the following order:
# 
# Type  Description
# 
# 1 string Should be always A
# 
# 2 string Contig/unitig name
# 
# 3 int  Contig/unitig start coordinate of subregion constructed by read
# 
# 4 char Read strand: “+” or “-”
# 
# 5 string Read name
# 
# 6 int Read start coordinate of subregion which is used to construct contig/unitig
# 
# 7 int Read end coordinate of subregion which is used to construct contig/unitig
# 
# 8 id:i:int Read ID
# 
# 9  HG:A:char Haplotype status of read. HG:A:a, HG:A:p, HG:A:m indicate read is non-binnable, father/hap1-specific and mother/hap2-specific, respectively.

cat ("\n\n --> importing GFA file         \n\n")


gfa = fread(gfa.file, header=FALSE,  fill=T)
names(gfa) = c("V1","contig","start","read.strand","read","read.start" ,"read.end","read.id","hap.status")
print(gfa)

gfN = gfa[ gfa$"V1"== "S", 1:5]
gfN$"coverage"= as.numeric( gsub("rd:i:","",gfN$"read"))
gfN$"ctg.len" = as.numeric( gsub("LN:i:","",gfN$"read.strand"))

cat ("\n\n --> GFA contigs statistics         \n\n")
print(gfN)


gfn = setDT(data.frame("contig"=gfN$"contig", "coverage"=gfN$"coverage",stringsAsFactors=FALSE))


cat ("\n\n --> Chromosome size         \n\n")

if ( is.na(agp.file)==FALSE  ) 
{
cat ("      --> from AGP file         \n\n")
    agp = fread( agp.file, verbose=FALSE )
    names(agp)=c("ref.chr","ref.start","ref.end","part_number","component_type","contig","q.start","q.end","strand")
    print(agp)
    
    
    agp$"ref.chr" = gsub("_RagTag","",agp$"ref.chr"  )
    gap = agp[ agp$"strand" =="align_genus",]
    agp = agp[ agp$"strand" !="align_genus",]
    
    
    chr = split(agp, agp$"ref.chr")
    chr = lapply(chr, function(x) data.frame("ref.chr"=x$"ref.chr"[1], "chr.len"=max(x$"ref.end", stringsAsFactors=FALSE)))
    chr = setDT( do.call(rbind,chr) )
    } else {
        
        cat ("      --> no AGP file  --> using contig list       \n\n")  
        agp = setDT(data.frame("ref.chr"=gfN$"contig", "ref.start"=1,"ref.end"=gfN$"ctg.len","part_number"=1,"component_type"="W","contig"=gfN$"contig","q.start"=1,"q.end"=gfN$"ctg.len","strand"="+",stringsAsFactors=FALSE))
        gap = agp[ agp$"strand" =="align_genus",]
        chr = setDT(data.frame("ref.chr"=gfN$"contig", "chr.len"=gfN$"ctg.len",stringsAsFactors=FALSE))
        print(agp)
        }
        
        
if ( is.na(ref.file)==FALSE & is.na(agp.file)==FALSE ) 
{
cat ("\n\n --> Reference fasta file         \n\n")

# REFERENCE FILE
ref.fai =paste0(ref.file,".fai")
if (file.exists(ref.fai)==FALSE ) { system( paste0( "samtools faidx ",ref.file )) }

ref=fread(ref.fai)[,1:2]
names(ref)=c("chr","ref.len")
ref= as.data.frame(ref)
#revert order
ref = ref[ order(rev(ref$"chr")) ,]
ref$"chr" = gsub("_RagTag","",ref$"chr")

# sort if the have progressive name
ref$"prefix" = substr(ref$"chr",1,2)
ref$"num"= as.numeric( gsub("CHR|Chr|LG|CM0|C|l|c|ptg|PTG|h1g|h2g|","",ref$"chr" ))
#print(ref)
if ( length(unique(ref$"prefix"))==1 & sum(is.na(ref$"num"))==0 & any(order(ref$"ref.len") != 1:nrow(ref)) & any(order(ref$"ref.len") != nrow(ref):1  ) )    ref = ref[ rev(order(ref$"num")) ,]
ref$"prefix" = NULL
ref$"num" = NULL

ref = ref[ ref$"ref.len" >= chr.filt ,]

chr_prefix = unique(substring(ref$"chr",1,3))

ref.name = system( paste(' grep ">" ',ref.file), inter=TRUE) 
ref.df= data.frame( name=ref.name, stringsAsFactors=FALSE     )


if ( (chr_prefix %in% c("Chr","chr")==FALSE)  & sum(grepl("Chr|chr|LG|chl|Chl|Mit", ref.name))>1 )
    {
    cat("   --> Reference complete names:\n")
    print(ref.df)
    
    ref.words = strsplit(ref.name, " ")
    ref.all.words = sort(unique(unlist(ref.words)))
    
    
     REF.ALL.WORDS = list()
     for (i in ref.all.words) REF.ALL.WORDS[[i]] = sapply(ref.words, function(x) { sum(x==i) })
     
     REF.ALL.WORDS.sin = REF.ALL.WORDS [ which(sapply(REF.ALL.WORDS, function(x) ! all(x==1))) ] 
     REF.ALL.WORDS.dup = REF.ALL.WORDS [ which(sapply(REF.ALL.WORDS, function(x)   all(x==1))) ] 
     ref.dup = names(REF.ALL.WORDS.dup)
     ref.dup = grep("Chr|chr|LG",ref.dup, invert=T, value = TRUE)
     ref.dup = paste(ref.dup, collapse= "|")
    
    ref.df$"simp_name" = ref.df$"name"
    ref.df$"simp_name" = gsub(ref.dup, "", ref.df$"simp_name" )
    ref.df$"simp_name" = gsub("       |     |   |  |,|;", "  ", ref.df$"simp_name" )
    
    ref.edit = data.frame(do.call(rbind, strsplit(ref.df$"simp_name","  " )))
    for ( ii in names(ref.edit)) if ( all( as.character(ref.edit[,ii])==" ") | all(as.character(ref.edit[,ii])=="")) ref.edit[ii] = NULL
    
    w1.col = names( which(sapply(ref.edit,function(x) length(grep(chr_prefix,x)))>0))
    w2.col = names( which(sapply(ref.edit,function(x) length(grep(chr_prefix,x)))==0))
    
    name.ref.edit = names(ref.edit)
    name.ref.edit[ name.ref.edit == w1.col] = "chr"
    name.ref.edit[ name.ref.edit == w2.col] = "chr.rename"
    names(ref.edit) = name.ref.edit
    ref.edit$"chr"        = gsub(">","", ref.edit$"chr")
    ref.edit$"chr.rename" = gsub("Chromosome |Chromosome|chromosome |chromosome","Chr", ref.edit$"chr.rename")
    
    ref.edit$"nchar" = nchar(ref.edit$"chr.rename")
    min_char = min(ref.edit$"nchar" )
    max_char = max(ref.edit$"nchar" )
    if( min_char == max_char-1  ) ref.edit[ ref.edit$"nchar" == min_char, ]$"chr.rename" = gsub("Chr","Chr0",ref.edit[ ref.edit$"nchar" == min_char, ]$"chr.rename")
    
    ref.edit$"nchar" = NULL
          
    cat('\n  --> Chromosomes renamed as below ! \n\n')
    print(ref.edit)
    
    ref.edi2 = ref.edit
    names(ref.edi2) = c("ref.chr","chr.rename")
    
    chr = merge(chr,ref.edi2, by="ref.chr" , all.x=T) ; w_rep = which( ! is.na(chr$"chr.rename")) ; chr$"ref.chr" [w_rep]     =chr$"chr.rename"[w_rep] ; chr$"chr.rename"=NULL
    agp = merge(agp,ref.edi2, by="ref.chr" , all.x=T) ; w_rep = which( ! is.na(agp$"chr.rename")) ; agp$"ref.chr" [w_rep]     =agp$"chr.rename"[w_rep] ; agp$"chr.rename"=NULL
    gap = merge(gap,ref.edi2, by="ref.chr" , all.x=T) ; w_rep = which( ! is.na(gap$"chr.rename")) ; gap$"ref.chr" [w_rep]     =gap$"chr.rename"[w_rep] ; gap$"chr.rename"=NULL

    
    }
}



chr = chr[ chr$"chr.len" >= chr.filt ,]
agp = agp[ agp$"ref.chr" %in% chr$"ref.chr" ,]
gap = gap[ gap$"ref.chr" %in% chr$"ref.chr"  ,]

# Sorting first chromosomes by name, then contings by size, if present
if( length(grep("CHR|Chr|LG|CM0|C|RagTag",chr$"ref.chr"))>0)
{
    ctg = chr[ -grep("CHR|Chr|LG|CM0|C|RagTag",chr$"ref.chr") , ] # Identify Chromosomes or Linkage Groups or 
    chh = chr[  grep("CHR|Chr|LG|CM0|C|RagTag",chr$"ref.chr") , ]
    ctg = ctg[ rev(order(ctg$"chr.len" )) , ]
    
    # Filt by CHR/contig size
    chr = rbind(chh,ctg)
    }



      
mat = chr$"chr.len"
names(mat) = chr$"ref.chr"
mat = mat [ rev(names(mat)) ]

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

cat(' \n --->   PLOT 1a) Contig coverage profile as reported on the GFA file  \n\n')


pdf( paste0("Coverage_Plot_01a.Contig_coverage_profile.",chr.filt.char,OUT,"from_gfa.pdf" )  ,width=plot.w, height=plot.h) ; 
par(mfcol=c(1    ,1    ), mar=c(2,6,4,1) , oma=c(1,3,1,1) )



par(mfcol=c(1    ,1    ), mar=c(2,7,4,1) , oma=c(1,4,1,1) )
bar = barplot(mat, horiz=T, beside=T, las=1, col="white",xlim=c(-max(mat)*0.07,max(mat)*1.15), main="Contig Coverage profile from .gfa assembly",cex.names=1.5, cex.main=2.5)
bbb = data.frame( "ref.chr"=names(mat), "chr.len"=mat, pos = bar[,1])
bbb$start = 1


AGP = merge(agp, gfN, by="contig",all.x=TRUE)
GAP = merge(agp, gfN, by="contig",all.x=TRUE)

AGP = merge(AGP, bbb, by="ref.chr",all.x=TRUE)
GAP = merge(GAP, bbb, by="ref.chr",all.x=TRUE)


max_cov = max(AGP$coverage)
q90_cov = quantile(AGP$coverage,probs=seq(0,1,by=0.01))["90%"]
q95_cov = quantile(AGP$coverage,probs=seq(0,1,by=0.01))["95%"]
q99_cov = quantile(AGP$coverage,probs=seq(0,1,by=0.01))["95%"]


COL_1 = "cornflowerblue"
COL_2 = "black"


high_cov = max_cov
names(high_cov) = "100%"
AGP$"col" = COL_1
AGP$"COV" = AGP$"coverage"

if( max_cov/q99_cov  > 2 )  { high_cov = q99_cov ; AGP[ AGP$"coverage" > high_cov, ]$"col" = COL_2;  AGP[ AGP$"coverage" > high_cov, ]$"COV" = high_cov }
if( q99_cov/q95_cov  > 2 )  { high_cov = q95_cov ; AGP[ AGP$"coverage" > high_cov, ]$"col" = COL_2 ;  AGP[ AGP$"coverage" > high_cov, ]$"COV" = high_cov }

AGP$"cov" = AGP$"COV" / high_cov
AGP$"p0" = AGP$"pos" -0.5
AGP$"p1" =AGP$"pos" -0.5 + AGP$"cov" 


cat('   -> Contig Coverage Statistics \n')
summary(AGP$coverage)
cat('  \n')

segments(bbb$"start", bbb$"pos"+0.3, bbb$"chr.len",bbb$"pos"+0.3, col="gray95", lwd=1)
segments(bbb$"start", bbb$"pos"+0.1, bbb$"chr.len",bbb$"pos"+0.1, col="gray95", lwd=1)
segments(bbb$"start", bbb$"pos"-0.1, bbb$"chr.len",bbb$"pos"-0.1, col="gray95", lwd=1)
segments(bbb$"start", bbb$"pos"-0.3, bbb$"chr.len",bbb$"pos"-0.3, col="gray95", lwd=1)



for (j in 1:nrow(AGP)) { polygon( c(AGP$"ref.start"[j], AGP$"ref.end"[j],AGP$"ref.end"[j],AGP$"ref.start"[j]) ,c( AGP$"p1"[j] , AGP$"p1"[j], AGP$"p0"[j],  AGP$"p0"[j]) , lwd=1 , col=AGP$"col"[j] , border=NA)    }


if ( nrow(GAP)> 0)  segments(GAP$"ref.end", GAP$"pos"-0.5, GAP$"ref.end",GAP$"pos"+0.5, col="gray", lwd=1)


# redraw a box
par(new = T,lwd = 2)
bar = barplot(mat, horiz=T, beside=T, las=1, col=NA,xlim=c(-max(mat)*0.07,max(mat)*1.15), cex.names=1.5, cex.main=2.5)


# legeng

cov_range = names(high_cov)
cov_range2 =paste0("Coverage <",cov_range,"tile (1-",high_cov,")"  )
cov_range3 =paste0("Coverage >",cov_range,"tile (",high_cov+1,"-",max_cov,")"  )

if (cov_range == "100%") LEGEND = data.frame( "name"=c( paste0("Coverage (1-",max_cov,")") ,"GAP")    , pch=c(15,3)    , col=c(COL_1,"gray"),        pt.bg=c(COL_1,"gray"),       cex=2.0, pt.cex=2.2, stringsAsFactors=FALSE)
if (cov_range != "100%") LEGEND = data.frame( "name"=c( cov_range2,cov_range3,"GAP") , pch=c(15,15,3) , col=c(COL_1,COL_2,"gray"),  pt.bg=c(COL_1,COL_2,"gray"), cex=1.5, pt.cex=1.7, stringsAsFactors=FALSE)


legend("bottomright", LEGEND$"name", pch=LEGEND$"pch",  col=LEGEND$col , pt.bg=LEGEND$"pt.bg", cex=2, bg="white",pt.cex=2.2)

dev.off()




cat(' \n --->   PLOT 1b) Contig coverage profile are reported on the GFA file with data \n\n')


pdf( paste0("Coverage_Plot_01b.Contig_coverage_profile.",chr.filt.char,OUT,"with_names.from_gfa.pdf" )  ,width=plot.w, height=plot.h) ; 
par(mfcol=c(1    ,1    ), mar=c(2,6,4,1) , oma=c(1,3,1,1) )



par(mfcol=c(1    ,1    ), mar=c(2,7,4,1) , oma=c(1,4,1,1) )
bar = barplot(mat, horiz=T, beside=T, las=1, col="white",xlim=c(-max(mat)*0.07,max(mat)*1.15), main="Contig Coverage profile from .gfa assembly",cex.names=1.5, cex.main=2.5)
bbb = data.frame( "ref.chr"=names(mat), "chr.len"=mat, pos = bar[,1])
bbb$start = 1

segments(bbb$"start", bbb$"pos"+0.3, bbb$"chr.len",bbb$"pos"+0.3, col="gray95", lwd=1)
segments(bbb$"start", bbb$"pos"+0.1, bbb$"chr.len",bbb$"pos"+0.1, col="gray95", lwd=1)
segments(bbb$"start", bbb$"pos"-0.1, bbb$"chr.len",bbb$"pos"-0.1, col="gray95", lwd=1)
segments(bbb$"start", bbb$"pos"-0.3, bbb$"chr.len",bbb$"pos"-0.3, col="gray95", lwd=1)


for (j in 1:nrow(AGP)) { polygon( c(AGP$"ref.start"[j], AGP$"ref.end"[j],AGP$"ref.end"[j],AGP$"ref.start"[j]) ,c( AGP$"p1"[j] , AGP$"p1"[j], AGP$"p0"[j],  AGP$"p0"[j]) , lwd=1 , col=AGP$"col"[j] , border=NA)    }


AGP$"ref.mid" = (AGP$"ref.start"+AGP$"ref.end") / 2
AGP$"ctg_name_2" = AGP$"contig" 
AGP$"cex" =1

ctg_nchar = unique(nchar(AGP$"contig")) 
ctg_substr_5 = unique(substr(AGP$"contig",1,ctg_nchar-5) )
ctg_substr_6 = unique(substr(AGP$"contig",1,ctg_nchar-6) )
if ( length(ctg_substr_5) ==1 ) { ctg_substr = ctg_substr_5} else  ctg_substr = ctg_substr_6
# max chr > 250Mb  -> min ctg = chr 8.8MB --> ratio 1/28th

substr_thrshld =  max(mat)*1.15 / 28
ctg_below_thrshld =  which( AGP$"ctg.len" < substr_thrshld )

if (any (  AGP$"ctg.len" < substr_thrshld )) { AGP$"ctg_name_2"[ ctg_below_thrshld ] = gsub(ctg_substr,".", AGP$"ctg_name_2"[ ctg_below_thrshld ]) ; AGP$"cex"[ ctg_below_thrshld ] = 0.6   }
#if (any (  AGP$"ctg.len" < substr_thrshld )) { AGP$"ctg_name_2"[ ctg_below_thrshld ] = gsub(ctg_substr,".", AGP$"ctg_name_2"[ ctg_below_thrshld ]) ; AGP$"cex"[ ctg_below_thrshld ] = 0.6   }



if (any (  AGP$"coverage" != AGP$"COV"  )) AGP[AGP$"coverage" != AGP$"COV" , ]$"p1" =AGP[AGP$"coverage" != AGP$"COV" , ]$"pos"



text( AGP$"ref.mid", AGP$"p1" , AGP$"ctg_name_2", cex=AGP$"cex" )   


if ( nrow(GAP)> 0)  segments(GAP$"ref.end", GAP$"pos"-0.5, GAP$"ref.end",GAP$"pos"+0.5, col="gray", lwd=1)


# redraw a box
par(new = T,lwd = 2)
bar = barplot(mat, horiz=T, beside=T, las=1, col=NA,xlim=c(-max(mat)*0.07,max(mat)*1.15), cex.names=1.5, cex.main=2.5)


# legeng

cov_range = names(high_cov)
cov_range2 =paste0("Coverage <",cov_range,"tile (1-",high_cov,")"  )
cov_range3 =paste0("Coverage >",cov_range,"tile (",high_cov+1,"-",max_cov,")"  )

if (cov_range == "100%") LEGEND = data.frame( "name"=c( paste0("Coverage (1-",max_cov,")") ,"GAP")    , pch=c(15,3)    , col=c(COL_1,"gray"),        pt.bg=c(COL_1,"gray"),       cex=2.0, pt.cex=2.2, stringsAsFactors=FALSE)
if (cov_range != "100%") LEGEND = data.frame( "name"=c( cov_range2,cov_range3,"GAP") , pch=c(15,15,3) , col=c(COL_1,COL_2,"gray"),  pt.bg=c(COL_1,COL_2,"gray"), cex=1.5, pt.cex=1.7, stringsAsFactors=FALSE)


legend("bottomright", LEGEND$"name", pch=LEGEND$"pch",  col=LEGEND$col , pt.bg=LEGEND$"pt.bg", cex=2, bg="white",pt.cex=2.2)

dev.off()



## high resolution plot

cat(' \n --->   PLOT 2) Hifi read density from .gfa assembly (',wind.char,' windows) \n\n')


gfL = setDT(data.frame("contig"=gfN$"contig", start=1, "ctg.len"=gfN$"ctg.len",stringsAsFactors=FALSE))


gfr = gfa [ gfa$"V1" == "A" , ]
#gfL = gfa [ gfa$"V1" == "L" , ]

gfR = setDT(data.frame("contig"=gfr$"contig", "start"=as.numeric(gfr$"start"), "end"=as.numeric(gfr$"start")+as.numeric(gfr$"read.end")-1, stringsAsFactors=FALSE))

write.table( gfL, file = "Contig_length.txt" , col.names=FALSE, row.names=F, quote=F, sep="\t")    
write.table( gfR, file = "Read_coordinates_on_contigs.bed" , col.names=FALSE, row.names=F, quote=F, sep="\t")    



     #   system( paste0( " bedtools makewindows -b Contig_length.txt -w 1000   -s 1000  > Contig_windows.1Kb.bed"))
     #   system( paste0( " bedtools makewindows -b Contig_length.txt -w 10000  -s 10000 > Contig_windows.10Kb.bed"))
        system( paste0( " bedtools makewindows -b Contig_length.txt -w ",wind," -s ",wind," > Contig_windows.",wind.char,".bed"))
     #  system( paste0( " bedtools intersect -c -a  Contig_windows.1Kb.bed   -b Read_coordinates_on_contigs.bed > Contig_vs_reads.1Kb_window.bed" ))
     #   system( paste0( " bedtools intersect -c -a  Contig_windows.10Kb.bed  -b Read_coordinates_on_contigs.bed > Contig_vs_reads.10Kb_window.bed" ))
        system( paste0( " bedtools intersect -c -a  Contig_windows.",wind.char,".bed -b Read_coordinates_on_contigs.bed > Contig_vs_reads.",wind.char,"_window.bed" ))



bed = fread(  paste0( "Contig_vs_reads.",wind.char,"_window.bed" ) )
names(bed) = c("contig","start","end","n")
bed$"mid" = (bed$"start"+bed$"end") / 2
bed$"start" = NULL
bed$"end" = NULL


# sum(AGP$"ref.start" < AGP$"ref.end" )
# sum(AGP$"q.start"   < AGP$"q.end" )


minicol = c("ref.chr","contig","cov","pos","ref.start","ref.end","q.start","strand")
miniAGP = AGP [, ..minicol]

bed = merge(miniAGP,bed, by="contig", all.x=TRUE)

# invert coordnates of "-" CTGs


w = which(bed$"strand" == "-" )
bed$"MID" = bed$"ref.start" + bed$"mid" 


if( length(w) > 0) bed$"MID"[w] =  bed$"ref.end"[w] -  bed$"mid"[w]


max_reads = max(bed$"n")
q90_cov = quantile(bed$"n",probs=seq(0,1,by=0.01))["90%"]
q95_cov = quantile(bed$"n",probs=seq(0,1,by=0.01))["95%"]
q99_cov = quantile(bed$"n",probs=seq(0,1,by=0.01))["95%"]

cat('   -> Hifi reads per ',wind.char,' Statistics  \n')
summary(bed$n)
cat('  \n')


high_reads = max_reads
names(high_reads) = "100%"
bed$"col" = COL_1
bed$"N" = bed$"n"

if( max_reads/q99_cov  > 2 )  { high_reads = q99_cov ; bed[ bed$"n" > high_reads, ]$"col" = COL_2 ;  bed[ bed$"n" > high_reads, ]$"N" = high_reads }
if( q99_cov/q95_cov    > 2 )  { high_reads = q95_cov ; bed[ bed$"n" > high_reads, ]$"col" = COL_2 ;  bed[ bed$"n" > high_reads, ]$"N" = high_reads }

bed$"nn" = bed$"N" / high_reads
bed$"p0" = bed$"pos" -0.5
bed$"p1" = bed$"pos" -0.5 + bed$"nn" 
bed$"pn" = bed$"pos" -0.5 + bed$"nn" 

bed = bed[ order(bed$"ref.chr",bed$"contig", bed$"MID"),]


pdf( paste0("Coverage_Plot_02.Reads_used_to_build_contigs_per_100Kb_window.",chr.filt.char,OUT,"from_gfa.pdf" )  ,width=plot.w, height=plot.h) ; 
par(mfcol=c(1    ,1    ), mar=c(2,6,4,1) , oma=c(1,3,1,1) )


par(mfcol=c(1    ,1    ), mar=c(2,7,4,1) , oma=c(1,4,1,1) )
bar = barplot(mat, horiz=T, beside=T, las=1, col="white",xlim=c(-max(mat)*0.07,max(mat)*1.15), main=paste0("Hifi read density from .gfa assembly (",wind.char,"windows)"),cex.names=1.5, cex.main=2.5)
bbb = data.frame( "ref.chr"=names(mat), "chr.len"=mat, pos = bar[,1])
bbb$start = 1


bs = split(bed, bed$"ref.chr")
#for (j in 1:length(bs)) { lines( bs[[j]]$"MID" , bs[[j]]$"p1", , lwd=1 , col="black" )  }

segments(bed$"MID"  , bed$"pos"-0.5, bed$"MID"  ,bed$"pn", col= COL_1  , lwd=4 , lend=1, border=NULL)

# redraw a box
par(new = T,lwd = 2)
bar = barplot(mat, horiz=T, beside=T, las=1, col=NA,xlim=c(-max(mat)*0.07,max(mat)*1.15),cex.names=1.5, cex.main=2.5 )
if ( nrow(GAP)> 0)  segments(GAP$"ref.end", GAP$"pos"-0.5, GAP$"ref.end",GAP$"pos"+0.5, col="black", lwd=1)


read_range = names(high_reads)
read_range2 =paste0("Coverage <",read_range,"tile (1-",high_reads,")"  )
read_range3 =paste0("Coverage >",read_range,"tile (",high_reads+1,"-",max_reads,")"  )

if (read_range == "100%") LEGEND = data.frame( "name"=c( paste0("Coverage (1-",max_cov,")") ,"GAP")    , pch=c(15,3)    , col=c(COL_1,"gray"),        pt.bg=c(COL_1,"gray"),       cex=2.0, pt.cex=2.2, stringsAsFactors=FALSE)
if (read_range != "100%") LEGEND = data.frame( "name"=c( read_range2,read_range3,"GAP") , pch=c(15,15,3) , col=c(COL_1,COL_2,"gray"),  pt.bg=c(COL_1,COL_2,"gray"), cex=1.5, pt.cex=1.7, stringsAsFactors=FALSE)


legend("bottomright", LEGEND$"name", pch=LEGEND$"pch",  col=LEGEND$col , pt.bg=LEGEND$"pt.bg", cex=2, bg="white",pt.cex=2.2)

dev.off()



END.time =Sys.time() 
TOT.time =difftime(END.time, START.time, units="mins")


cat( paste( "==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n\n" ) )




