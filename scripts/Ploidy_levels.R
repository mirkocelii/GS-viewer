
#!/bin/R

# R

args <- commandArgs(TRUE)



   if("--help" %in% args) 
    {
        cat('\n')
        cat(
        " Plot ploidy level from a RagTag scaffold \n",
        " \n",
        "USAGE: Rscript Ploidy_levels.R --scaffold=ragtag.scaffold.fasta --ref=reference_genome.fasta \n",
        " \n",
        "Arguments:\n",
        "  --scaffold=         |  ragtag.scaffold.fasta file                      ","\n",
        "  --ref=              |  Reference FASTA file                           [ optional, default absent   used for importing real Chr names   ]","\n",
        "  --overlap=          |  Max contig overlap for beein at the same level [ 0-to-1, default 0.02                                           ] ","\n",
        "  --ploidy=           |  Set ploidy level                               [ default determined with 90-95% of reference covered            ]  ","\n",
        "  --chr.rename=       |  Replace GenBank identifier with Chr/Contig     [ default FALSE                                                  ]    ","\n",
        "  --chr.filt=         |  minimal Chr size plotted                       [ default = 1Mb/1000Kb/1000000; numeric,Kb,Mb character accepted ]  ","\n",
        "  --hap.fa=           |  output fasta file of each haplotype set draft  [ default absent   ] ","\n",
        "  --plot.h=           |  Height in inches of PDF output files           [ default 18       ]  ","\n",
        "  --plot.w=           |  Width in inches of PDF output files            [ default 35       ]  ","\n",
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
cat("###                             Ploidy_levels.r                            ###\n")
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
    
	fas.file        = as.character(ARG[ ARG$X1=="--scaffold" ,]$"X2")
	ref.file        = as.character(ARG[ ARG$X1=="--ref"      ,]$"X2")
	chr.filt        = as.character(ARG[ ARG$X1=="--chr.filt" ,]$"X2")
	chr.rename      = as.character(ARG[ ARG$X1=="--chr.rename" ,]$"X2")
	ploidy          = as.character(ARG[ ARG$X1=="--ploidy"   ,]$"X2")
	plot.h          = as.character(ARG[ ARG$X1=="--plot.h"   ,]$"X2")
	plot.w          = as.character(ARG[ ARG$X1=="--plot.w"   ,]$"X2")
	OVERLAP         = as.character(ARG[ ARG$X1=="--overlap"  ,]$"X2")
	hap.fa          = as.character(ARG[ ARG$X1=="--hap.fa"   ,]$"X2")
		
	plot.h = as.numeric(plot.h)
	plot.w = as.numeric(plot.w)
	OVERLAP = as.numeric(OVERLAP)


#  DEFAULT VALUES

    if ( length(chr.rename)==0) chr.rename = FALSE  else if ( chr.rename=="--chr.rename") chr.rename = TRUE
    if ( length(chr.filt)==0  ) chr.filt   = 1 * 1e6
    if ( length(plot.h)==0    ) plot.h      = 18
    if ( length(plot.w)==0    ) plot.w      = 35
    if ( length(OVERLAP)==0   ) OVERLAP     = 0.02
    if ( length(ploidy)==0    ) ploidy      = NA
    if ( length(hap.fa)==0    ) hap.fa      = FALSE else if ( hap.fa=="--hap.fa") hap.fa = TRUE
    
 
# TEST PARAMETERS for interactive session
if( 2 ==1 )
{

# Rscript "$Ploidy_levels" --ref=$ref --scaffold=RAGTAG/ragtag_scaffold.hifiasm.l2_hap1+hap2.100Kb.confusa/ragtag.scaffold.fasta  --chr.rename --chr.filt=1Mb --hap.fa=TRUE
# ref="Acacia_confusa.ASM3707449v1/data/GCA_037074495.1/GCA_037074495.1_ASM3707449v1_genomic.fna"; 

    fas.file="ragtag.scaffold.fasta"
    ref.file="/ibex/scratch/projects/c2067/celiim/Salvadora_persica/REF.Dobera_glabra/Dobera_genome_11CHR.fasta"
    chr.rename=FALSE
    chr.filt="1Mb"
    ploidy = 4
    plot.h = 18
    plot.w = 35    
    out    = NA
    OVERLAP = 0.02
    hap.fa = FALSE

  }


cat('\n')
cat("############################################################################\n")
cat('------------------------ Summary of Parameters -----------------------------\n')
cat("############################################################################\n")
cat('\n')
    cat("   scaffold   = ", fas.file,"\n")
    cat("   ref.file   = ", ref.file,"\n")
    cat("   chr.rename = ", chr.rename,"\n")
    cat("   chr.filt   = ", chr.filt,"\n")
    cat("   ploidy     = ", ploidy,"\n")
    cat("   plot.h     = ", plot.h,"\n")
    cat("   plot.w     = ", plot.w,"\n")
    cat("   overlap    = ", OVERLAP,"\n")
    cat("   hap.fa     = ", hap.fa,"\n")

cat('\n')

library(data.table)
options(width=250)
options(scipen=999)

user.ploidy = ploidy

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


# cd  ragtag_scaffold.hifiasm.hap1.100Kb.jujuba
# cd /ibex/scratch/projects/c2067/celiim/Ziziphus_spina-christi/ragtag_scaffold.hifiasm.both_haps.100Kb.jujuba

ref.fai  =  paste0(ref.file,".fai")
if (file.exists(ref.fai)==FALSE ) { system( paste0( "samtools faidx ",ref.file )) }


fai.file = gsub(".fasta",".fasta.fai",fas.file)
agp.file = gsub(".fasta",".agp",fas.file)
asm.file = gsub(".fasta",".asm.paf",fas.file)
log.file = gsub(".fasta",".asm.paf.log",fas.file)


REF = fread(ref.fai , fill=T)[,1:2]
FAI = fread(fai.file , fill=T)[,1:2]
AGP = fread(agp.file , fill=F)
ASM = fread(asm.file , fill=T)
LOG = system( paste( "cat ", log.file ), inter=T)


# select assembly file from LOG
LOG =  tail( unlist( strsplit( grep("fasta",LOG, value=TRUE), " ")),1)



agp.name = c("ref.chr","ref.start","ref.end","part_number","component_type","query","q.start","q.end","strand")
asm.name = c("query","q.len","q.start","q.end","strand","ref.chr","ref.len","ref.start","ref.end","l1","l2","V12","V13","V14","V15","V16","V17")
ref.name = c("ref.chr","ref.len")
fai.name = c("ref.chr","len")


n_paf = ncol(ASM)
if( all(n_paf==18)) { asm.name = c(asm.name,"V18")}



names(AGP) = agp.name
names(ASM) = asm.name
names(REF) = ref.name
names(FAI) = fai.name

AGP$"ref.chr" = gsub("_RagTag","",AGP$"ref.chr")
FAI$"ref.chr" = gsub("_RagTag","",FAI$"ref.chr")

GAP = AGP[ AGP$"strand" =="align_genus",]
AGP = AGP[ AGP$"strand" !="align_genus",]

ASM$"ctg.len"  =  ASM$"q.end"-ASM$"q.start"+1; 
ASM$"scaf.perc"= round( ASM$"ctg.len"/ASM$"q.len",2);
ASM$"ctg.ref"  =  ASM$"ref.end"-ASM$"ref.start"+1

AGP$"q.len" = as.numeric(AGP$"q.end")- as.numeric(AGP$"q.start")+1;


ASM = ASM[ order(ASM$"ref.chr"),] 



# sort if the have progressive name
REF$"prefix" = substr(REF$"ref.chr",1,2)
REF$"num"= as.numeric( gsub("CHR|Chr|LG|CM0|C|l|c|ptg|PTG|h1g|h2g|","",REF$"ref.chr" ))
#print(ref)
if ( length(unique(REF$"prefix"))==1 & sum(is.na(REF$"num"))==0 & any(order(REF$"ref.len") != 1:nrow(REF)) & any(order(REF$"ref.len") != nrow(REF):1  ) )    REF = REF[ rev(order(REF$"num")) ,]
REF$"prefix" = NULL
REF$"num" = NULL

REF = REF[ REF$"ref.len" >= chr.filt ,]


chr_prefix = unique(substring(REF$"ref.chr",1,3))

ref.name = system( paste(' grep ">" ',ref.file), inter=TRUE) 
ref.df= data.frame( name=ref.name, stringsAsFactors=FALSE     )


if ( all(chr_prefix %in% c("Chr","chr")==FALSE)  & sum(grepl("Chr|chr|LG|chl|Chl|Mit", ref.name))>1 )
    {
    cat("\n\n   --> Reference complete names:\n")
    print(ref.df)
    
    if ( chr.rename == TRUE  )
    {
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
          
        cat('\n  --> Chromosomes renamed ! \n\n')
        print(ref.edit)
        
        ref.edi2 = ref.edit
        names(ref.edi2) = c("ref.chr","chr.rename")
        
        REF = merge(REF,ref.edi2, by="ref.chr" , all.x=T) ; w_rep = which( ! is.na(REF$"chr.rename")) ; REF$"ref.chr" [w_rep] =REF$"chr.rename"[w_rep] ; REF$"chr.rename"=NULL
        FAI = merge(FAI,ref.edi2, by="ref.chr" , all.x=T) ; w_rep = which( ! is.na(FAI$"chr.rename")) ; FAI$"ref.chr" [w_rep] =FAI$"chr.rename"[w_rep] ; FAI$"chr.rename"=NULL
        AGP = merge(AGP,ref.edi2, by="ref.chr" , all.x=T) ; w_rep = which( ! is.na(AGP$"chr.rename")) ; AGP$"ref.chr" [w_rep] =AGP$"chr.rename"[w_rep] ; AGP$"chr.rename"=NULL
        GAP = merge(GAP,ref.edi2, by="ref.chr" , all.x=T) ; w_rep = which( ! is.na(GAP$"chr.rename")) ; GAP$"ref.chr" [w_rep] =GAP$"chr.rename"[w_rep] ; GAP$"chr.rename"=NULL
        ASM = merge(ASM,ref.edi2, by="ref.chr" , all.x=T) ; w_rep = which( ! is.na(ASM$"chr.rename")) ; ASM$"ref.chr" [w_rep] =ASM$"chr.rename"[w_rep] ; ASM$"chr.rename"=NULL
        } else { cat('\n  --> If you want to rename chromosomes  ===>> USE OPTION  --chr.rename !! \n\n') }
        

        
        }
        

 cat('\n  --> Chromosome set: \n\n')

print(REF)

ASM = ASM[ ASM$"ref.chr" %in% REF$"ref.chr", ]
AGP = AGP[ AGP$"ref.chr" %in% REF$"ref.chr", ]
GAP = GAP[ GAP$"ref.chr" %in% REF$"ref.chr", ]
FAI = FAI[ FAI$"ref.chr" %in% REF$"ref.chr", ]



ASMs = split(ASM,ASM$"ref.chr")
AGPs = split(AGP,AGP$"ref.chr")
GAPs = split(GAP,GAP$"ref.chr")
REFs = split(REF,REF$"ref.chr")
FAIs = split(FAI,FAI$"ref.chr")


RR.df = list()
AS.df = list()
CT.df = list()


pdf( paste0("Ploidy_Levels_01.Chromosome_detail_ploidy_levels.pdf"), width=plot.w, height=plot.h )

for ( i in names(ASMs))
{
    
    # i = names(ASMs)[1]
    
    asm = ASMs[[i]]
    agp = AGPs[[i]]
    gap = GAPs[[i]]
    ref = REFs[[i]]
    fai = FAIs[[i]]
    
    pla = asm[   asm$"query" %in% agp$"query",]
    unp = asm[ ! asm$"query" %in% agp$"query",]   
    
    
      cat("\n",i,"\t")
    a = asm
    g = agp
    u = unp
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
        a = a[ order(a$"ref.start"),]          

        # sort by biggest contig
        as_max = sapply(as, function(x) sum(x$"ctg.len") )
        as = as [ rev(order(as_max)) ]
        
      #  OVERLAP = 0.02
       # All position occupied in yhr chromosomes
        a2 = data.frame(a)
        ALL_pos_occupied = lapply( 1:nrow(a2), function(i) { seq(a2[i, "ref.start"], a2[i, "ref.end"])   } )   
        ALL_pos_occupied = sort(unique(unlist(ALL_pos_occupied)))      
        rm(a2)
        
        length(ALL_pos_occupied) /  max(ALL_pos_occupied)
          
        RR = setDT(data.frame( pos=1:max(a$"ref.len"), stringsAsFactors=FALSE))
      #  RR = setDT(data.frame( pos=ALL_pos_occupied  , stringsAsFactors=FALSE))
        RR$"n1"= 0
        # n = ploidy level 1
        
        log_ctg = list()
        
        cat(length(as),"\n")
        k=0
                
          length(ALL_pos_occupied) /  max(a$"ref.end")
        
       # NNN=length(as)
        for ( ctg in names(as))
      # NNN = 5; for ( ctg in names(as)[1:NNN]) # for debugging
        {
        #  ctg = names(as)[NNN+1]
        #  ctg = names(as)[1]
          k=k+1
          cat( k,"" )
          
          aa = as[[ctg]]
          aa = aa[ aa$"ctg.len" > 5000, ]
          aa = aa[ order(aa$"ref.start"),]    
          
          if( nrow(aa)> 0)
          {
          
              ctg_len = a$"q.len"[1]
              
              N_lev = ncol(RR)-1
              
              bb=as.data.frame(aa[, 8:9])
              
              # count REF coordinates occupied by this contig
              aa_pos = lapply( 1:nrow(bb), function(i) { seq(bb[i, "ref.start"], bb[i, "ref.end"])   } ) # verified       
              aa_pos = sort( unique( unlist(aa_pos) ))
              qq_ref = quantile(aa$"ctg.ref", 0:20/20)[c(2,20)] 
              qq_pos = quantile(aa_pos, 0:20/20)[c(2,20)]   
              qq_len = qq_pos[ "95%"] - qq_pos[ "5%"]

#               {
#                   print(summary(aa_pos))
#                   print(quantile(aa_pos, 0:20/20))    
#                   cat( "Max - min = ", summary(aa_pos)[6] - summary(aa_pos)[1], "\n") 
#                   cat( "Q95 - Q05 = ", qq_len, "\n") 
#                   cat( "CTG bp/ref= ",length(aa_pos), "\n") 
#                   cat( "CTG len   = ",ctg_len, "\n") 
#                   cat( "% pos/len = ",length(aa_pos)/ctg_len, "\n") 
#                   }
              
              # Loop untill you can place the contig  on  free reference genome level      
              for ( i_lev in 1:length(as)  )
              {
               
                  ni_lev = paste0("n", i_lev)
                  if ( is.null(RR[[ ni_lev ]] )==TRUE) RR[[ ni_lev ]]=0 # create the level if not present
                  
                  # free regions in the level   = 0
                  # previously occupied regions = 1
                  
                  # check pos occipied
                  ref_occup.bp = sum( RR[[ ni_lev ]] [ aa_pos ])                  # bp occupied
                  ref_occup.pc = round( ref_occup.bp/length(aa_pos),3 )
                  
                  # check whole regione
                  ref_region.bp = sum( RR[[ ni_lev ]] [ (qq_pos[ "5%"] : qq_pos[ "95%"])])
                  ref_region.pc = round( ref_region.bp/length(aa_pos),3 )                 
                  
                  
                #  cat( ref_occup.pc, "")
                  
                  if ( ref_occup.pc <= OVERLAP)
                  {
                     RR[[ ni_lev ]] [ aa_pos ] =1
                     as[[ctg]]$"level" = ni_lev
                     as[[ctg]]$"ploidy" = i_lev
                     as[[ctg]]$"cum.sum.ref" = sum(as[[ctg]]$"ref.width")
                     as[[ctg]]$"ploidy_bp"   = length(aa_pos)

                     log_ctg [[ ctg ]] = data.frame("ctg"=ctg, "level"=ni_lev, min=min(aa_pos), max=max(aa_pos), "len.bp"=length(aa_pos), width=length(aa_pos) ,"bp.overlap"=ref_occup.bp, "perc.overlap"=ref_occup.pc,  stringsAsFactors=FALSE)
                    cat( i_lev, "; ")
                    #   cat("\n")
                     break 
                     
                     }
                  }
                  }
                  }
          log_ctg = do.call(rbind,log_ctg)
          
          # FILT OUTPUT 
          
          AS = lapply(as, function(x) { x= x[ x$"ctg.len" > 20000,];  x$"ref.start.min" = min(x$"ref.start") ; x$"ref.end.max" = max(x$"ref.end")  ; return(x) } )
          AS = AS[  which( sapply(AS, function(x) { nrow( x[ x$"ctg.len" > 20000,] ) } ) >0 ) ]
          AS=do.call(rbind,AS)

          AS$"ploidy" = as.numeric(gsub("n","",AS$"level"))
          AS$"ref.mid" = ( AS$"ref.start.min" + AS$"ref.end.max" ) /2
          AS$"ref.len.filt" = AS$"ref.end.max" - AS$"ref.start.min" 

        
        transparent_blue <- rgb (0, 0, 1, alpha = 0.3)
        transparent_yell <- rgb (1, 1, 0, alpha = 0.4)
          mini.col = c("query","q.len","ref.start.min","ref.end.max","ploidy","level")
          pol = AS[ , ..mini.col]  
          pol = unique(pol)  
          pol = pol[ order(pol$"ref.start.min"),]    
          pol = pol[ order(pol$"q.len"),]    
          pol$"scaffold" = "-"
          pol$"col" = "-"
          pol[  which(   pol$"query" %in% agp$"query") , ]$"scaffold" = "placed"
          pol[  which( ! pol$"query" %in% agp$"query") , ]$"scaffold" = "unplaced"
          pol[  which(   pol$"query" %in% agp$"query") , ]$"col" = transparent_blue
          pol[  which( ! pol$"query" %in% agp$"query") , ]$"col" = transparent_yell
          pol$"ref.mid" = ( pol$"ref.start.min" + pol$"ref.end.max" ) /2

          Ymax = max(AS$"ploidy")*1.1

    #    AS = AS[ AS$"query" %in% names(as)[1:NNN],] ; pol = pol[ pol$"query" %in% names(as)[1:NNN],] ; # for debugging
       
        par(mfcol=c(1    ,1    ), mar=c(5,7,4,1) , oma=c(1,4,1,1) )
        plot(c(0, chr.max ), c(0, 0), type="l", las=1, xlim=c(-chr.max/10,chr.max), ylim=c(-Ymax/10, Ymax), cex.axis=2, main=i, cex.main=3, ylab="Ploidy level", xlab= paste(i, "coordinates"),cex.lab=2)

        bottom5.pc = quantile(AS$"cum.sum.ref",seq(0,1,0.05))["5%"]


        AS = AS[ AS$"ctg.len"  > 20000, ]
        A2 = AS[ AS$"ref.len.filt"  > bottom5.pc , ]
        A3 = AS[ AS$"ref.len.filt" <= bottom5.pc , ]

 
        segments( AS$"ref.start" , AS$"ploidy", AS$"ref.start"  , AS$"ploidy", lwd=10 , col=AS$"col") 
            text( A2$"ref.mid"   , A2$"ploidy"-0.2, A2$"query"  ) 
            text( A3$"ref.mid"   , A3$"ploidy"-0.2, "*" , cex=2 ) 
       
   for(n in 1:nrow(pol)) polygon( c(pol$"ref.start.min"[n], pol$"ref.end.max"[n], pol$"ref.end.max"[n],pol$"ref.start.min"[n]) ,c(pol$"ploidy"[n]-0.3, pol$"ploidy"[n]-0.3, pol$"ploidy"[n]+0.2, pol$"ploidy"[n]+0.2 ), lwd=1 ,col=pol$"col"[n])
         
# legeng
LEGEND = data.frame( "name"=c("Scaffolded contigs","Unplaced contigs","short contigs","FW orientation","REV orientation") , pch=c(22,22,8,16,16) , col=c("black","black","black","darkred","blue4"),  pt.bg=c(transparent_blue,transparent_yell,"black","darkred","blue4"), stringsAsFactors=FALSE)

legend("topright", LEGEND$"name", pch=LEGEND$"pch",  col=LEGEND$col , pt.bg=LEGEND$"pt.bg", cex=2, bg="white",pt.cex=2.2)
          
 
          cat("\n")
          RR.df[[ i]]  = RR
          AS.df[[ i]]  = AS
          CT.df[[ i]]  = log_ctg
                   }

}

dev.off()

# correction step

# cat('\n  --> CORRECTION \n\n')

      
 cat('\n  --> PLOIDY STATS \n\n')


mat = REF$"ref.len"
names(mat) = REF$"ref.chr"
names(mat) = gsub("_RagTag","",names(mat) )
mat = mat [ rev(names(mat)) ]


bar = barplot(mat, horiz=T, beside=T, las=1, col=gray(0.90) ,xlim=c(-max(mat)*0.07,max(mat)*1.15),lend=1, cex.axis=3, main=i, cex.main=3, ylab="Ploidy level", xlab= paste(i, "coordinates"),cex.lab=3)
bbb = data.frame( "ref.chr"=names(mat), pos = bar[,1])

CT.DF = do.call(rbind,CT.df)
AS.DF = do.call(rbind,AS.df)
AS.DF = merge(AS.DF,bbb, by="ref.chr")
max.ploidy = max(AS.DF$"ploidy")

AS.DF$"ploidy_level" = AS.DF$"ploidy"/max.ploidy
AS.DF = AS.DF[ AS.DF$"ctg.len" > 5000, ]

# check highest ploify

AS.pl = split(AS.DF, AS.DF$"ploidy")

cc=c("ref.chr","query","ploidy_bp")
AS.pl_stat = lapply(AS.pl, function(x) { y=x[,..cc]; y=y[!duplicated(y),]; z=data.frame("level"=x$"level"[1], "ploidy"=x$"ploidy"[1],  "n.chr"=length(unique(x$"ref.chr")), "n.ctg"=length(unique(x$"query")), "tot.ctg.bp"=sum(y$"ploidy_bp"), stringsAsFactors=FALSE) ; return(z) } ) 
AS.pl_stat = do.call(rbind,AS.pl_stat)


RR.stat = lapply( RR.df, function(x) { pl_col= grep("n",names(x),value=T) ; y= colSums(x[,..pl_col]) ; return(y) }     )
max_col = max(sapply(RR.stat,length))
RR.stat = lapply( RR.stat, function(x) { new_col= paste0("n",length(x):max_col) ; if( length(new_col)>0 ) for( nn in new_col ) x[[nn]]=0; return(x) }     )
RR.stat = do.call(rbind,RR.stat)

print(RR.stat)
print(colSums(RR.stat))
print(rowSums(RR.stat))



RR.score = sum(RR.stat)
RR.perc= RR.stat/RR.score


cat('\n  --> Top 10 levels \n\n')


RR.perc_t = t(RR.perc)
Lev_sums=rowSums(RR.perc_t)
Lev_cumsums=cumsum(Lev_sums)


RR.perc_t2 = cbind( RR.perc_t, Lev_sums, Lev_cumsums)
RR.perc_t3 = round(RR.perc_t2,3)
print( head(RR.perc_t3, 10))

Lev_cumsums90 = Lev_cumsums[ Lev_cumsums >= 0.90 ]

name_lev90        = names(Lev_cumsums90[1])
name_lev90_plus_1 = names(Lev_cumsums90[2])

Lev_cumsums90_95 = Lev_cumsums90 [ Lev_cumsums90 < 0.95 ]

diff90_95 = max(Lev_cumsums90_95) - min(Lev_cumsums90_95)

RR.stat.full = RR.stat

#  ploidy level
ploidy = as.numeric(gsub("n","",name_lev90))


AS.DF = AS.DF[ AS.DF$"ploidy" <= ploidy +2 , ]
max.ploidy = max(AS.DF$"ploidy")


 
cat('\n  --> Max  Ploidy: ',max.ploidy,' \n\n')
cat('\n  --> User Ploidy: ',user.ploidy,' \n\n')
cat('\n  --> Ploidy set to :',ploidy,', determined with 90-95% of reference covered   \n\n')





# if ( is.na(ploidy) == FALSE) max.ploidy =  ploidy +2
# if ( is.na(ploidy) == FALSE) AS.DF$"ploidy_level" = AS.DF$"ploidy"/max.ploidy

RR.stat = RR.stat.full [ , 1:(ploidy+2) ]
RR.stat

col.ploidy = data.frame( "ploidy" = 1:max.ploidy   , col = colorRampPalette(c("black","yellow","blue","red"))(max.ploidy), stringsAsFactors=FALSE)
col.chromm = data.frame( "ploidy" = 1:nrow(RR.stat), col = colorRampPalette(c("black","yellow","blue","red"))(nrow(RR.stat)), stringsAsFactors=FALSE)

 cat('\n  -->  - Plot 2a) Ploidy Levels Cumulated \n\n')

pdf( paste0("Ploidy_Levels_02a.Chromosome_ploidy_levels.statistics.pdf"), width=plot.w, height=plot.h )
 
    par(mfcol=c(1    ,2    ), mar=c(7,7,4,1) , oma=c(1,4,1,1) )
    barplot(RR.stat   , las=1, col=col.chromm$"col", legend=T, xlim=c(0,ncol(RR.stat)*1.4),args.legend=list(cex=2), cex.names=1.5)
    barplot(t(RR.stat), las=2, col=col.ploidy$"col" , legend=T, xlim=c(0,nrow(RR.stat)*1.4),args.legend=list(cex=2), cex.names=1)
dev.off()


 cat('\n  -->  - Plot 2b) Ploidy Levels Separated \n\n')

pdf( paste0("Ploidy_Levels_02b.Chromosome_ploidy_levels.statistics_v2.pdf"), width=plot.w, height=plot.h )
    par(mfcol=c(2    ,1    ), mar=c(6,4,0,1) , oma=c(1,4,1,1) )
    barplot(RR.stat   , las=1, col=col.chromm$"col" , legend=T, xlim=c(0,ncol(RR.stat)*nrow(RR.stat)*1.4), beside=TRUE,args.legend=list(cex=2), cex.names=1.5)
    barplot(t(RR.stat), las=2, col=col.ploidy$"col" , legend=T, xlim=c(0,ncol(RR.stat)*nrow(RR.stat)*1.4), beside=TRUE,args.legend=list(cex=2), cex.names=1.5)
dev.off()


 cat('\n  -->  - Plot 3a) Ploidy Levels: scaffolded vs unscaffolded, with transparency \n\n')


pdf( paste0("Ploidy_Levels_03a.Chromosome_ploidy_levels.scaffold_vs_unplaced.pdf"), width=plot.w, height=plot.h )

par(mfcol=c(1    ,1    ), mar=c(2,7,4,1) , oma=c(1,4,1,1) )

      
mat = REF$"ref.len"
names(mat) = REF$"ref.chr"
names(mat) = gsub("_RagTag","",names(mat) )
mat = mat [ rev(names(mat)) ]



bar = barplot(mat, horiz=T, beside=T, las=1, col=gray(0.90) ,xlim=c(-max(mat)*0.07,max(mat)*1.15), main=" Ploidy level profile", border="white", cex.main=2.5, cex.axis=2.0,cex.names=2.0)
bbb = data.frame( "ref.chr"=names(mat), pos = bar[,1])
bb2 = data.frame( "ref.chr"=names(mat), "ref.len"=as.numeric(mat), pos = bar[,1], start=0)


          mini.col = c("ref.chr","query","q.len","ref.start.min","ref.end.max","ploidy","level","ref.len","ref.mid","ref.len.filt","pos","ploidy_level")
          pol = AS.DF[ , ..mini.col]  
          pol = unique(pol)  
          pol$"scaffold" = "-"
          pol$"col" = "-"
          pol[  which( ! pol$"query" %in% AGP$"ref.chr") , ]$"scaffold" = "placed"
          pol[  which(   pol$"query" %in% AGP$"ref.chr") , ]$"scaffold" = "unplaced"
          pol[  which( ! pol$"query" %in% AGP$"ref.chr") , ]$"col" = transparent_blue
          pol[  which(   pol$"query" %in% AGP$"ref.chr") , ]$"col" = transparent_yell 
          pol$"ref.mid" = ( pol$"ref.start.min" + pol$"ref.end.max" ) /2
          pol$"POS"  = pol$"pos" -0.5 +  pol$"ploidy_level" - 1/max.ploidy
          pol$"POS2" = pol$"pos" -0.5 +  pol$"ploidy_level" 
          
          
          chr.col = c("ref.chr","query","q.len","ref.start.min","ref.end.max","ploidy","level","ref.len","ref.mid","ref.len.filt","pos","ploidy_level")

for(n in 1:nrow(pol)) polygon( c(pol$"ref.start.min"[n], pol$"ref.end.max"[n], pol$"ref.end.max"[n],pol$"ref.start.min"[n]) ,c(pol$"POS"[n], pol$"POS"[n], pol$"POS2"[n], pol$"POS2"[n] ), lwd=1 ,col=pol$"col"[n] , border=NA)
LEGEND = data.frame( "name"=c("Scaffolded contigs","Unplaced contigs") , pch=c(22,22) , col=c("black","black"),  pt.bg=c(transparent_blue,transparent_yell), stringsAsFactors=FALSE)
legend("bottomright", LEGEND$"name", pch=LEGEND$"pch",  col=LEGEND$col , pt.bg=LEGEND$"pt.bg", cex=2.5, bg="white",pt.cex=2.2)
for (n in 1:max.ploidy)    text( rep(-max(mat)*0.005,max.ploidy), bb2$"pos"-0.5 + n/max.ploidy -0.5/max.ploidy , n , col="black",cex=0.5 , lend=1)   
dev.off()

    #    segments( bb2$"start" , bb2$"pos"-0.5, bb2$"ref.len"  , bb2$"pos"-0.5 , col="black",lwd=2 , lend=1) 
  
  # v2 

# Crea una funzione per la palette sfumata

col.ploidy = data.frame( "ploidy" = 1:max.ploidy, col = colorRampPalette(c("black","yellow","blue","red"))(max.ploidy), stringsAsFactors=FALSE)
pol$"col" = NULL
pol=merge(pol,col.ploidy,by="ploidy")            
            
 cat('\n  -->  - Plot 3b) Ploidy Levels: Each level with different color \n\n')
       
pdf( paste0("Ploidy_Levels_03b.Chromosome_ploidy_levels.one_level_one_color.pdf"), width=plot.w, height=plot.h )

par(mfcol=c(1    ,1    ), mar=c(2,7,4,1) , oma=c(1,4,1,1) )
bar = barplot(mat, horiz=T, beside=T, las=1, col=gray(0.90) ,xlim=c(-max(mat)*0.07,max(mat)*1.15), main=" Ploidy level profile", border="white", cex.main=2.5, cex.axis=2.0,cex.names=2.0)


for(n in 1:nrow(pol)) polygon( c(pol$"ref.start.min"[n], pol$"ref.end.max"[n], pol$"ref.end.max"[n],pol$"ref.start.min"[n]) ,c(pol$"POS"[n], pol$"POS"[n], pol$"POS2"[n], pol$"POS2"[n] ), lwd=1 ,col=pol$"col"[n] , border=NA)
for (n in 1:max.ploidy)    text( rep(-max(mat)*0.005,max.ploidy), bb2$"pos"-0.5 + n/max.ploidy -0.5/max.ploidy , n , col="black",cex=0.5 , lend=1)   

LEGEND = data.frame( "name"= paste("n =", col.ploidy$"ploidy") , pch=22 , col=rep("black",max.ploidy),  pt.bg=col.ploidy$"col", stringsAsFactors=FALSE)
legend("bottomright", as.character(LEGEND$"name"), pch=LEGEND$"pch",  col=LEGEND$col , pt.bg=LEGEND$"pt.bg", cex=2, bg="white",pt.cex=2.2)
dev.off()



chr_sums = rowSums(RR.stat)
plo_sums = colSums(RR.stat)

plo_frac = plo_sums/sum(plo_sums)
plo_cumf = cumsum(plo_frac) 

plo_95 = head(names(plo_cumf[plo_cumf>0.95]),1)
plo_90 = head(names(plo_cumf[plo_cumf>0.90]),1)
plo_85 = head(names(plo_cumf[plo_cumf>0.85]),1)

plo_95_num = as.numeric(gsub("n","",plo_95))
plo_90_num = as.numeric(gsub("n","",plo_90))
plo_85_num = as.numeric(gsub("n","",plo_85))

plo_95_num_even = 2*ceiling(plo_95_num/2)
plo_90_num_even = 2*ceiling(plo_90_num/2)
plo_85_num_even = 2*ceiling(plo_85_num/2)

# if ( plo_95_num_even == plo_90_num_even) PLOIDY_MAX = plo_90_num_even
# if ( plo_95_num_even  > plo_90_num_even) PLOIDY_MAX = plo_95_num_even
# if ( plo_85_num_even  > plo_85_num_even) PLOIDY_MAX = plo_85_num_even

PLOIDY_MAX = plo_90_num_even

#if ( is.na(ploidy)== FALSE ) PLOIDY_MAX = ploidy

PLOIDY_MAX.char  = paste0("n",PLOIDY_MAX)
PLOIDY_OVER.char  = paste0("n>",PLOIDY_MAX+1)

cat('\n  -->  DETERMINING PLOIDY LEVEL \n\n')
cat('\n    --> higher ploidy level :',PLOIDY_MAX ,'\n\n')





max.ploidy = PLOIDY_MAX

RR.stat = RR.stat.full[ , 1:max.ploidy ]
RR.stat

AS.DF = AS.DF[ AS.DF$"ploidy" <= ploidy +2 , ]


col.ploidy = data.frame( "ploidy" = 1:max.ploidy   , col = colorRampPalette(c("black","yellow","blue","red"))(max.ploidy), stringsAsFactors=FALSE)
col.chromm = data.frame( "ploidy" = 1:nrow(RR.stat), col = colorRampPalette(c("black","yellow","blue","red"))(nrow(RR.stat)), stringsAsFactors=FALSE)

 cat('\n  -->  - Plot 2a) Ploidy Levels Cumulated \n\n')

pdf( paste0("Ploidy_Levels_02a.Chromosome_ploidy_levels.statistics.ploidy_max_",PLOIDY_MAX,".pdf"), width=plot.w, height=plot.h )
 
    par(mfcol=c(1    ,2    ), mar=c(7,7,4,1) , oma=c(1,4,1,1) )
    barplot(RR.stat   , las=1, col=col.chromm$"col", legend=T, xlim=c(0,ncol(RR.stat)*1.4),args.legend=list(cex=2), cex.names=2)
    barplot(t(RR.stat), las=2, col=col.ploidy$"col" , legend=T, xlim=c(0,nrow(RR.stat)*1.4),args.legend=list(cex=2), cex.names=1.5)

dev.off()


 cat('\n  -->  - Plot 2b) Ploidy Levels Separated \n\n')

pdf( paste0("Ploidy_Levels_02b.Chromosome_ploidy_levels.statistics_v2.ploidy_max_",PLOIDY_MAX,".pdf"), width=plot.w, height=plot.h )
    par(mfcol=c(2    ,1    ), mar=c(6,4,0,1) , oma=c(1,4,1,1) )
    barplot(RR.stat   , las=1, col=col.chromm$"col" , legend=T, xlim=c(0,ncol(RR.stat)*nrow(RR.stat)*1.4), beside=TRUE,args.legend=list(cex=2), cex.names=2)
    barplot(t(RR.stat), las=2, col=col.ploidy$"col" , legend=T, xlim=c(0,ncol(RR.stat)*nrow(RR.stat)*1.4), beside=TRUE,args.legend=list(cex=2), cex.names=1.5)
dev.off()






pdf( paste0("Ploidy_Levels_03c.Chromosome_ploidy_levels.scaffold_vs_unplaced.ploidy_max_",PLOIDY_MAX,".pdf"), width=plot.w, height=plot.h )

par(mfcol=c(1    ,1    ), mar=c(2,7,4,1) , oma=c(1,4,1,1) )

      
mat = REF$"ref.len"
names(mat) = REF$"ref.chr"
names(mat) = gsub("_RagTag","",names(mat) )
mat = mat [ rev(names(mat)) ]



bar = barplot(mat, horiz=T, beside=T, las=1, col=gray(0.90) ,xlim=c(-max(mat)*0.07,max(mat)*1.15), main=" Ploidy level profile", border="white", cex.main=2.5, cex.axis=2.0,cex.names=2.0)
bbb = data.frame( "ref.chr"=names(mat), pos = bar[,1])
bb2 = data.frame( "ref.chr"=names(mat), "ref.len"=as.numeric(mat), pos = bar[,1], start=0)


    mini.col = c("ref.chr","query","q.len","ref.start.min","ref.end.max","ploidy","level","ref.len","ref.mid","ref.len.filt","pos","ploidy_level")
    pol = AS.DF[ , ..mini.col]  
    pol = unique(pol)  
    pol[  which( pol$"ploidy" >PLOIDY_MAX ) , ]$"ploidy" = PLOIDY_MAX+1
    pol[  which( pol$"ploidy" >PLOIDY_MAX ) , ]$"level" = PLOIDY_OVER.char
    pol$"ploidy_level" = pol$"ploidy"/max.ploidy
    pol$"scaffold" = "-"
    pol$"col" = "-"
    pol[  which( ! pol$"query" %in% AGP$"ref.chr") , ]$"scaffold" = "placed"
    pol[  which(   pol$"query" %in% AGP$"ref.chr") , ]$"scaffold" = "unplaced"
    pol[  which( ! pol$"query" %in% AGP$"ref.chr") , ]$"col" = transparent_blue
    pol[  which(   pol$"query" %in% AGP$"ref.chr") , ]$"col" = transparent_yell 
    pol$"ref.mid" = ( pol$"ref.start.min" + pol$"ref.end.max" ) /2
    pol$"POS"  = pol$"pos" -0.5 +  pol$"ploidy_level" - 1/max.ploidy
    pol$"POS2" = pol$"pos" -0.5 +  pol$"ploidy_level" 
    
    
    chr.col = c("ref.chr","query","q.len","ref.start.min","ref.end.max","ploidy","level","ref.len","ref.mid","ref.len.filt","pos","ploidy_level")

for(n in 1:nrow(pol)) polygon( c(pol$"ref.start.min"[n], pol$"ref.end.max"[n], pol$"ref.end.max"[n],pol$"ref.start.min"[n]) ,c(pol$"POS"[n], pol$"POS"[n], pol$"POS2"[n], pol$"POS2"[n] ), lwd=1 ,col=pol$"col"[n] , border=NA)
LEGEND = data.frame( "name"=c("Scaffolded contigs","Unplaced contigs") , pch=c(22,22) , col=c("black","black"),  pt.bg=c(transparent_blue,transparent_yell), stringsAsFactors=FALSE)
legend("topright", LEGEND$"name", pch=LEGEND$"pch",  col=LEGEND$col , pt.bg=LEGEND$"pt.bg", cex=2.5, bg="white",pt.cex=2.2)
for (n in 1:max.ploidy)    text( rep(-max(mat)*0.005,max.ploidy), bb2$"pos"-0.5 + n/max.ploidy -0.5/max.ploidy , n , col="black",cex=0.5 , lend=1)   
dev.off()

        segments( bb2$"start" , bb2$"pos"-0.5, bb2$"ref.len"  , bb2$"pos"-0.5 , col="black",lwd=2 , lend=1) 
  
  # v2 

# Crea una funzione per la palette sfumata

col.ploidy = data.frame( "ploidy" = 1:max.ploidy, col = c(colorRampPalette(c("black","yellow","blue","red"))(max.ploidy-1),"black"), stringsAsFactors=FALSE)

pol$"col" = NULL
pol=merge(pol,col.ploidy,by="ploidy")            
            
            
pdf( paste0("Ploidy_Levels_03d.Chromosome_ploidy_levels.one_level_one_color.ploidy_max_",PLOIDY_MAX,".pdf"), width=plot.w, height=plot.h )

par(mfcol=c(1    ,1    ), mar=c(2,7,4,1) , oma=c(1,4,1,1) )
bar = barplot(mat, horiz=T, beside=T, las=1, col=gray(0.90) ,xlim=c(-max(mat)*0.07,max(mat)*1.15), main=" Ploidy level profile", border="white", cex.main=2.5, cex.axis=2.0,cex.names=2.0)


for(n in 1:nrow(pol)) polygon( c(pol$"ref.start.min"[n], pol$"ref.end.max"[n], pol$"ref.end.max"[n],pol$"ref.start.min"[n]) ,c(pol$"POS"[n], pol$"POS"[n], pol$"POS2"[n], pol$"POS2"[n] ), lwd=1 ,col=pol$"col"[n] , border=NA)
for (n in 1:max.ploidy)    text( rep(-max(mat)*0.005,max.ploidy), bb2$"pos"-0.5 + n/max.ploidy -0.5/max.ploidy , n , col="black",cex=0.5 , lend=1)   

LEGEND = data.frame( "name"= paste("n =", col.ploidy$"ploidy") , pch=22 , col=rep("black",max.ploidy),  pt.bg=col.ploidy$"col", stringsAsFactors=FALSE)
LEGEND$"name"[ max.ploidy ] = paste0("n >",PLOIDY_MAX)

legend("topright", as.character(LEGEND$"name"), pch=LEGEND$"pch",  col=LEGEND$col , pt.bg=LEGEND$"pt.bg", cex=2, bg="white",pt.cex=2.2)
dev.off()





sub.leveles = split(AS.DF$"query",AS.DF$"level")
sub.leveles = lapply(sub.leveles, function(x) data.frame(ctg=unique(x)))

for (nn in names(sub.leveles)) write.table( sub.leveles[[nn]], file = paste0("Contigs_",nn,".list"), col.names=FALSE, row.names=F, quote=F, sep="\t")    

if( hap.fa == TRUE)
{
for (nn in names(sub.leveles)) { cat(nn,"") ; system( paste0("  seqtk subseq ", LOG, " Contigs_",nn,".list  > Contigs_",nn,".fasta") )}
}


END.time =Sys.time() 
TOT.time =difftime(END.time, START.time, units="mins")


cat( paste( "==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n\n" ) )

## FILTERS


if ( is.na(user.ploidy) == FALSE) { ploidy = user.ploidy ; max.ploidy =user.ploidy +2}
if ( is.na(user.ploidy) == TRUE ) { ploidy = 4           ; max.ploidy =ploidy +2}


cat('\n  --> User Ploidy   : ',user.ploidy,' \n\n')
cat('\n  --> New Ploidy    : ',ploidy,' \n\n')
cat('\n  --> New Max Ploidy: ',max.ploidy,' \n\n')



RR.stat = RR.stat.full[ , 1:(ploidy+2) ]
RR.stat

AS.DF = AS.DF[ AS.DF$"ploidy" <= ploidy +2 , ]


col.ploidy = data.frame( "ploidy" = 1:max.ploidy   , col = colorRampPalette(c("black","yellow","blue","red"))(max.ploidy), stringsAsFactors=FALSE)
col.chromm = data.frame( "ploidy" = 1:nrow(RR.stat), col = colorRampPalette(c("black","yellow","blue","red"))(nrow(RR.stat)), stringsAsFactors=FALSE)

 cat('\n  -->  - Plot 2a) Ploidy Levels Cumulated \n\n')

pdf( paste0("Ploidy_Levels_02a.Chromosome_ploidy_levels.statistics..ploidy_max_",ploidy,".pdf"), width=plot.w, height=plot.h )
 
    par(mfcol=c(1    ,2    ), mar=c(7,7,4,1) , oma=c(1,4,1,1) )
    barplot(RR.stat   , las=1, col=col.chromm$"col", legend=T, xlim=c(0,ncol(RR.stat)*1.4),args.legend=list(cex=2), cex.names=2)
    barplot(t(RR.stat), las=2, col=col.ploidy$"col" , legend=T, xlim=c(0,nrow(RR.stat)*1.4),args.legend=list(cex=2), cex.names=1)
dev.off()


 cat('\n  -->  - Plot 2b) Ploidy Levels Separated \n\n')

pdf( paste0("Ploidy_Levels_02b.Chromosome_ploidy_levels.statistics_v2..ploidy_max_",ploidy,".pdf"), width=plot.w, height=plot.h )
    par(mfcol=c(2    ,1    ), mar=c(6,4,0,1) , oma=c(1,4,1,1) )
    barplot(RR.stat   , las=1, col=col.chromm$"col" , legend=T, xlim=c(0,ncol(RR.stat)*nrow(RR.stat)*1.4), beside=TRUE,args.legend=list(cex=2), cex.names=2)
    barplot(t(RR.stat), las=1, col=col.ploidy$"col" , legend=T, xlim=c(0,ncol(RR.stat)*nrow(RR.stat)*1.4), beside=TRUE,args.legend=list(cex=2), cex.names=2)
dev.off()


 cat('\n  -->  - Plot 3a) Ploidy Levels: scaffolded vs unscaffolded, with transparency \n\n')


pdf( paste0("Ploidy_Levels_03a.Chromosome_ploidy_levels.scaffold_vs_unplaced..ploidy_max_",ploidy,".pdf"), width=plot.w, height=plot.h )

par(mfcol=c(1    ,1    ), mar=c(2,7,4,1) , oma=c(1,4,1,1) )

      
mat = REF$"ref.len"
names(mat) = REF$"ref.chr"
names(mat) = gsub("_RagTag","",names(mat) )
mat = mat [ rev(names(mat)) ]



bar = barplot(mat, horiz=T, beside=T, las=1, col=gray(0.90) ,xlim=c(-max(mat)*0.07,max(mat)*1.15), main=" Ploidy level profile", border="white", cex.main=2.5, cex.axis=2.0,cex.names=2.0)
bbb = data.frame( "ref.chr"=names(mat), pos = bar[,1])
bb2 = data.frame( "ref.chr"=names(mat), "ref.len"=as.numeric(mat), pos = bar[,1], start=0)


          mini.col = c("ref.chr","query","q.len","ref.start.min","ref.end.max","ploidy","level","ref.len","ref.mid","ref.len.filt","pos","ploidy_level")
          pol = AS.DF[ , ..mini.col]  
          pol = unique(pol)  
          pol$"scaffold" = "-"
          pol$"col" = "-"
          pol[  which( ! pol$"query" %in% AGP$"ref.chr") , ]$"scaffold" = "placed"
          pol[  which(   pol$"query" %in% AGP$"ref.chr") , ]$"scaffold" = "unplaced"
          pol[  which( ! pol$"query" %in% AGP$"ref.chr") , ]$"col" = transparent_blue
          pol[  which(   pol$"query" %in% AGP$"ref.chr") , ]$"col" = transparent_yell 
          pol$"ref.mid" = ( pol$"ref.start.min" + pol$"ref.end.max" ) /2
          pol$"POS"  = pol$"pos" -0.5 +  pol$"ploidy_level" - 1/max.ploidy
          pol$"POS2" = pol$"pos" -0.5 +  pol$"ploidy_level" 
          
          
          chr.col = c("ref.chr","query","q.len","ref.start.min","ref.end.max","ploidy","level","ref.len","ref.mid","ref.len.filt","pos","ploidy_level")

for(n in 1:nrow(pol)) polygon( c(pol$"ref.start.min"[n], pol$"ref.end.max"[n], pol$"ref.end.max"[n],pol$"ref.start.min"[n]) ,c(pol$"POS"[n], pol$"POS"[n], pol$"POS2"[n], pol$"POS2"[n] ), lwd=1 ,col=pol$"col"[n] , border=NA)
LEGEND = data.frame( "name"=c("Scaffolded contigs","Unplaced contigs") , pch=c(22,22) , col=c("black","black"),  pt.bg=c(transparent_blue,transparent_yell), stringsAsFactors=FALSE)
legend("bottomright", LEGEND$"name", pch=LEGEND$"pch",  col=LEGEND$col , pt.bg=LEGEND$"pt.bg", cex=2.5, bg="white",pt.cex=2.2)
for (n in 1:max.ploidy)    text( rep(-max(mat)*0.005,max.ploidy), bb2$"pos"-0.5 + n/max.ploidy -0.5/max.ploidy , n , col="black",cex=0.5 , lend=1)   
dev.off()

    #    segments( bb2$"start" , bb2$"pos"-0.5, bb2$"ref.len"  , bb2$"pos"-0.5 , col="black",lwd=2 , lend=1) 
  
  # v2 

# Crea una funzione per la palette sfumata

col.ploidy = data.frame( "ploidy" = 1:max.ploidy, col = colorRampPalette(c("black","yellow","blue","red"))(max.ploidy), stringsAsFactors=FALSE)
pol$"col" = NULL
pol=merge(pol,col.ploidy,by="ploidy")            
            
 cat('\n  -->  - Plot 3b) Ploidy Levels: Each level with different color \n\n')
       
pdf( paste0("Ploidy_Levels_03b.Chromosome_ploidy_levels.one_level_one_color.ploidy_max_",ploidy,".pdf"), width=plot.w, height=plot.h )

par(mfcol=c(1    ,1    ), mar=c(2,7,4,1) , oma=c(1,4,1,1) )
bar = barplot(mat, horiz=T, beside=T, las=1, col=gray(0.90) ,xlim=c(-max(mat)*0.07,max(mat)*1.15), main=" Ploidy level profile", border="white", cex.main=2.5, cex.axis=2.0,cex.names=2.0)


for(n in 1:nrow(pol)) polygon( c(pol$"ref.start.min"[n], pol$"ref.end.max"[n], pol$"ref.end.max"[n],pol$"ref.start.min"[n]) ,c(pol$"POS"[n], pol$"POS"[n], pol$"POS2"[n], pol$"POS2"[n] ), lwd=1 ,col=pol$"col"[n] , border=NA)
for (n in 1:max.ploidy)    text( rep(-max(mat)*0.005,max.ploidy), bb2$"pos"-0.5 + n/max.ploidy -0.5/max.ploidy , n , col="black",cex=0.5 , lend=1)   

LEGEND = data.frame( "name"= paste("n =", col.ploidy$"ploidy") , pch=22 , col=rep("black",max.ploidy),  pt.bg=col.ploidy$"col", stringsAsFactors=FALSE)
legend("bottomright", as.character(LEGEND$"name"), pch=LEGEND$"pch",  col=LEGEND$col , pt.bg=LEGEND$"pt.bg", cex=2, bg="white",pt.cex=2.2)
dev.off()



cat( paste( "==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n\n" ) )

