
#!/bin/R

# R

args <- commandArgs(TRUE)



   if("--help" %in% args) 
    {
        cat('\n')
        cat(
        " Identify rDNA sequences from a given genome using rDNA sequences from a close species \n",
        " \n",
        "USAGE: Rscript rDNA_finder.R --ref.rDNA=plant_rDNA.fasta. --genome=my_assembly.fasta  --out=my_species_rDNA.fasta \n",
        " \n",
        "Arguments:\n",
        "  --ref.rDNA=            | rDNA fasta file of a close species ] ","\n",
        "  --genome=              | genome  ","\n",
        "  --out=                 | output file name [ default = genome.rDNA.fasta ] ","\n",
        " \n",
        "blast parameters:\n",
        "  --num_threads=         |  num_threads   [ default = 32      ]  ","\n",
        "  --perc_identity=       |  perc_identity [ default = 70      ]  ","\n",
        "  --evalue=              |  evalue        [ default = 0.000001]  ","\n",
        " \n",
        "Modules required: \n",
        " - R \n",
        " - R packages: data.table  \n",
        " - R packages: seqinr \n",
        " - blast \n",
        " - samtools \n",
        " - bedtools \n",
        "\n",
        " you can find rDNA sequences on NCBI or SILVA [ -> https://www.arb-silva.de/ ]",
        "\n\n",
        "make sure you rDNA sequences have comparable values with the following:",
        "\n     --> rDNA 18S  1800 bp ",
        "\n     --> rDNA 28S  3500 bp ",
        "\n     --> rDNA 5.8S  155 bp ",
        "\n     --> rDNA 5S    120 bp \n",
         "\n"
        )

        cat('\n')
        q(save="no")
    }



cat('\n')
cat('\n')
cat("############################################################################\n")
cat("############################################################################\n")
cat("###                             rDNA.finder.r                            ###\n")
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
    
	ref.rDNA        = as.character(ARG[ ARG$X1=="--ref.rDNA"     ,]$"X2")
	genome          = as.character(ARG[ ARG$X1=="--genome"       ,]$"X2")
	num_threads     = as.character(ARG[ ARG$X1=="--num_threads"  ,]$"X2")
	perc_identity   = as.character(ARG[ ARG$X1=="--perc_identity",]$"X2")
	evalue          = as.character(ARG[ ARG$X1=="--evalue"       ,]$"X2")
	out             = as.character(ARG[ ARG$X1=="--out"          ,]$"X2")

#  DEFAULT VALUES

    if ( length(num_threads      )==0  ) num_threads = 16  
    if ( length(perc_identity    )==0  ) perc_identity = 70 
    if ( length(evalue           )==0  ) evalue = 0.000001 
    if ( length(out              )==0  ) out = gsub(".fasta|.fa",".rDNA.fasta",genome) 
    
# TEST PARAMETERS for interactive session
if( 2 ==1 )
{
    ref.rDNA="ragtag.scaffold.fasta"
    genome="Barhee_BC4_genes_on_Ajwa_220x_18seq_nogap.gff3"
    num_threads = 16  
    perc_identity = 70 
    evalue = 0.000001 
    out = gsub(".fasta|.fa",".rDNA.fasta",genome) 

 }


cat('\n')
cat("############################################################################\n")
cat('------------------------ Summary of Parameters -----------------------------\n')
cat("############################################################################\n")
cat('\n')
    cat("   ref.rDNA      = ", ref.rDNA,"\n")
    cat("   genome        = ", genome,"\n")
    cat("   num_threads   = ", num_threads,"\n")
    cat("   perc_identity = ", perc_identity,"\n")
    cat("   evalue        = ", evalue,"\n")
    cat("   out           = ", out,"\n")
cat('\n')


##check database

genome.ndb = paste0(genome,".ndb")

if ( file.exists(genome.ndb)== FALSE ) { cat ("\n\n --> makeblastdb \n") ; system( paste( " makeblastdb  -dbtype nucl -in ", genome ))  } else { cat ("\n\n --> blastdb existing: ",genome.ndb," \n") }

    

cat ("\n\n --> blastn        ")
system( paste( "blastn -task blastn -num_threads ",num_threads," -perc_identity ",perc_identity," -evalue ",evalue," -db ",genome,"  -query ",ref.rDNA," -outfmt 6 -out Blast_ref.rDNA_vs_genome.txt" ))


library(data.table)
library(seqinr)
options(width=300)

cat ("\n\n --> blastn output processing         \n\n")

ref.fai = paste0(ref.rDNA,".fai") 
out.fai = paste0(out     ,".fai") 

if ( file.exists(ref.fai)== FALSE ) { cat (" --> generating rDNA  TElibrary         \n");  system( paste("samtools faidx", ref.rDNA ) ) } else { cat (" --> FAIDX rDNA file present \n") }
ref.len = fread( ref.fai )[,1:2]
names(ref.len) = c("gene","len")


blast = fread("Blast_ref.rDNA_vs_genome.txt")
names(blast) = c("gene","chr","id","align","mism","gaps" ,"g.start","g.end","start","end","evalue","bitscore")

blast = merge(blast,ref.len,by="gene", all.x=TRUE)
blast$"perc.rep" = round(100*blast$"align" / blast$"len",2)


bs  = split(blast,blast$"gene")
bs2 = lapply( bs, function(x) { x[ x$"bitscore" == max(x$"bitscore"),]   })

med_id =  sapply( bs, function(x) { summary(x$"id")["Median"] })
qua_id =  sapply( bs, function(x) { summary(x$"id")["3rd Qu."] })
med_pc =  sapply( bs, function(x) { summary(x$"perc.rep")["Median"] })
qua_pc =  sapply( bs, function(x) { summary(x$"perc.rep")["3rd Qu."] })
med_bs =  sapply( bs, function(x) { summary(x$"bitscore")["Median"] })
qua_bs =  sapply( bs, function(x) { summary(x$"bitscore")["3rd Qu."] })
qua_min =  min(qua_id)
qua_mi2 =  min(qua_pc)

#filter top 25%
bs3 = lapply( bs, function(x) { x[ x$"bitscore" >= summary(x$"bitscore")["3rd Qu."],]   })
bs3 = do.call(rbind, bs3)


bs3$"strand"="---"
bs3[ which(bs3$"start" > bs3$"end") , ]$"strand" = "-"
bs3[ which(bs3$"start" < bs3$"end") , ]$"strand" = "+"

bs3$"x" = paste0( "ID=", bs3$"id",",%Ref=",bs3$"perc.rep")

bs4 = setDT( data.frame( chr=bs3$"chr", start=bs3$"start", end=bs3$"end", name=bs3$"gene", x=bs3$"x", strand=bs3$"strand"))

aa = which(bs4$"start" > bs4$"end") 

bs4[ aa, ]$"start" -> x
bs4[ aa , ]$"end"  -> bs4[ aa , ]$"start" 
x                  -> bs4[ aa, ]$"end" 

write.table(bs4, file="Blast_ref.rDNA_vs_genome.best.bed",  sep="\t", row.names = FALSE, col.names=FALSE ,quote=FALSE)

cat ("\n --> extracting top 25% fasta entries        \n\n")


system( paste( " bedtools getfasta -fi ",genome," -bed Blast_ref.rDNA_vs_genome.best.bed -name -fo Blast_ref.rDNA_vs_genome.best.fasta" ))

fasta = read.fasta("Blast_ref.rDNA_vs_genome.best.fasta" )

# collapse nt in a single string
fa = sapply(fasta, paste, collapse="")
fa = sapply(fa, toupper)

# make a dataframe
fa = setDT( data.frame( name=names(fa), seq=as.character(fa), stringsAsFactors=FALSE))
# according to seqinr versiones, fasta names can contain also gebnome coordinates ,so remove them
fa$"name" = as.character(do.call(rbind,strsplit(fa$"name","::"))[,1])
fa = split(fa$"seq", fa$"name")

cat ("\n --> Names:", paste(names(fa)),"        \n\n")

fa = lapply(fa,function(x) rev(sort(table(x))))

top1 = lapply(fa, head,1)
top1 = sapply(top1, names)
names(top1) = paste0(">rDNA_",names(top1))

top1 = t( data.frame( NAME=names(top1), SEQ=as.character(top1)))
top1 = as.vector(top1)
cat ("\n\n --> saving best candidates   \n")


# modify names
top1 = gsub("_ribosomal_RNA","#rDNA",top1)
top1 = gsub("rDNA_rDNA","rDNA",top1)


write.table(top1, file=out,  sep="\t", row.names = FALSE, col.names=FALSE ,quote=FALSE)


cat ("\n--> QUALITY CHECK:        \n")
cat ("\n    make sure you rDNA sequences have comparable values with the following:")
cat ("\n     --> rDNA 18S  1800 bp ")
cat ("\n     --> rDNA 28S  3500 bp ")
cat ("\n     --> rDNA 5.8S  155 bp ")
cat ("\n     --> rDNA 5S    120 bp \n \n ")

cat ("\n ==> rDNA_finder has found these sequences :\n\n")
#system( paste( " head ",out," -c20000"    ))
system( paste( " samtools faidx ",out  ))
system( paste( " paste <(cut -f1-2 ", out.fai,") <( echo bp$'\n'bp$'\n'bp$'\n'bp )" ) )
cat ("\n")

END.time =Sys.time() 
TOT.time =difftime(END.time, START.time, units="mins")

system( paste( " rm Blast_ref.rDNA_vs_genome.txt" ))
system( paste( " rm Blast_ref.rDNA_vs_genome.best.bed" ))
system( paste( " rm Blast_ref.rDNA_vs_genome.best.fasta" ))


cat( paste( "==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n\n" ) )


