#!/bin/R

# R

args <- commandArgs(TRUE)



   if("--help" %in% args) 
    {
        cat('\n')
        cat(
        " Diplace busco scorere and exact duplication percentages \n",
        " \n",
        "USAGE: Rscript Plot_Busco_dup.R --grep=Busco_folder \n",
        " \n",
        "Arguments:\n",
        "  --grep=             | grep folders containing mysummary.txt          [ default all folders having full_table.tsv inside ]   linux regular expression allowed,  ","\n",
        "  --skip=             | skip folders containing these strings          [ default absent  ]  ","\n",
        "  --main.dir=         | main directory containing busco folders        [ default current  ]  ","\n",
        "  --thr.plot=         | min. n° of dupl. genes to be writtin o barplot [ default 30      ]  ","\n",
        "  --out=              | output file suffix                             [ default Busco_plot ]         ","\n",
        "  --plot.h=           | Hight in inches of PDF output files            [ default 18      ]  ","\n",
        "  --plot.w=           | Width in inches of PDF output files            [ default 35      ]  ","\n",
        "Modules required: \n",
        " - R packages: data.table  \n",
         "\n"
        )

        cat('\n')
        q(save="no")
    }



cat('\n')
cat('\n')
cat("############################################################################\n")
cat("############################################################################\n")
cat("###                           Plot_Busco_dup.r                           ###\n")
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
    
	grep         = as.character(ARG[ ARG$X1=="--grep"        ,]$"X2")
	skip         = as.character(ARG[ ARG$X1=="--skip"        ,]$"X2")
	main.dir     = as.character(ARG[ ARG$X1=="--main.dir"    ,]$"X2")
	plot.h       = as.character(ARG[ ARG$X1=="--plot.h"      ,]$"X2")
	plot.w       = as.character(ARG[ ARG$X1=="--plot.w"      ,]$"X2")
	thr.plot     = as.character(ARG[ ARG$X1=="--thr.plot"    ,]$"X2")
	out          = as.character(ARG[ ARG$X1=="--out"      ,]$"X2")
		
	plot.h = as.numeric(plot.h)
	plot.w = as.numeric(plot.w)
	thr.plot = as.numeric(thr.plot)


#  DEFAULT VALUES

    if ( length(grep)==0      ) grep = NA
    if ( length(skip)==0      ) skip = NA
    if ( length(main.dir)==0  ) main.dir = getwd()
    if ( length(thr.plot)==0  ) thr.plot = 30
    if ( length(plot.h)==0    ) plot.h   = 18
    if ( length(plot.w)==0    ) plot.w   = 35
    if ( length(out)==0       ) out   = "Busco_plot" else  out   = paste0("Busco_plot.",out)


# TEST PARAMETERS for interactive session
if( 2 ==1 )
{
   grep   =  NA
   skip   =  NA
   plot.h = 18
   plot.w = 35
   main.dir = getwd()
   thr.plot = 30
   out   = "Busco_plot"
}


# R
library(data.table)
options(width=300)
setwd(main.dir)


files = system(" ls */run_*/full_table.tsv", inter=T)
LIST = list()
#for (i in files) { j=gsub("Busco.|.fasta/summary.txt","",i); print(j) ; LIST[[j]]=fread(i)  }
#for (i in files) { j=strsplit(i,"/")[[1]][1] ; LIST[[j]]=as.data.frame(fread(i, sep="\t", fill=TRUE))  }



# grep
if ( is.na(grep)==TRUE  ) grep="*"
grep.split = unlist(strsplit(grep,","))
if( length(grep.split)==1) grep.words = grep
if( length(grep.split)>=2) grep.words = paste0("{",grep,"}")
DIR_SELECTION = paste0(grep.words ) 
if ( is.na(grep)==FALSE ) { files = grep( DIR_SELECTION, files, value=TRUE) }

# skip
SKIP=gsub(",","|", skip)
if ( is.na(skip)==FALSE ) { files =  grep( SKIP, files, value=TRUE, invert=TRUE)}


# LOAD FILES
cat("\n   --> Loading files in R \n")

for (i in files) { j=strsplit(i,"/")[[1]][1] ; LIST[[j]]=as.data.frame(fread(   cmd=paste('grep "Complete\\|Missing\\|Duplicated\\|Fragmented\\|Status" ',i )     , sep="\t", fill=TRUE))  }

# replace "" with NA
LIST = lapply(LIST, function(x) { if(any(x=="",na.rm=T)==TRUE) { x[ x=="" ]=NA } ; return(x) } )



list_ncol = sapply(LIST,ncol)
list_names= sapply(LIST,names)

all_names= names(LIST[[ which.max(list_ncol) ]] )
les_names= names(LIST[[ which.min(list_ncol) ]] )

# files complete?
cat("\n   --> full_table.tsv with all columns? \n")

print(as.data.frame(list_ncol))

if( any(list_ncol<10) )
    {
    w_miss = which(list_ncol<10)
    cat("\n   --> Adding columns to ", names(w_miss)," \n")
    for (  ii in w_miss ) { add_col = setdiff( all_names,  names(LIST[[ ii ]] ) ) ;   for (jj in add_col ) { LIST[[ ii ]] [[jj]]=NA  } }
    
    }

# Simplify names 
rm.str=""

file.words = gsub(",|_|-|/",".", names(LIST))
file.words = lapply(file.words, function(x) unlist(strsplit(x,".", fixed=TRUE)))
fi.u.words = lapply(file.words, unique)
all.words = unique(unlist(fi.u.words))

# if one word is repeat twice in only 1 or few names (e.g. same name for reference.suffix and sample), remove it before all the others
ALL.WORDS = list()
for (i in all.words) ALL.WORDS[[i]] = sapply(file.words, function(x) { sum(x==i) })


ALL.WORDS = data.frame( do.call(cbind,ALL.WORDS))


DUP.WORDS = apply(ALL.WORDS,2,function(x) all(x)==1)
DUP.WORDS = DUP.WORDS[ DUP.WORDS==TRUE ]
DUP.NAMES = names(DUP.WORDS)

dedup_file.words = lapply(file.words, function(x) x[ ! x%in% DUP.NAMES])
dedup_file.words = sapply(dedup_file.words, paste, collapse=".")



cat("\n   --> File Names \n")


NAMES = data.frame( full = files, samples = names(LIST), short =dedup_file.words, stringsAsFactors=FALSE)
names(LIST)= NAMES$"short"

print(NAMES)



cat("\n   --> Calculating  \n")

colname = as.character ( names(LIST[[1]]) )
colname = gsub("# ","",colname)
colname = gsub(" ",".",colname)


stat = sapply(LIST, head,1 )

# L2 = lapply(LIST, function(x) x[ -grep("BUSCOs|Busco", x[,1]) ,] )
# L3 = lapply(L2, function(x) { do.call(rbind, strsplit(x,"\t")) } )
#L4 = lapply(LIST, function(x) { x=setDT(as.data.frame(x)); names(x) =colname ; return(x) } )

# add factors level
L4 = lapply(LIST, function(x) { x=setDT(as.data.frame(x)); names(x) =colname ; x$"Status" <- factor(x$"Status", levels = c("Missing", "Complete", "Duplicated", "Fragmented")); return(x) } )

L5 = lapply(L4, function(x) { split(x,x$"Busco.id") } )
L6 = lapply(L5, lapply, function(x) { y=x[1,] ; y$"Sequence" = paste(x$"Sequence",collapse=",") ; y$"Gene.Start" = paste(x$"Gene.Start",collapse=","); y$"Gene.End" = paste(x$"Gene.End",collapse=","); y$"Strand" = paste(x$"Strand",collapse=","); y$"Strand" = paste(x$"Strand",collapse=",");  y$"Score" = paste(x$"Score",collapse=",");y$"Length" = paste(x$"Length",collapse=",") ;  y$"dup"=nrow(x); return(y)} )

L7 = lapply(L6,function(x) { do.call(rbind, x) } )


L7.stat = as.data.frame(t(sapply(L7,function(x) {table(x$"Status")})))
L7.dup= lapply(L7,function(x) {table(x[ x$Status=="Duplicated", ]$"dup")})
dup_n0 = sapply(L7.dup,length)==0
if ( any(dup_n0==0) )  L7.dup[dup_n0]=NA
L7.top= lapply(L7.dup,function(x) { w = which(as.numeric(names(x))<=10) ; x=x[w] ; return(x)})
L7.n10= lapply(L7.dup,function(x) { w = which(as.numeric(names(x))>10) ; x=x[w] ; return(x)})

dup_n0 = sapply(L7.top,length)==0 ; if ( any(dup_n0==0) )  L7.top[dup_n0]=NA
dup_n0 = sapply(L7.n10,length)==0 ; if ( any(dup_n0==0) )  L7.n10[dup_n0]=NA


L7.top= lapply(L7.top,function(x) { names(x)=paste0(names(x),"x"); return(x)})



dup = data.frame(n = paste0(c(2:10,">10"),"x"), stringsAsFactors=FALSE)


for (i in names(L7.top))
{
    a = L7.top[[i]] 
    b = L7.n10[[i]] 
    
    if(length(b)==0) b=0
    
    a = data.frame( "n" = names(a), genes =     as.numeric(a) , stringsAsFactors=FALSE)
    b = data.frame( "n" = ">10x"  , genes = sum(as.numeric(b)), stringsAsFactors=FALSE)
    
    
    ab=rbind(a,b)
    names(ab) = c("n", i)
    
    dup = merge(dup, ab, by="n", all.x=TRUE)
}

rownames(dup)= dup$"n"
dup = dup[ paste0(c(2:10,">10"),"x"), ]

dup$"n" = NULL
dup = t(dup)
dup [ which(is.na(dup)) ] = 0


#as.data.frame( L7.stat$"Duplicated" == rowSums(dup) )



#rownames(L7.stat) == rownames(dup)
samples = rownames(L7.stat) 

cat("\n   --> Sum of Busco entries [ it must be a single value ]: ", unique(rowSums(L7.stat)) ,"   \n")


tot.busco =  as.numeric( rowSums(L7.stat)[1])


cat("\n   --> Duplicated frequencies Matching with Duplicated values: ",all(rowSums(dup) == L7.stat$Duplicated),"   \n")




L7_dup = cbind(L7.stat,dup)
L7_per = round(100*L7_dup/tot.busco,1)

L7_dup = cbind(samples,L7_dup)
L7_per = cbind(samples,L7_per)

rownames(L7_dup) = NULL
rownames(L7_per) = NULL


cat("\n   --> Raw Results:  \n")
print(L7_dup)
cat("\n   --> Results in %:  \n")
print(L7_per)



cat("\n   --> Plotting  \n")


M2 = L7.stat


M2= t(as.matrix(M2))
# rownames(M2) =c("Complete","Duplicated","Fragmented","I","Missing")
# 
M2 = M2[ c("Complete","Duplicated","Fragmented","Missing") , ]

my_colors = c("#56B4E9", "#3492C7", "#F0E442", "#F04442")
my_colors2 = c("#56B4E9", rep("#3492C7", 10) , "#F0E442", "#F04442")
my_colors3 = c(colorRampPalette(c("#56B4E9", "black"))(11), "#F0E442", "#F04442")



#tot.busco = max(colSums(M2))


P2 = round(100*M2/tot.busco,2)

ss = data.frame( "Single" = L7.stat$"Complete", stringsAsFactors=FALSE)
fm = data.frame( "Fragmented" = L7.stat$"Fragmented", "Missing" = L7.stat$"Missing", stringsAsFactors=FALSE)

D2 = t( cbind(ss,dup,fm))



par(mar=c(3,12,3,3))

barplot(M2,horiz=T, las=1, legend=T, col=my_colors , xlim=c(0,tot.busco*1.195), main="Busco assestment results", xlab=" Busco Genes")
barplot(P2,horiz=T, las=1, legend=T, col=my_colors , xlim=c(0,100      *1.195), main="Busco assestment results", xlab=" % Busco ")
barplot(D2,horiz=T, las=1, legend=T, col=my_colors3, xlim=c(0,tot.busco *1.195), main="Busco Duplication assestment results", xlab=" Busco Genes") -> D2_box


D3 = t(D2)
D4 = D3
for( i in 1:nrow(D3)) D4[i,] = cumsum(D3[i,])


D5 = D4-(D3/2)

#threshold = 30
if( any(D3 < thr.plot )) D5 [ D3 < thr.plot ] = NA
D5[,"Single"] = D5[,"Fragmented"] = D5[,"Missing"] = NA
for( i in 1:nrow(D3)) {  d5 = D5[i,];text( d5 , rep(D2_box[i], length(d5)), names(d5) ) }


PP2 = round(100*D2/tot.busco,2)
barplot(PP2,horiz=T, las=1, legend=T, col=my_colors3 , xlim=c(0,100      *1.195), main="Busco assestment results", xlab=" % Busco ") -> PP2_box
PP5= round(100*D5/tot.busco,2)
for( i in 1:nrow(PP5)) {  pp5 = PP5[i,];text( pp5 , rep(PP2_box[i], length(pp5)), names(pp5) ) }

write.table( L7_dup, file = "Busco_duplications_stats.number_of_genes.txt" , col.names=TRUE, row.names=F, quote=F, sep="\t")    
write.table( L7_per, file = "Busco_duplications_stats.percentage_of_genes.txt" , col.names=TRUE, row.names=F, quote=F, sep="\t")    


Plot_N  = paste0(out,".number_of_genes.pdf")
Plot_P  = paste0(out,".percentage_of_genes.pdf")
Plot_ND = paste0(out,".number_of_genes.duplicated_fraction.pdf")
Plot_NP = paste0(out,".percentage_of_genes.duplicated_fraction.pdf")

pdf ( Plot_N , width=12, height= nrow(M2)+3 ) ; par(mar=c(3,9,3,3)); barplot(M2,horiz=T, las=1, legend=T, col=my_colors , xlim=c(0,tot.busco*1.195), main="Busco assestment results", xlab=" Busco Genes") ; dev.off()
pdf ( Plot_P , width=12, height= nrow(M2)+3 ) ; par(mar=c(3,9,3,3)); barplot(P2,horiz=T, las=1, legend=T, col=my_colors , xlim=c(0,100      *1.195), main="Busco assestment results", xlab=" % Busco ")    ; dev.off()
pdf ( Plot_ND , width=12, height= nrow(M2)+3 ) ; par(mar=c(3,9,3,3)); barplot(D2,horiz=T, las=1, legend=T, col=my_colors3, xlim=c(0,tot.busco *1.195), main="Busco Duplication assestment results", xlab=" Busco Genes") -> D2_box
for( i in 1:nrow(D3)) {  d5 = D5[i,];text( d5 , rep(D2_box[i], length(d5)), names(d5) ) }
dev.off()

pdf ( Plot_NP , width=12, height= nrow(M2)+3 ); par(mar=c(3,9,3,3)); barplot(PP2,horiz=T, las=1, legend=T, col=my_colors3 , xlim=c(0,100      *1.195), main="Busco assestment results", xlab=" % Busco ") -> PP2_box
for( i in 1:nrow(PP5)) {  pp5 = PP5[i,];text( pp5 , rep(PP2_box[i], length(pp5)), names(pp5) ) }
dev.off()

cat("\n   --> Saving Gene Tables  \n")


NAMES$"out" = gsub("run_embryophyta_odb10/full_table.tsv","Busco_genes_with_duplication.tsv",NAMES$"full" )
for (i in names(L7)) { j=NAMES[ NAMES$"short"==i, ]$"out" ; write.table( L7[[i]], file = j, col.names=TRUE, row.names=F, quote=F, sep="\t")    }



END.time =Sys.time() 
TOT.time =difftime(END.time, START.time, units="mins")


cat( paste( "==>> Elapsed time:  ", gsub("Time difference of"," ", round(as.numeric(TOT.time),1)  ), "minutes \n\n\n" ) )
cat( paste( "==>> Output folder: ", main.dir, " \n" ) )
cat( paste( "==>>     file path:  ",main.dir,"/",out,"*pdf \n", sep=""))
cat("\n")
cat("\n")


