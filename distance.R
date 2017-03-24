## for calculating and extracting pairwise dist and Ka/Ks values from alignments

library(seqinr)
library(reshape2)

## colours:
aric<-"#cd7058"
avag<-"#f9a65a"
rmac<-"#599ad3"
rmag<-"#79c36a"
rtar<-"#6e7f8f"

## glob all nt alignments
files<-Sys.glob("nt_aln/*nt_ali.fasta")
cat('Number of nt aln files:',length(files),'\n')

if (file.exists("distances_nt.tab")) file.remove("distances_nt.tab")

## iterate thru and print all pairwise distances
for (file in files){
  ali<-read.alignment(file,format="fasta")
  dist<-dist.alignment(ali)

	## transform to identity; replace upper tri with NA
	identity<-1-as.matrix(dist)
	identity[upper.tri(identity,diag=T)]<-NA

	## print only upper/lower triangle of matrix
	df<-subset(melt(identity))
	df<-subset(df,df$value != "NA")

	## write to file
	write.table(df,file="distances_nt.tab",append=T,sep="\t",quote=F,row.names=F,col.names=F)
  
}

distances.nt<-read.table("distances_nt.tab",head=F)
hist(distances.nt$V3,breaks=100)

## do the same for aa alignments
files<-Sys.glob("aa_aln/*fasta")
cat('Number of aa aln files:',length(files),'\n')

if (file.exists("distances_aa.tab")) file.remove("distances_aa.tab")

for (file in files){
	ali<-read.alignment(file,format="fasta")
	dist<-dist.alignment(ali,matrix=c("identity"))
	
	## transform to identity
	identity<-1-as.matrix(dist)
	
	## print only upper/lower triangle of matrix
	df<-subset(melt(identity))
	df<-subset(df,df$Var1 != df$Var2)
	
	## write to file
	write.table(df,file="distances_aa.tab",append=T,sep="\t",quote=F,row.names=F,col.names=F)
	
}
