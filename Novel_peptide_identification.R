library(dplyr)
library(plyr)
library(tidyr)
library(data.table)
f1<-list.files(pattern="*Peptide Report.txt",recursive=TRUE)
f2<-list.files(pattern="*PSM Report.txt",recursive=TRUE)
for (k in 1:length(f1))
{
r1<-read.delim(f1[k])
print(f1[k])
r2<-subset(r1, !grepl("^NP_", r1$Protein.s.))
r3<-r2[!duplicated(r2$Sequence),]
r3<-subset(r3,Validation=="Confident")
r4<-read.delim(f2[k])
print(f2[k])
m<-merge(r4,r3,by="Sequence")
m1<-m[,c("Sequence","Protein.s..x","Modified.Sequence.x","Variable.Modifications.x","Fixed.Modifications.x","Spectrum.Title","Confidence.....x","Validation.x")]
m2<-subset(m1,Validation.x=="Confident")
list <- strsplit(as.character(f1[k]), split="/")
df <- ldply(list)
write.csv(m2,file=paste("Novel_Confident_Peptides_PSMs",df$V1,'csv',sep='.'))
}
f1<-list.files(pattern="*.Novel_Confident_Peptides_PSMs",recursive=TRUE)
for (k in 1:length(f1))
{
r1<-read.csv(f1[k])
print(f1[k])
r2<-r1[!duplicated(r1$Sequence),]
list <- strsplit(as.character(f1[k]), split="/")
df1 <- ldply(list)
r3<-r2$Sequence
r3<-as.data.frame(r3)
write.table(r3,row.names=F,col.names = F,quote=FALSE,file=paste("novel_peptides_for_ACTG",df1$V1,'txt',sep='.'))
}
f1<-list.files(pattern="*novel_peptides_for_ACTG",recursive=TRUE)
for (k in 1:length(f1))
{
data <- xmlParse("mapping_params.xml")
invisible(replaceNodes(data[["//Input/text()"]], newXMLTextNode(f1[k])))
list <- strsplit(as.character(f1[k]), split="/")
df1 <- ldply(list)
saveXML(data,file=paste(df1$V1,'xml',sep='.'))
}
f1<-list.files(pattern="*txt.xml",recursive=TRUE)
for (k in 1:length(f1))
{
i<-f1[k]
print(i)
cmd<-paste("java -Xmx8G -Xss2m -jar ACTG_mapping.jar ",i)
system(cmd)
}
f1<-list.files(pattern="*.flat$",recursive=TRUE)
f2<-list.files(pattern="*.gff$",recursive=TRUE)
for (k in 1:length(f1))
{
r1<-read.delim(f1[k])
print(f1[k])
single<-names(which(table(r1$Peptide)==1))
r2<-r1[r1$Peptide %in% single,]
write.csv(r2,file=paste(f1[k],'csv',sep='.'))
r3<-read.delim(f2[k],header=F)
list <- strsplit(as.character(r3$V9), split="=")
df <- ldply(list)
library(plyr)
df <- ldply(list)
r3$GFFID<-df$V2
r5<-merge(r2,r3,by="GFFID")
write.csv(r5,paste(f1[k],"ACTG_peptides_gff_combined",'csv',sep='.'))
}
f1<-list.files(pattern="*Novel_Confident_Peptides_PSMs",recursive=TRUE)
f2<-list.files(pattern="*ACTG_peptides_gff_combined",recursive=TRUE)
for(k in 1:length(f1))
{
r1<-read.csv(f1[k])
r2<-read.csv(paste(f2[k],sep=''))
colnames(r2)[which(names(r2) == "Peptide")] <- "Sequence"
m<-merge(r2,r1,by="Sequence")
write.csv(m,paste(f1[k],"combined_ACTG_peptides_PSMs",'csv',sep='.'))
}
system(find . -name "*.combined_ACTG_peptides_PSMs.csv" -exec cat {} \;> ALL_NOVEL_PEPTIDES_ATCG.csv)

r1<-read.csv("ALL_NOVEL_PEPTIDES_ATCG.csv")
r2<-r1[!duplicated(r1$Spectrum.Title),]
r3<-count(r2,'Sequence')
m<-merge(r3,r1,by="Sequence")
r4<-r2[!duplicated(r2$Sequence),]
r5<-count(r4,'GeneID')
m1<-merge(m,r5,by="GeneID")
write.csv(m1,"GENE_PEPTIDE_FILE_NUMBER_OF_PSMs_PEPTIDE_All.csv")

