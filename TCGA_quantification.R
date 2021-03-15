r1<-read.csv("GENE_PEPTIDE_FILE_NUMBER_OF_PSMs_PEPTIDE.csv") #PEPTIDE PSM FILE AFTER ACTG
r2<-r1[,c(2,3,8,11,12)]
r2$id<-paste(r2$V1,r2$V4,r2$V5,r2$GeneID,r2$Sequence,sep="_") #Extract chromosomal locations
r3<-r2[!duplicated(r2$id),]
mat<-r3
r3$id<-NULL
r3$GeneID<-NULL
r3$Sequence<-NULL
write.table(r3,row.names=F,col.names = F,quote=FALSE,file=paste("Novel_peptides_chromosomal_locations.bed"),sep="\t")

#TCGA BAM file quantification

system(find /media/DSRG4new/Hari/Bamfiles -type f -name \*.bam | parallel -j 8 'java -jar dist/bamstats04.jar -B Novel_peptides_chromosomal_locations.bed {} > {.}_coverage.txt';) # Quantification of peptide loci in each bam file  
system(find /media/DSRG4new/Hari/Bamfiles -type f -name \*covergae.txt -exec cat {} + >mergedfile.txt)  # merging all the coverage files.
system(grep -vwE "(#chrom)" mergedfile.txt > merged_file_HEADER_REMOVED.txt) # remove headers
r1<-read.delim("merged_file_HEADER_REMOVED.txt",header=FALSE)
r2<-r1[,c(5,8)] #extract average expression of each peptide
r2[r2==0]<-NA #converting 0 to NA
write.table(r2,"BRCA_PROTEOGENOMICS_mergedfile_header_removed_NA_added.txt",sep="\t",col.names=F)
awk '1 {if (a[$2]) {a[$2] = a[$2]" "$3} else {a[$2] = $3}} END {for (i in a) { print i,a[i]}}' BRCA_PROTEOGENOMICS_mergedfile_header_removed_NA_added.txt >BRCA_PROTEOGENOMICS_ALL_TRANSPOSED.txt
r1<-read.delim("BRCA_PROTEOGENOMICS_ALL_TRANSPOSED.txt",header=FALSE,sep=" ")
rownames(r1)<-r1[,1]
r1<-r1[,-1]
colnames(r1)<-mat$id
r1[is.na(r1)] <- 0
r1$id<-rownames(r1)
write.csv(r1,"BRCA_PROTEOGENOMICS_all_transposed_with_na_header_added.csv")
r1$id <- substr(r1$id, 0, 16)
r3<-read.csv("brca_mapped_reads.csv") # Each Sample total mapped reads
m<-merge(r3,r1,by="id") 
m1<-m[!duplicated(m$id),]
rownames(m1)<-m1[,1]
m1<-m1[,-1]
m2<-(m1/(m1$readcount))*1000000 #CPM normalisation
m2$readcount<-NULL
write.csv(m2,"GENE_PEPTIDE_FILE_NUMBER_OF_PEPTIDES_NORMALISED_COUNT.csv")
