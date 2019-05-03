library(tidyr)
library(readr)
library(dplyr)
library(coloc)
library(qvalue)

overlap = function(chr, start, end){
#       print(chr)
#       print(start)
#       print(end)
        toreturn = sentinel %>% filter(X2 == chr) %>% filter((X3 >= start) & (X3 <= end))
        return(nrow(toreturn))
}

args = commandArgs(trailing=TRUE)
tissue = args[1]
p2 = as.numeric(args[2])


#LD blocks from Joe Pickrell
ld = read_delim("/home/shwetar/t2d/snps/ld_hg38.bed", col_names=FALSE,delim="\t")

#List of all GWAS SNPs with Pvalue < 5e-08
sentinel = read_delim("/home/shwetar/glgc/analysis_april2/ldl/sentinelsnps.txt",delim="\t",col_names=FALSE)

block_in_sentinel = apply(ld, 1, function(x){overlap(x[1],as.numeric(x[2]),as.numeric(x[3]));})

p1 = length(which(block_in_sentinel > 0)) / length(block_in_sentinel)
print(p1)

blocks_with_sentinel = which(block_in_sentinel > 0)

#find sentinel snp in each block
chrs = c()
starts = c()
sentinelsperblock = c()


for(i in blocks_with_sentinel){
	ldblock = ld[i,]
	sentinels = sentinel %>% filter(X2 == as.character(ldblock[1])) %>% filter(X3 >= as.numeric(ldblock[2])) %>% filter(X3 <= as.numeric(ldblock[3])) %>%  arrange(X9)
	chrs = c(chrs, as.character(sentinels[1,2]))
	starts = c(starts, as.numeric(sentinels[1,3]))
	sentinelsperblock = c(sentinelsperblock, paste0(as.character(sentinels[1,2]), "_", as.character(sentinels[1,3]),"_"))
}

write.table(sentinelsperblock, file="/home/shwetar/glgc/analysis_april2/ldl/sentinelsperblock.txt",row.names=F, col.names=F,quote=F)

#p12
#for each sentinel snp, read in all tests, and compute q value

#gene coordinates
geneids = read_delim("~/t2d/snps/gencode_geneids.txt",comment="#",delim="\t",col_names=FALSE)


#From all GTEx results, extract only those SNP_gene pairs of which the SNP is found in the sentinels per block
command = paste0("zcat /project/chrbrolab/gtex/results/GTEx_Analysis_v8_eQTL/GTEx_Analysis_v8_eQTL_all_associations/",tissue,".allpairs.txt.gz | grep -f /home/shwetar/glgc/analysis_april2/ldl/sentinelsperblock.txt > /home/shwetar/glgc/analysis_april2/ldl/",tissue,"/t1.txt")
system(command)
allps = read.table(paste0("/home/shwetar/glgc/analysis_april2/ldl/",tissue,"/t1.txt"),stringsAsFactors=F)
qs = qvalue(allps[,7])
p12 = 1 - qs$pi0
print(p12)


#read in gwasdata: obsolete
#gwas = read_delim("/home/shwetar/glgc/data/hg38/forcoloc_metal_ldl.txt",delim="\t")
#gwas$SNP = gsub(":","_",gwas$SNP)

#now submit coloc analyses for each ld block
j = 1

coloc_results = tibble(chr="", start="", end="", snp = "", gene="", pval1 = "", pval2 = "", pval3="",pval4="")
for(i in blocks_with_sentinel){
	print(paste0("Starting analysis for LD block ", i))
	ldblock = ld[i,]
	sentinel = sentinelsperblock[j]
	chr = strsplit(sentinel, "_")[[1]][1]
	pos = as.numeric(strsplit(sentinel, "_")[[1]][2])

	command = paste0("bsub -M 40000 -J ", chr, ":", pos, " -o /home/shwetar/glgc/analysis_april2/ldl/",tissue,"/error",chr,"_",pos, ".o -e /home/shwetar/glgc/analysis_april2/ldl/",tissue,"/error",chr,"_",pos, ".e Rscript /home/shwetar/glgc/analysis_april2/ldl/",tissue,"/coloc_per_block_pertissue.R ", chr, " ", pos, " ", sentinel, " ", p1, " ", p2, " ", tissue, " ", as.character(ldblock[1,1]), " ", as.numeric(ldblock[1,2]), " ", as.numeric(ldblock[1,3]))
	print(command)
	system(command)
	j = j + 1
}
