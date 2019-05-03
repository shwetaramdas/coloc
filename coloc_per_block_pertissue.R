library(coloc)
library(qvalue)
library(dplyr)
library(readr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)
print(args)
chr = args[1]
pos = args[2]
sentinel = args[3]

p1 = as.numeric(args[4])
p2 = as.numeric(args[5])
tissue = args[6]
ldblock = c(args[7], as.numeric(args[8]), as.numeric(args[9]))

print(ldblock)
print(chr)
pos = as.numeric(pos)

allps = read.table(paste0("/home/shwetar/glgc/analysis_april2/ldl/",tissue,"/t1.txt"),stringsAsFactors=F)
qs = qvalue(allps[,7])
p12 = 1 - qs$pi0
print(p12)


#gene coordinates
geneids = read_delim("/home/shwetar/t2d/snps/gencode_geneids.txt",comment="#",delim="\t",col_names=FALSE)

gwas = read_delim("/home/shwetar/glgc/data/hg38/forcoloc_metal_ldl.txt",delim="\t")
gwas$SNP = gsub(":","_",gwas$SNP)

#find genes within 1Mb
cisgenes = geneids %>% filter(X1 == chr) %>% filter(((X2 > pos - 1000000) & (X2 < pos + 1000000)) | ((X3 > pos - 1000000) & (X3 < pos + 1000000)) | ((X2 <  pos - 1000000) & (X3 > pos + 1000000))) %>% select(X4) %>% distinct()

print(cisgenes)

#remove the version number from the ensembl gene id
cisgenes$X4 = gsub("\\..*","",cisgenes$X4)

#keep only those egenes that are also cis to the sentinel snp
egenes = read.table(gzfile(paste0("/project/chrbrolab/gtex/results/GTEx_Analysis_v8_eQTL/GTEx_Analysis_v8_eQTL/",tissue,".v8.egenes.txt.gz")), sep="\t",header=T)
egenes[,1] = as.character(egenes[,1])
egenes = as.tbl(egenes)
egenes = egenes %>% filter(qval < 0.05) %>% select(gene_id) %>% distinct()
egenes$gene_id = gsub("\\..*","", egenes$gene_id)
egenes = inner_join(egenes, cisgenes, by=c('gene_id'='X4'))

print(egenes)


#subset gwas data to within ld block

coloc_results = tibble(chr="", start="", end="", snp = "", gene="", pval1 = "", pval2 = "", pval3="",pval4="")

#get sample size for eqtl
command = paste0("cat /project/chrbrolab/gtex/results/GTEx_Analysis_v8_eQTL//GTEx_Analysis_v8_eQTL_covariates/",tissue,".v8.covariates.txt | awk '{print NF-1}' | head -1")

num_samples_expn = as.numeric(system(command, intern=TRUE))


for(gene in egenes$gene_id){
        print(gene)
        #read in eqtl data for egenes
        files = list.files(paste0("/project/chrbrolab/gtex/results//GTEx_Analysis_v8_eQTL/GTEx_Analysis_v8_eQTL_all_associations_per_gene/",tissue))
        infile = paste0("/project/chrbrolab/gtex/results//GTEx_Analysis_v8_eQTL/GTEx_Analysis_v8_eQTL_all_associations_per_gene/",tissue,"/" ,grep(gene, files,value=TRUE)[2])
        eqtl = read_delim(infile,delim="\t")
        eqtl$variant = gsub("_[A-Z]*_[A-Z]*_b38","",eqtl$variant_id)

        #merge with gwas data
        m = inner_join(gwas, eqtl, by=c("SNP" = "variant"))
        r1 = which(m$maf == 0)
        r2 = which(m$EAF == 0)

        if(length(union(r1, r2)) > 0){
                m = m[-union(r1,r2),]
        }

        print(nrow(m))
        dataset1 = list()
        dataset2 = list()

        dataset1$pvalues = m$Pvalue
        dataset1$N = m$N[1]

        dataset1$MAF = pmin(m$EAF, 1 - m$EAF)
        dataset1$beta = m$Beta
        dataset1$varbeta = m$SE * m$SE
        dataset1$type = "quant"
#        dataset1$s = 74124/824006

        dataset2$pvalues = m$pval_nominal
        dataset2$N = num_samples_expn
        dataset2$MAF = m$maf
        dataset2$beta = m$slope
        dataset2$varbeta = m$slope_se * m$slope_se
        dataset2$type = "quant"

        #run coloc

        if(nrow(m) == 0){
                next;
        }
#       COLOC = coloc.abf(dataset1, dataset2, p1=p1, p2=p2,p12=p12)
        COLOC = coloc.abf(dataset1, dataset2)
        coloc_results = add_row(coloc_results, chr=chr, start=as.numeric(ldblock[2]), end = as.numeric(ldblock[3]), snp = sentinel, gene = gene, pval1 = COLOC$summary[3], pval2 = COLOC$summary[4], pval3 = COLOC$summary[5],pval4 = COLOC$summary[6])

}

write.table(coloc_results, file=paste0("/home/shwetar/glgc/analysis_april2/ldl/",tissue,"/egenes_noparameters_coloc_",chr,"_",pos,".txt"),row.names=F,col.names=F,sep="\t",quote=F)
