args = commandArgs(trailingOnly = TRUE)
mutation_file <- args[1]
signature_file <- args[2]
outpath <- args[3]
cores <- 4


# TESTING FILE
#mutation_file <- '/workspace/users/opich/bg/mutfootprints_2020/data/hartwig/signatures/mSigAct/Colon-Rectum_muts.tsv'
#signature_file <- '/workspace/users/opich/bg/mutfootprints_2020/data/hartwig/signatures/mSigAct/Colon-Rectum_SIGS_Capecitabine.tsv'
#outpath <- '/workspace/projects/reverse_calling/data/mSigAct/'

name_file <- basename(mutation_file)

outpath_name <- paste(outpath, name_file, sep ='/')
dir.create(outpath_name, showWarnings = FALSE)
cores <- as.integer("4")


# load mSigTools and mSigAct
source("./mSigAct.v0.10.R")
source("./mSigTools.v0.13.R")

library(pbmcapply)

sigs <- as.matrix(read.table(signature_file, sep = '\t', header = T, row.names = 1))
muts <- as.matrix(read.table(mutation_file, sep = '\t', header = T, row.names = 1))

names_sigs <- colnames(sigs)
target_signature <- names_sigs[length(names_sigs)]

mSigAct <- process.one.group(muts, sigs,
                             target.sig.name = target_signature,
                             path.root = outpath_name,
                             obj.fun = obj.fun.nbinom.maxlh,
                             nbinom.size=10, ## = dispersion parameter
                             mc.cores=cores) ## = number of cores to use

pval<-mSigAct$pval

apval<-p.adjust(pval,"bonferroni")
exposure<-mSigAct$exposure
df<-t(rbind(pval,apval,exposure))
df<-df[order(df[,1],decreasing = F),]
print(name_file)
name_outfile <- paste("results", name_file, "mSigAct", target_signature, "tsv", sep ='.')
write.table(df, file = paste(outpath_name, name_outfile, sep ="/"), sep ='\t')

name_outfile <- paste("results", name_file, "mSigAct", target_signature, "rds", sep ='.')
saveRDS(mSigAct, file = paste(outpath_name, name_outfile, sep ="\t"))
