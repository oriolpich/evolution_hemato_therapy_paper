library(deconstructSigs)

# read the cancer type name
input_file='decons_input.tsv'
output_file='decons_input.out.tsv'

# define which type of count we should be using!
counts_method <- 'genome'

# Convert to deconstructSigs input
# this step generates the matrix suitable for the program
sigs.input <- read.table(input_file, header=T, sep="\t",  row.names = 1, 
                         stringsAsFactors = FALSE,  check.names = FALSE)
# now we run deconstructSigs for each sample in our input list
flag = 0

new_sigs <- read.csv('signature_decons.tsv',
                     sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)



for (sample in unique(rownames(sigs.input)))
{
  test = whichSignatures(tumor.ref = sigs.input,
                         signatures.ref = new_sigs,
                         sample.id = sample,
                         contexts.needed = TRUE,
                         tri.counts.method = counts_method,
                         signature.cutoff = 0.06,
                         signatures.limit = 6,
  )
  
  a = test$weights  # save the weights for each signature.
  a['SSE']  = round(sqrt(sum(test$diff * test$diff)), digits = 3)  # compute the error rate
  a['mutation_count'] = nrow(x[which(x$Sample==sample),])  # number of mutations
  # append the results of each sample in to dataframe
  if (flag == 0){total = a; flag=1}
  else{total <- rbind(total, a)}
}

# prepare CSV file
myDF <- cbind(sample_id = rownames(total), total)  # assign row names
rownames(myDF) <- NULL

# write the output to a file
write.table(myDF, file=output_file, sep="\t", col.names = TRUE, row.names=FALSE)

