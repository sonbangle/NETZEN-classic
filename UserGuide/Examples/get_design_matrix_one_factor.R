
library(readr)
infile = "sample_group_table.csv"
sample_groups <- read_delim(infile, 
"\t", escape_double = FALSE, trim_ws = TRUE)

fmla <- as.formula(paste("~ 0 + group"))
print(fmla)
design <- model.matrix(fmla, data=sample_groups)
rownames(design) = as.vector(sample_groups$sample)
write.table(design, file ="design.csv", row.names = TRUE, col.names = TRUE, quote=FALSE, sep="\t")
