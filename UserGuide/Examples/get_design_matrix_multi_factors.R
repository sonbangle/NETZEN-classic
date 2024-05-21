library(readr)
library(limma)
sample_groups_file = "sample_group_table.csv"

get_design_matrix_multi_factors = function(sample_groups_file, outdir =
                                             ".")
{
  sample_groups <- read_delim(sample_groups_file,
                              "\t",
                              escape_double = FALSE,
                              trim_ws = TRUE)
  sample_groups = as.data.frame(sample_groups)

  levels(sample_groups$Treatment) 
  sample_groups$Treatment = factor(sample_groups$Treatment, levels = c("SC",   "A",     "Me",    "N",     "My",    "P",     "MeP",   "MePA",  "MePN", "ANMy"))
  sample_groups$Day = factor(sample_groups$Day, levels = c("two", "three") )

  xnam <- colnames(sample_groups)[2:ncol(sample_groups)]
  #xnam = xnam[2:4]
  fmla <-
    as.formula(paste("~ 0 +  ", paste(xnam, collapse = "+")))
  print(fmla)
  
  design <- model.matrix(fmla, data = sample_groups)
  rownames(design) = as.vector(sample_groups[,1])
  print(design)
  print(outdir)
  design_file =  paste0(outdir, "/design.csv")
  full_rank = is.fullrank(design)
  
  if (full_rank)
  {
    print("design matrix has full rank")
    write.table(
      design,
      file = design_file,
      sep = "\t",
      row.names = TRUE,
      quote = FALSE
    )
    return(design)
  }else
  {
    print("design matrix is not full rank. Redesign design matrix")
  }

}

design  =get_design_matrix_multi_factors(sample_groups_file)