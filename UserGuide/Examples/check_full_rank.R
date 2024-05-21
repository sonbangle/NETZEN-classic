library(readr)
library(limma)

check_full_rank = function(design_file = "design.csv")
{

design <- read.delim(design_file)


full_rank = is.fullrank(design)

if (full_rank)
{
  print("design matrix has full rank")
  
}else
{
  print("design matrix is not full rank. Redesign design matrix")
}
return (full_rank)
}

check_full_rank()