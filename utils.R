
#function to trim white space of string. x: string, can be a list, or vector
library(data.table)
library(readr)
trim <- function (x) 
{
  result<-sapply(x, function(t) gsub("^\\s+|\\s+$", "", t))
  return(as.vector(result))              
}


import_file = function(input)
{
  if (class(input) =="character")
  {
    out =read_delim(input,"\t", escape_double = FALSE, trim_ws = TRUE)
  }
  else
  {
    out = input
  }
  return(out)
}

fast_merge = function(x,y, by.x, by.y, all.x=TRUE, all.y=FALSE, sort=FALSE, allow.cartesian=FALSE)
{
  
  x$row_ID =1:nrow(x)
  x = data.table(x, key = by.x)
  y = data.table(y, key= by.y)
  out = merge(x, y, by.x =by.x, by.y =by.y, all.x=all.x, all.y=all.y, sort=sort, allow.cartesian = allow.cartesian)
  out = out[order(out$row_ID), ]
  
  out = data.frame(out, check.names = FALSE)
  cols = colnames(out)
  new_cols = c()
  for (col in cols)
  {
    if (col != "row_ID")
      new_cols = c(new_cols, col)
  } 
  out = out[, new_cols]
  return(out)
  
}

# remove rows with duplicated values in one colulmn in a data frame
rm_row_with_duplicated_value_in_col = function(data, col)
{
  duplicated_values = duplicated(data[,col])
  out = data[ duplicated_values==FALSE, ]
  return(out)
}

write_table = function(x,file)
{
write.table(x, file=file , sep = "\t", quote=FALSE, row.names = FALSE )
}

consolidate_column_data <-
  function(file_pattern="*/*",
           column_index=1, outfile="consolidated_table.csv",
           header = NULL)
    # Consolidate data from a particular column from tables that share the same file patterns
    #Args: file_pattern: pattern to search for file using glob
    # column_index : index of column to extract data from
    #header: the column names of output files. If not provided, then the names of input files will be used as the header of output file
    #outfile: name of ouput file
    #Output: a table with values of each columns come from individual column in extracted table.
  {
    
    files <- Sys.glob(file_pattern)
    i = 0
    for (f in files)
    {
      print(f)
      possibleError <-
        tryCatch({
          f_data <- fread(f, data.table = F, stringsAsFactors = FALSE)
        }, error = function(e)
          e)
      if (inherits(possibleError, "error")) {
        print(possibleError)
        next
      }
      f_data <- data.frame(extracted_data = f_data[, column_index])
      i = i + 1
      colnames(f_data) <- c(f)
      if (i == 1)
        f_table <- f_data
      else
      {
        if (nrow(f_table) > nrow(f_data))
        {
          f_data  = data.frame(c(as.vector(f_data[,1]), rep(NA, nrow(f_table) - nrow(f_data))))
        }
        if (nrow(f_table) < nrow(f_data))
        {
          patch = matrix(nrow = nrow(f_data) - nrow(f_table), ncol= ncol(f_table))
          colnames(patch) = colnames(f_table)
          f_table = rbind(f_table, patch)
        }
        
        f_table <-
          data.frame(cbind(f_table, as.data.frame(f_data)))
        colnames(f_table)[ncol(f_table)] <- f
      }
    }
    if (! is.null(header))
    {
      if (length(header) == ncol(f_table))
      {
        colnames(f_table) = header
      }else
      {
        stop("Number of input files is different from the provided header. Please check the header parameter")
      }
    }
    #print(f_table)
    write.table(
      f_table,
      file = outfile,
      row.names = F,
      sep = "\t",
      quote = F
    )
    return(f_table)
  }



merge_nSCORE_rank_tables <-
  function(file_pattern="*/*", outfile="consolidated_table.csv",
           header = NULL, is_score_table=FALSE, ntop=NULL)
    # Consolidate rank data from nSCORE ranks.csv tables or scores.csv, also create top genes table similar to the consolidated table but only include the top genes
    #Args: file_pattern: pattern to search for rank file using glob
    #header: the column names of output files. If not provided, then the names of input files will be used as the header of output file
    #outfile: name of ouput file
    #Output: a table  with firt column is the name of genes and following columns are ranks values for each sample.
    # is_score_table: input files can be score table or rank table
    # ntop: number of top genes to be extracted
  {
    
    files <- Sys.glob(file_pattern)
    i = 0
    print(144)
    print(files)
    print(146)
    for (f in files)
    {
      print(f)
      possibleError <-
        tryCatch({
          f_data <- fread(f, data.table = F, stringsAsFactors = FALSE)
        }, error = function(e)
          e)
      if (inherits(possibleError, "error")) {
        print(possibleError)
        next
      }
      i = i + 1
      print(is_score_table)
      if (is_score_table)
      {
        f_data= f_data[,c("gene","score")]
        f_data[,"score"] = round(f_data[,"score"],digits = 1)
      }
      comparison = substring(f, unlist(regexpr("comparisons",f)[1] + 12) , unlist(regexpr("nSCORE",f)[1] -2 ))
      colnames(f_data)[] <- c("Gene",comparison)
      
      if (!is.null(ntop))
      {
        if (is_score_table)
        {
          sort_decreasing=TRUE
          
         
        }else
        {
          sort_decreasing=FALSE
        }
        f_data = f_data[order(f_data[,comparison],decreasing = sort_decreasing),]
        f_data = f_data[c(1:ntop), ]
      }
      if (i == 1)
      {
        f_table <- f_data}
      
      else
      {
        
        f_table <-
          fast_merge(f_table, f_data, by.x = "Gene", by.y = "Gene", all.x = TRUE, all.y = TRUE )
      }
    }
    if (! is.null(header))
    {
      if (length(header) == ncol(f_table))
      {
        colnames(f_table) = header
      }else
      {
        stop("Number of input files is different from the provided header. Please check the header parameter")
      }
    }
    print(f_table)
    if (!is_score_table)
    {
    f_table[is.na(f_table)] = nrow(f_table) 
    f_table$Total = apply(f_table[2:ncol(f_table)], 1, sum)
    f_table = f_table[order(f_table$Total, decreasing = FALSE),]
    f_table[f_table==nrow(f_table)] = NA
    }else
    {
      f_table[is.na(f_table)] = 0 
      f_table$Total= apply(f_table[2:ncol(f_table)], 1, sum)
      f_table = f_table[order(f_table$Total, decreasing = TRUE),]
      f_table[f_table==0] = NA
    }
    write.table(
      f_table,
      file = outfile,
      row.names = F,
      sep = "\t",
      quote = F,
      na=""
    )
    return(f_table)
  }
    