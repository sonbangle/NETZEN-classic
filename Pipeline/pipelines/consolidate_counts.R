library(readr)
# This script consolidate counts from different CONSOLIDATED_COUNTS folders

consolidate_counts = function(folders=c("CONSOLIDATED_COUNTS", "Kasumi_AML_CONS_COUNTS"), outfolder="Combined_Immune_AML_cell_lines")
{

  outfolder = paste0(outfolder, "/CONSOLIDATED_COUNTS")
  dir.create(outfolder)
  
  count_types = c("count", "cpm", "rpkm", "tpm")
  collapsed_levels =c("gene_level", "transcript_level")
  translations=c("", "_translated")
  for (count_type in count_types)
  {
    dir.create(paste0(outfolder, "/", count_type))
    for(level in collapsed_levels)
    {
      dir.create(paste0(outfolder, "/", count_type, "/", level))
      if (level=="gene_level")
      {
        translations=c("", "_translated")
      }else
      {
        translations=c("")
      }
      
      for (translation in translations)
      {
        
        filename = paste0(count_type, "/",level, "/consolidated_", count_type, "_table", translation, ".csv")
        outfilename = paste0(outfolder, "/", filename)
        outfile = NULL
        for(folder in folders)
        {
          infile = paste0(folder, "/", filename)
          print(paste("processing file:", infile))
          indata <- read_delim(infile,"\t", escape_double = FALSE, trim_ws = TRUE)
           
          if (is.null(outfile))
          {
            outfile = indata
          }else
          {
            outfile = merge(outfile, indata, by = colnames(indata)[1], all.x=TRUE, all.y=TRUE, sort=FALSE)
          }
  
        }
        write.table(outfile, file=outfilename, quote = FALSE, sep="\t", row.names = FALSE)
        
      }
    }
  }
  
}

#consolidate_counts(folders=c("CONSOLIDATED_COUNTS", "Kasumi_AML_CONS_COUNTS"), outfolder="Combined_Immune_AML_cell_lines")

#This script extract predefined samples from different CONSOLIDATED_COUNTS folders and put into a new CONSOLIDATED_COUNTS folder

consolidate_selective_samples = function(sample_group_table = "/run/user/1835174/gvfs/sftp:host=hpg.rc.ufl.edu/orange/dtran/Collaborator/David/MoonShot_Conversion/Hepatocytes/sample_group_table.csv", 
                                         folders = c("/orange/dtran/Human_Immune_Cells/CONSOLIDATED_COUNTS","/orange/dtran/Collaborator/David/MoonShot_Conversion/Hepatocytes/Hepatocytes_Only/CONSOLIDATED_COUNTS"),
                                         outfolder ="/run/user/1835174/gvfs/sftp:host=hpg.rc.ufl.edu/orange/dtran/Collaborator/David/MoonShot_Conversion/Hepatocytes/" )
{
  
  #sample_group_table = "/run/user/1835174/gvfs/sftp:host=hpg.rc.ufl.edu/orange/dtran/Collaborator/David/MoonShot_Conversion/Hepatocytes/sample_group_table.csv"
  #folders = c("/run/user/1835174/gvfs/sftp:host=hpg.rc.ufl.edu/orange/dtran/Human_Immune_Cells/CONSOLIDATED_COUNTS","/run/user/1835174/gvfs/sftp:host=hpg.rc.ufl.edu/orange/dtran/Collaborator/David/MoonShot_Conversion/Hepatocytes/Hepatocytes_Only/CONSOLIDATED_COUNTS")
  #outfolder ="/run/user/1835174/gvfs/sftp:host=hpg.rc.ufl.edu/orange/dtran/Collaborator/David/MoonShot_Conversion/Hepatocytes/"
  
  sample_group_table = read.delim(sample_group_table,stringsAsFactors=FALSE)
  samples = sample_group_table$sample
  
  outfolder = paste0(outfolder, "/CONSOLIDATED_COUNTS")
  dir.create(outfolder)
  
  count_types = c("count", "cpm", "rpkm", "tpm")
  collapsed_levels =c("gene_level", "transcript_level")
  translations=c("", "_translated")
  for (count_type in count_types)
  {
    dir.create(paste0(outfolder, "/", count_type))
    for(level in collapsed_levels)
    {
      dir.create(paste0(outfolder, "/", count_type, "/", level))
      if (level=="gene_level")
      {
        translations=c("", "_translated")
      }else
      {
        translations=c("")
      }
      
      for (translation in translations)
      {
        
        filename = paste0(count_type, "/",level, "/consolidated_", count_type, "_table", translation, ".csv")
        outfilename = paste0(outfolder, "/", filename)
        outfile = NULL
        for(folder in folders)
        {
          infile = paste0(folder, "/", filename)
          print(paste("processing file:", infile))
          indata <- read_delim(infile,"\t", escape_double = FALSE, trim_ws = TRUE)
          print(paste("samples:", samples))
          selected_samples = intersect(colnames(indata)[2:ncol(indata)], samples)
          print(paste("selected sample:", selected_samples))
          indata_selected = indata[, selected_samples]
          indata = indata[,1]
          indata = cbind(indata, indata_selected)
          
          if (is.null(outfile))
          {
            outfile = indata
          }else
          {
            outfile = merge(outfile, indata, by = colnames(indata)[1], all.x=FALSE, all.y=FALSE, sort=FALSE)
          }
          
        }
        write.table(outfile, file=outfilename, quote = FALSE, sep="\t", row.names = FALSE)
        
      }
    }
  }
  
  
  
}

