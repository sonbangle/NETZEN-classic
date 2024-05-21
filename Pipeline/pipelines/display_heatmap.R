library(optparse)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(heatmaply)

source_folder = "$SOURCE"
source_path = system(paste("echo ", source_folder), intern = TRUE)

cluster_config_file= paste0(source_path, "/pipelines/config_templates/cluster_configs/cluster_config_15Gb.conf")

display_heatmap = function(data,
                           # input data frame in format columns are samples, rows are gene expression values with row names are gene names, no "SYMBOL" column. Names of columns are names of samples
                           file = "heatmap_data.csv",
                           # file name to save to hard disk
                           title = "all genes",
                           # chart type name
                           width = 9,
                           # width of plot in inches
                           height = 6,
                           # height of plot in inches
                           min = -5,
                           # minimal value of gene expression value for color range display
                           max = 5,
                           # maximal value of gene expression value for color range display
                           cluster_row = FALSE,
                           # drawing row clustering
                           # if doing row clustering and show in the heatmap
                           cluster_col = TRUE,
                           # drawing column clustering
                           # if doing column clustering and show in the heatmap
                           mat_col = NULL,
                           # A data frame of multifactor of groups to draw as a sample group annotation bar ndicating the groups that the samples in data belong to.The number of row of mat_col should equal the number of columns of data and in the same order.
                           # for pheatmap: data frame that specifies the annotations shown on top of the heatmap. Each row defines the features for a specific row. The rows in the data and in the annotation are matched using corresponding row names. Note that color schemes takes into account if variable is continuous or discrete.
                           mat_colors = NULL,
                           # A  named list of colours to color each sample group.
                           #A named list for specifying annotation_row and annotation_col track colors manually. It is possible to define the colors for only some of the features.
                           outdir = ".",
                           use_cluster=FALSE
                        )

# colour for the sampe group bar. Number of colors should be the same as the number of unique groups in mat_col
{
  # replace () with < > in title to avoid error in pandoc conversion
  title = gsub("\\(", "_", title)
  title = gsub("\\)", "_", title)
  title = gsub("\\,", "_", title)
  title = gsub("\\/", "_", title)
  
  title_no_space = gsub(" ", "_", title)
  
  # draw pheatmap and heatmaply
  if (use_cluster == TRUE)
  {
    if (!exists("job_count"))
    {
      job_count  <<- 0
    }
    if (!exists("job_log_dir"))
    {
      print("job_log_dir does not exist, createing new job_log_dir")
      job_log_dir <<- "./jobs"
      dir.create(job_log_dir)
    }else
    {
      print(paste("59 display_heatmap.R job_log_dir:", job_log_dir))
    }
    job_count <<- job_count+ 1
    print(paste("62 display_heatmap.R job count", job_count))
    inputs = list(data=data,
                  file = file,
                  title=title,
                  width = width,
                  height = height,
                  min = min,
                  max = max,
                  cluster_row = cluster_row,
                  cluster_col = cluster_col,
                  mat_col = mat_col,
                  mat_colors = mat_colors,
                  outdir = outdir)
    
    RDS_dir = paste0(job_log_dir, "/RDS")
    dir.create(RDS_dir)
    outfile = paste0(RDS_dir,"/", job_count,".RDS")
    saveRDS(inputs, file=outfile)
    print(paste("64 display_heatmap_R", job_log_dir))
    if (is.null(job_log_dir))
    {
      job_log_dir = outdir
    }
   print(paste("69 display_heatmap_R", job_log_dir))
    cmdline = paste0("$SOURCE/pipelines/submit_cluster.py --job_log_dir " , job_log_dir," --simulated False --jobname heatmap_", job_count, ".sh", " --cluster_config ", cluster_config_file,
                     " --modules R/3.5.1 Rscript $SOURCE/display_heatmap_cmdline.R --heatmap_RDS_file ",  outfile)
    print(cmdline)
    system(cmdline)
    return(0)
  }
  print(paste("heatmap data dim:", dim(data)))
  print(paste("heatmap title:", title))
  if (nrow(data) > 0)
    # not empty data
  {
    data_folder = paste0(outdir,"/data")
    dir.create(data_folder)
    write.table(
      data,
      file = paste0(data_folder,"/",file),
      quote = FALSE,
      sep = "\t",
      row.names = TRUE
    )
    data[data > max] = max
    data[data < min] = min
    data_max = max(data)
    data_min = min(data)
    data_lim = ceiling(max(abs(data_max), abs(data_min)))
    
    # manually create color range
    #myColors = c("green", "black", "red")
    myColors = c("blue", "black", "red")
    # expand the color range
    myColors = colorRampPalette(myColors)(100)
    fontsize_col = width * 50 / ncol(data)
    if (fontsize_col > 18)
    {
      fontsize_col = 18
    }
    #interactive heatmap
    if (nrow(data) < 3)
    {
      cluster_row = FALSE  # too few samples to do clustering
    }
    if (ncol(data) < 3)
    {
      cluster_col = FALSE  # too few samples to do clustering
    }
    
    dendrogram = "none"
    if (cluster_col & cluster_row)
    {
      dendrogram = "both"
    }
    if (cluster_col & !cluster_row)
    {
      dendrogram = "column"
    }
    if (!cluster_col &  cluster_row)
    {
      dendrogram = "row"
    }
    if (!cluster_col &  !cluster_row)
    {
      dendrogram = "none"
    }
    print(paste("112 display_heatmap.R cluster_col", cluster_col, "ncol:", ncol(data)))
    print(paste("113 display_heatmap.R cluster_row", cluster_row, "nrow:", nrow(data)))
    print(dendrogram)
    #main = paste("Heatmap of subnetwork", subnet_name, ";", chart_type)


    html_folder = paste0(outdir,"/interactive")
    dir.create(html_folder)
    png_folder = paste0(outdir,"/static")
    dir.create(png_folder)
    filename = substring(file,1,nchar(file)-4)
    
    pheatmap_outfile = paste(png_folder,
                             "/",
                             filename,
                             ".png",
                             sep = "")
    heatmaply_outfile = paste(html_folder,
                              "/",
                              filename,
                              ".html",
                              sep = "")

    
    if (!is.null(mat_col) & !is.null(mat_colors))
    {
      mat_col = data.frame(group = mat_col)
      rownames(mat_col) <- colnames(data)
      
      heatmaply(
        data,
        file = heatmaply_outfile,
        fontsize_row = 8,
        fontsize_col = 8,
        colors = myColors,
        k_col = NA,
        k_row = NA,
        dendrogram = dendrogram,
        col_side_colors = mat_col,
        col_side_palette = mat_colors$group,
        main = title,
        limits=c(-data_lim, data_lim)
      )
      breaksList = seq(-data_lim, data_lim, length.out = 101)
      print("mat colors:")
      print(mat_colors)
      pheatmap(
        mat               = data,
        color = myColors,
        border_color      = NA,
        show_colnames     = TRUE,
        show_rownames     = TRUE,
        annotation_col    = mat_col,
        annotation_colors = mat_colors,
        drop_levels       = TRUE,
        fontsize          = 12,
        cluster_rows = cluster_row,
        cluster_cols = cluster_col,
        fontsize_row = min(height * 50 / nrow(data), height  * 50 / 20),
        fontsize_col = fontsize_col,
        main = title ,
        filename = pheatmap_outfile ,
        width = width,
        height = height,
        treeheight_row = 5,
        breaks = breaksList
      )
      
      
      
    } else
    {
      print(4754)
      heatmaply(
        data,
        file = heatmaply_outifle,
        fontsize_row = 8,
        fontsize_col = 8,
        colors = myColors,
        k_col = NA,
        k_row = NA,
        dendrogram = dendrogram,
        main = title
      )
      print("4738")
      
      pheatmap(
        mat               = data,
        color = myColors,
        border_color      = NA,
        show_colnames     = TRUE,
        show_rownames     = TRUE,
        drop_levels       = TRUE,
        fontsize          = 12,
        cluster_rows = cluster_row,
        cluster_cols = cluster_col,
        fontsize_row = min(height * 50 / nrow(data), height  * 50 / 20),
        fontsize_col = fontsize_col,
        main = title,
        filename = pheatmap_outfile,
        width = width,
        height = height,
        treeheight_row = 5
      )
      
      
    }
    
  }else
  {
    print("empty data, not drawing heatmap")
  }
  
}


