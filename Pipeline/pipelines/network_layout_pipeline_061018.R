

library(data.table)
seed = 0
set.seed(seed)


get_full_source = function(source_file, source_folder = "$SOURCE")
{
  source_path = system(paste("echo ", source_folder), intern = TRUE)
  script = paste0(source_path, "/", source_file)
  return(script)
}
network_layout_script = "network_layout_Son_053018.R"
utilities_script = "utils.R"
master_regulator_scoring_script = "Master_Regulator_Scoring_v12.R"
differential_expression_analysis_script = "differential_expression_analysis_06102018.R"

network_layout_script = get_full_source(network_layout_script)
utilities_script = get_full_source(utilities_script)
print(utilities_script)
master_regulator_scoring_script = get_full_source(master_regulator_scoring_script)
differential_expression_analysis_script = get_full_source(differential_expression_analysis_script)

source(network_layout_script)
source(utilities_script)
source(differential_expression_analysis_script)


# This pipeline takes network file in ncol format (Gene1, Gene2, connection weight), gene expression file , sample group file and other to output: subnetwork that  change in different pair comparison
pipeline_after_MR_analysis = function(network_file = "../data/Math_dataset/mouse_embryonic_stem_cell_net.csv" ,
                                      graph_format = "ncol",
                                      node_attributes_file = "../data/Math_dataset/Net_Zene/Fibroblast2ESC/nSCORE/Fibroblast2ESC_dea/scores.csv",
                                      # result from Master Regulator Score run,
                                      differential_expression_analysis_file = "../data/Math_dataset/Net_Zene/Fibroblast2ESC/diff/Fibroblast2ESC_dea.txt"   ,
                                      # input file for Master Regulator Score analysis script, containing columns: gene, logFC, pvalue, fdr, LR
                                      gene_expression_table = "../data/Math_dataset/Net_Zene/cpm_table.csv",
                                      sample_group_table = "../data/Math_dataset/sample_groups.csv",
                                      displayed_groups = c("Fibroblast", "ESC", "iPSC", "GL261_Control"),
                                      selected_subnets = NULL,
                                      #c("1.7", "1.5.5.3"),
                                      selected_genes = NULL,
                                      #c("NKX6-2","ASCL1","MYCN","BASP1"),
                                      clustering_method = "louvain",
                                      max_cluster_size = 100 ,
                                      ndim = 2,
                                      layout_method = "fr",
                                      space_ratios = c(1, 4),
                                      node_score_field = "logFC",
                                      enrichment_score_field = c("node_smooth_score"),
                                      sampling_size = 1000,
                                      beta = 0.5,
                                      n_top_annotation = 3,
                                      importance_field = "score",
                                      n_top_genes = 3,
                                      outdir = "../data/Math_dataset/Net_Zene/Fibroblast2ESC",
                                      gene_description_file = "../data/all_gene_summary.csv",
                                      return_result = FALSE,
                                      draw_heatmap = FALSE,
                                      n_top_genes_heatmap = 50,
                                      selected_fields = c(
                                        "node_name",
                                        "score",
                                        node_score_field,
                                        "description",
                                        "summary",
                                        "subnet_name",
                                        "degree",
                                        "pvalue",
                                        "fdr",
                                        "annotation"
                                      ),
                                      is_log_gene_expression_table = FALSE,
                                      do_extract_subnet_data = FALSE,
                                      baseline_samples = NULL,
                                      use_cluster_to_draw_heatmap = FALSE,
                                      export_top_master_genes_data = FALSE,
                                      baseline_table=NULL)
{
  # n_top_genes_heatmap=50
  print(paste("node attributes file:", node_attributes_file))
  if (draw_heatmap)
  {gene_expression_table = import_file(gene_expression_table)} else
  {gene_expression_table =NULL}
  
  differential_expression = import_file(differential_expression_analysis_file)
  colnames(differential_expression)[1] = "SYMBOL"
  if (!is.null(node_attributes_file))
  {
    node_attributes = import_file(node_attributes_file)
    colnames(node_attributes)[1] = "SYMBOL"
    node_attributes_cols = colnames(node_attributes)
    cols = c()
    # copy node_attributes , excluding columns fdr, logFC as it willl merge with differential_expression having the same columns and therefore forming cols like logFC.x, logFC.y
    for (col in node_attributes_cols)
    {
      if (!is.element(col, c("fdr", "logFC")))
      {
        cols = c(cols, col)
      }
    }
    node_attributes$SYMBOL = toupper(node_attributes$SYMBOL)
    differential_expression$SYMBOL = toupper(differential_expression$SYMBOL) # need to capitalise before merging in case one of the dataframes has gene in lowercase format
    node_attributes = fast_merge(
      differential_expression,
      node_attributes[, cols],
      by.x = "SYMBOL" ,
      by.y = "SYMBOL",
      all.x = TRUE,
      all.y = TRUE
    )
    
  }
  else
  {
    node_attributes = differential_expression
  }
  # Run the network layout script

  full_pipeline(
    network_file = network_file,
    node_attributes = node_attributes,
    graph_format = graph_format,
    clustering_method = clustering_method,
    max_cluster_size = max_cluster_size,
    ndim = ndim,
    layout_method = layout_method,
    space_ratios = space_ratios,
    node_score_field = node_score_field,
    enrichment_score_field = enrichment_score_field,
    sampling_size = sampling_size,
    beta = beta,
    n_top_annotation = n_top_annotation,
    importance_field = importance_field,
    n_top_genes = n_top_genes,
    outdir = outdir,
    selected_subnets = selected_subnets,
    selected_genes = selected_genes,
    gene_description_file = gene_description_file,
    return_result = return_result,
    draw_heatmap = draw_heatmap,
    gene_expression_table = gene_expression_table,
    sample_group_table = sample_group_table,
    displayed_groups = displayed_groups,
    n_top_genes_heatmap = n_top_genes_heatmap,
    selected_fields = selected_fields,
    is_log_gene_expression_table = is_log_gene_expression_table,
    do_extract_subnet_data = do_extract_subnet_data,
    baseline_samples = baseline_samples,
    use_cluster = use_cluster_to_draw_heatmap,
    export_top_master_genes_data = export_top_master_genes_data,
    baseline_table=baseline_table
  )
  
  
}


# This pipeline takes network file in ncol format (Gene1, Gene2, connection weight), gene expression file , differential expression analysis file,  sample group file and other to output: subnetwork that  change in different pair comparison
pipeline_before_MR = function(network_file =  "../data/GBM-rnaseq.final.tr.ncol",
                              graph_format = "ncol",
                              gene_expression_table = "~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_log.txt",
                              sample_group_table = "~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_sample_description.csv",
                              pair_comparison_groups = c("GSC", "NHA"),
                              displayed_groups = c("NHA", "GSC", "NSC"),
                              outdir = "../data/test_before_MR_analysis_pipeline",
                              selected_subnets = c("1.7", "1.5.5.3"),
                              selected_genes = c("NKX6-2", "ASCL1", "MYCN", "BASP1"),
                              differential_expression_analysis_file = "../data/NHA_vs_GSC.csv"   ,
                              # input file for Master Regulator Score analysis script, containing columns: gene, logFC, pvalue, fdr, LR
                              
                              clustering_method = "louvain",
                              max_cluster_size = 100 ,
                              ndim = 2,
                              layout_method = "fr",
                              space_ratios = c(1, 6),
                              node_score_field = "logFC",
                              enrichment_score_field = c("node_smooth_score"),
                              sampling_size = 1000,
                              beta = 0.5,
                              n_top_annotation = 3,
                              importance_field = "score",
                              n_top_genes = 3,
                              
                              gene_description_file = "../data/all_gene_summary.csv",
                              return_result = FALSE,
                              draw_heatmap = FALSE,
                              n_top_genes_heatmap = 50,
                              
                              consider_positive_values_only = TRUE,
                              gene_statistics_list = "logFC,degree,betweenness",
                              fdr_to_confidence = T,
                              nround = 2,
                              steps_combined = F,
                              top_genes = 0.3,
                              source_node_inclusion = "p",
                              neighbor_aggregation_method = "a",
                              
                              is_log_gene_expression_table = FALSE,
                              do_extract_subnet_data = FALSE,
                              baseline_samples = NULL,
                              use_cluster_to_draw_heatmap = FALSE,
                              export_top_master_genes_data = FALSE,
                              get_nScore_only = FALSE,
                              baseline_table = NULL)
#' @param get_nScore_only bool only calculate nScore, not doing network layout

{
  # Run Master Regulator Scoring
  nSCORE_outdir = paste0(outdir, "/nSCORE")
  dir.create(nSCORE_outdir)
  
  node_attributes = get_nScore(
    network = network_file,
    gene_expression_statistics = differential_expression_analysis_file,
    consider_positive_values_only = consider_positive_values_only,
    master_genes = selected_genes,
    outdir = nSCORE_outdir ,
    gene_statistics_list = gene_statistics_list,
    fdr_to_confidence = fdr_to_confidence,
    nround = nround,
    steps_combined = steps_combined,
    top_genes = top_genes,
    source_node_inclusion = source_node_inclusion,
    neighbor_aggregation_method = neighbor_aggregation_method
  )

  if (get_nScore_only == FALSE)
  {
  pipeline_after_MR_analysis(
    network_file = network_file ,
    graph_format = graph_format,
    node_attributes_file = node_attributes,
    # result from Master Regulator Score run,
    differential_expression_analysis_file = differential_expression_analysis_file  ,
    # input file for Master Regulator Score analysis script, containing columns: gene, logFC, pvalue, fdr, LR
    gene_expression_table = gene_expression_table,
    sample_group_table = sample_group_table,
    displayed_groups = displayed_groups,
    selected_subnets = selected_subnets,
    selected_genes = selected_genes,
    clustering_method = clustering_method,
    max_cluster_size = max_cluster_size ,
    ndim = ndim,
    layout_method = layout_method,
    space_ratios = space_ratios,
    node_score_field = node_score_field,
    enrichment_score_field = enrichment_score_field,
    sampling_size = sampling_size,
    beta = beta,
    n_top_annotation = n_top_annotation,
    importance_field = importance_field,
    n_top_genes = n_top_genes,
    outdir = outdir,
    gene_description_file = gene_description_file,
    return_result = return_result,
    draw_heatmap = draw_heatmap,
    n_top_genes_heatmap = n_top_genes_heatmap,
    selected_fields = c(
      "node_name",
      "score",
      "logFC",
      "description",
      "summary",
      "subnet_name",
      "degree",
      "pvalue",
      "fdr",
      "annotation"
    ),
    is_log_gene_expression_table = is_log_gene_expression_table,
    do_extract_subnet_data = do_extract_subnet_data,
    baseline_samples = baseline_samples,
    use_cluster_to_draw_heatmap = use_cluster_to_draw_heatmap,
    export_top_master_genes_data = export_top_master_genes_data,
    baseline_table = baseline_table
  )
  
  }
  
}

# This pipeline takes network file in ncol format (Gene1, Gene2, connection weight), gene expression file  or gene count, sample group file and other to output: subnetwork that  change in different pair comparison
main_pipeline = function(network_file =  "../data/GBM-rnaseq.final.tr.ncol",
                         graph_format = "ncol",
                         gene_expression_table = "~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_log.txt",
                         sample_group_table = "~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_sample_description.csv",
                         pair_comparison_groups = c("GSC", "NHA"),
                         displayed_groups = c("NHA", "GSC", "NSC"),
                         outdir = "../data/test_before_MR_analysis_pipeline",
                         selected_subnets = c("1.7", "1.5.5.3"),
                         selected_genes = c("NKX6-2", "ASCL1", "MYCN", "BASP1"),
                         gene_count_table = NULL ,
                         
                         
                         clustering_method = "louvain",
                         max_cluster_size = 100 ,
                         ndim = 2,
                         layout_method = "fr",
                         space_ratios = c(1, 6),
                         node_score_field = "logFC",
                         enrichment_score_field = c("node_smooth_score"),
                         sampling_size = 1000,
                         beta = 0.5,
                         n_top_annotation = 3,
                         importance_field = "score",
                         n_top_genes = 3,
                         
                         gene_description_file = "../data/all_gene_summary.csv",
                         return_result = FALSE,
                         draw_heatmap = TRUE,
                         n_top_genes_heatmap = 50,
                         
                         
                         consider_positive_values_only = TRUE,
                         gene_statistics_list = "logFC,degree,betweenness",
                         fdr_to_confidence = T,
                         nround = 2,
                         steps_combined = F,
                         top_genes = 0.3,
                         source_node_inclusion = "p",
                         neighbor_aggregation_method = "a",
                         is_log_gene_expression_table = FAlSE,
                         do_extract_subnet_data = FALSE)

{
  # Get differential expression analysis for Master Regulator Scoring
  
  if (is.null(gene_count_table))
  {
    is_gene_expression_count_data = FALSE
  } else
  {
    is_gene_expression_count_data = TRUE
  }
  
  differential_expression = differential_expression_analysis(
    gene_count_table,
    sample_group_table,
    pair_comparison_groups,
    is_gene_expression_count_data = is_gene_expression_count_data,
    outdir = outdir
  )
  differential_expression_analysis_file = differential_expression$filename
  
  pipeline_before_MR(
    network_file =  network_file,
    graph_format = graph_format,
    gene_expression_table = gene_expression_table,
    sample_group_table = sample_group_table,
    pair_comparison_groups = pair_comparison_groups,
    displayed_groups = displayed_groups,
    outdir = outdir,
    selected_subnets = selected_subnets,
    selected_genes = selected_genes,
    differential_expression_analysis_file = differential_expression_analysis_file,
    
    clustering_method = clustering_method,
    max_cluster_size = max_cluster_size  ,
    ndim = ndim,
    layout_method = layout_method,
    space_ratios = space_ratios,
    node_score_field = node_score_field,
    enrichment_score_field = enrichment_score_field,
    sampling_size = sampling_size,
    beta = beta,
    n_top_annotation = n_top_annotation,
    importance_field = importance_field,
    n_top_genes = n_top_genes ,
    
    gene_description_file = gene_description_file,
    return_result = return_result,
    draw_heatmap = draw_heatmap,
    n_top_genes_heatmap = n_top_genes_heatmap,
    
    consider_positive_values_only = consider_positive_values_only,
    gene_statistics_list = gene_statistics_list,
    fdr_to_confidence =  fdr_to_confidence,
    nround = nround,
    steps_combined = steps_combined,
    top_genes = top_genes,
    source_node_inclusion = source_node_inclusion,
    neighbor_aggregation_method = neighbor_aggregation_method,
    
    is_log_gene_expression_table = is_log_gene_expression_table,
    do_extract_subnet_data = do_extract_subnet_data
  )
  
}








get_nScore = function(network = "../data/GBM-rnaseq.final.tr.ncol",
                      gene_expression_statistics = "../data/NHA_vs_GSC.csv",
                      consider_positive_values_only = T,
                      master_genes = c("NKX6-2", "ASCL1", "MYCN", "BASP1"),
                      outdir = "../data/nSCORE_results",
                      gene_statistics_list = "logFC,degree,betweenness",
                      fdr_to_confidence = T,
                      nround = 2,
                      steps_combined = F,
                      top_genes = 0.3,
                      source_node_inclusion = "p",
                      neighbor_aggregation_method = "a")
{
  cmdline = paste(
    "Rscript ",
    master_regulator_scoring_script,
    " --network=",
    network,
    " --cmd=run",
    " --gene_ex_stat=",
    gene_expression_statistics,
    " --outdir=",
    outdir,
    " --gene_statistics_list=",
    gene_statistics_list,
    " --fdr_to_confidence=",
    fdr_to_confidence,
    " --nround=",
    nround,
    " --steps_combined=",
    steps_combined,
    " --top_genes=",
    top_genes,
    " --source_node_inclusion=",
    source_node_inclusion,
    " --neighbor_aggregation_method=",
    neighbor_aggregation_method,
    sep = ""
  )
  if (!is.null(master_genes))
  {
    cmdline = paste0(cmdline,
                     " --master_genes=",
                     paste0(master_genes, collapse = ","))
  }
  if (consider_positive_values_only)
  {
    cmdline = paste0(cmdline, " --consider_positive_values_only")
  }
  print(paste("nSCORE command:", cmdline))
  submit_response = system(cmdline, intern = TRUE)
  print(submit_response)
  gep_filename <-
    tools::file_path_sans_ext(basename(gene_expression_statistics))
  
  result_file <- paste0(outdir, "/", gep_filename, "/scores.csv")
  return(result_file)
  
}

todo = function()
{
  print("to do")
}




# This pipeline take fit data, pair comparisons group as main inputs. It will create differential analysis file and then call pipeline_before_MR
pipeline_from_fit_data = function(fit_data = "../data/Math_dataset/Net_Zene/edgeR_fit_result/RUVs_k50edgeR_original.RDS",
                                  network_file =  "../data/Math_dataset/mouse_embryonic_stem_cell_net.csv",
                                  graph_format = "ncol",
                                  gene_expression_table = "../data/Math_dataset/Net_Zene/cpm_table.csv",
                                  sample_group_table = "../data/Math_dataset/sample_groups.csv",
                                  
                                  pair_comparison_groups = c("Fibroblast", "ESC"),
                                  displayed_groups = c("Fibroblast", "ESC", "iPSC", "GL261_Control"),
                                  outdir = "../data/Math_dataset/Net_Zene/Fibroblast2ESC",
                                  
                                  selected_subnets = NULL,
                                  selected_genes = NULL,
                                  
                                  
                                  clustering_method = "louvain",
                                  max_cluster_size = 100 ,
                                  ndim = 2,
                                  layout_method = "fr",
                                  space_ratios = c(1, 6),
                                  node_score_field = "logFC",
                                  enrichment_score_field = c("node_smooth_score"),
                                  sampling_size = 1000,
                                  beta = 0.5,
                                  n_top_annotation = 3,
                                  importance_field = "score",
                                  n_top_genes = 3,
                                  
                                  gene_description_file = "../data/all_gene_summary.csv",
                                  return_result = FALSE,
                                  draw_heatmap = FALSE,
                                  n_top_genes_heatmap = 50,
                                  
                                  consider_positive_values_only = TRUE,
                                  gene_statistics_list = "logFC,fdr,betweenness",
                                  fdr_to_confidence = F,
                                  nround = 1,
                                  steps_combined = T,
                                  top_genes = 0.2,
                                  source_node_inclusion = "p",
                                  neighbor_aggregation_method = "s",
                                  
                                  is_log_gene_expression_table = FALSE,
                                  do_extract_subnet_data = FALSE,
                                  contrast=NULL,
                                  baseline_samples = NULL,
                                  use_cluster_to_draw_heatmap = FALSE,
                                  export_top_master_genes_data = FALSE,
                                  get_nScore_only = FALSE,
                                  baseline_table = NULL
                                  )  
# "contrast for edgeR algorithm as a string of comma separated list of values for design factors. If not provided, then contrast will be made from pair_comparison groups
#    To know which factors the contrast codes fore, need to look at design table. For example --contrast=-0.5,-0.5,1. Number of elements must equal to the number of columns of design. Contrast can be made by using makeContrasts function in R limma package. For exmaple:#
#     BvsA <- makeContrasts(B-A, levels=design)

{

  dir.create(outdir)
  #save.image(file=paste0(outdir,"/pipeline_from_fit_data.Rdata"))
  source_group = pair_comparison_groups[1]
  target_group = pair_comparison_groups[2]
  dea_outdir = paste0(outdir, "/diff")

  if (! is.null(contrast) & class(contrast) == "character")
  {
    # convert string of comma separated value into a vector
    contrast = unlist(strsplit(contrast, split=","))
    contrast = trimws(contrast)
    contrast = as.numeric(contrast)
    contrast = as.vector(contrast)

  }
  
  dea(
    fit_data_name = fit_data,
    source_group = source_group,
    target_group = target_group,
    outdir = dea_outdir,
    RDS_out = T,
    contrast = contrast
  )
  if (is.null(contrast))
  {
  #differential_expression_analysis_file = paste0(dea_outdir, "/", source_group, 2, target_group, "_dea.txt")    # input file for Master Regulator Score analysis script, containing columns: gene, logFC, pvalue, fdr, LR
    differential_expression_analysis_file = paste0(dea_outdir, "/dea.txt") 
  }else
  {
    differential_expression_analysis_file = paste0(dea_outdir, "/dea.txt")  
  }
  print("524 network layout pipeline 061018")
  pipeline_before_MR (
    differential_expression_analysis_file = differential_expression_analysis_file,
    network_file = network_file,
    graph_format =  graph_format ,
    gene_expression_table = gene_expression_table,
    sample_group_table = sample_group_table,
    pair_comparison_groups = pair_comparison_groups,
    displayed_groups = displayed_groups,
    outdir = outdir,
    selected_subnets = selected_subnets,
    selected_genes = selected_genes,
    
    
    
    clustering_method = clustering_method,
    max_cluster_size = max_cluster_size  ,
    ndim = ndim,
    layout_method = layout_method,
    space_ratios = space_ratios,
    node_score_field = node_score_field,
    enrichment_score_field = enrichment_score_field,
    sampling_size = sampling_size,
    beta = beta,
    n_top_annotation = n_top_annotation,
    importance_field = importance_field ,
    n_top_genes = n_top_genes,
    
    gene_description_file = gene_description_file ,
    return_result = return_result,
    draw_heatmap = draw_heatmap,
    n_top_genes_heatmap = n_top_genes_heatmap,
    
    consider_positive_values_only = consider_positive_values_only,
    gene_statistics_list = gene_statistics_list,
    fdr_to_confidence = fdr_to_confidence,
    nround = nround,
    steps_combined = steps_combined,
    top_genes = top_genes,
    source_node_inclusion = source_node_inclusion,
    neighbor_aggregation_method = neighbor_aggregation_method,
    
    is_log_gene_expression_table = is_log_gene_expression_table,
    do_extract_subnet_data = do_extract_subnet_data,
    baseline_samples = baseline_samples,
    use_cluster_to_draw_heatmap = use_cluster_to_draw_heatmap,
    export_top_master_genes_data = export_top_master_genes_data,
    get_nScore_only = get_nScore_only,
    baseline_table=baseline_table
  )
  
  
  
  
}



# This function break count table and sample_group_table into separated smaller count tables and sample group tables (for example for each cancer type) so that can process with get_fit_data. (If intact, it will be very slow for big data)
break_count_data = function(count_table_file= "consolidated_count_table_translated_with_groups.csv",
                            sample_group_table="actual_sample_groups.csv",
                            sample_big_group_table="actual_sample_big_groups.csv",
                            outdir = "cancers")
{
  
  # Break data into each group (for example ech cancer subtype)
  sample_groups= fread(sample_group_table, data.table = FALSE)
  colnames(sample_groups)[c(1,2)]=c("Sample", "Group")
  
  count_table = fread(count_table_file, data.table = FALSE)
  colnames(count_table)[1] = "SYMBOL"
  
  sample_big_groups = fread(sample_big_group_table, data.table = FALSE)
  colnames(sample_big_groups)[c(1,2)]=c("Sample", "Big_Group")
  
  big_groups = unique(sample_big_groups$Big_Group)
  dir.create(outdir)
  for (big_group in big_groups)
  {
    big_group_dir = paste0(outdir,"/", big_group)
    dir.create(big_group_dir)
    big_group_samples = sample_big_groups[sample_big_groups$Big_Group == big_group, "Sample"]
    count_table_big_group = count_table[,c("SYMBOL", big_group_samples)]
    sample_group_big_group = sample_groups[sample_groups$Sample %in% big_group_samples, ]
    write.table(count_table_big_group, file=paste0(big_group_dir,"/count_table.csv"), quote=FALSE, row.names=FALSE, sep="\t")
    write.table(sample_group_big_group, file=paste0(big_group_dir,"/sample_groups.csv"), quote=FALSE, row.names=FALSE, sep="\t")            
  }
  write.table(big_groups, file=paste0(outdir,"/big_groups_list.csv"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep="t")
}






suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))
option_list <- list(
  make_option(
    "--cmd",
    default = "help",
    help = "run network layout pipeline . Accepted values:score_de_attributes_in_network, score_de_attributes_not_in_network, main, pipeline_from_fit_data, get_fit_data_from_count",
    type = "character"
  ),
  
  
  make_option(
    "--network",
    default = "../data/Math_dataset/mouse_embryonic_stem_cell_net.csv",
    help = "The gene regulatory network in the tab delimited format gene1 gene 2 MI from generep pipeline",
    type = "character"
  ),
  
  make_option(
    c("--node_attributes"),
    default = NULL,
    help = "node attribute file, ",
    type = "character"
  )  ,

  make_option(
    c("--graph_format"),
    default = "ncol",
    help = "format of input graph file",
    type = "character"
  ),

  make_option(
    "--clustering_method",
    default = "louvain",
    help = "clustering method, choose from louvain",
    type = "character"
  ),

  make_option(
    "--max_cluster_size",
    default = 100,
    help = "maximum cluster size to breakdown from large initial network",
    type = "integer"
  ),

  make_option(
    "--ndim",
    default = 2,
    help = "number of dimmensions for network visualization in Paraview",
    type = "integer"
  ),

  make_option(
    "--layout_method",
    default = "fr",
    help = "graph layout method",
    type = "character"
  ),

  make_option(
    "--space_ratios",
    default = "1,6",
    help = "ratio of distance between two nodes/sum(radius both nodes). The more value the ratio, the more sparse the network",
    type = "character"
  ),

  make_option(
    "--node_score_field",
    default = "logFC",
    help = "field to calculated node score",
    type = "character"
  ),

  make_option(
    "--sampling_size",
    default = 1000,
    help = "number of sample to generate random permutation for null hypothesis",
    type = "integer"
  ),

  make_option(
    "--beta",
    default = 0.5,
    help = "beta network smoothing ratio",
    type = "double"
  ),

  make_option(
    "--n_top_annotation",
    default = 3,
    help = "number of annotations for each subnet",
    type = "integer"
  ),

  make_option(
    "--importance_field",
    default = "score",
    help = "field of importance score",
    type = "character"
  ),

  make_option(
    "--n_top_genes",
    default = 3,
    help = "number of top genes to display for each subnet",
    type = "integer"
  ),

  make_option(
    "--outdir",
    default = "../data/Math_dataset/Net_Zene/Fibroblast2ESC",
    help = "output directory",
    type = "character"
  ),

  make_option(
    "--selected_genes",
    default = NULL,
    help = "comma separated list of selected genes to extract data",
    type = "character"
  ),

  make_option(
    "--selected_subnets",
    default = NULL,
    help = "comma separated list of selected subnets to extract data",
    type = "character"
  ),

  make_option(
    "--gene_description_file",
    default = "../data/all_gene_summary.csv",
    help = "gene description file name",
    type = "character"
  ),
  make_option(
    "--draw_heatmap",
    default = FALSE,
    help = "draw heatmap or not",
    action = "store_true"
  ),

  make_option(
    "--gene_expression_table",
    default = "../data/Math_dataset/Net_Zene/cpm_table.csv",
    help = "gene expression table with column SYMBOL, Sample1, Sample2, etc, rows are gene names, values are gene expression value",
    type = "character"
  ),

  make_option(
    "--sample_group_table",
    default = "../data/Math_dataset/sample_groups.csv",
    help = "sample group table with columns sample, group",
    type = "character"
  ),
  make_option(
    "--sample_big_group_table",
    default = NULL,
    help = "sample big group table with columns Sample, Big_Group. For example, Big Group is a cancer type, wheree group is cell type ",
    type = "character"
  ),
  make_option(
    "--displayed_groups",
    default = "Fibroblast,ESC,iPSC,GL261_Control",
    help = "comma separated list of sample groups to draw heatmap",
    type = "character"
  ),
  make_option(
    "--pair_comparison_groups",
    default = "Fibroblast,ESC",
    help = "list of two comma separated groups to compared each against other like EXP vs CTL"
  ),
  make_option(
    "--n_top_genes_heatmap",
    default = 50,
    help = "number of top most important genes, for drawing heatmap",
    type = "integer"
  ),

  make_option(
    "--selected_fields",
    default = "node_name, score, diff, description, summary, subnet_name, degree, pval, pfdr, annotation",
    help = "comma separated list of sample groups to draw heatmap",
    type = "character"
  ),
  make_option(
    c("--diff_expr_file"),
    default = "../data/NHA_vs_GSC.csv",
    help = "differential expression analysis file, ",
    type = "character"
  ),

  make_option("--gene_count_table", default = NULL, help = "name of gene count table "),

  make_option(
    c("-u", "--consider_positive_values_only"),
    action = "store_true",
    default = FALSE,
    help = "Only consider positive values of input statistic (like logFC) only",
    type = "logical"
  ),
  make_option(
    c("-f", "--fdr_to_confidence"),
    action = "store_true",
    default = F,
    help = "Convert fdr value to confidence score",
    type = "logical"
  ),

  make_option(
    c("-s", "--step"),
    action = "store",
    default = 2,
    help = "The number of steps from the source node to calculate neigborhood genes",
    type = "integer"
  ),
  make_option(
    c("-g", "--top_gene_statistics"),
    action = "store",
    default = "logFC",
    help = "The gene expression statistics used to substract top genes subnetwork. Choose from: logFC,pvalue,fdr,LR ",
    type = "character"
  ),
  make_option(
    c("-l", "--gene_statistics_list"),
    action = "store",
    default = "logFC,fdr,betweenness",
    #"logFC,fdr,betweenness",
    help = "The comma separated list of gene statistics used to calculate master score. Choose from: logFC,pvalue,fdr,LR,degree,betweenness,coreness,pagerank,eigen",
    type = "character"
  ),
  make_option(
    c("-r", "--nround"),
    action = "store",
    default = 1,
    help = "number of rounds in iterative calculation, the score from previous run serves as the input for the next run",
    type = "integer"
  ),
  make_option(
    c("-z", "--steps_combined"),
    action = "store_true",
    default = T,
    help = "combine scores at different steps or only use the specified step score",
    type = "logical"
  ),
  make_option(
    c("-t", "--top_genes"),
    action = "store",
    default = 0.2,
    help = "the proportion of top genes to extract subnetwork from whole network. Using top nodes subnetwork will decrease the computing time and may increase robustness but loosing sensitivity",
    type = "double"
  ),
  make_option(
    c("-d", "--source_node_inclusion"),
    action = "store",
    default = "p",
    help = "Include source node statistics in the calculation or not. Choose from n (no), s(sum) or p (product), m: use source node statistics only as inputs",
    type = "character"
  ),
  make_option(
    c("-a", "--neighbor_aggregation_method"),
    action = "store",
    default = "s",
    help = "Method to aggregate neighborhood score: sum, average or median (s,a,m)",
    type = "character"
  ),
  make_option(
    c("--fit_data"),
    action = "store",
    default = "../data/Math_dataset/Net_Zene/edgeR_fit_result/RUVs_k50edgeR_original.RDS",
    help = "differential expression fit data filename from function get_fit_data call",
    type = "character"
  ),
  make_option(
    c("--design_table"),
    action = "store",
    default = NULL,
    help = "external provided design table file. If not provided (design=Null) then the function will autogenerate design table based on sample_factor_info using simple linear model where each factor is independent, no interactions.
rows are sample IDs, columns are factors of experiments. The table should have row names (sampleIDs), no column for sample IDs. Example of design table in Examples/design.csv",
    type = "character"
  ),
  make_option(
    c("--contrast"),
    action = "store",
    default = NULL,
    help = "contrast for edgeR algorithm as a string of comma separated list of values for design factors.
    To know which factors the contrast codes fore, need to look at design table. For example --contrast=-0.5,-0.5,1. Number of elements must equal to the number of columns of design. Contrast can be made by using makeContrasts function in R limma package. For exmaple:
     BvsA <- makeContrasts(B-A, levels=design)
    ",
    type = "character"
  ),

  make_option(
    "--is_log_gene_expression_table",
    action = "store_true",
    default = FALSE,
    help = "if gene expression table contains log value or not",
    type = "logical"
  ),
  make_option(
    "--use_edgeR_GLM_robust",
    action = "store_true",
    default = FALSE,
    help = "if use GLM robust option when calculating dispersion in GLM. Enable it would give more accurate , robust result in trade off of long calculating time",
    type = "logical"
  ),
  make_option(
    "--k",
    action = "store",
    default = NULL,
    help = "the k parameter for batch effect removal algorithm, representing the number of unwanted factors. The more k the more normalized datasets but the cons is dataset is overnormalized that can decrease the sensitivity of catching master regulators in nScore. If k=0: not doing batch effect removal",
    type = "integer"
  ),
  make_option(
    "--do_extract_subnet_data",
    action = "store_true",
    default = FALSE,
    help = "extract data for each subnet or not. The root data will always be extracted",
    type = "logical"
 ),
 make_option(
   "--baseline_samples",
   default = NULL,
   help = "comma separated list of baseline samples to calculaate baseline level of expression when  drawing heatmap",
   type = "character"
 ),
 make_option(
   "--baseline_table",
   default = NULL,
   help = "name of table describing for each sample which other samples should be used as baseline to rescale for heatmap drawing.  
   Contain two colulmn: SampleID, Baseline_samples, where in the column Baseline_samples is the comma separated list of basseline samples for corresponding SampleID.
   This table takes priority over baseline_samples parameter",
   type = "character"
 ),
 make_option(
   "--use_cluster_to_draw_heatmap",
   action = "store_true",
   default = FALSE,
   help = "submit job to slurm cluster to draw heatmap",
   type = "logical"
 ),
 make_option(
   "--export_top_master_genes_data",
   action = "store_true",
   default = FALSE,
   help = "export heatmap, subnetworks of top master genes",
   type = "logical"
 ),
 make_option(
   "--get_nscore_only" ,
   action = "store_true",
   default = FALSE,
   help = "only get nScore, not doing network layout",
   type = "logical"
 )
 
)

#get command line options, if help option encountered print help and exit,
#otherwise if options not found on command line then set defaults,
print(947)
opt <- parse_args(OptionParser(option_list = option_list))
print(paste(580, opt))
if (opt$cmd != "help")
{

  print(opt)
  space_ratios = unlist(lapply(strsplit(opt$space_ratios, split = ","), as.integer))
  print(space_ratios)
  if (!is.null(opt$selected_subnets))
  {
    selected_subnets <- trim(strsplit(opt$selected_subnets, ","))
  }
  else
  {
    selected_subnets = NULL
  }
  if (!is.null(opt$selected_genes))
  {
    selected_genes = trim(strsplit(opt$selected_genes, ","))
  }
  else
  {
    selected_genes = NULL
  }
  displayed_groups = trim(strsplit(opt$displayed_groups, ","))
  #print(606)
  #print(displayed_groups)
  selected_fields  = trim(strsplit(opt$selected_fields, ","))
  pair_comparison_groups = trim(strsplit(opt$pair_comparison_groups, ","))
  #print(paste(442, pair_comparison_groups))
  switch(
    opt$cmd,
    "score_de_attributes_in_network" =
    {
      full_pipeline(
        network_file = opt$network ,
        node_attributes = opt$node_attributes,
        graph_format = opt$graph_format,
        clustering_method = opt$clustering_method,
        max_cluster_size = opt$max_cluster_size ,
        ndim = opt$ndim,
        layout_method = opt$layout_method,
        space_ratios = space_ratios,
        node_score_field = opt$node_score_field,
        enrichment_score_field = c("node_smooth_score"),
        sampling_size = opt$sampling_size,
        beta = opt$beta,
        n_top_annotation = opt$n_top_annotation,
        importance_field = opt$importance_field,
        n_top_genes = opt$n_top_genes,
        outdir = opt$outdir,
        selected_subnets = selected_subnets ,
        selected_genes = selected_genes,
        gene_description_file = opt$gene_description_file,
        draw_heatmap = opt$draw_heatmap,
        gene_expression_table = opt$gene_expression_table,
        sample_group_table = opt$sample_group_table,
        displayed_groups = displayed_groups,
        n_top_genes_heatmap = opt$n_top_genes_heatmap,
        selected_fields = selected_fields,
        is_log_gene_expression_table = opt$is_log_gene_expression_table,
        do_extract_subnet_data = opt$do_extract_subnet_data
      )


    },
    "score_de_attributes_not_in_network" =
    {
      pipeline_after_MR_analysis(
        network_file = opt$network ,
        graph_format = opt$graph_format,
        node_attributes_file = opt$node_attributes,
        differential_expression_analysis_file = opt$diff_expr_file   ,
        # input file for Master Regulator Score analysis script, containing columns: gene, logFC, pvalue, fdr, LR
        gene_expression_table = opt$gene_expression_table,
        sample_group_table = opt$sample_group_table,
        displayed_groups = displayed_groups,
        selected_subnets = selected_subnets,
        selected_genes = selected_genes,
        clustering_method = opt$clustering_method,
        max_cluster_size = opt$max_cluster_size  ,
        ndim = opt$ndim,
        layout_method = opt$layout_method,
        space_ratios = space_ratios,
        node_score_field = opt$node_score_field,
        enrichment_score_field = c("node_smooth_score"),
        sampling_size = opt$sampling_size,
        beta = opt$beta,
        n_top_annotation = opt$n_top_annotation,
        importance_field = opt$importance_field,
        n_top_genes = opt$n_top_genes,
        outdir = opt$outdir,
        gene_description_file = opt$gene_description_file,
        return_result = FALSE,
        draw_heatmap = opt$draw_heatmap,
        n_top_genes_heatmap = opt$n_top_genes_heatmap,
        selected_fields = selected_fields,
        is_log_gene_expression_table = opt$is_log_gene_expression_table,
        do_extract_subnet_data = opt$do_extract_subnet_data
      )
    },
    "main" =
    {
      main_pipeline(
        network_file =  opt$network,
        graph_format = opt$graph_format,
        gene_expression_table = opt$gene_expression_table,
        sample_group_table = opt$sample_group_table,
        pair_comparison_groups = pair_comparison_groups,
        displayed_groups = displayed_groups,
        outdir = opt$outdir,
        selected_subnets = selected_subnets,
        selected_genes = selected_genes,


        clustering_method = opt$clustering_method,
        max_cluster_size = opt$max_cluster_size,
        ndim = opt$ndim,
        layout_method = opt$layout_method,
        space_ratios = space_ratios,
        node_score_field = opt$node_score_field,
        enrichment_score_field = c("node_smooth_score"),
        sampling_size = opt$sampling_size,
        beta =  opt$beta,
        n_top_annotation = opt$n_top_annotation,
        importance_field = opt$importance_field,
        n_top_genes = opt$n_top_genes,

        gene_description_file = opt$gene_description_file,
        return_result = FALSE,
        draw_heatmap = opt$draw_heatmap,
        n_top_genes_heatmap = opt$n_top_genes_heatmap,

        gene_count_table = opt$gene_count_table,

        consider_positive_values_only = opt$consider_positive_values_only,
        gene_statistics_list = opt$gene_statistics_list,
        fdr_to_confidence = opt$fdr_to_confidence,
        nround = opt$nround,
        steps_combined = opt$steps_combined,
        top_genes = opt$top_genes,
        source_node_inclusion = opt$source_node_inclusion,
        neighbor_aggregation_method = opt$neighbor_aggregation_method,
        is_log_gene_expression_table = opt$is_log_gene_expression_table,
        do_extract_subnet_data = opt$do_extract_subnet_data
      )



    },
    "pipeline_from_fit_data" =
    {
      pipeline_from_fit_data(
        fit_data = opt$fit_data,
        network_file =   opt$network,
        graph_format = opt$graph_format,
        gene_expression_table = opt$gene_expression_table,
        sample_group_table = opt$sample_group_table,

        pair_comparison_groups = pair_comparison_groups,
        displayed_groups = displayed_groups,
        outdir = opt$outdir,

        selected_subnets = selected_subnets,
        selected_genes = selected_genes,


        clustering_method = opt$clustering_method,
        max_cluster_size = opt$max_cluster_size,
        ndim = opt$ndim,
        layout_method = opt$layout_method,
        space_ratios = space_ratios,
        node_score_field = opt$node_score_field,
        enrichment_score_field = c("node_smooth_score"),
        sampling_size = opt$sampling_size,
        beta =  opt$beta,
        n_top_annotation = opt$n_top_annotation,
        importance_field = opt$importance_field,
        n_top_genes = opt$n_top_genes,

        gene_description_file = opt$gene_description_file,
        return_result = FALSE,
        draw_heatmap = opt$draw_heatmap,
        n_top_genes_heatmap = opt$n_top_genes_heatmap,

        consider_positive_values_only = opt$consider_positive_values_only,
        gene_statistics_list = opt$gene_statistics_list,
        fdr_to_confidence = opt$fdr_to_confidence,
        nround = opt$nround,
        steps_combined = opt$steps_combined,
        top_genes = opt$top_genes,
        source_node_inclusion = opt$source_node_inclusion,
        neighbor_aggregation_method = opt$neighbor_aggregation_method,
        is_log_gene_expression_table = opt$is_log_gene_expression_table,
        do_extract_subnet_data = opt$do_extract_subnet_data,
        contrast=opt$contrast,
        baseline_samples = opt$baseline_samples,
        use_cluster_to_draw_heatmap = opt$use_cluster_to_draw_heatmap,
        export_top_master_genes_data = opt$export_top_master_genes_data,
        get_nScore_only = opt$get_nscore_only ,
        baseline_table = opt$baseline_table
      )

    },
    "get_fit_data_from_count" =
    {
      print(paste(772, "use_edgeR_GLM_robust", opt$use_edgeR_GLM_robust))
      get_fit_data_from_count(
        count_table_file = opt$gene_count_table,
        sample_factor_info  = opt$sample_group_table,
        outdir = opt$outdir,
        use_edgeR_GLM_robust = opt$use_edgeR_GLM_robust,
        k = opt$k,
        design_table = opt$design_table
      )
    },
    "break_count_data" =
    {
      break_count_data(count_table_file= opt$gene_count_table,
                                  sample_group_table=opt$sample_group_table,
                                  sample_big_group_table=opt$sample_big_group_table,
                                  outdir = opt$outdir)

    },
    
    "pipeline_before_MR" =
    {
      pipeline_before_MR(
        network_file =   opt$network,
        graph_format = opt$graph_format,
        gene_expression_table = opt$gene_expression_table,
        sample_group_table = opt$sample_group_table,
        
        pair_comparison_groups = pair_comparison_groups,
        displayed_groups = displayed_groups,
        outdir = opt$outdir,
        
        selected_subnets = selected_subnets,
        selected_genes = selected_genes,
        
        differential_expression_analysis_file = opt$diff_expr_file,
        # input file for Master Regulator Score analysis script, containing columns: gene, logFC, pvalue, fdr, LR
        
        
        clustering_method = opt$clustering_method,
        max_cluster_size = opt$max_cluster_size,
        ndim = opt$ndim,
        layout_method = opt$layout_method,
        space_ratios = space_ratios,
        node_score_field = opt$node_score_field,
        enrichment_score_field = c("node_smooth_score"),
        sampling_size = opt$sampling_size,
        beta =  opt$beta,
        n_top_annotation = opt$n_top_annotation,
        importance_field = opt$importance_field,
        n_top_genes = opt$n_top_genes,
        
        gene_description_file = opt$gene_description_file,
        return_result = FALSE,
        draw_heatmap = opt$draw_heatmap,
        n_top_genes_heatmap = opt$n_top_genes_heatmap,
        
        consider_positive_values_only = opt$consider_positive_values_only,
        gene_statistics_list = opt$gene_statistics_list,
        fdr_to_confidence = opt$fdr_to_confidence,
        nround = opt$nround,
        steps_combined = opt$steps_combined,
        top_genes = opt$top_genes,
        source_node_inclusion = opt$source_node_inclusion,
        neighbor_aggregation_method = opt$neighbor_aggregation_method,
        is_log_gene_expression_table = opt$is_log_gene_expression_table,
        do_extract_subnet_data = opt$do_extract_subnet_data,
        baseline_samples = opt$baseline_samples,
        use_cluster_to_draw_heatmap = opt$use_cluster_to_draw_heatmap,
        export_top_master_genes_data = opt$export_top_master_genes_data,
        get_nScore_only = opt$get_nscore_only ,
        baseline_table = opt$baseline_table
      )
      
    }

  )

} else
{
  print(
    " This program will convert from network into graphml network with node coordinate and hierarchial grouping"
  )
  print_help(OptionParser(option_list = option_list))
}






