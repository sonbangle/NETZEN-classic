#pipe_post_processing_script = "pipe_post_processing_no_cmdline_053118.R"
pipe_post_processing_script = "pipe_post_processing_no_cmdline_011519.R"
source_path = system("echo $SOURCE", intern = TRUE)
#source_path = "/run/user/1000/gvfs/sftp:host=sftp.rc.ufl.edu/ufrc/dtran/son.le/Blender_Workout/BlenderExercises/vtk_test/source"
#source_path = "/run/user/1835174/gvfs/sftp:host=sftp.rc.ufl.edu/ufrc/dtran/son.le/Blender_Workout/BlenderExercises/vtk_test/source"
pipe_post_processing_script = paste0(source_path, "/", pipe_post_processing_script)
source(pipe_post_processing_script, chdir = TRUE)
library(data.table)


# Differential expression analysis , return the output with format for Master Regulator Score: gene, logFC, pvalue, fdr, LR

differential_expression_analysis = function(gene_expression_table="~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_log.txt",
sample_group_table="~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_sample_description.csv",
pair_comparison=c("NHA", "GSC"),
is_gene_expression_count_data=FALSE,
outdir ="../data/test_de")
{
    gene_expression = import_file(gene_expression_table)
    colnames(gene_expression)[1] = "SYMBOL"
    sample_group = import_file(sample_group_table)
    colnames(sample_group)[1] = "Sample"
    colnames(sample_group)[2] = "Group"

    gene_expression = rm_row_with_duplicated_value_in_col(gene_expression, "SYMBOL")
    group1 = pair_comparison[1]
    group2 = pair_comparison[2]
    group1_samples = as.character(unlist(sample_group[sample_group$Group == group1, "Sample"]))
    group2_samples = as.character(unlist(sample_group[sample_group$Group == group2, "Sample"]))
    group1_gene_expression = gene_expression[, c("SYMBOL", group1_samples)]
    group2_gene_expression = gene_expression[, c("SYMBOL", group2_samples)]

    if (is_gene_expression_count_data)
    {
        print("analysis RNAseq data using edgeR")
        out = differential_expression_analysis_RNAseq(gene_count_table = gene_expression_table,
        sample_group_table = sample_group_table,
        pair_comparison = pair_comparison, outdir = outdir)
        return(out)
    }else
    {
        print("t test")
        ngenes = nrow(gene_expression)
        #gene, logFC, pvalue, fdr, LR
        differential_expression = data.frame(gene = rep(NA, ngenes), logFC = rep(NA, ngenes), pvalue = rep(NA, ngenes), fdr = rep(NA, ngenes), LR = rep(NA, ngenes))
        for (i in c(1 : ngenes))
        {
            gene = as.character(gene_expression[i, "SYMBOL"])
            print(paste(i, "calculating t test for :", gene))
            gr1_expr = as.numeric(group1_gene_expression[i, group1_samples])
            gr2_expr = as.numeric(group2_gene_expression[i, group2_samples])
            logFC = mean(gr2_expr) - mean(gr1_expr)
            tryCatch(
            {test = t.test(x = gr1_expr, y = gr2_expr, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
                pvalue = test$p.value
                LR = test$statistic
            },
            error =
            {
                pvalue = 1.
                LR = 0.
            })
            differential_expression[i, "gene"] = gene
            differential_expression[i, "logFC"] = logFC
            differential_expression[i, "pvalue"] = pvalue
            differential_expression[i, "LR"] = LR
        }
        differential_expression$fdr = p.adjust(differential_expression$pvalue, method = "BH")
        if (! is.null(outdir))
        {
            dir.create(outdir)
            outfile = paste(outdir, "/differential_expression_analysis_", group1, "vs", group2, ".txt", sep = "")
            write_table(differential_expression, file = outfile)
        }

        return(list(filename = outfile, differential_expression = differential_expression))
    }
}


# Get differential expression fit data from count table, using edgeR. The pipeline includes follwing steps: input data visualization, batch effect removal, robust fitting 
# The output is RDS file
get_fit_data_from_count = function(count_table_file="consolidated_count_table_translated.csv", # RNA seq count table with column (SYMBOLS, sample1, sample2, ...)
sample_factor_info="sample_groups.csv", # Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ...
outdir ="fit_data",
use_edgeR_GLM_robust=TRUE, # this option allow good result but a little bit slow,
k=NULL, k_max=50, break_data_into_big_groups=FALSE, sample_big_group_table=NULL,
use_batch_effect_removal=TRUE,
design_table=NULL # external provided design table file. If not provided (design=Null) then the function will autogenerate design table based on sample_factor_info using simple linear model where each factor is independent, no interactions. 
# design table :rows are sample IDs, columns are factors of experiments. The table should have row names (sampleIDs), no column for sample IDs. Example of design table in Examples/design.csv
)
{

    dir.create(outdir)
    print("Importing count table")
    count_table <- fread(count_table_file, data.table = F) # read count table
    library_size <- colSums(count_table[, 2 : ncol(count_table)], na.rm = T)

    # make cpm table for visualization
    for (i in 2 : ncol(count_table))
    {
        if (library_size[i - 1] == 0)count_table[, i] <- 0 else count_table[, i] <- count_table[, i] * 10 ^ 6 / library_size[i - 1]
    }
    cpm_file = paste0(outdir, "/cpm_table.csv")
    colnames(count_table)[1] = "SYMBOL"
    count_table$SYMBOL = toupper(count_table$SYMBOL)
    write.table(count_table, file = cpm_file, row.names = F, sep = "\t", quote = F)

    sample_group_table <- fread(sample_factor_info, data.table = F) # read sample group table
    sample_info_table = data.frame(Group = apply(sample_group_table,1, function(x) paste(as.character(x[2:length(x)]),collapse="_")), sample = sample_group_table[, 1])  # Reverse the order of colums to fit the pipe post processing script
    n_group = length(unique(sample_info_table$Group))
    sample_info_table_file = paste0(outdir, "/sample_info_table.csv")
    write.table(sample_info_table, file = sample_info_table_file, row.names = F, sep = "\t", quote = F)
    # Visualize data before batch effect removal
    sampleVisualization(count_table = cpm_file,
    sample_factor_info = sample_factor_info,
    method = "edgeR",
    distance = "dist",
    outdir = paste0(outdir, "/SampleVisualisation_before_batch_effect_removal")
    )

    #Batch effect removal
    
    batch_effect_removal_outdir = paste0(outdir, "/batch_effect_removal")
 
    if (is.null(k))
    {
      print("k is Null, calculating k from count table")
        nsample = ncol(count_table) - 1
        k = ceiling(nsample / 4) - n_group - 4   # For each factor, there should be approximately 4 samples for the estimation to be correct. Therefore need to divide tottal samples to four. Total number of factors would be k + number of sample groups (n_group)
        if (k > k_max)
        {k = k_max}
        if (k <2)
        {k = 2 }
      print(paste("calculated value of k:", k))
    }

    
    if (use_batch_effect_removal & k > 0 & is.null(design_table))
    {batch_removal_method="RUVs"} else
    {
      batch_removal_method="original"
      k = 0
    }
    print(paste("143 differential_expression_analysis.R k:", k, "beginning batch effect removal", "use_batch_effect_removal:", use_batch_effect_removal))
    
    batch_effect_removal(count_table = count_table_file,
    sample_factor_info = sample_factor_info,
    outdir = batch_effect_removal_outdir,
    method = batch_removal_method, k = k)



    # Robust fit using edgeR_GLM_robust
    fit_RDS_outdir = paste0(outdir, "/edgeR_fit_result")
    if (use_edgeR_GLM_robust)
    {
        dea_method = "edgeR_GLM_robust"
    }else
    {
        dea_method = "edgeR_original"
    }
    
    if (is.null(design_table))
    {design=NULL}else
    {
      design=  read.table(design_table, header = TRUE,sep="\t")
    }
    fit_RDS(batch_effect_removal_RDS = paste0(batch_effect_removal_outdir, "/RDS_out/batch_effect_removal.RDS"),
            dea_method = dea_method, 
            outdir = fit_RDS_outdir,
            break_data_into_big_groups = break_data_into_big_groups,
            sample_big_group_table = sample_big_group_table, design = design)

    #Robust differential comparison
    fit_data = paste0(fit_RDS_outdir, "/fit_result.RDS")

    return(fit_data)
}




# multiple combinations of comparisons derived from sample groups.
#  If target_group  is not none, then only process for combinations of comparisons containing target group as second group in pair comparison. Otherwise, all possible comparisons (n^2 -n)
differential_expression_analysis_RNAseq_multipair_comparison = function(count_table_file="../data/Math_dataset/consolidated_count_table_translated.csv", # RNA seq count table with column (SYMBOLS, sample1, sample2, ...)
sample_group_table="../data/Math_dataset/sample_groups.csv", # sample description table , with columns (sample, group)
target_group="ESC",
outdir ="../data/Math_dataset/Net_Zene",
use_cluster=FALSE,
use_edgeR_GLM_robust=FALSE # this option allow good result but a little bit slow
)
{


    fit_data = get_fit_data_from_count(count_table_file = count_table_file,
    sample_group_table = sample_group_table,
    outdir = outdir,
    use_edgeR_GLM_robust = use_edgeR_GLM_robust) # this option allow good result but a little bit slow)

    sample_groups = unique(sample_group_table[, 2])
    print(paste(146, "sample groups:", sample_groups))
    dea_robust_outdir = paste0(outdir, "/robust_differenential_comparison")
    dir.create((dea_robust_outdir))
    if (! is.null(target_group))
    {
        for (source_group in sample_groups)
        {
            if (source_group != target_group)
            {
                pair = paste0(source_group, "_vs_", target_group)
                print(paste("comparing pair:", pair))
                dea(fit_data_name = fit_data, source_group = source_group, target_group = target_group, outdir = paste0(dea_robust_outdir , "/", pair), RDS_out = T)}
        }
    }else
    {
        for (target_group in sample_groups)
        {
            for (source_group in sample_groups)
            {
                if (source_group != target_group)
                {
                    pair = paste0(source_group, "_vs_", target_group)
                    print(paste("comparing pair:", pair))
                    dea(fit_data_name = fit_data, source_group = group, target_group = target_group, outdir = paste0(outdir, "/robust_differenential_comparison/", group, "_vs_", target_group), RDS_out = T)}
            }
        }
    }
}



differential_expression_analysis_RNAseq = function(gene_count_table="~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_log.txt",
sample_group_table="~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_sample_description.csv",
pair_comparison=c("NHA", "GSC"),
outdir ="../data/test_de")
{
    print("todo")
}




