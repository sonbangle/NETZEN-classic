library(optparse)
source_folder = "$SOURCE"
source_path = system(paste("echo ", source_folder), intern = TRUE)
#source_path = "/run/user/1000/gvfs/sftp:host=sftp.rc.ufl.edu/ufrc/dtran/son.le/Blender_Workout/BlenderExercises/vtk_test/source"
#source_path = "/run/user/1835174/gvfs/sftp:host=sftp.rc.ufl.edu/ufrc/dtran/son.le/Blender_Workout/BlenderExercises/vtk_test/source/"
source(paste0(source_path, "/display_heatmap.R"))

option_list <- list(
  make_option(
    "--heatmap_RDS_file",
    default = NULL,
    help = "file name of heatmap RDS to draw heatmap",
    type = "character"
  )
)


print(paste("4723 network_layout_Son"))
opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

inputs = readRDS(opt$heatmap_RDS_file)

    
display_heatmap(data= inputs$data,
                    file = inputs$file,
                    title = inputs$title,
                    width = inputs$width,
                    height = inputs$height,
                    min = inputs$min,
                    max = inputs$max,
                    cluster_row = inputs$cluster_row,
                    cluster_col = inputs$cluster_col,
                    mat_col = inputs$mat_col,
                    mat_colors = inputs$mat_colors,
                    outdir = inputs$outdir,
                    use_cluster = FALSE)

system(paste("rm -f ",opt$heatmap_RDS_file ))
    
