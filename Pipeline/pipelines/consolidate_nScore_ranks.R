library(optparse)

source_path = system("echo $SOURCE", intern = TRUE)
utilities_script = "utils.R"
source(paste0(source_path, "/",utilities_script))


option_list <- list(
make_option(
"--file_pattern",
default = "comparisons/*/nSCORE/dea/rank.csv",
help = "pattern of nScore rank file",
type = "character"
),

make_option(
  "--outfile",
  default = "nScore_genes_rank_consolidated.csv",
  help = "consolidated rank file name",
  type = "character"
),

make_option(
  "--header",
  default = NULL,
  help = "header of consolidated rank file; comma separed list of colum names for consolidated rank file",
  type = "character"
),

make_option(
  "--datatype",
  default = "rank",
  help = "header of consolidated rank file; comma separed list of colum names for consolidated rank file",
  type = "character"
),

make_option(
  "--ntop",
  default = NULL,
  help = "number of top genes",
  type = "integer"
)
)
print(30)
opt <- parse_args(OptionParser(option_list = option_list))
print(32)
header = opt$header

if (!is.null(header))
{
  header = strsplit(header, ",")
  header = trimws(header)
}
if (opt$datatype=="rank")
{
  is_score_table = FALSE
}else
{
  is_score_table = TRUE
}
merge_nSCORE_rank_tables(file_pattern=opt$file_pattern, outfile=opt$outfile,
           header = header, is_score_table = is_score_table,  ntop=opt$ntop)