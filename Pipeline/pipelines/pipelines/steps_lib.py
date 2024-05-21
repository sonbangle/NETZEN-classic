from pipelines.utils import *
import pickle
from pipelines.base_layer import Layer, pipe_object
from pipelines.submit_layer import submit_layer

source_path = get_script_path()
from collections import OrderedDict


class get_fastq_from_sra(submit_layer):
    modules = None #""#"sra" #"sra/2.8.0"
    create_individual_data_outdir = True
    outdir = "FASTQS"
    cluster_config_template = "cluster_config_3Gb_48h.conf"
    @property
    def inputs_dataset_attributes(self):
        return OrderedDict([("accession", "accession number of SRR or the SRR filename")])

    @property
    def output_dataset_attributes(self):
        return OrderedDict([("fastq_dir", "directory containing fastq files for each sample")])

    @property
    def allowed_inputs_kwargs(self):
        return OrderedDict([('sra_list', 'list of sra accession numbers')])

    def get_data_outdir(self, input, index):
        return input['accession']

    def set_output_data(self, inputs, pipe_jobId=None, index=None, outdir=None, **kwargs):
        input = inputs.dataset[index]
        data = {}
        data['fastq_dir'] = outdir
        data['sampleID'] = input['sampleID']
        return data

    def get_cmdline(self, outdir=None, accession=None, **kwargs):
        cmdline = """
echo "beginning sra to fastq at $(date)"
retry=0
while [ $retry -le 5 ]
do
{{

    module load sra && fastq-dump --split-e -L 5  -O {}  {} && break
    
}} ||  {{ 

   echo error loading sra
   retry=$(( $retry + 1 ))
   sleep 120
   
}}
done 


echo "done sra to fastq at $(date "+%Y-%m-%d %H:%M:%S") "
""".format(outdir, accession)
        #print(cmdline)
        return cmdline

    @classmethod
    def get_inputs(cls, **kwargs):
        srr_list = kwargs['sra_list']
        out = pipe_object()
        with open(srr_list, "r") as f:
            for line in f:
                accession = line.rstrip()
                data = {"accession": accession}
                data["sampleID"] = accession
                out.add(data)
        return out


class get_count_from_sra(submit_layer):
    # This program convert from sra accession number to count table
    modules = ["sra", "python"] #"sra/2.8.0"
    create_individual_data_outdir = True
    cluster_config_template = "cluster_config_90Gb.conf"
    outdir = "COUNTS"

    @property
    def inputs_keys(self):
        inputs_keys = [
            "accession",  # "accession number of SRR",
            "organism",  # "organism to align to reference genome. Either human or mouse",
        ]
        return inputs_keys

    def get_cmdline(self, accession, organism="human", ncpu=30, **kwargs):
        cmdline = """
echo "beginning sra to count at $(date)"

sra_accession_to_count.sh  {}  {} {} {}

echo "done sra to count at $(date "+%Y-%m-%d %H:%M:%S") "
""".format(accession, organism, self.outdir + "/" + accession, ncpu)
        print(cmdline)
        return cmdline

    def get_data_outdir(self, input, index):
        outdir = input['accession']
        return outdir

    def set_output(self, output, inputs, **kwargs):
        output.organism = getattr(inputs, "organism", "human")
        return output

    @classmethod
    def get_inputs(cls, **kwargs):
        allowed_kwargs = ["sra_list", "organism"]
        cls.check_kw_error(kwargs, allowed_kwargs)
        srr_list = kwargs['sra_list']
        organism = kwargs.get('organism', "human")

        out = pipe_object()
        with open(srr_list, "r") as f:
            for line in f:
                accession = line.rstrip()
                data = {"accession": accession, "organism": organism}
                out.add(data)
        out.organism = organism
        return out


class get_count_from_fastq(submit_layer):
    """This program convert from fastq files  to count data
    """
    #modules = ["sra", "python"] #"sra/2.8.0"
    modules = ["gcc/8.3.0", "python/3.9.2"]
    create_individual_data_outdir = True
    cluster_config_template = "cluster_config_90Gb_cpu_64.conf"
    outdir = "COUNTS"
    run_config_template = "get_count_from_fastq_run_config.conf"

    @property
    def configs_keys(self):
        return OrderedDict([("ncpu", "Number of cpu used for star alignment"),
                            ("remove_bam", "remove bam files after alignment to save storage")])

    @property
    def inputs_dataset_attributes(self):
        return OrderedDict([("sampleID", "sample Identification"),
                            ("organism", "Species of sample (human or mouse)"),
                            ("fastq_folder", "parent folder containing all fastq files"),
                            ("fastq_dir", "directory containing fastq files for each sample (used in case of pulling fastq from SRA) in contrast of fastq_folder that contains fastq files for all samples"),
                            ("fastq1",
                             "fastq file of sample, first file in case of pair end or the only file in case of single end"),
                            ("fastq2", "second fastq file in case of pair end sequencing. Ignore if not pair end")])

    @property
    def output_attributes(self):
        return OrderedDict([("organism", "Species of sample (human or mouse)")])

    def get_cmdline(self, sampleID, organism="human", outdir=".", ncpu=62, fastq1=None, fastq2=None, remove_bam=False, fastq_dir=None, **kwargs):
        remove_bam = str2bool(remove_bam)
        if remove_bam:
            remove_bam = "1"
        else:
            remove_bam = "0"
        if fastq_dir:  # case fastq from SRA, need to figure out fastq file names from fastq_dir

            cmdline = """
echo "beginning fastq to count at $(date)"
fastq_files=({}/*)

fastq_to_count.sh  {}  {} {} {} {} """.format(fastq_dir, sampleID, organism, outdir, ncpu, remove_bam)

            cmdline += """${fastq_files[0]}  ${fastq_files[1]}

echo "done fastq to count at $(date "+%Y-%m-%d %H:%M:%S") "
           """

        else:
            if fastq2 is not None:
                cmdline = """
echo "beginning fastq to count at $(date)"
    
fastq_to_count.sh  {}  {} {} {} {} {} {}
    
echo "done fastq to count at $(date "+%Y-%m-%d %H:%M:%S") "
            """.format(sampleID, organism, outdir, ncpu, remove_bam, fastq1, fastq2)
            else:
                cmdline = """
echo "beginning fastq to count at $(date)"

fastq_to_count.sh  {}  {} {} {} {} {}

echo "done fastq to count at $(date "+%Y-%m-%d %H:%M:%S") "
                """.format(sampleID, organism, outdir, ncpu, remove_bam, fastq1)
        return cmdline

    def get_data_outdir(self, input, index):
        outdir = input['sampleID']
        return outdir

    def set_output(self, output, inputs, **kwargs):
        output.organism = getattr(inputs, "organism", "human")
        return output

    @classmethod
    def get_inputs(cls, fastq_object=None, organism="human", sample_fastq_table=None, fastq_folder=".", **kwargs):
        import pandas as pd
        if fastq_object is None:
            out = pipe_object()
            out.organism = organism
            sample_fastq = pd.read_csv(sample_fastq_table, sep="\t")
            sample_fastq.columns = ["Sample", "Fastq"]
            samples = sample_fastq.Sample.unique().tolist()
            for sample in samples:
                #print(210, sample)
                fasq_files = sample_fastq[sample_fastq["Sample"] == sample]["Fastq"].tolist()
                data = {}
                data["sampleID"] = sample
                data["organism"] = organism
                #print(214, fastq_folder)
                #print(215, fasq_files)
                data["fastq1"] = fastq_folder + "/" + fasq_files[0]
                if len(fasq_files) == 2:
                    data["fastq2"] = fastq_folder + "/" + fasq_files[1]
                else:
                    data["fastq2"] = None
                out.add(data)
            return out
        else:
            fastq_object.organism = organism
            for sample in fastq_object.dataset:
                sample['organism'] = organism
            return fastq_object


    @property
    def allowed_inputs_kwargs(self):
        return OrderedDict([
            ("organism", "Species of sample (human or mouse)"),
            ("sample_fastq_table",
             "table with two columns (Sample, Fastq), describing where fastq files located for each sampleID"),
            ("fastq_object", "pipe object from get_fastq_from_sra step")
        ])


class check_count_and_consolidate_step(submit_layer):
    run_config_template = "check_counts_consolidate_config.conf"
    cluster_config_template = "cluster_config_15Gb.conf"

    @property
    def inputs_dataset_attributes(self):
        inputs_keys = OrderedDict([
            ("counts_dir", "Directory containing counts for all sample, default=COUNTS"),
            ("organism", "Species of sample"),
            ("outdir", "output directory")])
        return inputs_keys

    def get_cmdline(self, counts_dir="COUNTS", organism="human", outdir=".", **kwargs):
        cmdline = "check_count_and_consolidate.sh " + counts_dir + " " + organism + " " + outdir
        return cmdline

    def set_output_data(self, inputs, pipe_jobId=None, index=None, outdir=None, **kwargs):
        data = {}
        # data["count_table_file"] = outdir + "/consolidated_counts_out/consolidated_count_table_translated.csv"
        data["count_table_file"] = outdir + "/count/gene_level/consolidated_count_table_translated.csv"
        return data

    @property
    def output_dataset_attributes(self):
        return OrderedDict([("count_table_file",
                             "table containing count values of all samples, columns are (SYMBOL, sample1, sample2, etc), row are count value for each gene")])


class get_fit_data_from_count_step(submit_layer):
    """ Apply batch effect removal and calculate fit data using edgeR from count table\
     Also create images of PCA, MDS, tsne of samples before and after batch effect removal.

     """

    Rscript_file = "network_layout_pipeline_061018.R"
    #modules = ["gcc/11.3.0",   "openblas/0.3.20", "python/3.11.3", "R/4.2.3"]
    #modules = ["gcc/11.3.0",   "openblas/0.3.20", "python/3.11.3", "r/4.2.3"]
    cluster_config_template = "cluster_config_30Gb.conf"

    out_file_pattern = "edgeR_fit_result/RUVs_k*.RDS"
    conda_env = "check_count_consolidate"

    # out_file_pattern = "consolidated_counts_out/consolidated*"# specific output file pattern , use to check if the pipeline is really done or not
    # n_out_files = 4  # should exist files: consolidated_count_table.csv , consolidated_cpm_table.csv, consolidated_count_table_translated.csv,  consolidated_cpm_table_translated.csv
    @property
    def configs_keys(self):
        return OrderedDict([('use_edger_glm_robust',
                             "Use glm robust option in edgeR to increase accuracy"),
                            ("k",
                             "the k parameter for batch effect removal algorithm, representing the number of unwanted factors."
                             " The more k the more normalized datasets but the cons is dataset is overnormalized that can decrease the sensitivity of catching master regulators in nScore"
                             "if k =0, then not execute batch effect removal algorithm ")])

    @property
    def inputs_dataset_attributes(self):
        return OrderedDict([
            ("count_table_file",
             "table containing count values of all samples, columns are (SYMBOL, sample1, sample2, etc), row are count value for each gene"),
            ("sample_group_table", "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("design_table",
             "external provided experiment design factor table file (generated using R and sample_group_table using model.matrix function). If not provided (design=Null) then autogenerate design table based on sample_factor_info,"
             " using simple linear model where each factor is independent, no interactions. design table :rows are sample IDs, columns are factors of experiments."
             " The table should have row names (sampleIDs), no column for sample IDs. Example of design table in Examples/design.csv")

        ])

    @property
    def allowed_inputs_kwargs(self):
        return OrderedDict([
            ("count_table_file",
             "table containing count values of all samples, columns are (SYMBOL, sample1, sample2, etc), row are count value for each gene"),
            ("sample_group_table",
             "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("design_table",
             "external provided experiment design factor table file (generated using R and sample_group_table using model.matrix function). If not provided (design=Null) then autogenerate design table based on sample_factor_info,"
             " using simple linear model where each factor is independent, no interactions. design table :rows are sample IDs, columns are factors of experiments."
             " The table should have row names (sampleIDs), no column for sample IDs. Example of design table in Examples/design.csv")
        ])

    @property
    def output_dataset_attributes(self):
        return OrderedDict([
            ('fit_data', "RDS data file from R session containing fit data"),
            ('gene_expression_table', 'cpm gene expression table , will be used for heatmap in follwing step'),
            ('sample_group_table', "Table describing samples. Contains two columns: Sample, Group")

        ])

    def get_cmdline(self, count_table_file,
                    sample_group_table,
                    use_edger_glm_robust=True,
                    outdir="fit_data",
                    k=None,
                    design_table=None):
        if outdir is None:
            outdir = "fit_data"
        assert count_table_file is not None, "count table file is None"
        assert sample_group_table is not None, "sample group table is None"
        cmdline = "source ~/.bashrc; conda activate {}\n".format(self.conda_env)

        cmdline += "Rscript $SOURCE/" + self.Rscript_file + \
                  " --cmd=get_fit_data_from_count" + \
                  " --gene_count_table=" + count_table_file + \
                  " --sample_group_table=" + sample_group_table + \
                  " --outdir=" + outdir
        if use_edger_glm_robust:
            cmdline += " --use_edgeR_GLM_robust"
        if k is not None:
            cmdline += " --k=" + str(k)
        if design_table is not None:
            cmdline += " --design_table=" + design_table
        # cmdline = "echo test"
        print(cmdline)
        return cmdline

    def set_output_data(self, inputs, pipe_jobId=None, index=None, outdir=None, **kwargs):
        input = inputs.dataset[index]
        data = {}
        if outdir is None:
            outdir = "fit_data"
        data['fit_data'] = outdir + "/edgeR_fit_result/fit_result.RDS"
        data['gene_expression_table'] = outdir + "/cpm_table.csv"
        data['sample_group_table'] = input.get('sample_group_table')
        return data

    @classmethod
    def get_inputs(cls, count_table_file=None, sample_group_table=None):
        data = {'count_table_file'  : count_table_file,
                'sample_group_table': sample_group_table}
        output = pipe_object()
        output.add(data)
        return output


class break_count_data(submit_layer):
    Rscript_file = "network_layout_pipeline_061018.R"

    def get_cmdline(self, count_table_file=None,
                    sample_group_table=None,
                    sample_big_group_table=None,
                    outdir=None,
                    **kwargs):  # this option allow good result but a little bit slow)
        print('160 steps lib', locals())
        cmdline = "Rscript $SOURCE/" + self.Rscript_file + \
                  " --cmd=break_count_data" + \
                  " --gene_count_table=" + count_table_file + \
                  " --sample_group_table=" + sample_group_table + \
                  " --sample_big_group_table=" + sample_big_group_table + \
                  " --outdir=" + outdir

        print(cmdline)
        return cmdline

    @property
    def inputs_keys(self):
        return [
            "count_table_file",
            "sample_group_table",
            "sample_big_group_table"
        ]

    @classmethod
    def get_inputs(cls, count_table_file=None, sample_group_table=None, sample_big_group_table=None):
        data = {'count_table_file'      : count_table_file,
                'sample_group_table'    : sample_group_table,
                'sample_big_group_table': sample_big_group_table}
        output = pipe_object()
        output.add(data)
        print(184, vars(output))
        return output


class zip_file(submit_layer):
    """
    Compress files using zip
    """

    def get_cmdline(self, input_file):
        cmdline = "gzip {}".format(input_file)
        return cmdline

    @property
    def inputs_dataset_attributes(self):
        return OrderedDict([("input_file", "name of input file to be zipped")])


