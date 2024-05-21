import sys

from pipelines.base_layer import *
from  pipelines.steps_lib import *
# import fastq_to_count as fq
from pipelines.utils import *


# class get_fastq_from_sra_list(pm.Manager):
#     input_options = {"srr_list": "file containing list of SRR accession number or SRR file location",
#                      "step_use_cluster": " if use cluster for each data set",
#                      "step_cluster_config": "cluster config file for use in each step"}
#     run_config_template = "get_fastq_from_SRR_list_config.csv"
#
#     def run(self, srr_list, step_use_cluster, step_cluster_config, **kwargs):
#         # SRR_list = self.inputs["srr_list"]
#         # step_use_cluster = self.inputs["step_use_cluster"]
#         # step_cluster_config = self.inputs["step_cluster_config"]
#         # print("11 pipelines lib, step_use_cluster:", step_use_cluster,"outdir:", self.outdir)
#
#         with open(srr_list, "r") as f:
#             for line in f:
#                 accession = line.rstrip()
#                 get_fastq_from_sra(accession=accession,
#                                    use_cluster=step_use_cluster,
#                                    cluster_config=step_cluster_config,
#                                    outdir=self.outdir + "/" + accession,
#                                    name=accession, job_log_dir=self.outdir + "/job_log_dir").controlled_run()


# class get_count_from_fastq(pm.Manager):
#     run_config_template = "fastq_to_count_config.csv"
#
#     def run(self, **kwargs):
#         p = fq.Pipeline()
#         p.loadConfiguration(sys.argv[2])
#         p.createjobs()
#         p.submitjobs()
#
#     pass



# class get_count_from_sra_list(pm.Manager):
#     input_options = {"srr_list": "file containing list of SRR accession number or SRR file location",
#                      "organism": " organism of reference genome",
#                      "ncpu": "number of cpu used for each sample alignment"}
#     run_config_template = "get_count_from_SRR_list_config.csv"
#
#     def run(self, srr_list, organism, ncpu, **kwargs):
#         # print("11 pipelines lib, step_use_cluster:", step_use_cluster,"outdir:", self.outdir)
#         print(srr_list, organism, ncpu, kwargs)
#         with open(srr_list, "r") as f:
#             for line in f:
#                 accession = line.rstrip()
#                 accession_outdir = self.outdir + "/" + accession
#                 get_count_from_sra(accession=accession, organism=organism, ncpu=ncpu,
#                                    outdir=accession_outdir,
#                                    name=accession, job_log_dir=self.outdir + "/job_log_dir",
#                                    ).controlled_run()


class consolidate_counts(Layer):
    """
    Consolidate star_read_counts from individual count file into one big count table where \
    columns are samples, row are genes. First row is symbol. Also translate gene EnsembleID into HUGO gene ID.
    Output files: consolidated_count_table.csv  consolidated_count_table_translated.csv  consolidated_cpm_table.csv  consolidated_cpm_table_translated.csv
    """
    outdir = 'CONSOLIDATED_COUNTS'

    @property
    def inputs_attributes(self):
        return OrderedDict([("outdir", "directory of the inputs"),
                            ("organism", "species of sample. Accepted values: human, mouse")])

    def execute(self, inputs, **kwargs):
        counts_dir = getattr(inputs, 'outdir', '.')
        organism = getattr(inputs, 'organism', None)

        pipe_ob = pipe_object()
        pipe_jobId = []
        for input in inputs.dataset:
            jId = input.get('pipe_jobId', None)
            if jId is not None:
                pipe_jobId.append(jId)
        data = {"counts_dir": counts_dir,
                "organism": organism,
                "pipe_jobId": pipe_jobId
                }

        pipe_ob.add(data)
        print('86 diff exp pipelines')
        out_obj = check_count_and_consolidate_step(simulated=self.simulated, outdir=self.outdir)(pipe_ob)

        # out_obj = pipe_object()
        # out_obj = self.set_general_output(out_obj, inputs)


        return out_obj

    @property
    def output_dataset_attributes(self):
        return OrderedDict([("count_table_file", "table containing count values of all samples, columns are (SYMBOL, sample1, sample2, etc), row are count value for each gene")])




class get_fit_data(dataset_layer):
    """ Apply batch effect removal and calculate fit data using edgeR from count table\
    Also create images of PCA, MDS, tsne of samples before and after batch effect removal.
    For samples that belongs to big group, there is option to break down into several small datasets.\
     each dataset is one big group that use the same network. Call get_fit_data_from_count_step as internal engine (that does not do breakdown to small datasets).
     Fit data can be used in multifactor experiment , allowing to do mach pair comparison (for example pair of slow vs fast GSC from the same patient).
     The fit data determines the coefficent for each factor in simple GLM model where each factor is independent from each other, no interactions.
      In the follwowing analysis step, need to define contrast, which factor, groups in this factor to do comparison
     """
    run_config_template = "get_fit_data_run_config.csv"
    cluster_config_template = "cluster_config.csv"
    out_file_pattern = "edgeR_fit_result/RUVs_k*.RDS"
    outdir = "fit_data"



    @property
    def configs_keys(self):
        return OrderedDict([('use_edger_glm_robust', "Use glm robust option in edgeR to increase accuracy"),
                            ("k",
                             "the k parameter for batch effect removal algorithm, representing the number of unwanted factors."
                             " The more k the more normalized datasets but the cons is dataset is overnormalized that can decrease the sensitivity of catching master regulators in nScore "
                             "if k =0, then not execute batch effect removal algorithm")
                            ])

    @property
    def output_dataset_attributes(self):
        return OrderedDict([
            ('fit_data', "RDS data file from R session containing fit data"),
            ('gene_expression_table', 'cpm gene expression table , will be used for heatmap in follwing step'),
            ('sample_group_table', "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..." )

        ])



    def call(self, input, inputs, output, index, outdir, **kwargs):
        count_table_file = input.get('count_table_file')
        sample_group_table = input.get('sample_group_table')
        break_data_big_group = input.get('break_data_big_group')
        sample_big_group_table = input.get('sample_big_group_table')
        use_edger_glm_robust = self.kwargs.get('use_edger_glm_robust')
        use_edger_glm_robust = str2bool(use_edger_glm_robust)
        break_data_big_group = str2bool(break_data_big_group)
        design_table = input.get('design_table', None)
        out = None

        if break_data_big_group and sample_big_group_table:
            raise ValueError("120 diff exp pipelines need to  do: make inputs, check for dependency (jobId)")
            break_count_data(count_table_file=count_table_file,
                             sample_group_table=sample_group_table,
                             sample_big_group_table=sample_big_group_table, outdir=self.outdir)
            big_group_list_file = self.outdir + "/big_groups_list.csv"
            with open(big_group_list_file, "r") as f:
                big_groups = f.readlines()
                for big_group in big_groups:
                    big_group = big_group.strip()
                    out = get_fit_data_from_count_step(use_edger_glm_robust=use_edger_glm_robust,
                                                       outdir=self.outdir + "/" + big_group, simulated= self.simulated).with_inputs(
                        count_table_file="{}/{}/count_table.csv".format(self.outdir, big_group),
                        sample_group_table="{}/{}/sample_groups.csv".format(self.outdir,
                                                                            big_group),
                        design_table=design_table,
                        **kwargs)

        else:
            in_ob = pipe_object()
            in_ob.add(input)
            out = get_fit_data_from_count_step(use_edger_glm_robust=use_edger_glm_robust,
                                               outdir=self.outdir, simulated=self.simulated)(in_ob)
            for data in out.dataset:
                output.add(data)
        return output



    @classmethod
    def get_inputs(cls, **kwargs):
        out = pipe_object()

        if "consolidated_counts" in kwargs:
            out = kwargs.get('consolidated_counts')
            if out is None:
                return None
            for key, value in kwargs.items():
                if key != "consolidated_counts":
                    for data in out.dataset:
                        data[key] = value
        else:
            data = {}
            for key, value in kwargs.items():
                data[key] = value
            out.add(data)

        return out

    @property
    def allowed_inputs_kwargs(self):
        return OrderedDict([
            ("count_table_file", "table containing count values of all samples, columns are (SYMBOL, sample1, sample2, etc), row are count value for each gene"),
            ("sample_group_table", "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("break_data_big_group", "True/False option to break data into several big group for making fit data for each big group. \
             This would make fitting faster as the data is smaller for each big group"),
            ("sample_big_group_table", "table assigning each sample to a big group, containing two columns (Sample, Big_Group)"),
            ("consolidated_counts", "pipe object from consolidate_counts"),
            ("design_table",
             "external provided experiment design factor table file (generated using R and sample_group_table using model.matrix function).  If not provided (design=Null) then autogenerate design table based on sample_factor_info,"
             " using simple linear model where each factor is independent, no interactions. design table :rows are sample IDs, columns are factors of experiments."
             " The table should have row names (sampleIDs), no column for sample IDs. Example of design table in Examples/design.csv")

        ])

    @property
    def inputs_dataset_attributes(self):
        return OrderedDict([
            ("count_table_file",
             "table containing count values of all samples, columns are (SYMBOL, sample1, sample2, etc), row are count value for each gene"),
            ("sample_group_table", "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("break_data_big_group", "True/False option to break data into several big group for making fit data for each big group. \
             This would make fitting faster as the data is smaller for each big group"),
            ("sample_big_group_table",
             "table assigning each sample to a big group, containing two columns (Sample, Big_Group)"),
            ("design_table",
             "external provided experiment design factor table file (generated using R and sample_group_table using model.matrix function).  If not provided (design=Null) then autogenerate design table based on sample_factor_info,"
             " using simple linear model where each factor is independent, no interactions. design table :rows are sample IDs, columns are factors of experiments."
             " The table should have row names (sampleIDs), no column for sample IDs. Example of design table in Examples/design.csv")

        ])
