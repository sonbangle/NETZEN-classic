from NetZen_pipeline_dev import *
from pipelines.diff_exp_pipelines import *
from pipelines.pipe_layer import Pipe


class NetZen_from_fit_data(Pipe):
    def execute(self, inputs, **kwargs):
        y = get_network_layout_from_fit_data()(inputs)
        subnets = get_subnets_from_network_layout()(y)
        print('11 NZp', subnets.dataset[0])
        image_jobs = get_subnet_image_jobs()(subnets)
        print('13 NZp', image_jobs)
        network_data = get_network_data_step()(image_jobs)
        network_images = save_network_image()(network_data)

        # Making web app
        dataset_name = getattr(inputs, 'dataset_name', 'RNAseq_dataset')
        data_folder = getattr(inputs, 'outdir', '.')
        web_outdir = getattr(inputs, 'web_outdir', '..')
        if data_folder is None:
            data_folder = "."
        print(112, dataset_name, data_folder)
        print(113, inputs.__dict__.keys())
        web_app = make_web_app()(network_images,
                                dataset_name=dataset_name,
                                data_folder=data_folder,
                                 web_outdir = web_outdir)

        return network_images


class NetZen_from_sra(Pipe):
    run_config_template = "NetZen_from_sra_run_config.csv"  # run config_template

    def execute(self, inputs, **kwargs):
        #counts = get_count_from_sra()(inputs)
        fastq = get_fastq_from_sra()(inputs)

        counts = get_count_from_fastq()(fastq)

        consolidated_counts = consolidate_counts()(counts)
        fit = get_fit_data()(consolidated_counts=consolidated_counts,
                                         sample_group_table=getattr(inputs, 'sample_group_table', None))

        comparisons = get_comparisons()(fit=fit,
                                                    big_group_networks_table=getattr(inputs, 'big_group_network_table',
                                                                                     None),
                                                    sample_big_group_table=getattr(inputs, 'sample_big_group_table',
                                                                                   None),
                                                    network_file=getattr(inputs, 'network_file', None))
        network_layouts = get_network_layout_from_fit_data()(comparisons)
        subnets = get_subnets_from_network_layout()(network_layouts)
        image_jobs = get_subnet_image_jobs()(subnets)
        network_data = get_network_data_step()(image_jobs)
        network_images = save_network_image()(network_data)

        # Making web app
        dataset_name = getattr(inputs, 'dataset_name', 'RNAseq_dataset')
        data_folder = getattr(inputs, 'outdir', '.')
        web_outdir = getattr(inputs, 'web_outdir', '..')
        if data_folder is None:
            data_folder = "."
        print(112, dataset_name, data_folder)
        print(113, inputs.__dict__.keys())
        web_app = make_web_app()(network_images,
                                dataset_name=dataset_name,
                                data_folder=data_folder,
                                 web_outdir = web_outdir)


        return network_images

    @classmethod
    def get_inputs(cls, **kwargs):
        inputs = get_fastq_from_sra.get_inputs(sra_list=kwargs.get('sra_list'))

        inputs.sample_group_table = kwargs.get('sample_group_table')
        inputs.big_group_network_table = kwargs.get('big_group_network_table')
        inputs.sample_big_group_table = kwargs.get('sample_big_group_table')
        inputs.network_file = kwargs.get('network_file')
        inputs.organism = kwargs.get('organism', 'human')
        return inputs

    @property
    def allowed_inputs_kwargs(self):
        kwargs = OrderedDict([('sra_list', 'list of sra accession numbers'),
                              ("organism", "Species of sample (human or mouse)"),
                              ("sample_group_table", "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
                              ("network_file", "Generep network used for calculating nSCORE"),
                              ("big_group_networks_table",
                               "table for network allocation for each big group, containing two column (Big_Group, Network)"),
                              ("sample_big_group_table",
                               "table assigning each sample to a big group, containing two columns (Sample, Big_Group)")

                              ])
        return kwargs

    @property
    def inputs_attributes(self):
        return OrderedDict([
            ("sample_group_table", "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("network_file", "Generep network used for calculating nSCORE"),
            ("big_group_network_table",
             "table for network allocation for each big group, containing two column (Big_Group, Network)"),
            ("sample_big_group_table",
             "table assigning each sample to a big group, containing two columns (Sample, Big_Group)"),
            ("organism", "Species of sample (human or mouse)"),
        ])

    @property
    def inputs_dataset_attributes(self):
        return OrderedDict([
                            ("accession", "accession number of SRR or the SRR filename")])




class NetZen_from_consolidated_count(Pipe):
    run_config_template = "NetZen_from_count_run_config.csv"  # run config_template
    def execute(self, inputs, **kwargs):
        fit = get_fit_data()(inputs)
        comparisons = get_comparisons()(fit=fit,
                                        big_group_networks_table=getattr(inputs, 'big_group_network_table', None),
                                        sample_big_group_table=getattr(inputs, 'sample_big_group_table', None),
                                        network_file=getattr(inputs, 'network_file', None),
                                        contrast_network_table=getattr(inputs, 'contrast_network_table', None),
                                        comparisons=getattr(inputs, 'comparisons', None))
        y = get_network_layout_from_fit_data()(comparisons)
        rank_consolidation = consolidate_nScore_ranks()(y)
        subnets = get_subnets_from_network_layout()(y)
        image_jobs = get_subnet_image_jobs()(subnets)
        network_data = get_network_data_step()(image_jobs)
        network_images = save_network_image()(network_data)

        # Making web app
        dataset_name = getattr(inputs, 'dataset_name', 'RNAseq_dataset')
        data_folder = getattr(inputs, 'outdir', '.')
        web_outdir = getattr(inputs, 'web_outdir', '..')
        if data_folder is None:
            data_folder = "."
        print(112, dataset_name, data_folder)
        print(113, inputs.__dict__.keys())
        web_app = make_web_app()(network_images,
                                dataset_name=dataset_name,
                                data_folder=data_folder,
                                 web_outdir = web_outdir)


        return network_images

    @property
    def allowed_inputs_kwargs(self):
        return OrderedDict([
            ("count_table_file",
             "table containing count values of all samples, columns are (SYMBOL, sample1, sample2, etc), row are count value for each gene"),
            ("sample_group_table", "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("network_file", "Generep network used for calculating nSCORE"),
            ("big_group_networks_table",
             "table for network allocation for each big group, containing two column (Big_Group, Network)"),
            ("sample_big_group_table",
             "table assigning each sample to a big group, containing two columns (Sample, Big_Group)"),
            ("design_table",
             "external provided experiment design table file (generated using R and sample_group_table using model.matrix function). If not provided (design=Null) then autogenerate design table based on sample_factor_info,"
             " using simple linear model where each factor is independent, no interactions. design table :rows are sample IDs, columns are factors of experiments."
             " The table should have row names (sampleIDs), no column for sample IDs. Example of design table in Examples/design.csv"),
            ("contrast_network_table", " table with column (name, contrast, network, baseline_samples or baseline_table)\
                     where name is name of comparison, "
                                       "contrast is a comma separarated string of values of factors to make contrast, "
                                       "network is name of Aracne network,"
                                       " baseline_samples or baseline_table are comma separated list of baseline_samples or name of baseline table."
                                       " If use baseline_samples, then all samples are normalized to the same baseline samples. "
                                       "baseline_table gives more flexiblity, can define different baseline samples for each sample"
                                       "Example: Examples/contrast_network_table.csv or Examples/contrast_network_table_with_baseline_table.csv."
                                       "This file usually generated using R (Examples/making_contrast_table.R "),
            ("break_data_big_group", "True/False option to break data into several big group for making fit data for each big group. \
                     This would make fitting faster as the data is smaller for each big group"),
            ("comparisons",
             "Comparisons pairs to compare. A list of comma separated comparisons. For example Fibroblast^ESC,GL261_Control^ESC"),
            ("dataset_name", "name of dataset/experiment"),
            ("data_folder", "data folder, parent folder of COUNTS folder"),
            ("web_outdir", "output folder for Web App")
        ])

    @property
    def inputs_attributes(self):
        return OrderedDict([
            ("count_table_file",
             "table containing count values of all samples, columns are (SYMBOL, sample1, sample2, etc), row are count value for each gene"),
            ("sample_group_table",
             "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("network_file", "Generep network used for calculating nSCORE"),
            ("big_group_network_table",
             "table for network allocation for each big group, containing two column (Big_Group, Network)"),
            ("sample_big_group_table",
             "table assigning each sample to a big group, containing two columns (Sample, Big_Group)"),
            ("design_table",
             "external provided experiment design factor table file (generated using R and sample_group_table using model.matrix function). If not provided (design=Null) then autogenerate design table based on sample_factor_info,"
             " using simple linear model where each factor is independent, no interactions. design table :rows are sample IDs, columns are factors of experiments."
             " The table should have row names (sampleIDs), no column for sample IDs. Example of design table in Examples/design.csv"),
            ("contrast_network_table", " table with column (name, contrast, network, baseline_samples or baseline_table)\
                     where name is name of comparison, "
                                       "contrast is a comma separarated string of values of factors to make contrast, "
                                       "network is name of Aracne network,"
                                       " baseline_samples or baseline_table are comma separated list of baseline_samples or name of baseline table."
                                       " If use baseline_samples, then all samples are normalized to the same baseline samples. "
                                       "baseline_table gives more flexiblity, can define different baseline samples for each sample"
                                       "Example: Examples/contrast_network_table.csv or Examples/contrast_network_table_with_baseline_table.csv."
                                       "This file usually generated using R (Examples/making_contrast_table.R "),
            ("comparisons",
             "Comparisons pairs to compare. A list of comma separated comparisons. For example Fibroblast^ESC,GL261_Control^ESC")
        ])

    @classmethod
    def get_inputs(cls,
                   dataset_name=None,
                   data_folder=None,
                   web_outdir="..",
                   **kwargs):
        inputs = get_fit_data.get_inputs(count_table_file=kwargs.get('count_table_file'),
                                         sample_group_table=kwargs.get('sample_group_table'),
                                         break_data_big_group=kwargs.get('break_data_big_group'),
                                         sample_big_group_table=kwargs.get('sample_big_group_table'),
                                         design_table=kwargs.get('design_table'))
        print("155 Nzplines", inputs.dataset)
        inputs.sample_group_table = kwargs.get('sample_group_table')
        inputs.big_group_network_table = kwargs.get('big_group_network_table')
        inputs.sample_big_group_table = kwargs.get('sample_big_group_table')
        inputs.network_file = kwargs.get('network_file')
        inputs.contrast_network_table = kwargs.get('contrast_network_table')
        inputs.comparisons = kwargs.get('comparisons')
        inputs.dataset_name = dataset_name
        inputs.data_folder = data_folder
        inputs.web_outdir = web_outdir
        return inputs

class NetZen_from_fastq(Pipe):
    """ nScore and subnet analysis package from fastq files with optional multifactor design table analysis"""
    run_config_template = "NetZen_from_fastq_run_config.csv"  # run config_template

    def execute(self, inputs, **kwargs):
        counts = get_count_from_fastq()(inputs)
        consolidated_counts = consolidate_counts()(counts)
        fit = get_fit_data()(consolidated_counts=consolidated_counts,
                             sample_group_table=getattr(inputs, 'sample_group_table', None),
                             design_table=getattr(inputs, 'design_table', None))

        comparisons = get_comparisons()(fit=fit,
                                        big_group_networks_table=getattr(inputs, 'big_group_network_table', None),
                                        sample_big_group_table=getattr(inputs, 'sample_big_group_table', None),
                                        network_file=getattr(inputs, 'network_file', None),
                                        contrast_network_table=getattr(inputs,'contrast_network_table', None),
                                        comparisons=getattr(inputs, 'comparisons', None))
        print("199 nzpipelines",comparisons.dataset)
        network_layouts = get_network_layout_from_fit_data()(comparisons)
        rank_consolidation = consolidate_nScore_ranks()(network_layouts)
        subnets = get_subnets_from_network_layout()(network_layouts)
        image_jobs = get_subnet_image_jobs()(subnets)
        network_data = get_network_data_step()(image_jobs)
        network_images = save_network_image()(network_data)

        dataset_name = getattr(inputs, 'dataset_name', 'RNAseq_dataset')
        data_folder = getattr(inputs, 'outdir', '.')
        if data_folder is None:
            data_folder = "."
        print(208, dataset_name, data_folder)
        web_app = make_web_app()(network_images,
                                dataset_name=dataset_name,
                                data_folder=data_folder)

        return network_images

    @classmethod
    def get_inputs(cls,
                   fastq_folder=None,
                   sample_fastq_table=None,
                   organism=None,
                   sample_group_table=None,
                   network_file=None,
                   big_group_network_table=None,
                   sample_big_group_table=None,
                   design_table=None,
                   contrast_network_table=None,
                   dataset_name=None,
                   data_folder=None,
                   **kwargs):
        inputs = get_count_from_fastq.get_inputs(organism=organism, sample_fastq_table=sample_fastq_table,
                                                 fastq_folder=fastq_folder)
        inputs.sample_group_table = sample_group_table
        inputs.big_group_network_table = big_group_network_table
        inputs.sample_big_group_table = sample_big_group_table
        inputs.network_file = network_file
        inputs.design_table = design_table
        inputs.contrast_network_table = contrast_network_table
        inputs.dataset_name = dataset_name
        inputs.data_folder = data_folder
        return inputs

    @property
    def allowed_inputs_kwargs(self):
        return OrderedDict([
            ("fastq_folder", "Folder containing fastq files"),
            ("organism", "Species of sample (human or mouse)"),
            ("sample_fastq_table",
             "table with two columns (Sample, Fastq), describing where fastq files located for each sampleID"),
            ("sample_group_table", "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("network_file", "Generep network used for calculating nSCORE"),
            ("big_group_networks_table",
             "table for network allocation for each big group, containing two column (Big_Group, Network)"),
            ("sample_big_group_table",
             "table assigning each sample to a big group, containing two columns (Sample, Big_Group)"),
            ("design_table",
             "external provided experiment design table file (generated using R and sample_group_table using model.matrix function). If not provided (design=Null) then autogenerate design table based on sample_factor_info,"
             " using simple linear model where each factor is independent, no interactions. design table :rows are sample IDs, columns are factors of experiments."
             " The table should have row names (sampleIDs), no column for sample IDs. Example of design table in Examples/design.csv"),
            ("contrast_network_table", " table with column (name, contrast, network, baseline_samples or baseline_table)\
                     where name is name of comparison, "
                                       "contrast is a comma separarated string of values of factors to make contrast, "
                                       "network is name of Aracne network,"
                                       " baseline_samples or baseline_table are comma separated list of baseline_samples or name of baseline table."
                                       " If use baseline_samples, then all samples are normalized to the same baseline samples. "
                                       "baseline_table gives more flexiblity, can define different baseline samples for each sample"
                                       "Example: Examples/contrast_network_table.csv or Examples/contrast_network_table_with_baseline_table.csv."
                                       "This file usually generated using R (Examples/making_contrast_table.R "),
            ("break_data_big_group", "True/False option to break data into several big group for making fit data for each big group. \
                 This would make fitting faster as the data is smaller for each big group"),
            ("comparisons",
             "Comparisons pairs to compare. A list of comma separated comparisons. For example Fibroblast^ESC,GL261_Control^ESC"),
            ("dataset_name", "name of dataset/experiment"),
            ("data_folder", "data folder, parent folder of COUNTS folder")
        ])

    @property
    def inputs_attributes(self):
        return OrderedDict([
            ("sample_group_table", "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("network_file", "Generep network used for calculating nSCORE"),
            ("big_group_network_table",
             "table for network allocation for each big group, containing two column (Big_Group, Network)"),
            ("sample_big_group_table",
             "table assigning each sample to a big group, containing two columns (Sample, Big_Group)"),
            ("design_table",
             "external provided experiment design factor table file (generated using R and sample_group_table using model.matrix function). If not provided (design=Null) then autogenerate design table based on sample_factor_info,"
             " using simple linear model where each factor is independent, no interactions. design table :rows are sample IDs, columns are factors of experiments."
             " The table should have row names (sampleIDs), no column for sample IDs. Example of design table in Examples/design.csv"),
            ("contrast_network_table", " table with column (name, contrast, network, baseline_samples or baseline_table)\
                     where name is name of comparison, "
                                       "contrast is a comma separarated string of values of factors to make contrast, "
                                       "network is name of Aracne network,"
                                       " baseline_samples or baseline_table are comma separated list of baseline_samples or name of baseline table."
                                       " If use baseline_samples, then all samples are normalized to the same baseline samples. "
                                       "baseline_table gives more flexiblity, can define different baseline samples for each sample"
                                       "Example: Examples/contrast_network_table.csv or Examples/contrast_network_table_with_baseline_table.csv."
                                       "This file usually generated using R (Examples/making_contrast_table.R "),
            ("comparisons",
             "Comparisons pairs to compare. A list of comma separated comparisons. For example Fibroblast^ESC,GL261_Control^ESC")
        ])

    @property
    def inputs_dataset_attributes(self):
        return OrderedDict([("sampleID", "sample Identification"),
                            ("organism", "Species of sample (human or mouse)"),
                            ("fastq_folder", "folder containing fastq files"),
                            ("fastq1",
                             "fastq file of sample, first file in case of pair end or the only file in case of single end"),
                            ("fastq2", "second fastq file in case of pair end sequencing. Ignore if not pair end")])


class fastq_to_consolidated_counts(Pipe):
    """ convert fastq files into consolidated count table"""
    run_config_template = "fastq_to_consolidated_counts_run_config.csv"  # run config_template

    def execute(self, inputs, **kwargs):
        counts = get_count_from_fastq()(inputs)
        consolidated_counts = consolidate_counts()(counts)
        if consolidate_counts is not None:
             return consolidated_counts

    @classmethod
    def get_inputs(cls,
                   fastq_folder=None,
                   sample_fastq_table=None,
                   organism=None,
                   **kwargs):
        inputs = get_count_from_fastq.get_inputs(organism=organism, sample_fastq_table=sample_fastq_table,
                                                 fastq_folder=fastq_folder)
        return inputs

    @property
    def allowed_inputs_kwargs(self):
        return OrderedDict([
            ("fastq_folder", "Folder containing fastq files"),
            ("organism", "Species of sample (human or mouse)"),
            ("sample_fastq_table",
             "table with two columns (Sample, Fastq), describing where fastq files located for each sampleID")
        ])


    @property
    def inputs_dataset_attributes(self):
        return OrderedDict([("sampleID", "sample Identification"),
                            ("organism", "Species of sample (human or mouse)"),
                            ("fastq_folder", "folder containing fastq files"),
                            ("fastq1",
                             "fastq file of sample, first file in case of pair end or the only file in case of single end"),
                            ("fastq2", "second fastq file in case of pair end sequencing. Ignore if not pair end")])



class nScore_from_consolidated_count(Pipe):
    run_config_template = "NetZen_from_count_run_config.csv"  # run config_template
    def execute(self, inputs, **kwargs):
        fit = get_fit_data()(inputs)
        comparisons = get_comparisons()(fit=fit,
                                        big_group_networks_table=getattr(inputs, 'big_group_network_table', None),
                                        sample_big_group_table=getattr(inputs, 'sample_big_group_table', None),
                                        network_file=getattr(inputs, 'network_file', None),
                                        contrast_network_table=getattr(inputs, 'contrast_network_table', None),
                                        comparisons=getattr(inputs, 'comparisons', None))
        y = get_network_layout_from_fit_data(get_nscore_only=True)(comparisons)
        rank_consolidation = consolidate_nScore_ranks()(y)

    @property
    def allowed_inputs_kwargs(self):
        return OrderedDict([
            ("count_table_file",
             "table containing count values of all samples, columns are (SYMBOL, sample1, sample2, etc), row are count value for each gene"),
            ("sample_group_table", "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("network_file", "Generep network used for calculating nSCORE"),
            ("big_group_networks_table",
             "table for network allocation for each big group, containing two column (Big_Group, Network)"),
            ("sample_big_group_table",
             "table assigning each sample to a big group, containing two columns (Sample, Big_Group)"),
            ("design_table",
             "external provided experiment design table file (generated using R and sample_group_table using model.matrix function). If not provided (design=Null) then autogenerate design table based on sample_factor_info,"
             " using simple linear model where each factor is independent, no interactions. design table :rows are sample IDs, columns are factors of experiments."
             " The table should have row names (sampleIDs), no column for sample IDs. Example of design table in Examples/design.csv"),
            ("contrast_network_table", " table with column (name, contrast, network, baseline_samples or baseline_table)\
                     where name is name of comparison, "
                                       "contrast is a comma separarated string of values of factors to make contrast, "
                                       "network is name of Aracne network,"
                                       " baseline_samples or baseline_table are comma separated list of baseline_samples or name of baseline table."
                                       " If use baseline_samples, then all samples are normalized to the same baseline samples. "
                                       "baseline_table gives more flexiblity, can define different baseline samples for each sample"
                                       "Example: Examples/contrast_network_table.csv or Examples/contrast_network_table_with_baseline_table.csv."
                                       "This file usually generated using R (Examples/making_contrast_table.R "),
            ("break_data_big_group", "True/False option to break data into several big group for making fit data for each big group. \
                     This would make fitting faster as the data is smaller for each big group"),
            ("comparisons",
             "Comparisons pairs to compare. A list of comma separated comparisons. For example Fibroblast^ESC,GL261_Control^ESC")
        ])

    @property
    def inputs_attributes(self):
        return OrderedDict([
            ("count_table_file",
             "table containing count values of all samples, columns are (SYMBOL, sample1, sample2, etc), row are count value for each gene"),
            ("sample_group_table",
             "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("network_file", "Generep network used for calculating nSCORE"),
            ("big_group_network_table",
             "table for network allocation for each big group, containing two column (Big_Group, Network)"),
            ("sample_big_group_table",
             "table assigning each sample to a big group, containing two columns (Sample, Big_Group)"),
            ("design_table",
             "external provided experiment design factor table file (generated using R and sample_group_table using model.matrix function). If not provided (design=Null) then autogenerate design table based on sample_factor_info,"
             " using simple linear model where each factor is independent, no interactions. design table :rows are sample IDs, columns are factors of experiments."
             " The table should have row names (sampleIDs), no column for sample IDs. Example of design table in Examples/design.csv"),
            ("contrast_network_table", " table with column (name, contrast, network, baseline_samples or baseline_table)\
                     where name is name of comparison, "
                                       "contrast is a comma separarated string of values of factors to make contrast, "
                                       "network is name of Aracne network,"
                                       " baseline_samples or baseline_table are comma separated list of baseline_samples or name of baseline table."
                                       " If use baseline_samples, then all samples are normalized to the same baseline samples. "
                                       "baseline_table gives more flexiblity, can define different baseline samples for each sample"
                                       "Example: Examples/contrast_network_table.csv or Examples/contrast_network_table_with_baseline_table.csv."
                                       "This file usually generated using R (Examples/making_contrast_table.R "),
            ("comparisons",
             "Comparisons pairs to compare. A list of comma separated comparisons. For example Fibroblast^ESC,GL261_Control^ESC")
        ])

    @classmethod
    def get_inputs(cls, **kwargs):
        inputs = get_fit_data.get_inputs(count_table_file=kwargs.get('count_table_file'),
                                         sample_group_table=kwargs.get('sample_group_table'),
                                         break_data_big_group=kwargs.get('break_data_big_group'),
                                         sample_big_group_table=kwargs.get('sample_big_group_table'),
                                         design_table=kwargs.get('design_table'))
        print("155 Nzplines", inputs.dataset)
        inputs.sample_group_table = kwargs.get('sample_group_table')
        inputs.big_group_network_table = kwargs.get('big_group_network_table')
        inputs.sample_big_group_table = kwargs.get('sample_big_group_table')
        inputs.network_file = kwargs.get('network_file')
        inputs.contrast_network_table = kwargs.get('contrast_network_table')
        inputs.comparisons = kwargs.get('comparisons')
        return inputs


class consolidated_count_from_sra(Pipe):
    def execute(self, inputs, **kwargs):
        fastq = get_fastq_from_sra()(inputs)
        print(396, fastq.dataset)
        counts = get_count_from_fastq()(fastq_object=fastq, organism=getattr(inputs, 'organism', 'human'))
        consolidated_counts = consolidate_counts()(counts)

    @classmethod
    def get_inputs(cls, **kwargs):
        inputs = get_fastq_from_sra.get_inputs(sra_list=kwargs.get('sra_list'))

        inputs.organism = kwargs.get('organism', 'human')
        return inputs

    @property
    def allowed_inputs_kwargs(self):
        kwargs = OrderedDict([('sra_list', 'list of sra accession numbers'),
                              ("organism", "Species of sample (human or mouse)")
                              ])
        return kwargs

    @property
    def inputs_attributes(self):
        return OrderedDict([
            ("organism", "Species of sample (human or mouse)"),
        ])

    @property
    def inputs_dataset_attributes(self):
        return OrderedDict([
                            ("accession", "accession number of SRR or the SRR filename")])




class NetZen_from_dea_file(Pipe):
    """ nScore and subnet analysis package from manually created differential expression analysis file with columns gene,logFC, pvaule, fdr, LR.
    (for example from t test result because of the input expression file is not in count format(FPKM or microarray data """
    run_config_template = "NetZen_from_fastq_run_config.csv"  # run config_template

    def execute(self, inputs, **kwargs):
        network_layouts = get_network_layout_from_fit_data(cmd="pipeline_before_MR")(inputs)
        #rank_consolidation = consolidate_nScore_ranks()(network_layouts)
        subnets = get_subnets_from_network_layout()(network_layouts)
        image_jobs = get_subnet_image_jobs()(subnets)
        network_data = get_network_data_step()(image_jobs)
        network_images = save_network_image()(network_data)

        return network_images

    @classmethod
    def get_inputs(cls,
                   sample_group_table=None,
                   network_file=None,
                   gene_expression_table=None,
                   displayed_groups = None,
                   baseline_samples = None,
                   baseline_tables = None,
                   differential_expression_analysis_file=None,
                   **kwargs):

        out = pipe_object()
        out.sample_group_table = sample_group_table
        data = {}
        data["sample_group_table"] = sample_group_table
        data["network_file"] = network_file
        data["gene_expression_table"] = gene_expression_table
        data["displayed_groups"] = displayed_groups
        data["baseline_samples"] = baseline_samples
        data["baseline_table"] = baseline_tables
        data["differential_expression_analysis_file"] = differential_expression_analysis_file
        out.add(data)

        return out

    @property
    def allowed_inputs_kwargs(self):
        return OrderedDict([
            ("sample_group_table", "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("network_file", "Generep network used for calculating nSCORE"),
            ("gene_expression_table", "cpm gene expression table , will be used for heatmap in following step"),
            ('displayed_groups',
             'groups of samples to be displayed in heatmap of gene expression among different groups'),
            ('pair_comparison_groups', 'pair of two comparing groups source - target'),
            ("baseline_samples", "A comma separated list  of baseline samples for heatmap normalization."
                                 " The expresion value of each gene for each sample will be divided to the averge value of baseline_samples."
                                 "If NULL, then average of all samples will be used as baseline to normalize gene expression values"),
            ("baseline_table",
             "name of table describing for each sample which other samples should be used as baseline to rescale for heatmap drawing."
             "Contain two colulmn: SampleID, Baseline_samples, where in the column Baseline_samples is the comma separated list of basseline samples for corresponding SampleID."
             " This table takes priority over baseline_samples parameter"),
            ("differential_expression_analysis_file",
             " table produced by t test or edgeR, containing column :gene, logFC, pvalue, fdr, LR. used for feeding Master Regulator Score algorithm")

        ])

    @property
    def inputs_attributes(self):
        return OrderedDict([])

    @property
    def inputs_dataset_attributes(self):
        return OrderedDict([
             ("sample_group_table",
             "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("network_file", "Generep network used for calculating nSCORE"),
            ("gene_expression_table", "cpm gene expression table , will be used for heatmap in following step"),
            ('displayed_groups',
             'groups of samples to be displayed in heatmap of gene expression among different groups'),
            ('pair_comparison_groups', 'pair of two comparing groups source - target'),
            ("baseline_samples", "A comma separated list  of baseline samples for heatmap normalization."
                                 " The expresion value of each gene for each sample will be divided to the averge value of baseline_samples."
                                 "If NULL, then average of all samples will be used as baseline to normalize gene expression values"),
            ("baseline_table",
             "name of table describing for each sample which other samples should be used as baseline to rescale for heatmap drawing."
             "Contain two colulmn: SampleID, Baseline_samples, where in the column Baseline_samples is the comma separated list of basseline samples for corresponding SampleID."
             " This table takes priority over baseline_samples parameter"),
            ("differential_expression_analysis_file",
             " table produced by t test or edgeR, containing column :gene, logFC, pvalue, fdr, LR. used for feeding Master Regulator Score algorithm")
    ])
