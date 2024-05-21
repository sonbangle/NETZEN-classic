# from network_parser_vtk_053118 import *
import os.path
import pickle

# import pipelines.pipeline_manager as pm
import pandas as pd

source_folder = "$SOURCE"
Rscript_file = "network_layout_pipeline_061018.R"
from pipelines.submit_layer import submit_layer
from pipelines.base_layer import *
import random
from pipelines.submit_cluster import submit

# from pipelines.diff_exp_pipelines import get_fit_data
source_path = get_script_path()
from pipelines import backend as B


# import pipeline
# print(8, os.environ["PATH"])


class get_network_data_step(submit_layer):
    """ Get vtu file from network layout"""
    cluster_config_template = "cluster_config_18Gb.conf"
    cluster_config = source_path + "/config_templates/cluster_configs" + "/" + cluster_config_template
    #modules = "vtk"
    conda_env = "vtk"
    python_script = "network_parser_vtk_020119.py"
    use_cluster = True
    create_individual_data_outdir = True

    @property
    def inputs_dataset_attributes(self):
        return OrderedDict([
            ("net_file", "network layout file"),
            ("size_field", "field to determine size of gene node"),
            ("filter_condition", "filter to select subnet"),
            ("node_size", "unit node size"),
            ("resolution", "node sphere resolution in Paraview"),
            ("max_node_size", "max node size for paraview"),
            ("scale", "layout coordinate scale in vtu"),
            ("get_edges", "True/False show edges in subnet image"),
            ("n_top_gene", " number of gene nodes in a simplified subnet"),
            ("outfile_vtu", "output subnet vtu file for importing into Paraview")
        ])

    @property
    def configs_keys(self):
        return OrderedDict([
            ("node_col_up", "color of upregulated subnet"),
            ("node_col_down", "color of downregulated subnet"),
            ("normalize_node_location", "True/False do normalize node location so center of system is at (0,0)"),
            ("text_display_option", "list of text to attached to subnet node. choose from (label, annotation, masters"),
            ("n_top_gene", "Number of master regulators to display for each subnet node")
        ])

    @property
    def output_dataset_attributes(self):
        return OrderedDict([
            ("net_file", "network layout file"),
            ("size_field", "field to determine size of gene node"),
            ("filter_condition", "filter to select subnet"),
            ("node_size", "unit node size"),
            ("resolution", "node sphere resolution in Paraview"),
            ("max_node_size", "max node size for paraview"),
            ("scale", "layout coordinate scale in vtu"),
            ("get_edges", "True/False show edges in subnet image"),
            ("n_top_gene", " number of gene nodes in a simplified subnet"),
            ("outfile_vtu", "output subnet vtu file for importing into Paraview"),
            ("vtu_file", "vtu file for importing into paraview")])

    def get_cmdline(self,
                    net_file,
                    intensity_field=None,
                    group_field=None,
                    intensity_scale_factor=2.,
                    get_edges=False,
                    outfile_vtu="network_data.vtu",
                    node_col_up=(1., 0.,
                                 0.,
                                 1.),
                    node_col_down=(
                            0.,
                            0.,
                            1.,
                            1.),
                    up_limit=3,
                    down_limit=3,
                    size_field="radius",
                    node_size=2,
                    resolution=8,
                    fc_limit=1,
                    # fold change limit, only include the nodes that outside the range(-fc_limit, fc_limit)
                    max_node_size=None,
                    scale=1000,
                    normalize_node_location=False,
                    export_to_vtk=True,
                    filter_condition=[
                        {
                            "field": "subnet_name",
                            "op": "startswith",
                            "compared_value": "1.1.1.1"}],
                    text_display_options=[
                        "label",
                        "annotation",
                        "masters"],
                    n_top_gene=3,
                    outdir=".",

                    **kwargs

                    ):
        pickle_file = outfile_vtu + "_filter_condition.pickle"
        with open(pickle_file, "wb") as f:
            pickle.dump(filter_condition, f, protocol=2)
            # Protocol is 2 for compatibility with python 2 in vtk
        cmdline = "source ~/.bashrc\n"
        cmdline += "conda activate {}\n".format(self.conda_env)
        cmdline += "python " + source_folder + "/" + self.python_script + \
                  " --net_file=" + net_file
        if intensity_field is not None:
            cmdline += " --intensity_field=" + intensity_field
        if group_field is not None:
            cmdline += "--group_field=" + group_field
        cmdline += " --intensity_scale_factor=" + str(intensity_scale_factor) + \
                   " --get_edges=" + str(get_edges) + \
                   " --outfile=" + outfile_vtu
        cmdline += " --node_col_up "
        for col in node_col_up:
            cmdline += " " + str(col)
        cmdline += " --node_col_down "
        for col in node_col_down:
            cmdline += " " + str(col)

        cmdline += " --up_limit=" + str(up_limit) + \
                   " --down_limit=" + str(down_limit) + \
                   " --size_field=" + size_field + \
                   " --node_size=" + str(node_size) + \
                   " --resolution=" + str(resolution) + \
                   " --fc_limit=" + str(fc_limit)
        if max_node_size is not None:
            cmdline += " --max_node_size=" + str(max_node_size)
        cmdline += " --scale=" + str(scale) + \
                   " --export_to_vtk=" + str(export_to_vtk) + \
                   " --filter_condition_pickle=" + pickle_file

        cmdline += " --text_display_options "
        for opt in text_display_options:
            cmdline += " " + opt
        if n_top_gene is not None:
            cmdline += " --n_top_gene " + str(n_top_gene)
        # print(cmdline)

        sys.stdout.flush()
        # self.outfile = outfile_vtu
        return cmdline

    def set_output_data(self, inputs, pipe_jobId=None, index=None, outdir=None, **kwargs):

        input = inputs.dataset[index]
        data = input
        vtu_file = input.get('outfile_vtu', None)
        if vtu_file is None:
            vtu_file = self.get_defaults().get('outfile_vtu', "default.vtu")
        data['vtu_file'] = vtu_file

        # To implement in childeren class
        return data

    def get_data_outdir(self, input, index):
        return input['data_outdir']


class save_network_image(submit_layer):
    """ Save paraview network image from vtu vtk file"""
    cluster_config_template = "cluster_config_9Gb.conf"

    run_config_template = "save_network_image_run_config.conf"
    #modules = ["vtk", "paraview/5.4.1"]
    #modules = ["gcc/8.3.0", "openmpi/4.0.1", "vtk/8.1.2","paraview/5.6.2"]
    check_ready_to_run = False
    DISPLAY = 1 + 2000 * random.randint(0, 10)  # Display id  value

    # for GUI partition to run paraview on Hipergator server. Random starting value so that not collisons between different batch jobs

    @property
    def inputs_dataset_attributes(self):
        inputs_keys = OrderedDict([
            ("vtu_file", "vtu file for importing into paraview"),
            ("outfile", "name of outfile png file"),
            ("coloring_field", "field used for coloring node"),
            ("min_val", "minimal value"),
            ("max_val", "max value"),
            ("opacity", "node opacity"),
            ("SetScalarBarVisibility", "True/False to show scalar bar in paraview rendered image"),
            ("ImageResolution_w", "Setting for paraview image resolution width in pixel"),
            ("ImageResolution_h", "Setting for paraview image resolution height in pixel"),
        ])

        return inputs_keys

    @property
    def output_dataset_attributes(self):
        return OrderedDict([("network_image_file", "name of subnet image file")])

    def get_cmdline(self,
                    vtu_file=None,
                    outfile=None,
                    ImageResolution_w=2000,
                    ImageResolution_h=1500,
                    coloring_field="radius",
                    min_val=-1.,
                    max_val=1.,
                    opacity=1.,
                    SetScalarBarVisibility=False,
                    **kwargs
                    ):
        if outfile is None:
            outfile = vtu_file + ".png"
        # cmdline = """
        #
        # export DISPLAY=:{}
        # Xvfb :{} -screen 0 1024x768x16 &
        #
        # """.format(self.__class__.DISPLAY, self.__class__.DISPLAY)

        cmdline = """

DISPLAY_ID=%s
export DISPLAY=:${DISPLAY_ID}
Xvfb :${DISPLAY_ID} -screen 0 1024x768x16 &

        """ % self.__class__.DISPLAY

        self.__class__.DISPLAY += 1
        if self.__class__.DISPLAY == 30000:
            self.__class__.DISPLAY = 1  # Resetting  DISPLAY ID

        # cmdline = """
        #
        # DISPLAY_ID=$RANDOM
        # export DISPLAY=:${DISPLAY_ID}
        # Xvfb :${DISPLAY_ID} -screen 0 1024x768x16 &
        #
        # """
        cmdline = ""
        cmdline += "export PATH=${MY_BIN}/ParaView-5.11.1-osmesa-MPI-Linux-Python3.9-x86_64/bin:$PATH\n"
        cmdline += "export LD_LIBRARY_PATH=${MY_BIN}/ParaView-5.11.1-osmesa-MPI-Linux-Python3.9-x86_64/lib:$LD_LIBRARY_PATH"
        cmdline += "echo $PATH\n"
        cmdline += "echo $LD_LIBRARY_PATH\n"
        cmdline += "pvpython --force-offscreen-rendering " + source_folder + "/paraview_save_network.py --infile=" + vtu_file + \
                   " --outfile=" + outfile + \
                   " --ImageResolution_w=" + str(ImageResolution_w) + \
                   " --ImageResolution_h=" + str(ImageResolution_h) + \
                   " --coloring_field=" + coloring_field + \
                   " --min_val=" + str(min_val) + \
                   " --max_val=" + str(max_val) + \
                   " --opacity=" + str(opacity)
        if SetScalarBarVisibility:
            cmdline += " --SetScalarBarVisibility=" + str(SetScalarBarVisibility)
        return cmdline

    def set_output_data(self, inputs, pipe_jobId=None, index=None, outdir=None, **kwargs):
        data = {}
        input = inputs.dataset[index]
        data['network_image_file'] = input['vtu_file'] + ".png"
        return data


# get subnet network where only genes are displayed
class get_subnet_image_jobs(dataset_layer):
    cluster_config_template = "cluster_config_30Gb.conf"
    cluster_config = source_path + "/config_templates/cluster_configs" + "/" + cluster_config_template

    @property
    def inputs_dataset_attributes(self):
        return OrderedDict([
            ("subnet", "subnetwork in a list of significant subnetworks from sign_subnet_table"),
            ("net_file_small_space_ratio", "Gene network layout with small space ratio, nodes are close to each other"),
            ("net_file_large_space_ratio", "Gene network layout with large space ratio, nodes are far to each other"),
            ("data_outdir", "data output directory")
        ])

    @property
    def configs_keys(self):
        return OrderedDict([
            ("size_field_gene", "field to determine size of gene node "),
            ("size_field_subnet", " field to determine size of subnet node"),
            ("n_top_gene", " number of gene nodes in a simplified subnet"),
            ("resolution", "node sphere resolution in Paraview"),
            ("scale", "layout coordinate scale in vtu"),
            ("node_size_ratio", "node size ratio"),
            ("SetScalarBarVisibility", "True/False to show scalar bar in paraview rendered image"),
            ("coloring_field", "field used for coloring node"),
            ("subnet_coloring_with_sign", "subnet colored with different hue for dow regulation and upregulation")
        ])

    @property
    def output_dataset_attributes(self):
        return OrderedDict([
            ("net_file", "network layout file"),
            ("size_field", "field to determine size of gene node"),
            ("filter_condition", "filter to select subnet"),
            ("node_size", "unit node size"),
            ("resolution", "node sphere resolution in Paraview"),
            ("max_node_size", "max node size for paraview"),
            ("scale", "layout coordinate scale in vtu"),
            ("get_edges", "True/False show edges in subnet image"),
            ("n_top_gene", " number of gene nodes in a simplified subnet"),
            ("check_done", "True/False need to check if this layer has been done"),
            ("outfile_png", "name of outfile png file"),
            ("coloring_field", "field used for coloring node"),
            ("min_val", "minimal value"),
            ("max_val", "max value"),
            ("opacity", "node opacity"),
            ("SetScalarBarVisibility", "True/False to show scalar bar in paraview rendered image"),
            ("data_outdir", "data output folder"),
            ("subnet", "subnet Id"),
            ("outfile_vtu", "output subnet vtu file for importing into Paraview"),
            ("ImageResolution_w", "Setting for paraview image resolution width in pixel"),
            ("ImageResolution_h", "Setting for paraview image resolution height in pixel")
        ])

    def run(self,
            subnet="1.1.1.1",
            net_file_small_space_ratio="layout/1/out_graph.graphml",
            net_file_large_space_ratio="layout/4/out_graph.graphml",
            size_field_gene="score",  # "score" is another option,
            size_field_subnet="radius",  # "score" is another option,
            n_top_gene=20,
            resolution=24,
            scale=1,
            node_size_ratio=1,
            SetScalarBarVisibility=False,
            coloring_field="logFC",
            subnet_coloring_with_sign=True,
            outdir="network/subnets/1.1.1.1/3D_view",
            **kwargs):

        # Convert into appropriate formats so that  no error when parsing from config file t
        print("get subnet image step")
        n_top_gene = int(n_top_gene)
        resolution = int(resolution)
        scale = float(scale)
        node_size_ratio = float(node_size_ratio)

        filter_con = [{"field": "subnet_name", "op": "startswith", "compared_value": subnet}]
        filter_con_node_only = [{"field": "subnet_name", "op": "startswith", "compared_value": subnet},
                                {"field": "net_type", "op": "==", "compared_value": 3}]
        filter_con_subnet_only = [{"field": "name", "op": "startswith", "compared_value": subnet}]
        mkdir(outdir)



        get_edges = [False, True]

        network_types_dic = {"gene_only": {"n_top_gene": None, "max_node_size": 10,
                                           "filter": filter_con_node_only,
                                           "size_field": size_field_gene,
                                           "node_size": (
                                               0.003 * node_size_ratio, 0.003 * node_size_ratio / 10),
                                           "opacity": 1,
                                           "ImageResolution_w": 20000,
                                           "ImageResolution_h": 15000,
                                           },
                             "gene_only_top_gene": {"n_top_gene": n_top_gene,
                                                    "max_node_size": 2,
                                                    "filter": filter_con_node_only,
                                                    "size_field": size_field_gene,
                                                    "node_size": (
                                                        0.015 * node_size_ratio, 0.015 * node_size_ratio / 5),
                                                    "opacity": 1,
                                                    "ImageResolution_w": 4000,
                                                    "ImageResolution_h": 3000,

                                                    },
                             "subnet_and_gene": {"n_top_gene": None,
                                                 "max_node_size": None,
                                                 "filter": filter_con,
                                                 "size_field": size_field_subnet,
                                                 "node_size": (0.8 * node_size_ratio, 2 * node_size_ratio),
                                                 "opacity": 0.65,
                                                 "ImageResolution_w": 20000,
                                                 "ImageResolution_h": 15000,
                                                 },
                             "subnet_only": {"n_top_gene": None,
                                             "max_node_size": None,
                                             "filter": filter_con_subnet_only,
                                             "size_field": size_field_subnet,
                                             "node_size": (0.8 * node_size_ratio, 2 * node_size_ratio),
                                             "opacity": 0.4,
                                             "ImageResolution_w": 20000,
                                             "ImageResolution_h": 15000,
                                             }
                             }
        # network_types = ["gene_only", "subnet_and_gene"]
        network_types = ["gene_only", "gene_only_top_gene", "subnet_and_gene", "subnet_only"]
        # network_types = ["gene_only_top_gene"]
        job_list = []

        for edges in get_edges:
            for network_type in network_types:
                network_config = network_types_dic[network_type]
                net_files = [net_file_small_space_ratio, net_file_large_space_ratio]
                for i in range(len(net_files)):
                    # net_file = network_config["net_file"]
                    net_file = net_files[i]
                    filter = network_config["filter"]
                    size_field = network_config["size_field"]
                    node_size = network_config["node_size"][i]
                    top_gene = network_config["n_top_gene"]
                    max_node_size = network_config["max_node_size"]
                    ImageResolution_w = network_config["ImageResolution_w"]
                    ImageResolution_h = network_config["ImageResolution_h"]

                    # outfile_vtu = subnet + "_edge" + str(edges) + "_" + network_type + str(i) + ".vtu"
                    outfile_vtu = outdir + "/e" + str(edges) + "_" + network_type + str(i) + ".vtu"

                    sys.stdout.flush()
                    outfile_png = outdir + "/" + subnet + "_edge" + str(edges) + "_" + network_type + str(i) + ".png"
                    # save network image
                    if (network_type == "subnet_and_gene") or (network_type == "subnet_only"):

                        max_val = 1.0
                        if subnet_coloring_with_sign:
                            min_val = -1.0
                            coloring_field_paraview = "enrichment_subnet_score_with_sign"
                        else:
                            min_val = 0.0
                            coloring_field_paraview = "enrichment_subnet_score"
                    else:

                        min_val = -1.0
                        max_val = 1.0
                        coloring_field_paraview = coloring_field
                    opacity = network_config["opacity"]
                    job = {
                        "net_file": net_file,
                        "size_field": size_field,
                        "filter_condition": filter,
                        "node_size": node_size,
                        "resolution": resolution,
                        "max_node_size": max_node_size,
                        "scale": scale,
                        "get_edges": edges,
                        "n_top_gene": top_gene,
                        "check_done": False,
                        "outfile_png": outfile_png,
                        "coloring_field": coloring_field_paraview,
                        "min_val": min_val,
                        "max_val": max_val,
                        "opacity": opacity,
                        "SetScalarBarVisibility": SetScalarBarVisibility,
                        "data_outdir": outdir,
                        "subnet": subnet,
                        "outfile_vtu": outfile_vtu,
                        "ImageResolution_w": ImageResolution_w,
                        "ImageResolution_h": ImageResolution_h,
                    }
                    job_list.append(job)

                    # print("282 NZpd", job['outfile_vtu'], job['coloring_field'])

        # self.PM.data["get_network_data_step"] = job_list
        # self.PM.outdir = outdir
        return job_list

    def call(self, input, inputs, output, index, outdir, subnet_dir="/network/subnets", **kwargs):
        # converting input kwargs into kwargs for run method
        if kwargs is None:
            kwargs = {}

        kwargs.update(self.kwargs)
        kwargs.update(input)

        # setting the output directory for this step:
        subnets_dir = input['data_outdir'] + subnet_dir
        mkdir(subnets_dir)
        subnet_dir = subnets_dir + "/" + input['subnet']
        mkdir(subnet_dir)
        view_dir = subnet_dir + "/3Dview"
        mkdir(view_dir)
        kwargs['outdir'] = view_dir

        job_list = self.run(**kwargs)
        for job in job_list:
            output.add(job)
        return output


def mkdir(dir):
    if not os.path.exists(dir):
        try:
            os.makedirs(dir)
        except Exception as e:
            print("918, cannot make dir, probably folder exists as another process is creating it")
            print(e)
            pass


class get_network_layout_from_fit_data(submit_layer):
    """
    From fit data calculate nScore. Based on nScore, make subnet analysis, network layout. Output network layout
    """
    # modules = "R/3.4.3, python, pandoc" ["gcc/11.3.0",   "openblas/0.3.20", "python/3.11.3", "R/4.2.3"]
    #modules = "R/3.5.1, python, pandoc"
    #modules = "gcc/11.3.0, openblas/0.3.20, python/3.11.3 r/4.2.3, pandoc"
    conda_env  = "check_count_consolidate"
    use_cluster = True
    cluster_config_template = "cluster_config_30Gb.conf"
    create_individual_data_outdir = True
    data_outdir_key = "pair_comparison_group"

    @property
    def inputs_dataset_attributes(self):
        return OrderedDict([
            ("network_file", "Generep network used for calculating nSCORE"),
            ("gene_expression_table", "cpm gene expression table , will be used for heatmap in following step"),
            ("sample_group_table", "Table describing samples. Contains two columns: Sample, Group"),
            ('displayed_groups',
             'groups of samples to be displayed in heatmap of gene expression among different groups'),
            ('pair_comparison_groups', 'pair of two comparing groups source - target'),
            ('fit_data', "RDS data file from R session containing fit data"),
            ("contrast", "contrast for edgeR algorithm as a string of comma separated list of values for design factors. If not provided, then contrast will be made from pair_comparison groups \
To know which factors the contrast codes fore, need to look at design table. Contrast example: '-0.5,-0.5,1.' \
 Number of elements must equal to the number of columns of design. Contrast can be made by using makeContrasts function in R limma package. For exmaple:\
BvsA <- makeContrasts(B-A, levels=design)"),
            ("baseline_samples", "A comma separated list  of baseline samples for heatmap normalization."
                                 " The expresion value of each gene for each sample will be divided to the averge value of baseline_samples."
                                 "If NULL, then average of all samples will be used as baseline to normalize gene expression values"),
            ("baseline_table",
             "name of table describing for each sample which other samples should be used as baseline to rescale for heatmap drawing."
             "Contain two colulmn: SampleID, Baseline_samples, where in the column Baseline_samples is the comma separated list of basseline samples for corresponding SampleID."
             " This table takes priority over baseline_samples parameter"),
            ("differential_expression_analysis_file", " table produced by t test or edgeR, containing column :gene, logFC, pvalue, fdr, LR. used for feeding Master Regulator Score algorithm")

        ])

    @property
    def output_dataset_attributes(self):
        return OrderedDict([('sign_subnet_file', "table about significantly disturbed subnetwork"),
                            ('network_layout_file',
                             'network file in graphml format where nodes are genes with layout coordinates by fr algorithm and node score ')])

    @property
    def configs_keys(self):
        return OrderedDict([
            ("cmd", "pipeline command. Accepted values: pipeline_from_fit_data, pipeline_before_MR"),
            ("graph_format", "graph format for igraph. Accepted values: ncol, graphml"),
            ("clustering_method", "clusetering method : louvain"),
            ("max_cluster_size", "maximum number of nodes in a subcluster. Use for breakdown big network"),
            ("ndim", " number of dimmension of network layout:[2, 3]"),
            ("layout_method", "network layout method from igraph layout methods:[fr]"),
            ("node_score_field",
             "field to calculate diffused (propagated) node score in subnetwork analysis algorithm:[logFC, FC, fdr]"),
            ("sampling_size", " sampling size in subnetwork analysis algorithm"),
            ("beta",
             "beta parameter, a leakage proportion from parent subnet to a child subnet in subnetwork analysis algorithm"),
            ("n_top_annotation", "number of top annotaion in subnet analysis"),
            ("importance_field",
             "importance field in subnet analysis algorithm to select Master Regulators to show in heatmap and pathway annotations"),
            ("n_top_genes", "number of top local master genes for a subnet in subnet analysis - annotation algorithm"),
            ("space_ratios",
             "space ratio for calculation network layout, the more space ratio, the fare away the nodes and subnet"),
            ("n_top_gene_network_drawing", "number of gene nodes in a simplified subnet"),
            ("resolution", "node sphere resolution in Paraview"),
            ("scale", "layout coordinate scale in vtu"),
            ("node_size_ratio", "node size ratio"),

            ("draw_heatmap", "drawing heatmap between comparison groups"),
            ("selected_fields", " selected fields"),
            ("subnet_coloring_with_sign", "subnet colored with different hue for dow regulation and upregulation"),

            ("consider_positive_values_only", " in nScore,only calculate score for overexpressed genes"),
            ("gene_statistics_list",
             "comma separted list of fields to consider when calculating nScore, for example:'logFC,fdr,betweenness' "),
            ("fdr_to_confidence", "convert fdr to confidence score before nScore analysis"),
            ("nround", "nround in nScore analysis"),
            ("steps_combined", "combine steps in nScore analysis"),
            ("top_genes", "proportion of top differentialy expressed genes to consider in nScore "),
            ("source_node_inclusion", "method of including sourced node in nScore analysis"),
            ("neighbor_aggregation_method", "method of aggregation onf neighbor score in nScore analysis"),
            ("is_log_gene_expression_table", " is gene expression table in log scale"),
            # ("do_get_subnet_image", "True/False make 3D subnet image"),
            ("do_extract_subnet_data", "True/False, extract subnet data for each subnet, "
                                       "including heatmap and all expression tables, gene layout."
                                       " If not True, then only data for the root subnet is extracted"),
            ("use_cluster_to_draw_heatmap",
             "True/False submit draw heatmap function to cluster to parallelize this function, making it run faster"),
            ("gene_description_file",
             "file describing gene functions from RefSeq. If not exists, then it will be autogenereted but ensure that the folder containing file exists"),
            ("export_top_master_genes_data", "True/False export heatmap, subnetworks of top master genes"),
            ("get_nscore_only", "True/False only get nScore, not calculating network layout")
        ])

    def get_cmdline(self,
                    cmd="pipeline_from_fit_data",
                    network_file=source_path + "/../../data/Math_dataset/mouse_embryonic_stem_cell_net.csv",
                    graph_format="ncol",
                    clustering_method="louvain",
                    max_cluster_size=100,
                    ndim=2,
                    layout_method="fr",
                    node_score_field="logFC",
                    sampling_size=1000,
                    beta=0.5,
                    n_top_annotation=5,
                    importance_field="score",
                    n_top_genes=5,
                    gene_description_file="/ufrc/dtran/son.le/MY_BIN/data//all_gene_summary.csv",
                    space_ratios=[1, 4],
                    n_top_gene_network_drawing=20,
                    resolution=24,
                    scale=1,
                    node_size_ratio=1,
                    draw_heatmap=False,
                    gene_expression_table=source_path + "/../../data//Math_dataset/Net_Zene/cpm_table.csv",
                    sample_group_table=source_path + "/../../data/Math_dataset/sample_groups.csv",
                    displayed_groups=None,  # ["Fibroblast", "ESC", "iPSC", "GL261_Control"],

                    selected_fields="node_name,score,logFC,description,summary,subnet_name,degree,pvalue,fdr,annotation",
                    subnet_coloring_with_sign=True,

                    pair_comparison_groups=None,  # ["Fibroblast", "ESC"],
                    gene_count_table=None,
                    consider_positive_values_only=True,
                    gene_statistics_list="logFC,fdr,betweenness",
                    fdr_to_confidence=False,
                    nround=1,
                    steps_combined=True,
                    top_genes=0.2,
                    source_node_inclusion="s",
                    neighbor_aggregation_method="s",

                    fit_data=source_path + "/../../data/Math_dataset/Net_Zene/edgeR_fit_result/RUVs_k50edgeR_original.RDS",
                    is_log_gene_expression_table=False,
                    # do_get_subnet_image=False,
                    do_extract_subnet_data=True,
                    outdir=".",
                    contrast=None,

                    baseline_samples=None,
                    use_cluster_to_draw_heatmap=False,
                    export_top_master_genes_data=False,
                    get_nscore_only=False,
                    baseline_table=None,
                    differential_expression_analysis_file=None,
                    **kwargs):
        # final_kwargs = locals()
        # self.update_final_kwargs(final_kwargs=final_kwargs)

        self.prt("beginning single pair comparison for comparison:", pair_comparison_groups)
        use_cluster_to_draw_heatmap = str2bool(use_cluster_to_draw_heatmap)
        max_cluster_size = int(max_cluster_size)
        ndim = int(ndim)
        sampling_size = int(sampling_size)
        beta = float(beta)
        n_top_annotation = int(n_top_annotation)
        n_top_genes = int(n_top_genes)
        if type(space_ratios).__name__ == 'str':
            space_ratios = str2int_list(space_ratios)
        n_top_gene_network_drawing = int(n_top_gene_network_drawing)
        resolution = int(resolution)
        scale = float(scale)
        node_size_ratio = float(node_size_ratio)
        draw_heatmap = str2bool(draw_heatmap)
        displayed_groups = str2list(displayed_groups)
        subnet_coloring_with_sign = str2bool(subnet_coloring_with_sign)
        pair_comparison_groups = str2list(pair_comparison_groups)

        consider_positive_values_only = str2bool(consider_positive_values_only)
        fdr_to_confidence = str2bool(fdr_to_confidence)
        nround = int(nround)
        steps_combined = str2bool(steps_combined)
        top_genes = float(top_genes)
        is_log_gene_expression_table = str2bool(is_log_gene_expression_table)
        # do_get_subnet_image = str2bool(do_get_subnet_image)
        do_extract_subnet_data = str2bool(do_extract_subnet_data)

        space_ratios_cmdl = list2str(space_ratios)
        export_top_master_genes_data = str2bool(export_top_master_genes_data)

        # print("674 NZpd, gene desc:", gene_description_file)
        options = " --cmd=" + cmd + \
                  " --network=" + network_file + \
                  " --graph_format=" + graph_format + \
                  " --clustering_method=" + clustering_method + \
                  " --max_cluster_size=" + str(max_cluster_size) + \
                  " --ndim=" + str(ndim) + \
                  " --layout_method=" + layout_method + \
                  " --space_ratios=" + space_ratios_cmdl + \
                  " --node_score_field=" + node_score_field + \
                  " --sampling_size=" + str(sampling_size) + \
                  " --beta=" + str(beta) + \
                  " --n_top_annotation=" + str(n_top_annotation) + \
                  " --importance_field=" + importance_field + \
                  " --n_top_genes=" + str(n_top_genes) + \
                  " --outdir=" + outdir + \
                  " --gene_description_file=" + gene_description_file
        if is_log_gene_expression_table:
            options += " --is_log_gene_expression_table"

        # if cmd == "score_de_attributes_not_in_network":  # need to provide differential expression analysis file with column gene, logFC, pvalue, fdr, LR
        #     options += " --diff_expr_file=" + diff_expr_file + \
        #                " --selected_fields=" + selected_fields
        elif cmd == "score_de_attributes_in_network":  # for example, node score, logfC, etc are imported in graphml network via Cytoscape
            options += " --selected_fields=" + selected_fields

        elif cmd == "main" or cmd == "pipeline_from_fit_data" or cmd == "pipeline_before_MR":

            options += \
                " --gene_statistics_list=" + gene_statistics_list + \
                " --nround=" + str(nround) + \
                " --top_genes=" + str(top_genes) + \
                " --source_node_inclusion=" + source_node_inclusion + \
                " --neighbor_aggregation_method=" + neighbor_aggregation_method
            if pair_comparison_groups is not None:
                pair_comparison_groups_cml = ",".join(pair_comparison_groups)
                options += " --pair_comparison_groups=" + pair_comparison_groups_cml

            if steps_combined:
                options += " --steps_combined"
            if consider_positive_values_only:
                options += " --consider_positive_values_only"
            if fdr_to_confidence:
                options += " --fdr_to_confidence"
            if cmd == "main":
                if gene_count_table is not None:
                    options += " --gene_count_table=" + gene_count_table
            elif cmd == "pipeline_from_fit_data":
                options += " --fit_data=" + fit_data
            elif cmd == "pipeline_before_MR":
                options += " --diff_expr_file=" + differential_expression_analysis_file

        if draw_heatmap and (displayed_groups is not None):
            displayed_groups_cml = ",".join(displayed_groups)
            options += " --draw_heatmap " + \
                       " --gene_expression_table=" + gene_expression_table + \
                       " --sample_group_table=" + sample_group_table + \
                       " --displayed_groups=" + displayed_groups_cml + \
                       " --n_top_genes_heatmap=" + str(n_top_gene_network_drawing)


        if do_extract_subnet_data:
            options += " --do_extract_subnet_data"
        if contrast is not None:
            options += " --contrast=" + contrast
        print("724 Nzpd", baseline_table)
        if baseline_table is not None:
            options += " --baseline_table=" + baseline_table
        else:
            if baseline_samples is not None:
                options += " --baseline_samples=" + baseline_samples

        if export_top_master_genes_data:
            options += " --export_top_master_genes_data"
        if use_cluster_to_draw_heatmap:
            options += " --use_cluster_to_draw_heatmap"
        if get_nscore_only:
            options += " --get_nscore_only"
        cmdline = "source ~/.bashrc \n"
        cmdline += "conda activate {} \n ".format(self.conda_env)
        cmdline += "Rscript " + source_folder + "/" + Rscript_file + options
        #         if draw_heatmap and use_cluster_to_draw_heatmap:
        #             job_dir = outdir + "/jobs/jobs"
        #             submit_jobs_cmdline = \
        #             """
        #
        # for job in $(ls %s/*)
        # do sbatch %s/${job}
        # done
        #             """ % (job_dir, job_dir)
        #             cmdline += submit_jobs_cmdline
        print(cmdline)

        return cmdline

    def set_output_data(self, inputs, pipe_jobId=None, index=None, outdir=None, **kwargs):
        if outdir is None:
            outdir = self.outdir
        if 'space_ratios' in self.kwargs:
            space_ratios = self.kwargs.get('space_ratios')
        else:
            space_ratios = self.get_defaults().get('space_ratios')

        sign_subnet_file = outdir + "/network/sign_subnets.txt"
        network_layout_file = outdir + "/network/layout/" + str(space_ratios[0]) + "/" + "out_graph.graphml"
        data = {}
        data["sign_subnet_file"] = sign_subnet_file
        data["network_layout_file"] = network_layout_file

        return data

    def get_data_outdir(self, input, index):
        comparison = input['pair_comparison_groups']
        input_outdir = input.get('outdir', ".")
        if type(comparison).__name__ == "str":
            outdir = input_outdir  # when pair comparison is actually manually entered in contrast_network table
        else:
            outdir = "{}/{}_{}".format(input_outdir, comparison[0], comparison[1])
        return outdir

    def set_output(self, output, inputs, **kwargs):
        output.rank_file_pattern = "comparisons\/\*\/nSCORE\/dea\/ranks.csv"
        output.score_file_pattern = "comparisons\/\*\/nSCORE\/dea\/scores.csv"
        return output


class get_subnets_from_network_layout(dataset_layer):
    """
    Create a list of subnet and network layot for significant subnets
    """
    use_cluster = False
    check_ready_to_run = True

    @property
    def inputs_dataset_attributes(self):
        input_keys = OrderedDict([
            ("sign_subnet_file", "table about significantly disturbed subnetwork")])
        return input_keys

    @property
    def configs_keys(self):
        return OrderedDict([("do_get_subnet_image", "True/False flag to get subnet image from vtu file")])

    def call(self, inputs, layout_dir="network/layout", **kwargs):
        output = pipe_object()
        if 'space_ratios' in inputs.layer_config:
            space_ratios = inputs.layer_config["space_ratios"]
        else:
            space_ratios = inputs.cmdline_defaults['space_ratios']

        for i, data in enumerate(inputs.dataset):

            sign_subnet_file = data["sign_subnet_file"]
            do_get_subnet_image = self.kwargs.get('do_get_subnet_image', True)
            print("814 NZpd sign_subnet_file:", sign_subnet_file)
            # sign_subnet_file = outdir + "/layout/" + str(space_ratios[0]) + "/" + "sign_subnets.txt"
            subnets = ["1"]
            net_file_small_space_ratio = data['data_outdir'] + "/" + layout_dir + "/" +\
                                         str(space_ratios[0]) + "/out_graph.graphml"
            net_file_large_space_ratio = data['data_outdir'] + "/" + layout_dir + "/" +\
                                         str(space_ratios[1]) + "/out_graph.graphml"
            if do_get_subnet_image:
                self.wait_for_file(sign_subnet_file)
                self.wait_for_file(net_file_small_space_ratio)
                self.wait_for_file(net_file_large_space_ratio)
                with open(sign_subnet_file, "r") as f:
                    line = f.readline()  # skip first line
                    while line:
                        line = f.readline()
                        #print("line", line)
                        parsed_line = line.split("\t")
                        if len(parsed_line) > 1:
                            subnet = parsed_line[0]
                            level = int(parsed_line[2])
                            # print(subnet, level)
                            if level <= 2:
                                subnets.append(subnet)

            # view_outdir = self.PM.outdir + "/layout/subnet_views"
            # self.PM.data["view_outdir"] = view_outdir

            # comparison = data['kwargs'].get("pair_comparison_groups")
            # comparison = "{}_{}".format(comparison[0], comparison[1])
            for subnet in subnets:
                out_data = {"subnet": subnet,
                            "net_file_small_space_ratio": net_file_small_space_ratio,
                            "net_file_large_space_ratio": net_file_large_space_ratio,
                            "data_outdir": data['data_outdir']}
                output.add(out_data)
        return output

    @property
    def output_dataset_attributes(self):
        return OrderedDict([
            ("subnet", "subnetwork in a list of significant subnetworks from sign_subnet_table"),
            ("net_file_small_space_ratio", "Gene network layout with small space ratio, nodes are close to each other"),
            ("net_file_large_space_ratio", "Gene network layout with large space ratio, nodes are far to each other"),
            ("data_outdir", "data output directory")
        ])


class get_comparisons(data_layer):
    """ This class makes contrasts (source - target )for following steps.
    In case samples do have big groups (brain cancer that contains both GBM and LGG) the comparisons are only among groups within big groups.
    """
    network_folder = B.source_path + "/../../data/networks"
    outdir = "comparisons"

    @property
    def inputs_attributes(self):
        inputs_attributes = OrderedDict([
            ("big_group_networks_table",
             "table for network allocation for each big group, containing two column (Big_Group, Network)"),
            ("gene_expression_table", "cpm gene expression table , will be used for heatmap in follwing step"),
            ("sample_group_table",
             "Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ..."),
            ("sample_big_group_table",
             "table assigning each sample to a big group, containing two columns (Sample, Big_Group)"),
            ('fit_data', "RDS data file from R session containing fit data"),
            ("network_file", "Generep network used for calculating nSCORE"),
            ("pipe_jobId", "internal jobId for a submit request"),
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
        return inputs_attributes

    @property
    def allowed_inputs_kwargs(self):
        return OrderedDict([
            ("fit", "output pipe object from get_fit_data Layer"),
            ("big_group_networks_table",
             "table for network allocation for each big group, containing two column (Big_Group, Network)"),
            ("sample_big_group_table",
             "table assigning each sample to a big group, containing two columns (Sample, Big_Group)"),
            ("network_file", "Generep network used for calculating nSCORE"),
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
    def output_dataset_attributes(self):
        return OrderedDict([
            ("network_file", "Generep network used for calculating nSCORE"),
            ("gene_expression_table", "cpm gene expression table , will be used for heatmap in follwing step"),
            ("sample_group_table", "Table describing samples. Contains two columns: Sample, Group"),
            ('pair_comparison_groups', 'pair of two comparing groups source - target'),
            ('fit_data', "RDS data file from R session containing fit data"),
            ('displayed_groups',
             'groups of samples to be displayed in heatmap of gene expression among different groups'),
            ('outdir', 'output directory'),
            ('pipe_jobId', 'pipe_jobId'),
            ("baseline_samples", "A comma separated list  of baseline samples for heatmap normalization."
                                 " The expresion value of each gene for each sample will be divided to the averge value of baseline_samples."
                                 "If NULL, then average of all samples will be used as baseline to normalize gene expression values")
        ])

    def run(self,
            gene_expression_table=None,
            sample_group_table=None,
            sample_big_groups_table=None,
            big_group_networks_table=None,
            network_folder=None,
            fit_data=None,
            network_file=None,
            pipe_jobId=None,
            contrast_network_table=None,
            comparisons=None,
            **kwargs):
        output = pipe_object()
        sample_groups = pd.read_csv(sample_group_table, sep="\t")
        # print("877 Nzpd", sample_group_table, sample_groups)
        # all_columns = list(sample_groups.columns)[1:]
        # sample_groups["Group"] = ""
        # for index, column_name in enumerate(all_columns):
        #     print("current column {column_name}".format(column_name=column_name))
        #     sample_groups["Group"] += "_" + sample_groups[column_name].astype(str)
        sample_groups["Group"] = ['_'.join(row.astype(str)) for row in sample_groups[sample_groups.columns[1:]].values]
        sample_groups.rename(columns={sample_groups.columns[0]: "Sample"}, inplace=True)
        # print("882 Nzpd", sample_groups)
        groups = sample_groups["Group"].unique()
        if contrast_network_table is None:
            if network_folder is None:
                network_folder = self.network_folder
            c = 0

            def add_data(groups, network, outdir, output, c, comparisons=None):
                groups = groups.astype(str).tolist()

                # print("Selected group:", groups)
                def _add_data(pair_comparison_groups, _outdir="."):
                    sys.stdout.flush()
                    # c += 1
                    data = {
                        'network_file': network,
                        'gene_expression_table': gene_expression_table,
                        'sample_group_table': sample_group_table,
                        'pair_comparison_groups': pair_comparison_groups,
                        'fit_data': fit_data,
                        'displayed_groups': groups,
                        'outdir': _outdir,
                        'pipe_jobId': pipe_jobId
                    }
                    output.add(data)

                    # print("counter:", c)
                    # print(pair_comparison_groups)
                    # print("kwargs:", kwargs)

                if comparisons is None:
                    for target in groups:
                        # if target.endswith("cancer_stem_cell"):
                        if True:
                            # print("target:", target)
                            target_dir = outdir + "/to_" + target
                            mkdir(target_dir)
                            for source in groups:
                                if source != target:
                                    # print("source:", source)
                                    pair_comparison_groups = [source, target]
                                    _add_data(pair_comparison_groups=pair_comparison_groups, _outdir=target_dir)
                else:
                    if type(comparisons).__name__ == 'str':
                        comparisons_new = []
                        input_list = str2list(comparisons)
                        for input in input_list:
                            comparison = input.split("^")
                            comparisons_new.append(comparison)
                        comparisons = comparisons_new
                    for comparison in comparisons:
                        _outdir = outdir + "/" + "^".join(comparison)
                        _add_data(pair_comparison_groups=comparison, _outdir=_outdir)

            if sample_big_groups_table is not None and big_group_networks_table is not None:
                sample_big_groups = pd.read_csv(sample_big_groups_table, sep="\t")
                sample_big_groups.columns = ["Sample", "Big_Group"]
                # print(sample_big_groups.loc[1:3, :])
                big_groups = sample_big_groups.Big_Group.unique()
                networks = pd.read_csv(big_group_networks_table, sep="\t")
                networks.columns = ["Big_Group", "Network"]

                # Select comparison


                for big_group in big_groups:  # For each cancer type, select samples belong to this cancer type to identify comparison only within groups
                    # print("processing big group:", big_group)

                    # Finding network for this particular big_group:
                    networks['Network'] = networks.Network.astype(str)
                    network = networks[networks.Big_Group == big_group]["Network"]
                    network = network.to_string(header=False, index=False)
                    network = network_folder + "/" + network
                    # print("204 NZ big_group network:", network)
                    # Select samples in this big groups:
                    big_group_outdir = self.outdir + "/" + big_group
                    mkdir(big_group_outdir)
                    samples = sample_big_groups[sample_big_groups.Big_Group == big_group]["Sample"]
                    # print(samples)
                    # select rows in sample groups table that contain those samples
                    groups = sample_groups[sample_groups.Sample.isin(samples)]["Group"].unique()
                    add_data(groups=groups, network=network, outdir=big_group_outdir, output=output, c=c)
            else:

                network = network_file
                outdir = self.outdir
                add_data(groups=groups, network=network, outdir=outdir, output=output, c=c, comparisons=comparisons)
        else:
            outdir = self.outdir
            with open(contrast_network_table) as contrast_networks:
                li = 0  # line index
                baseline_type = None  # two baseline_type: baseline_samples and baseline_table , defined in contract_network_table, last column
                for line in contrast_networks:
                    line = line.strip().split("\t")
                    if li == 0:
                        if len(line) == 4:
                            baseline_type = line[3]
                    if li > 0:  # skip heading
                        self.prt(945, line)
                        # name, contrast, network, baseline_samples = line.strip().split("\t")
                        name = line[0]
                        contrast = line[1]
                        network = line[2]
                        if len(line) == 4:
                            baseline = line[3]
                        else:
                            baseline = None
                        data = {
                            'network_file': network,
                            'gene_expression_table': gene_expression_table,
                            'sample_group_table': sample_group_table,
                            'fit_data': fit_data,
                            'outdir': outdir + "/" + name,
                            'pipe_jobId': pipe_jobId,
                            'contrast': contrast,
                            'pair_comparison_groups': name,
                            'displayed_groups': groups
                        }
                        if baseline_type is not None:
                            baseline = line[3]
                            data[baseline_type] = baseline
                        #print(1077, baseline_type, data)
                        output.add(data)
                    li += 1

        self.prt(output)
        return output

    @classmethod
    def get_inputs(cls, fit=None,
                   big_group_networks_table=None,
                   sample_big_group_table=None,
                   network_file=None,
                   contrast_network_table=None,
                   comparisons=None,
                   **kwargs):
        ob = pipe_object
        ob.big_group_network_table = big_group_networks_table
        ob.sample_big_group_table = sample_big_group_table
        ob.network_file = network_file
        ob.contrast_network_table = contrast_network_table
        ob.comparisons = comparisons
        if fit is not None:
            data = fit.dataset[0]  # only one data in fit dataset
            ob.fit_data = data.get('fit_data')
            ob.gene_expression_table = data.get('gene_expression_table')
            ob.sample_group_table = data.get('sample_group_table')
            ob.pipe_jobId = data.get('pipe_jobId')
        return ob

class consolidate_nScore_ranks(Layer):
    script = "consolidate_nScore_ranks.R"
    outfile = "nScore_genes_consolidated.csv"
    ntop = 10
    def execute(self, inputs, **kwargs):
        outfile = self.outdir + "/rank_" + self.outfile
        rank_file_pattern = inputs.rank_file_pattern
        score_file_pattern = inputs.score_file_pattern
        after = []
        if inputs.dataset:
            for data in inputs.dataset:
                try:
                    pipe_jobId = data['pipe_jobId']
                    slurm_jobId = B.cluster_jobId_dic[pipe_jobId]
                    if slurm_jobId:
                        after.append(slurm_jobId)
                except Exception as e:
                    print(e)
        cmdline = "Rscript {}/{} --file_pattern {} --outfile {} --datatype=rank".format(source_folder, self.script, rank_file_pattern,
                                                                        outfile)
        cmdline += "\n" + "Rscript {}/{} --file_pattern {} --outfile {} --datatype=rank --ntop {}".format(source_folder, self.script, rank_file_pattern,
        self.outdir +"/top" + str(self.ntop) + "gene_" + self.outfile, self.ntop)

        cmdline += "\n" + "Rscript {}/{} --file_pattern {} --outfile {} --datatype=score".format(source_folder, self.script, score_file_pattern,
                                                                        self.outdir +"/score_" + self.outfile)
        cmdline += "\n" + "Rscript {}/{} --file_pattern {} --outfile {} --datatype=score --ntop {}".format(source_folder, self.script,
                                                                                          score_file_pattern,
                                                                                          self.outdir + "/top" + str(self.ntop) + "gene_score_" + self.outfile, self.ntop)




        print("1156 nzpd", cmdline)
	
        submit(cmdline, modules=["gcc/11.3.0",   "openblas/0.3.20", "python/3.11.3", "r/4.2.3"], after=after)


class make_web_app(Layer):
    def execute(self, inputs, dataset_name="RNAseq", data_folder=".", web_outdir="../", **kwargs):
        after = []
        if inputs.dataset:
            for data in inputs.dataset:
                try:
                    pipe_jobId = data['pipe_jobId']
                    slurm_jobId = B.cluster_jobId_dic[pipe_jobId]
                    if slurm_jobId:
                        after.append(slurm_jobId)
                except Exception as e:
                    print(e)
        cmdline = "$NETZEN/make_Web_app {} {} {} ".format(dataset_name, data_folder, web_outdir)
        print(cmdline)
        submit(cmdline, modules=["gcc/11.3.0",   "openblas/0.3.20", "python/3.11.3"], after=after)

    @property
    def configs_keys(self):
        return OrderedDict([
            ("dataset_name", "name of dataset/experiment"),
            ("data_folder", "data folder, parent folder of COUNTS folder") ,
            ("web_outdir", "Folder for NETZEN web app")
        ])