#!/usr/bin/env python
import os, sys
import argparse
import ConfigParser
import inspect
import importlib
from pipelines.utils import *
#from  pipelines_lib import *
from pipelines.pipelines_lib import *

#from  NetZen_pipelines import *


def usage():

    available_pipelines = get_class_names_from_parent_class("pipelines.pipelines_lib", "Layer")
    usage_msg ='''
    Run different pipelines on cluster. \n
    Current available pipelines: {}
        '''.format(available_pipelines)

    return(usage_msg)


# class run_pipe(object):
#
#     def __init__(self, pipe):
#         self.pipe = pipe
#
#     def __call__(self,**kwargs):


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=usage())
    parser.add_argument("-r", "--run_conf", help="submit to cluster pipeline described in run configuration file", default=None)
    parser.add_argument("-l", "--list_pipes", help="configuration file for the pipeline", action="store_true", default=False)
    parser.add_argument("-n", "--get_conf_temp", help="get configuration templates file for the pipeline", default=None)
    parser.add_argument("-s", "--get_step_conf_temp", help="get configuration templates file for a step", default=None)
    parser.add_argument("-c", "--use_cluster", help="use cluster for dispatching this pipeline or not", default="True")

    args = parser.parse_args()
    if args.run_conf is not None:
        config =args.run_conf
        Conf = ConfigParser.ConfigParser()
        Conf.read(config)
        print(config)
        pipe_name = getConf(Conf,"pipeline_name", "Pipeline general")
        pipe = locals()[pipe_name]
        use_cluster = getConf(Conf,"use_cluster", "Pipeline general", "True")
        use_cluster = str2bool(use_cluster)
        #print("37 run pipe", use_cluster)  # That is from run config
        # now priority is given to command line so no recursive submitting
        if args.use_cluster.lower() in ('yes', 'true', 't', 'y', '1') and use_cluster: # should satisfied both conditions
            use_cluster = True
        else:
            use_cluster = False
        if use_cluster:
            cluster_config = getConf(Conf, "cluster_config", "Pipeline general")
            cluster_outdir = getConf(Conf, "outdir", "General", "RESULT")
            if cluster_config is None:
                cluster_config = pipe.cluster_config_template_dir + "/" + pipe.cluster_config_template
            #cmdline = "module load python/2.7.6\n"
            #cmdline = "submit_cluster --cluster_config={} --cluster_outdir={} --modules={} run_pipe --run_conf {} --use_cluster=False --jobname={}"\
            #    .format(cluster_config, cluster_outdir, modules, args.run_conf, "pipe_" + pipe_name)

            cmdline = "submit_cluster --cluster_config={} --cluster_outdir={} --modules={} python -u $SOURCE/run_pipe --run_conf {} --use_cluster=False --jobname={}" \
                    .format(cluster_config, cluster_outdir,  "python", args.run_conf, "pipe_" + pipe_name)

            print(cmdline)
            jobId = system_call(cmdline)
            print("jobID:", jobId)
        else:

            pipe(run_config=config)
    elif args.list_pipes:
        available_pipelines = get_class_names_from_parent_class("pipelines_lib", "Pipe")
        print(available_pipelines)
    elif args.get_conf_temp:
        pipe = locals()[args.get_conf_temp]
        pipe.get_config_templates()



    else:
        parser.print_help()