#!/usr/bin/env python
import os, sys
from pipelines.utils import *
import configparser
import time
import argparse
import pickle
import pipelines.backend as B
script_dir = get_script_path()


B.use_cluster = True

def submit(cmdline,
           cluster_config=None,
           outdir="./RESULT",
           jobname=None,
           modules="Python/2.7.6",
           after=None,
           simulated=False,
           debug=True,
           wait=False,
           job_log_dir=None):
    # after: this job dispatch only after all the jobs in the comma separated list
    if cluster_config is None:
        cluster_config = script_dir + "/config_templates/cluster_configs/cluster_config.conf"

    if jobname is None:
        jobname ="generic_job" + str(int(time.time()))
    if job_log_dir is None:
        job_log_dir = B.job_log_dir
    mkdir(job_log_dir)
    jobdir = job_log_dir + "/jobs"
    logdir = job_log_dir + "/logs"
    B.cluster_jobId_dic_file = job_log_dir +"/cluster_jobId_dic.csv"
    mkdir(jobdir)
    mkdir(logdir)
    #print(12, cluster_config)
    # Create jobfile
    job = write_job(cmdline, cluster_config, jobdir=jobdir, logdir=logdir, jobname=jobname, modules=modules)
    print("cmdline:", cmdline)
    pipe_jobId, cluster_jobId = submit_one_job(job, after, simulated, wait)
    cluster_jobId = str(cluster_jobId).strip()
    if cluster_jobId is not None:
        sys.stdout.write("jobID:" + cluster_jobId+ "\n")
    if not str2bool(debug):
        os.remove(job)
    return pipe_jobId, cluster_jobId


def submit_one_job(job, after=None, simulated=False, wait=False):
    cmdstr = "unset TMPDIR; sbatch --parsable -D `pwd` "
    if wait:
        cmdstr += " --wait"
    jobIDs=""
    if after is not None:
        if type(after).__name__ == 'list':
            IDs = after
        elif type(after).__name__ == 'str':
            IDs = after.strip()
            IDs =IDs.split(",")
        else:
            IDs = [str(after)]
        if len(IDs) > 0:
            for ID in IDs:
                ID = str(ID) # make sure ID is string
                ID = ID.strip()
                ID = ":" + ID
                jobIDs += ID
            cmdstr += " --dependency=afterok{}".format(jobIDs)
    cmdstr += " " + job
    sys.stderr.write("submitting command:" + cmdstr + "\n")

    if not B.use_cluster:
        cmdstr = "chmod +x {} ; {}".format(job, job)
        print(" Not submitting to cluster, instead using desktop CPU")
    if not str2bool(simulated):
        jobID = system_call(cmdstr)
        jobID = jobID.decode()
    else:
        #print("Simulated job!!!")
        #jobID = str(B.job_counter)
        jobID = None

    B.cluster_jobId_dic[B.job_counter] = jobID
    with open(B.cluster_jobId_dic_file, "a") as f:
        f.write("{}\t{}\t{}\t{}\n".format(jobID, B.job_counter,job,jobIDs))

    return B.job_counter, jobID



def write_job(cmdline, cluster_config, jobdir, logdir, jobname="generic_submit.sh", modules=None, requeue=False,
              premodules_string="source /etc/profile.d/modules.sh"):  # modules are modules to load in ufrc
    Conf = configparser.ConfigParser()
    Conf.read(cluster_config)
    #sys.stderr.write("35 submit, cluster config: " + cluster_config + "\n")
    #sys.stderr.write("submit modules:", modules)
    B.job_counter += 1
    jobfile = jobdir + "/" + str(B.job_counter) + "_" + jobname
    with open(jobfile, "w+") as job:
        job_script = """#!/bin/sh
#SBATCH --job-name={}_job # Job name
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
""".format(jobname)

        for slurm_option in Conf.items("General"):
            job_script += "#SBATCH --{}={}\n".format(slurm_option[0], slurm_option[1])
        job_script += "#SBATCH --output={}/{}_%j.out".format(logdir, str(B.job_counter) + "_" + jobname)  # Standard output and error log
#         job.write("""#!/bin/sh
# #SBATCH --job-name={}_job # Job name
#
# #SBATCH --mail-user={} # Where to send mail
# #SBATCH --nodes=1 # Use one node
# #SBATCH --ntasks=1                   # Run a single task
# #SBATCH --cpus-per-task={}            # Number of CPU cores per task
# #SBATCH --mem={}gb # Memory limit
# #SBATCH --time={} # Time limit hrs:min:sec
        # #SBATCH --output={}/{}_%j.out # Standard output and error log
# #SBATCH --qos={}
# #SBATH  --partition={}
# pwd; hostname; echo starting:`date` \n
# module purge
# """.format(jobname,
#            getConf(Conf, "email"),
#            getConf(Conf, "cpus-per-task"),
#            getConf(Conf, "mem"),
#            getConf(Conf, "time"),
#            getConf(Conf, "partition", default=""),
#            logdir, str(B.job_counter) + "_" + jobname, getConf(Conf, "qos")))

        job_script += "\n\n"
        job_script += "echo started:`date` \n"
        job_script += "module purge \n"
        #job_script += premodules_string + "\n"
        if modules is not None:
            if type(modules).__name__ == "str":  # if modules is a string, then split into list.
                modules = modules.split(",")
            for module in modules:
                job_script += "module load " + module + "\n"
                # job.write("echo module load :" + module + "\n")
        #job.write("\necho -e submitted cmdline:{}\n\n\n".format(cmdline))
        job_script += cmdline
        if requeue:
            job_script += "|| scontrol requeue $SLURM_JOB_ID \n"
        else:
            job_script += "\n\n"
        job_script += "echo terminated:`date` \n"
        job.write(job_script)
    #sys.stderr.write("jobfile: {} written \n".format(jobfile))
    return jobfile



if __name__ == "__main__":
    cmdline = ""

    #cluster_config = script_dir + "/config_templates/cluster_configs/cluster_config.csv"  # default cluster config
    # outdir = "./SUBMIT_RESULT"
    # cmdline = ""
    # jobname ="generic_job"
    modules = None
    parser =argparse.ArgumentParser(description="This progamm submit commandline to slurm cluster",
                                    epilog='''Copyright (c) 2018, Son Le (son.le@neurosurgery.ufl.edu)\n
                                    University of Florida
                                    ''')
    parser.add_argument("--cluster_config", help="path to cluster configuration file. Example: cluster_config.csv", default=None)
    parser.add_argument("--cluster_outdir",
                        help="Folder containings the job result. If folder does not exits, a new one will be created. Default is ./RESULT", default="./RESULT")
    parser.add_argument("--jobname", help="name of job submitting", default="generic_job")
    parser.add_argument("--modules", help="comma separated list of modules to be loaded before running the script", default=None)
    parser.add_argument("--after", help="comma separated list of jobIDs that should finish with OK status before this job runs", default=None)
    parser.add_argument("--simulated", help="simluated run, not excecute script, only print out script", default="False")
    parser.add_argument("--debug", help="debuging mode, keep the submitted job file", default="True")
    parser.add_argument("--wait", help=" Do not exit until the submitted job terminates. ", default="False")
    parser.add_argument("cmdline", help="command/script to submit to cluster ")
    parser.add_argument("--job_log_dir", help="folder for job and logs", default=None)
    args = parser.parse_known_args()

    args = list(args)
    submit_args = vars(args[0])
    command_args = " ".join(args[1])

    cmdline = submit_args["cmdline"]
    if command_args is not None:
        cmdline +=  " " + command_args
    #print(cmdline)

    #if str(submit_args['cluster_config']) != "None":
    cluster_config = submit_args['cluster_config']

    submit(cmdline, cluster_config=cluster_config,
           outdir=submit_args['cluster_outdir'],
           jobname=submit_args['jobname'],
           modules=submit_args['modules'],
           after=submit_args['after'],
           simulated=submit_args['simulated'],
           wait=str2bool(submit_args['wait']),
           job_log_dir=submit_args['job_log_dir'])
