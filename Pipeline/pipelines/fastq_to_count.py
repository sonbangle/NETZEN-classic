#!/usr/bin/env python
import os
import os.path
import sys
import glob
import ConfigParser
import subprocess
import shlex
import re  # regular expression


# sys.path.insert(0, '$MY_BIN')
# to do: read cluster config file



#### use this pipeline to create jobs for processing RNAseq raw fastq or SRA data. It will call create jobs for
### sra_to_count_star_pipeline.sh


# Main class

class Pipeline():
    indir = ""
    outdir = "OUT"
    ncpu = 24
    mem = 72
    configFile = None
    Conf = None
    pipe = None
    MY_BIN = os.environ['MY_BIN']
    cluster_conf = MY_BIN + "/cluster_conf.csv"
    print(cluster_conf)
    email = "son.le@neurosurgery.ufl.edu"
    runtime = "96:00:00"
    qos = "dtran-b"
    jobdir = None
    logdir = None

    # Internal methods (not meant to be called by user)

    def __init__(self):
        self.configFile = "sample_config.csv"

    # Support for configuration files

    def loadConfiguration(self, filename):
        print("Load configuration file `filename'")
        if os.path.isfile(filename):
            config_dir = os.path.dirname(os.path.abspath(filename))
            # print "config dir: {}".format(config_dir)
            # print filename
            self.configFile = filename
            self.Conf = ConfigParser.ConfigParser()
            self.Conf.read(filename)
            self.pipe = self.getConf("pipe")
            # print self.pipe
            # print self.getConf("indir")
            self.indir = os.path.abspath(self.getConf("indir"))
            # print self.indir
            self.outdir = os.path.abspath(self.getConf("outdir"))
            self.ncpu = self.getConf("nCPU_per_task")
            self.mem = int(self.ncpu) * 3
            print(self.outdir)
            self.jobdir = "{}/jobs".format(self.outdir)
            self.logdir = "{}/logs".format(self.outdir)
            self.job_pattern = self.getConf("job_pattern")
            self.nright_trim = int(self.getConf("nright_trim"))
            os.chdir(config_dir)
            print(69, self.job_pattern)
            return self.Conf
        else:
            print "Error: configuration file {} not found or not readable.".format(filename)
            exit(1)

    def getConf(self, entry, section="General", default=None):
        """Return the value for `entry' in `section', if present, or the
        value of default (None if unspecified). `section' defaults to General."""
        print section
        try:
            if self.Conf.has_section(section):
                return self.Conf.get(section, entry)
            else:
                return default
        except ConfigParser.NoOptionError:
            return default

    def getConfBoolean(self, entry, section="General", default=None):
        """Return the boolean value for `entry' in `section', if present, or the
value of default (None if unspecified). `section' defaults to General."""
        try:
            return self.Conf.getboolean(section, entry)
        except ConfigParser.NoOptionError:
            return default

    def getConfList(self, entry, section="General", default=[]):
        """Return the value for `entry' in `section' as a comma-delimited list.
        section' defaults to General."""
        try:
            a = self.Conf.get(section, entry)
            return [w.strip(" ") for w in a.split(",")]
        except ConfigParser.NoOptionError:
            return default

            # def createjobs(self, job_pattern=self.job_pattern, nright_trim=self.nright_trim):

    # def createjobs(self, job_pattern=self.job_pattern, nright_trim=self.nright_trim):
    #     print "Making  outdir ..."
    def createjobs(self):
        job_pattern = self.job_pattern
        nright_trim = self.nright_trim

        print self.outdir
        mkdir(self.outdir)
        mkdir(self.jobdir)
        mkdir(self.logdir)

        files = os.listdir(self.indir)
        ext = self.getConfList("ext")
        print ext
        ncpu = self.ncpu
        for f in files:
            if os.path.isdir(
                    f):  # for cases of pair_end fastq where each sample has its own directory containing two fastq files
                print(f)
                fastq_files = os.listdir(f)
                print(124, fastq_files)
                fastq_files.sort()
                base = fastq_files[0]
                print(base)
                f = f + "/" + base
            else:
                base = f
            f_ext = os.path.splitext(f)[1][1:]
            if f_ext in ext:
                # print f
                jobname = "job_" + base[0:re.search(job_pattern, base).end() - nright_trim] + ".sh"
                print(133, jobname)
                self.write_job(ncpu, f, jobname)

    def write_job(self, ncpu, f, jobname):
        print("write job for :" + f)
        jobfile = os.path.join(self.jobdir, jobname)
        mem = self.mem
        with open(jobfile, "w+") as job:
            job.write("""#!/bin/sh
#SBATCH --job-name={}_job # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user={} # Where to send mail
#SBATCH --nodes=1 # Use one node
#SBATCH --ntasks=1                   # Run a single task	
#SBATCH --cpus-per-task={}            # Number of CPU cores per task
#SBATCH --mem={}gb # Memory limit
#SBATCH --time={} # Time limit hrs:min:sec
#SBATCH --output={}/logs/{}_%j.out # Standard output and error log
#SBATCH --qos={}    
pwd; hostname; date
{} {}/{} {} {}
date
""".format(self.pipe, self.email, ncpu, mem, self.runtime, self.outdir, jobname, self.qos, self.pipe, self.indir, f,
           self.outdir, ncpu)
                      )
        print "jobfile {} written".format(jobfile)

    def submit_one_job(self, job):
        cmdstr = "sbatch " + job
        print("cmdstr :" + cmdstr)
        os.system(cmdstr)

    def submitjobs(self):

        # check for memory exceed and resubmit failed

        def checkLogs_resubmit(self, error_pattern="memory", ncpu=30):
            cmdstr = "grep -l {}  {}/logs/*  ".format(error_pattern, self.outdir)
            failed_jobs = os.popen(cmdstr).read()
            # print failed_jobs
            failed_job_list = "{}/failed_jobs_list.csv".format(self.logdir)
            with open(failed_job_list, "w+") as fl:
                fl.writelines(failed_jobs)
            failed_jobs = failed_jobs.splitlines()
            for job in failed_jobs:
                prefix = "{}/logs/".format(self.outdir)
                # print prefix
                job = job[len(prefix):]
                underscore_p = job.find("_")
                job = job[:underscore_p]
                print job
                self.write_job(ncpu, job)
                jobfile = os.path.join(self.jobdir, "job_{}_{}.sh".format(self.pipe, job))
                self.submit_one_job(jobfile)

        jobs = os.listdir(self.jobdir)
        i = 0
        for job in jobs:
            job = os.path.join(self.jobdir, job)
            self.submit_one_job(job)
            i = i + 1
        print "Total {} jobs submitted".format(i)


        # Create directory for fastq pair_end samples if all pair_end samples are all in one directory


def fastq_pair_end_make_directory(indir=".", glob_file_pattern="*.fastq", reg_sample_pattern="_EGAR\w*_\d*_1_"):  #
    # get file list:
    dirs = []
    print(indir)
    os.chdir(indir)
    files = glob.glob(glob_file_pattern)
    for f in files:
        print(f)
        if os.path.isfile(f):
            print("Moving file:{}".format(f))
            # extract directory name, if not exist, create new directory
            reg = re.match(reg_sample_pattern, f)
            # print(203, reg,f)
            subdir = f[0:reg.end() - 7]
            print(204, subdir)
            mkdir(subdir)
            os.rename(f, subdir + "/" + f)
    print("Done with moving files")


def usage():
    print """
usage: pipeline.py [command] [configfilename]
[command] :run :run the pipeline
[command] check : check and rerun in case of error
for example:
pipeline.py run SCLC_config.csv    

  Executes pipeline defined in configfilename 

Copyright (c) 2016, Son Le (son.le@neurosurgery.ufl.edu)
University of Florida
"""


def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)


if __name__ == "__main__":
    p = Pipeline()
    if len(sys.argv) == 1:
        usage()
        quit()
    p.loadConfiguration(sys.argv[2])
    if sys.argv[1] == "run":
        print "Creating job..."
        p.createjobs()
        p.submitjobs()
    if sys.argv[1] == "check":  # check and rerun
        p.checkLogs_resubmit()