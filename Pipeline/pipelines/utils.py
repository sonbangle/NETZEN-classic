import os, sys
import configparser
import inspect
import argparse
import subprocess
import time
import pandas as pd
import glob
def mkdir(dir):
    try:
        if not os.path.exists(dir):
            os.makedirs(dir)
    except:
        pass

def get_script_path():
    path = os.path.dirname(os.path.realpath(__file__))
    return path


# Get value for line in configuration file in Section, Entry
def getConf(Conf, entry, section="General", default=None):
    """Return the value for `entry' in `section', if present, or the
    value of default (None if unspecified). `section' defaults to General."""
    try:
        if Conf.has_section(section):
            return Conf.get(section, entry)
        else:
            return default
    except configparser.NoOptionError:
        return default


# Return Conf instance
def Conf(config):
    Conf = configparser.ConfigParser()
    Conf.read(config)
    return Conf

# find all subclasses that belong to this parent class
def inheritors(klass):
    subclasses = set()
    work = [klass]
    while work:
        parent = work.pop()
        for child in parent.__subclasses__():
            if child not in subclasses:
                subclasses.add(child)
                work.append(child)
    return subclasses

# get names of all classess in a module
def get_class_names(module):
    current_module = sys.modules[module]
    class_names = []
    for name, obj in inspect.getmembers(current_module):
        if inspect.isclass(obj):
            class_names.append(obj.__name__)
            print(obj.__bases__)
    return class_names

# get names of all children classess  that belong to a parent class in a module. Usefull for list children class in a library
def get_class_names_from_parent_class(module=None, class_name=None):
    current_module = sys.modules[module]
    class_names = []
    for name, obj in inspect.getmembers(current_module):
        if inspect.isclass(obj):
            parents = inspect.getmro(obj)
            for parent in parents:
                if parent.__name__ == class_name:
                    class_names.append(obj.__name__)
    return class_names

def str2bool(v):
    if type(v).__name__ == 'str':
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')
    else:
        return v


def system_call(cmdline, sleep_time=600):
    while True:
        try:
            print("Calling subprocess for this:", cmdline, flush=True)
            output = subprocess.check_output(cmdline, shell=True, stderr=subprocess.STDOUT)
            print("Done calling subprocess", output, flush=True)
            break
        except  subprocess.CalledProcessError as e:
            print("System call failed:" + cmdline, flush=True)
            error_output = e.output.decode("utf-8")
            print(error_output, flush=True)
            print("got error when submitting, most likely too many submission, resubmitting after " + str(sleep_time) + "s")
            #time.sleep(sleep_time)
            break

    return output

def get_localtime(): # return formatted local time
    localtime = time.asctime(time.localtime(time.time()))
    return localtime


#convert a string of  comma separated list of intergers in config file into list of integers.
# for example a string  '1,2' to [1,2]
def str2int_list(str_list):
    if type(str_list).__name__ == 'str':
        items = str_list.split(",")
        outlist =[]
        for item in items:
            item = int(item)
            outlist.append(item)
        return outlist
    else:
        return


def str2list(str_list):
    if type(str_list).__name__ == 'str':
        mylist = str_list.split(",")
        mylist = [item.strip() for item in mylist]

        return mylist
    else:
        return str_list

# convert a list to a comma separated string. for example [1,2] to "1,2"
def list2str(v):
    if type(v).__name__ == 'list':
        out =[]
        for item in v:
            item = str(item)
            out.append(item)
        out = ",".join(out)
        return out
    else:
        return v

# return unique groups from a column in a table
def get_unique_group(sample_group_table, column="Group"):
        sample_groups = pd.read_csv(sample_group_table, sep="\t")
        groups = sample_groups[column].unique()
        groups = groups.astype(str).tolist()
        return groups

def submit_jobs_by_pattern(job_folder="cluster_log/jobs",
               job_pattern ="heatmap_*.sh",
               submit_command ="sbatch"):
    '''
submit jobs from job folder based on job pattern.
Used to submit jobs that generated but have not been submitted (for example from R script
:param job_folder: folder containing job scripts
:param job_pattern: glob pattern to search for job file.
:return: None, just submit to cluster
'''
    jobs = glob.glob(job_folder + "/" + job_pattern)
    for job in jobs:
        os.system(submit_command + " " + job)
