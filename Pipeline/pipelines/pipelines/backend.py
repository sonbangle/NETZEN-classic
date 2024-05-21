from pipelines.utils import get_script_path
"""
This module keeps the inventory of all layers, dictionary of internal jobID - cluster_jobID
"""

uid = {}
source_path = get_script_path()
job_counter = 0
cluster_jobId_dic = {}  # A dictionary pipe_jobId - slurm_jobId
use_cluster = True
job_log_dir = "cluster_logs"
cluster_jobId_dic_file = job_log_dir + "/cluster_jobId_dic.csv"
config_registry = {}
''' registry of all Layer class  configs. Config priority with increasing priority:
 class default,
 default config,
 config_registry,
 user config,
 cmdline

'''
step_simulation_registry = {}  # Register if step simulation is True . Highest priority
step_list = [] # list of step names
step_configs = {}  #dictionary of step configs file for each step.{step: config}. Use to initialize steps. Define in Pipe class
step_detailed_configs = {}  # dictionary of {step:{item, value}}. Differ from step_configs in that here the config for each item is define.
#  This configs have higher priority than step_configs


def get_uid(prefix=''):
    """Associates a string prefix with an integer counter in pipeline
  Arguments:
    prefix: String prefix to index.
  Returns:
    Unique integer ID.
  Example:
  ```
    >>> get_uid('save_network_image')
    1
    >>> get_uid('save_network_image')
    2
  ```
    """

    layer_uid = uid.get(prefix, 0) + 1
    uid[prefix] = layer_uid

    return layer_uid