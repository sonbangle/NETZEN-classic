Step to run a pipeline

1) Get template by run_pipe -n {pipe_name}
This will copies {pipe_name}_run_config.csv and {pipe_name}_cluster_config.csv in the current folder.

2) Change the run_config.csv file as appropriate: The file contains following sections:

[Input data]
Contain information about samples

[General]
Information about outdirectory, name

[Pipeline general]
Items that exit in every pipeline

[Pipeline specific]
Items - configurations parameters uniques for the pipeline/Layer


[Step configs]
Configution files for steps in pipeline. Each row begin of the step name. Value is the config file for this step

Then the sections containing specific detailed configuration for each step in pipeline.
For example:
[get_fit_data]
glm_robust=True.



Priority of configuration in increasing priorities:
Class default
Configs file defined in pipeline config file
Detailed config for each step defined in pipeline config file
User_defined configs file
cmdline


