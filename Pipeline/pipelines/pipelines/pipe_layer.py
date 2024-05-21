from .base_layer import Layer
from .utils import *
from . import backend as B
class Pipe(Layer):
    """
    Use this class to build pipeline. Main logic define in execute method
    """


    cluster_config_template = "cluster_config_30Gb.conf"

    @property
    def properties_to_define(self):
        """
        Properties to be defined at a child class
        Returns:
            list
        """
        properties = [
            'configs_keys',
            'inputs_attributes',
            'inputs_dataset_attributes',
            'output_attributes',
            'output_dataset_attributes',
            'allowed_inputs_kwargs']
        return properties

    @property
    def methods_to_define(self):
        """
        Methods to be defined at a child class
        :return: list
        """
        methods = ['execute',
                   'get_inputs']
        return methods


    def __init__(self, steps=None, **kwargs):
        #system_call("backup.sh $SOURCE {}".format(self.outdir))  # backup whole source folder.
        #print("41 pipe_layer, Finished backup source code")
        super(Pipe, self).__init__(**kwargs)
        self.steplist = None
        self.parse_steps(steps)
        self.parse_step_configs()
        self.autocall()

    def parse_steps(self, steps=None):
        """Reading steps config from user defined config files or from steps cmdline.
        adding step to the backend step_simulation_registry with simulated status if (-) exists
        Args: steps: list of steps to run.
        """
        if steps is None:
            self.steplist = str2list(getConf(self.conf, "steps", "Pipeline general", None))
        else:
            self.steplist = steps

        if self.steplist is not None:
            for step in self.steplist:
                if step.startswith("-") or step.startswith("no"):
                    simulated = True
                    step = step[1:]
                else:
                    simulated = False
                B.step_simulation_registry[step] = simulated


    def parse_step_configs(self):
        # getting step configs from user provided pipeline config file and put into backend dictionary
        try:
            step_configs_data = self.conf.items("Step configs")
            for item in step_configs_data:
                B.step_configs[item[0]] = item[1]
        except Exception as e:
            print('40 pipe layer, ', e)

        # parse individual step config in the main pipeline config. The step config is in each section with the name of section is the name of step
        pipeline_sections  =["Input data", "General", "Pipeline general", "Pipeline specific", "Step configs"]
        for section in self.conf.sections():
            if section not in pipeline_sections:  # section is step config:
                step_item_configs = self.conf.items(section)
                if section not in B.step_detailed_configs:  # section detailed config has not defined yet, initialization is needed
                    B.step_detailed_configs[section] = {}
                for item in step_item_configs:
                    B.step_detailed_configs[section][item[0]] = item[1]
                    self.prt('85 pipe_layer', B.step_detailed_configs)





    @classmethod
    def get_steps(cls):
        """
        Print list of steps in the pipeline with step names.
        Use this to write steps parameter in user defined config
        """
        Layer.execute_layer = False
        out = cls().execute(inputs=None)
        for step in B.step_list:
            print(step)
        Layer.execute_layer = True
        return B.step_list


    def report(self, inputs, out, **kwargs):
        steps = ""
        for step in B.step_list:
            steps +=   "<li><a href='test.html'>%s</a></li>\n" % step

        page = '''
<!DOCTYPE html>
<html>
<head>
<style>
ul {
list-style-type: none;
margin: 0;
padding: 0;
width: 400px;
background-color: #f1f1f1;
border: 1px solid #555;
}

li a {
display: block;
color: #000;
padding: 8px 16px;
text-decoration: none;
}

li {
text-align: left;
border-bottom: 1px solid #555;
}

li:last-child {
border-bottom: none;
}

li a.active {
background-color: #2f509d;
color: white;
}

li a:hover:not(.active) {
background-color: #555;
color: white;
}
h1 {
text-align: center;
background-color: #2f509d;
color: white;
}
</style>
</head>
<body>
<h1>  %s </h1>
<p>

</p>
<ul>
<li><a class="active" href="#home" STEPS </a></li>
%s
</ul>
</body>
</html>''' % (self.name.upper(), steps)

        with open("report.html", "w+") as report:
            report.write(page)

    @classmethod
    def make_config_templates(cls, outdir=".", **kwargs):
        print(173)
        steps_string = ", ".join(cls.get_steps()[1:])
        print(steps_string)
        steps = """
# Steps to run. if '-' before step name, this step will be run in simulated mode, not submitting to cluster but printout the sbatch command
# steps = %s""" % steps_string

        step_config_string = """
# Configution files for steps in pipeline. Each row begin of the step name. Value is the config file for this step
[Step configs]

# Detailed configs for each step in a pipeline
# [Layer name]
# Layer_parameter = 
"""
        super(Pipe, cls).make_config_templates(outdir=outdir, steps=steps, step_config_string=step_config_string)




