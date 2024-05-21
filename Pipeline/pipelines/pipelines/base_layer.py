from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import re
import pipelines.backend as B
from collections import OrderedDict
from inspect import  getmro

source_path = B.source_path
from pipelines.utils import *


class Layer(object):
    """Abstract base layer class.
Functional style. take inputs and produce output. inputs and output are pipe object that has attributes and dataset attributes.

    """
    inputs_options = {"help": "This help"}

    cluster_config_template = "cluster_config.conf"  # template of cluster config
    run_config_template = "generic_run_config.conf"  # run config_template

    cluster_config_template_dir = source_path + "/config_templates/cluster_configs"
    run_config_template_dir = source_path + "/config_templates"
    out_file_pattern = "jobs_done/*_done"  # specific output file pattern , use to check if the pipeline is really done or not
    n_out_files = 1  # number of specific output file to count
    sleep_time = 5  # Time in second to poll the result to see if the pipeline is done or not

    # set default:
    use_cluster = False
    check_done = False
    check_ready_to_run = False
    simulated = False
    verbose = False

    modules = "python"  # module to load
    default_run_config = source_path + "/config_templates/" + run_config_template
    default_cluster_config = source_path + "/config_templates/cluster_configs/" + cluster_config_template



    outdir = "."
    execute_layer = True

    @classmethod
    def usage(cls):
        """ print class usage"""
        empty_instance = cls()
        def get_list_element_name(list, description):
            description += ":\n"
            for item in list:
                description += '\t' + item + "\n"
            return description
        properties_to_define = get_list_element_name(empty_instance.properties_to_define,"Properties to defined in a child class")

        methods_to_define = get_list_element_name(empty_instance.methods_to_define,"Methods to be defined in a child class")

        configs_keys = cls.dic_description(empty_instance.configs_keys,"config keys during initialization")
        inputs_attributes = cls.dic_description(empty_instance.inputs_attributes, "inputs attributes")
        inputs_dataset_attributes = cls.dic_description(empty_instance.inputs_dataset_attributes, "inputs dataset attributes")
        
        output_attributes = cls.dic_description(empty_instance.output_attributes, "output attributes")
        output_dataset_attributes = cls.dic_description(empty_instance.output_dataset_attributes, "output dataset attributes")

        
        usage = '''
{} class.
Description: {}
Usage:
output = {}(**config)(inputs, **kwargs)
        
**config: config keys, applied for all data in dataset (if dataset exists in inputs) 
inputs: inputs pipe object.
output: output pipe object



{}\n
'''.format(cls.__name__,cls.__doc__, cls.__name__,  configs_keys )
        usage += 'inputs:\n'
        if inputs_attributes is not None:
            usage += inputs_attributes + "\n"
        if inputs_dataset_attributes is not None:
            usage += inputs_dataset_attributes + "\n"

        usage += 'output:\n'

        if output_attributes is not None:
            usage += output_attributes + "\n"
        if output_dataset_attributes is not None:
            usage += output_dataset_attributes + "\n"

        allowed_kw_inputs = cls.dic_description(empty_instance.allowed_inputs_kwargs, "Allowed keyword arguments if calling directly without using pipe object")
        if allowed_kw_inputs is not None:
            usage += allowed_kw_inputs


        usage += "run config template to defined:\n\t " + cls.__name__ + "_run_config.conf" + "\n"
        if properties_to_define is not None:
            usage +=  properties_to_define
        if methods_to_define is not None:
            usage += methods_to_define
        print(usage)
        return usage

    @staticmethod
    def dic_description(dic, dic_name=None):
        """ 
        Make a description string for  dictionary , each item is on one line
        If use print(dic), then everything is one one line, difficult to read

        Args:
             dic : dictionary
             dic_name: name of dictionary to be displayed as a description title
        Returns:
             str: a description string for dictionary
        """
        description = ""
        if dic is None:
            return None
        if dic_name is not None:
            description = dic_name + ":\n"
        else:
            dic_name = "Dictionary"
        for key, item in dic.items():
            description += "\t{} : {}\n".format(key, item)
        return description

    @property
    def inputs_attributes(self):
        """
        attributes of inputs pipe objects

        Returns:
           OrderedDict
        """

        inputs_attributes = None
        return inputs_attributes

    @property
    def inputs_dataset_attributes(self):
        """
        attributes of data in the inputs.dataset (sample attributes)

        Returns:
            OrderedDict
        """
        inputs_dataset_attributes = None
        return inputs_dataset_attributes

    @property
    def output_attributes(self):
        """
        attributes of output pipe objects

        Returns:
           OrderedDict
        """
        output_attributes = None
        return output_attributes
        
    @property
    def output_dataset_attributes(self):
        """
        attributes of daa in output.dataset

        Returns:
           OrderedDict
        """

        output_dataset_attributes = None
        return output_dataset_attributes

    @property
    def allowed_kwargs(self):
        """
        Allowed configuration keyword args
        :return: OrderedDict
        """
        allowed_kwargs = self.base_kwargs.copy()
        if self.configs_keys is not None:
            for key, item in self.configs_keys.items():
                allowed_kwargs[key] = item
        return allowed_kwargs

    @property
    def allowed_inputs_kwargs(self):
        """
        allowed keyword arguments for use when calling object by with_inputs method
        :return: OrderedDict
        """
    @property
    def base_kwargs(self):
        """
        Base config keyword args that all children classes use
        :return: OrderedDict
        """
        base_kwargs = OrderedDict([
            ('name', 'Name of layer'),
            ('cluster_config', ' Cluster configuration file for this layer, for defining memory, number of CPU to submit to cluster'),
            ('run_config', ' user defined  configuration file for this layer, optionally including information about input data '),
            ('check_done', ' need to check if the job is finished)'),
             ('check_ready_to_run',' check if all conditions are met before running, usually a presence of a specific file'),
             ('simulated', ' dry run if True, then only write the job batch file, but not actually submit to cluster '),
             ('outdir', 'output directory'),
             ('steps', 'list of steps to run'),
             ('create_individual_data_outdir', 'create subfolder for the output result of each data in dataset'), # Create individual output folder for each data in dataset
             ('verbose', 'print more information about pipeline when running. Modify behaviour of self.prt method')
        ])
        return base_kwargs

    @property
    def configs_keys(self):
        """
        OrderedDict of configs keyword args specific for a child class
        Config arguments are the ones applied  to the whole dataset
        To be override in child class.
        """
        configs_keys = None
        return configs_keys

    @property
    def method_with_defaults(self):  # Method that setup default value, usually it's get_cmdline method and run method
        return self.__class__.run

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
        methods = ['execute',  # Main logic
                   'set_output',  # output attributes
                   'set_output_data',  # output data attributes
                   'get_inputs']   # method to generate inputs pipe object from keyword arguments
        return methods

    def run(self, **kwargs):
        pass

    @classmethod
    def check_kw_error(cls, kwargs, allowed_kwargs=None):
        """ check if the input kwargs are in allowed_kwargs. allowed_kwargs is a OrderedDict(key, description)
        """
        allowed_kwargs_dic_print = ""
        print(257, cls.__name__, cls, allowed_kwargs)
        for key, value in allowed_kwargs.items():
            allowed_kwargs_dic_print += key + ": " +  value + "\n"
        for kwarg in kwargs:
            if kwarg not in allowed_kwargs:
                raise TypeError('Keyword argument not understood for class {}:{}. \nAllowed kwargs:\n{}'.
                                format(cls.__name__, kwarg, allowed_kwargs_dic_print))

    def check_kwargs_error(self, kwargs):
        """ check if the input kwargs are in class allowed_kwargs. allowed_kwargs is a OrderedDict(key, description). Differ from check_kw_error is that this is instance method
        """
        self.__class__.check_kw_error(kwargs=kwargs, allowed_kwargs=self.allowed_kwargs)

    def __init__(self, **kwargs):
        self.inputs = None
        self.output = None
        # These properties should be set by the user via keyword arguments
        self.check_kwargs_error(kwargs)
        name = kwargs.get('name')

        if not name:
            prefix = self.__class__.__name__
            class_index = B.get_uid(prefix)
            if class_index > 1:
                name = _to_snake_case(prefix) + '@' + str(B.get_uid(prefix))
            else:
                name = _to_snake_case(prefix)
        self.name = name
        print("\n\n 86 base_layer, step_name:", self.name, "\n")
        B.step_list.append(self.name)
        #self.outdir = self.name

        # Priority in increasing order is: class, clas config, user_config, PM data, command line
        # Set cluster config and run config to default configs
        self.cluster_config = source_path + "/config_templates/cluster_configs/" + self.cluster_config_template
        self.run_config = source_path + "/config_templates/" + self.run_config_template
                # set config

        # set self.inputs ( a dictionary of arguments) from config file
        # parse default config first
        # self.config = self.__class__.default_run_config
        default_configs = self.parse_config(self.__class__.default_run_config)["run_config"]
        backend_config = B.step_configs.get(self.__class__.__name__, None)
        backend_run_configs = None
        if backend_config is not None:
            self.run_config = backend_config
            backend_run_configs = self.parse_config(backend_config)["run_config"]
            print('145 base layer, backend run config:', self.name, backend_config, backend_run_configs)

        run_config = kwargs.get("run_config", None)
        user_input_data = None
        user_run_configs = None
        if run_config is not None:
            do_parse = False
            if not os.path.isfile(run_config):
                # if not exist config file, then first check in the config library.
                run_config = source_path + "/config_templates/" + run_config
                if os.path.isfile(run_config):
                    # If run_config in config library, then do_parse
                    do_parse = True
            else:
                do_parse = True
            if do_parse:
                self.run_config = run_config
                user_configs = self.parse_config(run_config)
                user_run_configs = user_configs['run_config']
                user_input_data = user_configs.get('data', None)

        # set command line
        if 'check_done' in kwargs:
            self.check_done = kwargs.get('check_done')
        if 'outdir' in kwargs:
            self.outdir = kwargs.get('outdir')
        mkdir(self.outdir)

        if 'cluster_config' in kwargs:
            self.cluster_config = kwargs.get('cluster_config')
        if 'check_ready_to_run' in kwargs:
            self.check_ready_to_run = kwargs.get('check_ready_to_run')
        if 'verbose' in kwargs:
            self.verbose = kwargs.get('verbose')

        # Direct args inputs will overide settings in run_config file
        self.kwargs = default_configs  # got all variable from parsed default config
        if backend_run_configs is not None:
            self.kwargs.update(backend_run_configs) # update backend run configs (provided by Pipe class)
        # update step details configs:
        backend_detailed_configs = B.step_detailed_configs.get(self.__class__.__name__, None)
        if backend_detailed_configs is not None:
            self.kwargs.update(backend_detailed_configs)
        if user_run_configs is not None:
            self.kwargs.update(user_run_configs)   # update user provided run configs through config file
        if kwargs is not None:
            for key, item in kwargs.items():
                if key not in self.base_kwargs:
                    self.kwargs[key] = item       # update config through command line
        self.done_dir = "jobs_done"
        mkdir(self.done_dir)
        self.modules = list2str(self.modules)
        if len(B.step_simulation_registry) > 0:  # steps specified:
            if self.name not in B.step_simulation_registry:
                self. simulated = False
                self.execute_layer = False
            else:
                self.simulated = B.step_simulation_registry[self.name]
        if 'simulated' in kwargs:
            self.simulated = kwargs['simulated']
        self.user_input_data = user_input_data

    def __call__(self, inputs=None, **kwargs):
        """
        Main entry point for layer
        Wrapper around self.call(), for handling internal references.
        # Arguments
            inputs: Can be a tensor or list/tuple of tensors.
            **kwargs: Additional keyword arguments to be passed to `call()`.
        # Returns
            Output of the layer's `call` method.
        """
        if inputs is None:
            if kwargs is None or len(kwargs) == 0:
                return None
            else:
                out = self.with_inputs(**kwargs)
        else:
            self.inputs = inputs
            out = self.pre_execute(inputs, **kwargs)
            if out is False:
                return None
            out = self.execute(inputs, **kwargs)
            out = self.set_general_output(out, inputs)  # setting additional properties for output
            out = self.set_output(out, inputs)  # setting additional properties for output
            out = self.post_execute(inputs, out, **kwargs)
            self.output = out
        return out

    def pre_execute(self, inputs, **kwargs):
        """
        Check for keyword args error, execute layer or not (specified in steps parameters in initializations of pipeline)
        :param inputs:
        :param kwargs:
        :return: if inputs is None or execute_layer is False, return False, otherwise return True
        """
        self.check_kwargs_error(kwargs)
        if inputs is None:
            return False
        elif self.execute_layer == False:
            return False
        else:
            return True

    def execute(self, inputs, **kwargs):
        """
        Main logic of Layer
        :param inputs: pipe_object - inputs object
        :param kwargs: additional keyword arguments
        :return: output pipe object
        """
        pass

    def post_execute(self, inputs,out, **kwargs):
        """ Make report. Clean up """
        self.report(inputs, out)
        self.cleanup(inputs,out)
        return out

    def report(self, inputs, out, **kwargs):
        """ Make html and pdf report"""
        pass

    def cleanup(self, inputs, out, **kwargs):
        """ archive, remove unneccesary files"""
        pass

    @classmethod
    def from_config(cls, config):
        """Creates a layer from its config.

        This method is the reverse of `get_config`,
        capable of instantiating the same layer from the config
        dictionary. It does not handle layer connectivity
        (handled by Network), nor weights (handled by `set_weights`).

        # Arguments
            config: A Python dictionary, typically the
                output of get_config.

        # Returns
            A layer instance.
        """
        return cls(**config)

    def get_config(self):
        """Returns the config of the layer.
        A layer config is a Python dictionary (serializable)
        containing the configuration of a layer.
        The same layer can be reinstantiated later
        from this configuration.

        # Returns
            Python dictionary.
        """
        general_config = {'name': self.name,
                  'cluster_config': self.cluster_config,
                  'run_config': self.run_config,
                  'check_done': self.check_done,
                  'check_ready_to_run': self.check_ready_to_run,
                  'simulated': self.simulated,
                  'outdir': self.outdir}
        specific_config = self.kwargs
        general_config.update(specific_config)

        return general_config

    def parse_config(self, run_config):
        """
        Parsing user defined config.
        Set up kwargs for initializations: outdir, name, check done, check_ready_to_run, conf,
        :param run_config: user defined config file name
        :return: dictionary with items(run_config:configuration , data: input data)
        """

        run_config_inputs = {}
        conf = Conf(run_config)
        try:
            pipeline_specific_inputs = conf.items("Pipeline specific")
        except Exception as e:
            pipeline_specific_inputs = []

        try:
            input_data = conf.items("Input data")
        except:
            input_data = []

        for item in pipeline_specific_inputs:
            run_config_inputs[item[0]] = item[1]
        self.check_kwargs_error(run_config_inputs)
        data_inputs = {}

        for item in input_data:
            data_inputs[item[0]] = item[1]
        #self.check_kwargs_error(data_inputs)

        # print(46, self.inputs)
        self.outdir = getConf(conf, "outdir", "General", self.outdir)
        self.name = getConf(conf, "run_name", "General", self.name)
        self.verbose = getConf(conf, "verbose", "General", self.verbose)

        self.check_done = getConf(conf, "check_done", "Pipeline general", str(self.check_done))
        self.check_done = str2bool(self.check_done)
        self.check_ready_to_run = getConf(conf, "check_ready_to_run", "Pipeline general",
                                          str(self.check_ready_to_run))
        self.check_ready_to_run = str2bool(self.check_ready_to_run)

        self.conf = conf


        return {"run_config": run_config_inputs, "data": data_inputs}

    def set_general_output(self, output, inputs):

        """
        Set up common output pipe object attributes (inputs, layer_name, layer_config, output_outdir)
        :param output: output pipe object
        :param inputs: inputs pipe object
        :return:
        """
        if output is not None:
            output.inputs = inputs
            output.layer_name = self.name
            output.layer_config = self.get_config()
            output.outdir = self.outdir
        return  output

    def set_output_data(self, **kwargs):
        """set class specific kwargs to output data (sample). Usually it's output files"""
        pass

    def set_output(self, output, inputs, **kwargs):
        # setting properties for whole pipe_object , for example copying certains properties from inputs.
        return output

    def wait_for_file(self, checkfile):
        """ Waiting until checkfile appears. Use when wait for a job from previous step to finish
        Args: checkfile: name of file that need to exist before move forward
        """
        if self.check_ready_to_run:
            print("waiting for file", checkfile)
            ready = os.path.isfile(checkfile)
            while not ready:
                ready = os.path.isfile(checkfile)
                if ready:
                    break
                sys.stderr.write(".")
                time.sleep(self.sleep_time)
            print("file {} exists". format(checkfile))

    def get_defaults(self):
        """
        Get Dictionary of Default Args from method with defaults (get_cmdline).
        :return:
        """
        args = inspect.getargspec(self.method_with_defaults)

        args_default = {}
        defaults = args.defaults
        if defaults is None:
            return {}
        n_defaults = len(defaults)
        n_args = len(args.args)
        arg_names = args.args[n_args - n_defaults:]  # self arg is at index 0, skip
        for i, arg in enumerate(arg_names):
            args_default[arg] = defaults[i]
        return args_default

    @classmethod
    def get_inputs(cls, **kwargs):
        """
        Make inputs pipe object from keyword args
        :param kwargs:
        :return: pipe object
        """
        pass

    def with_inputs(self, **kwargs):
        """Make inputs pipe object from kwargs and call the __call__ method
        Use when the inputs is not a pipe object but keyword args
        """
        self.__class__.check_kw_error(kwargs, self.allowed_inputs_kwargs)
        inputs = self.__class__.get_inputs(**kwargs)
        out = self.__call__(inputs=inputs)
        #print('384 base layer out', out)
        return  out

    def autocall(self):
        """
        Auto call after intialization based on user provided config file if the config file [Input data] Section is not empty

        Usage
        :return: output pipe object
        """
        if self.user_input_data is not None and len(self.user_input_data) > 0:
            print("Running from user supplied run config {}".format(self.run_config))
            out = self.with_inputs(**self.user_input_data)
            return out

    def prt(self, *msg):
        """
        print out if verbose=True
        :param msg: message to print
        :return:
        """
        if self.verbose:
            print(*msg)

    @classmethod
    def get_config_templates(cls, outdir=".", steps="", step_config_string=""):
        """ Copy run_config and cluster_config templates into user defined directory"""
        cluster_config = cls.cluster_config_template_dir + "/" + cls.cluster_config_template
        run_config = cls.run_config_template_dir + "/" + cls.__name__ + "_run_config.conf"

        final_cluster_config = outdir + "/" + cls.__name__ + "_cluster_config.conf"
        final_run_config = outdir + "/" + cls.__name__ + "_run_config.conf"
        if not os.path.isfile(run_config):
            print("template does not exist, making new template")
            cls.make_config_templates(outdir=cls.run_config_template_dir)
        cmdline = "cp {} {} ; cp {} {}".format(cluster_config,
                                               final_cluster_config,
                                               run_config,
                                               final_run_config)
        print("401 base layer", run_config, final_run_config)

        subprocess.call(cmdline, shell=True)
        print("templates {}, {} have been copied into ".format(final_run_config, final_cluster_config) + outdir)

    @classmethod
    def make_config_templates(cls, outdir=".", steps="", step_config_string=""):
        """
        Making run config template
        :param outdir: output directory
        :param steps: steps string to specified. Used for pipe child class
        :param step_config_string: step config string to specified. used for pipe child class
        :return: write down the template into output directory with file name {cls.__name__}_run_config.conf
        """
        template = outdir + "/" + cls.__name__ + "_run_config.conf"
        empty = cls()
        input_data = ""
        if empty.allowed_inputs_kwargs:
           for key, item in empty.allowed_inputs_kwargs.items():
               input_data += """
#{}
#{} =
               """.format(item, key)

        pipeline_configs = ""
        if empty.configs_keys:
            default_values = empty.get_defaults()
            for key, item in empty.configs_keys.items():
                default_value = default_values.get(key, None)
                pipeline_configs += """
#{}
#{} = {}
""".format(item,key, default_value)

        outstring = """
[Input data]
# add input data (samples) description
%s

[General]
# General information about

#Output directory
#outdir = RESULT

#Name of the run
#run_name = run1

##################################################

[Pipeline general]
#Items that exit in every pipeline

# pipeline name
pipeline_name =   %s

#cluster_config = cluster.config.conf

%s

[Pipeline specific]
#Items - configurations parameters uniques for the pipeline/Layer
## pipeline specific intialization variables


""" % (input_data, cls.__name__, steps) + pipeline_configs + step_config_string

        with open(template, "w+") as out:
            out.write(outstring)


class dataset_layer(Layer):
    """ use this class when inputs contain datset, for example RNA seq dataset of samples and the layer does not call bash script (for example Python script). In case of bash script, then use submit_layer.
    When each data in dataset is sufficient to run the execution. The layer will loop for every data in dataset,
     No share common attributes are needed
    """
    create_individual_data_outdir = False
    data_outdir_key = None
    @property
    def methods_to_define(self):
        """
        Methods to be defined at a child class
        :return: list
        """
        methods = [
                   'set_output',
                   'set_output_data',
                   'get_inputs',
                   'call']  # Main logic for each individual data

        return methods
    @property
    def method_with_defaults(self):  # Method that setup default value, usually it's get_cmdline method and run method
        return self.__class__.call

    def get_data_outdir(self, input, index):
        """return individual output directory name for each input in inputs.dataset
        Use when class attribute     create_individual_data_outdir = True
        """
        outdir = self.outdir
        return outdir

    def execute(self, inputs, **kwargs):
        output = pipe_object()
        for i, input in enumerate(inputs.dataset):
            if self.create_individual_data_outdir:
                try:
                    data_outdir = self.get_data_outdir(input, index=i)
                    data_outdir = self.outdir + "/" + data_outdir
                    mkdir(data_outdir)
                except Exception as e:
                    print(e)
                    data_outdir = self.outdir
            else:
                data_outdir = self.outdir
            output = self.call(input=input, inputs=inputs, output=output, index=i, outdir=data_outdir, **kwargs)


        return output



    def call(self, input, inputs, output, index,outdir,  **kwargs):

        """This is where the layer's logic lives.
        call to individual input
        # Arguments
        input: individual data
        inputs: inputs pipe object, containing dataset for each data, or list/tuple.
        output: output pipe object
        index: index of input in inputs.dataset
        **kwargs: Additional keyword arguments.

        # Returns
        output pipe object with information about the result.
        """
        #output = None
        # to be implemented

        return output

class data_layer(Layer):
    """
    Use this class when inputs does not have dataset structure (for get_example comparisons),
    the kwargs is from inputs.__dic__
    if the inputs has dataset, but need to access shared inputs attribute, then use Layer class, not this class,  overide execute method in generic Layer.
    The get_inputs method of this class will make pipe object with pipe object outter attributes (not dataset attributes) defined by keyword args
"""
    def execute(self, inputs, **kwargs):
        #self.check_kwargs_error(kwargs)

        if kwargs is None:
            kwargs = {}
        kwargs.update(self.kwargs)
        inputs_dic = inputs.__dict__
        for key in self.inputs_attributes:
            if key in inputs_dic:
                input_key = inputs_dic[key]
                kwargs[key] = input_key
        output = self.run(**kwargs)

        if self.check_done:
            self.is_done()

        output = self.set_general_output(output, inputs)  # setting additional properties for output
        return output

    def run(self, **kwargs):
        # Main logic here. to be implemented at children class

        pass

    @property
    def methods_to_define(self):
        """
        Methods to be defined at a child class
        :return: list
        """
        methods = [
            'set_output',
            'set_output_data',
            'get_inputs',
            'run']  # Main logic for each individual data
        return methods

    @classmethod
    def get_inputs(cls, **kwargs):
        ob = pipe_object
        for key in cls().inputs_dataset_attributes:
            if key in kwargs:
                ob.__dict__[key] = kwargs.get(key, None)
        return ob



def _to_snake_case(name):
    intermediate = re.sub('(.)([A-Z][a-z0-9]+)', r'\1_\2', name)
    insecure = re.sub('([a-z])([A-Z])', r'\1_\2', intermediate).lower()
    # If the class is private the name starts with "_" which is not secure
    # for creating scopes. We prefix the name with "private" in this case.
    if insecure[0] != '_':
        return insecure
    return 'private' + insecure

class pipe_object(object):
    ''' This class represents object of the pipeline: inputs and outputs of each step
    It has attributes and dataset attribute. Dataset (sample collection) attribute is a list of data (sample).\n
     Each data is a dictionary of data attributes'''

    def __init__(self):
        self.dataset = []
        self.layer_name = None
        self.layer_config = None
        self.inputs = None
        self.outdir = None
    def add(self, data):
        """ add data to dataset attribute"""
        self.dataset.append(data)