from pipelines.base_layer import *
from pipelines.submit_cluster import submit
import subprocess
import os, sys
from pipelines.backend import cluster_jobId_dic
import inspect
from pipelines.utils import str2list
from abc import abstractmethod

class submit_layer(dataset_layer):
    '''
    This class for layer that need to submit to cluster. It will generate batch file
    and call submit_cluster.submit to submit to slurm cluster.
    For each children class, need to implement:
    get_cmdline
    set_output_data
    inputs keys
    get_data_outdir if create_individual_outdir = True


    '''
    check_done = False
    # to be implemented in each inherited class

    modules = ["python", "R"]  # module to load

    @property
    def method_with_defaults(self):  # Method that setup default value, usually it's get_cmdline method and run method
        return self.__class__.get_cmdline

    @property
    def methods_to_define(self):
        """
        Methods to be defined at a child class
        :return: list
        """
        methods = [
                   'set_output',  # optional
                   'set_output_data',
                   'get_inputs',  # optional
        'get_cmdline',
        'get_data_outdir']  # optional, when create_individual_data_outdir = True
        return methods

    def call(self, input, inputs, output, index, outdir, **kwargs):
        if kwargs is None:
            kwargs = {}
        kwargs.update(self.kwargs)
        for key in self.inputs_dataset_attributes.keys():
            if key in input:
                input_key = input[key]
                kwargs[key] = input_key

        input_jobId = input.get('pipe_jobId', None)
        after = None
        if input_jobId is not None:
            after = input_jobId
        cmdline = self.get_cmdline(outdir=outdir, **kwargs)
        pipe_jobId, slurm_jobId = self.run(cmdline=cmdline, index= index, after=after)
        if self.check_done:
            self.is_done()
        data = self.set_output_data(inputs=inputs, index=index, pipe_jobId=pipe_jobId, outdir=outdir)

        data['data_outdir'] = outdir
        data['kwargs'] = kwargs
        data['pipe_jobId'] = pipe_jobId
        data['index'] = index
        output.add(data)
        output.cmdline_defaults = self.get_defaults()

        return output

    def set_output_data(self, inputs,pipe_jobId=None, index= None, outdir=None, **kwargs):
        """
        Set attributs for each data in output.dataset
        :param inputs: inputs object
        :param pipe_jobId: pipe_jobID for this data piece
        :param index:  index of this output data in output.dataset
        :param outdir: output directory for this data
        :param kwargs: Addtional kwargs
        :return: data : dictionary of attributes of data
        """
        data = {}
        # To implement in childeren class
        return data

    @abstractmethod
    def get_cmdline(self,**kwargs):
        """
        Get cmdline for (third party) bash program to submit to cluster
        To be implemented for each individual step.
        Usually it's a name of result folder
        """
        cmdline = ""  # to be implemented based on the input variables
        return cmdline

    def run(self, cmdline, index=0, after=None, wait=False):

        # After is list of jobIds that should have finished before this job to be submitted. wait is a slurm command sbatch command
        # Wait: Do not exit from sbatch  until the submitted job terminates.
        #if after is None:
        #    after = self.after

        # Translate pipe_jobId into slurm_jobId

        if after is not None:
            after_slurm = []
            if type(after).__name__ != 'list':
                after = [after]
            for jobId in after:
                slurm_jobId = cluster_jobId_dic.get(jobId, None)
                if slurm_jobId is not None:
                    after_slurm.append(slurm_jobId)
            after = after_slurm

        pipe_jobId = None
        slurm_jobId = None
        sys.stdout.flush()
        pre_run_cmd = ""
        cmdline = cmdline + self.done(index)  # add touch done every time finished job in cluster.

        jobname = self.name + "_" + str(index) +".sh"
        if self.check_ready_to_run:
            checkfile = self.get_ready_to_run_check_file()
            pre_run_cmd = '''
while [ ! -f  {} ]; do\n
echo .
sleep {}\n
done\n
echo {} exists, beginning job
'''.format(checkfile, self.sleep_time, checkfile)
            final_cmdline = pre_run_cmd + cmdline

            if self.cluster_config is not None:  # real job dispatch from watcher in case that cluster_config is defined (big resource is going to consumed)
                run_string = "submit_cluster "
                if after is not None:
                    run_string += " --after=" + after
                if wait:
                    run_string += " --wait=True"
                run_string += " --cluster_config={} --cluster_outdir={} --jobname={} --modules={} {}". \
                    format(self.cluster_config, self.outdir, self.name + ".sh", self.modules, cmdline)
                final_cmdline = pre_run_cmd + " " + run_string
                pipe_jobId, slurm_jobId = submit(final_cmdline, jobname=self.name + "_dispatcher.sh", simulated=self.simulated)
            else:  # watcher run the job by itself
                pipe_jobId, slurm_jobId = submit(final_cmdline,  outdir=self.outdir,
                                      jobname=jobname, modules=self.modules, after=after, wait=wait, simulated= self.simulated)
        else:
            pipe_jobId, slurm_jobId = submit(cmdline, cluster_config=self.cluster_config, outdir=self.outdir,
                                      jobname=jobname , modules=self.modules, after=after, wait=wait, simulated=self.simulated)

        sys.stdout.flush()
        # print("return value:", returned_valu

        return pipe_jobId, slurm_jobId

    def done(self, index):
        """Return empry done file signalling that the step is finished
        :param index: index of data in inputs dataset
        :return: done string for including in submit job.
        """
        done_string = "\ntouch " + self.done_dir + "/" + self.name + "_" + str(index) + ".done"
        return done_string

    def get_ready_to_run_check_file(self):
        return None

    def is_done(self):
        pass

