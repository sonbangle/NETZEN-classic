#!/bin/sh
module purge
module load python
run_pipe -r NetZen_from_consolidated_count_run_config.conf
