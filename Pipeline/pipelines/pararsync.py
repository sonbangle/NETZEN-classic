import os
import sys

indir =""
outdir = ""
# Copy the whole indir into the outdir
ncpus = 100 # Number of cpu
nprocs = 5 # number of process per cpu
#crawl each directory and rsync files only, if find directory: crawling further.



cmd = "ls /srv/mail | xargs -n1 -P" + nprocs + " -I% rsync -azhuP % myserver.com:/srv/mail/"
