from utils import *


#def system_call(cmdline):
#    print("Calling subprocess for cmdline:", cmdline, flush=True)
#    output = subprocess.check_output(cmdline, shell=True, stderr=subprocess.STDOUT)
#   print("Done calling subprocess", output, flush=True)


#cmdline  ="ls -lht"
#system_call(cmdline)
cmdline ="submit_cluster ls"
system_call(cmdline)