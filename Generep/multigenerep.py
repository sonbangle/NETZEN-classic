#!/usr/bin/env python

import sys
import os.path
import time
import subprocess as sp

# multigenerep.py gbm.conf 1 10 100

def write_conf(infile, outfile, pfactor, pfiltration, origaracne, maxp):
    with open(outfile, "w") as out:
        with open(infile, "r") as f:
            for line in f:
                if line[0] == '[':
                    out.write(line)
                    continue
                key = line.split(" ")[0]
                if key == "pfactor":
                    out.write("pfactor = {}\n".format(pfactor))
                elif key == "steps":
                    if origaracne:
                        if maxp:
                            out.write("steps = init, randomize, bootstrap, aracne, aracnecons, consensus, filter1, filter2, filter3, filter4, cx, aracnepropermi\n")
                        else:
                            out.write("steps = init, randomize, bootstrap, aracne, aracnecons, consensus, filter1, filter2, filter3, filter4, cx\n")
                    else:
                        out.write("steps = init, randomize, bootstrap, aracne, consensus, filter1, filter2, filter3, filter4, cx\n")
                    out.write("pfiltration = {}\n".format(pfiltration))
                elif key == "pfiltration":
                    pass
                elif key == "title":
                    out.write("title = genereprun_pf{}\n".format(pfactor))
                else:
                    out.write(line)

def get_title(conffile):
    with open(conffile, "r") as f:
        for line in f:
            if line[0] == '[':
                continue
            parts = line.split("=")
            key = parts[0].strip()
            if key == "title":
                return parts[1].strip()

def waitFor(f, timestamp, wait=60):
    sys.stderr.write("Waiting for {}\n".format(f))
    while True:
        if os.path.isfile(f) and os.path.getmtime(f) > timestamp:
            return True
        time.sleep(wait)

class MGR(object):
    infile = ""
    pfactors = []
    mode = 3
    confs = []

    def __init__(self, args):
        self.pfactors = []
        self.confs = []
        for a in args:
            if a == "-1":
                self.mode = 1   # write conf files and run generep
            elif a == "-2":
                self.mode = 2   # write conf files and run aracne
            elif a == "-3":
                self.mode = 3   # do all
            elif a == "-4":
                self.mode = 4   # collect edges only
            elif a == "-0":
                self.mode = 0   # only write confs file
            elif self.infile:
                self.pfactors.append(int(a))
            else:
                self.infile = a

    def generate_confs(self):
        (base, ext) = os.path.splitext(self.infile)
        maxp = max(self.pfactors)
        for p in self.pfactors:
            outfile = "{}-pf{}{}".format(base, p, ext)
            self.confs.append(outfile)
            sys.stderr.write("Writing " + outfile + "\n")
            write_conf(self.infile, outfile, p, True, True, p==maxp)
        return True

    def run(self):
        conf0 = self.infile
        now = time.time()
        base0 = os.path.splitext(conf0)[0]
        dir0 = get_title(conf0)
        target = dir0 + "/ARACNE.done"
        cmdline = "/blue/dtran/share/alberto/TCGA/bin/subone.sh {}".format(base0)
        donefiles = []
        if self.mode & 1:
            sys.stderr.write("Executing " + cmdline + "\n")
            sp.call(cmdline, shell=True)
            donefiles.append(target)
        if self.mode & 2:
            for conf in self.confs:
                base = os.path.splitext(conf)[0]
                donefiles.append(base + ".done")
                cmdline = "/blue/dtran/share/alberto/TCGA/bin/subone.sh {}".format(base)
                sys.stderr.write("Executing " + cmdline + "\n")
                sp.call(cmdline, shell=True)
        for d in donefiles:
            waitFor(d, now)
        if self.mode:
            sp.call("/blue/dtran/share/alberto/TCGA/bin/collect_nedges.py {} {} > NEDGES".format(conf0, " ".join(self.confs)), shell=True)

def usage():
    sys.stdout.write("""Usage: multigenerep.py conffile pfactor...

Example: multigenerep.py gbm.conf 1 10 100

""")

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) == 0 or "-h" in args:
        usage()
    else:
        m = MGR(sys.argv[1:])
        m.generate_confs()
        if m.mode > 0:
            m.run()
