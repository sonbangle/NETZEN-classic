#!/usr/bin/env python

#####################################################
#                                                   #
# ndextools.py - Interface to the NDEx server       #
#                                                   #
# (c) 2017, Alberto Riva (ariva@ufl.edu)            #
#           Son Le (nuibang@gmail.com)              #
#           University of Florida                   #
#####################################################

import sys
import time
import os.path

import ndex.client
from requests.exceptions import HTTPError

CONFFILE = "~/.ndexrc"

class NdexClient():
    ndex = None
    host = "http://test.ndexbio.org"
    username = ""
    password = ""

    def readConfiguration(self):
        with open(os.path.expanduser(CONFFILE), "r") as f:
            self.host = f.readline().strip()
            self.username = f.readline().strip()
            self.password = f.readline().strip()

    def __enter__(self):
        self.readConfiguration()
        self.ndex = ndex.client.Ndex(host=self.host, username=self.username, password=self.password)
        return self.ndex

    def __exit__(self, a1, a2, a3):
        pass

def searchNetwork(term, account=None):
    with NdexClient() as n:
        nets = n.search_networks(search_string=term, account_name=account)
        return nets

def uploadNetwork(cxfile):
    with NdexClient() as n:
        with open(cxfile, "r") as f:
            url = n.save_cx_stream_as_new_network(f)
        sp = url.rfind("/")
        netid = url[sp+1:]
    sys.stderr.write("Network {} uploaded, ID={}\n".format(cxfile, netid))
    return netid

def verifyNetwork(n, netid, delay=5, retries=5):
    sys.stderr.write("Verifying network {}\n".format(netid))
    while True:
        try:
            summary = n.get_network_summary(netid)
        except HTTPError as e:
            return (False, e)
        valid = summary['isValid']
        if valid:
            return (True, netid)
        elif summary['errorMessage']:
            return (False, summary['errorMessage'])
        else:
            retries -= 1
            sys.stderr.write("Retrying {}\n".format(retries))
            if retries == 0:
                return (False, False)
        time.sleep(delay)

def verifyNetworks(netids, delay=5, retries=5):
    with NdexClient() as n:
        for netid in netids:
            (flag, msg) = verifyNetwork(n, netid, delay=delay, retries=retries)
            if flag:
                sys.stdout.write("{}\tuploaded\tOk\n".format(netid))
            elif msg:
                sys.stdout.write("{}\terror\t{}\n".format(netid, msg))
            else:
                sys.stdout.write("{}\terror\tunknown - retries exhausted\n".format(netid))

def setNetworkSample(netid, cxfile):
    with open(cxfile, "r") as f:
        cxstring = f.read()

    with NdexClient() as n:
        n.set_network_sample(netid, cxstring)
    return True # n.get_sample_network(netid)

def displaySummary(summ):
    for field in ['externalId', 'owner', 'name', 'description', 'nodeCount', 'edgeCount']:
        sys.stdout.write("{}: {}\n".format(field, summ[field]))
    props = summ['properties']
    sys.stdout.write("Properties:\n")
    for prop in props:
        sys.stdout.write("  {predicateString}: {value}\n".format(**prop))

def networkSummary(netid):
    with NdexClient() as n:
        try:
            result = n.get_network_summary(netid)
        except HTTPError as e:
            if e.response.status_code == 404:
                sys.stdout.write("Error: network {} does not exist.\n".format(netid))                
            else:
                sys.stdout.write("Error: {}\n".format(e))
            return None
        displaySummary(result)
        return result

def deleteNetwork(n, netid):
    try:
        n.delete_network(netid)
        return (True, netid)
    except HTTPError as e:
        return (False, e)

def deleteNetworks(netids):
    with NdexClient() as n:
        for netid in netids:
            (flag, msg) = deleteNetwork(n, netid)
            if flag:
                sys.stdout.write("{}\tdeleted\tOk\n".format(netid))
            else:
                sys.stdout.write("{}\terror\t{}\n".format(netid, msg))

### Top level

def usage(what=None):
    sys.stdout.write("""ndextools - Command-line interface to the ndexbio website

Usage: ndextools.py command [arguments]

where command is one of: upload, sample, verify, search, summary, delete.

""")

if __name__ == "__main__":
    args = sys.argv[1:]
    if not args or "-h" in args:
        usage()
    else:
        if args[0] == "upload":
            netids = []
            for cxfile in args[1:]:
                netids.append(uploadNetwork(cxfile))
            verifyNetworks(netids)
        elif args[0] == "verify":
            verifyNetworks(args[1:])
        elif args[0] == "search":
            term = args[1]
            if len(args) > 2:
                acc = args[2]
            else:
                acc = None
            result = searchNetwork(term, account=acc)
            nets = result['networks']
            for net in nets:
                sys.stdout.write("{}\t{}\n".format(net['externalId'], net['name']))
        elif args[0] == "summary":
            networkSummary(args[1])
        elif args[0] == "delete":
            deleteNetworks(args[1:])
        elif args[0] == "sample":
            netid = args[1]
            cxfile = args[2]
            print setNetworkSample(netid, cxfile)
        else:
            usage()
