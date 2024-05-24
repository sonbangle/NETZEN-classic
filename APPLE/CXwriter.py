#!/usr/bin/env python

import sys

# TODO: 
# Add visualProperties section
#  make node size proportional to degree
# Add metadata for edge attributes
# Add clustering (igraph)

class AttrReader():
    attnames = []
    attrs = {}

    def __init__(self, filename):
        print "reading " + filename
        with open(filename, "r") as f:
            self.attnames = f.readline().rstrip("\r\n").split("\t")[1:]
            self.attnames = [ "<b>Approved name</b>" if x == "Approved name" else x for x in self.attnames]
            for line in f:
                fields = line.rstrip("\r\n").split("\t")
                self.attrs[fields[0]] = fields[1:]

    def attributesFor(self, name):
        if name in self.attrs:
            attrlist = self.attrs[name]
            return [ (self.attnames[i], attrlist[i]) for i in range(len(self.attnames)) ]
        else:
            return None

class CXwriter():
    adjfile = None
    nodesdb = {}
    nedges = {}                 # Number of edges for each node (degree)
    nodecnt = 0                 # Counter for @id fields in nodes
    edgecnt = 0                 # Counter for @id fields in edges
    visibleNodes = {}           # Used for filtering
    mindegree = None            # If set, show only nodes with this degree or higher
    maxnodes = None             # If set, show this number of genes in order of decreasing degree
    nodesWritten = 0            # Number of nodes actually written
    edgesWritten = 0            # Number of edges actually written
    nEdgeAttributes = 1         # Number of attributes for each edge
    maxNodeSize = 100.0
    minNodeSize = 10.0
    minFontSize = 12.0
    maxFontSize = 24.0
    minEdgeSize = 150.0
    maxEdgeSize = 255.0
    mis    = []                 # Data associated with each edge
    FPRcol = False              # Column containing FPR in mis

    def __init__(self, adjfile):
        self.adjfile = adjfile
        self.nodesdb = {}
        self.nedges = {}
        self.nodecnt = 0
        self.edgecnt = 0
        self.visibleNodes = {}
        self.nodesWritten = 0
        self.edgesWritten = 0
        self.mis = []

    def writePreamble(self, stream):
        stream.write("""
[ {
  "numberVerification" : [ {
    "longNumber" : 281474976710655
  } ]
}, """)

    def writeNetworkAttributes(self, stream, fields):
        comma = False
        stream.write("""{
  "networkAttributes" : [ """)
        for name, value in fields.iteritems():
            if comma:
                stream.write(", ")
            stream.write('{\n    "n" : "' + name + '",\n    "v" : "' + value + '"\n  }')
            comma = True
        stream.write('  ]\n}')

    def writeMetadata(self, stream):
        stream.write(""", {
  "metaData" : [ {
    "consistencyGroup" : 1,
    "elementCount" : """ + str(self.nodesWritten) + """,
    "idCounter" : 1,
    "name" : "nodes",
    "properties" : [ ],
    "version" : "1.0"
  }, {
    "consistencyGroup" : 1,
    "elementCount" : """ + str(self.edgesWritten) + """,
    "idCounter" : 1,
    "name" : "edges",
    "properties" : [ ],
    "version" : "1.0"
  }, {
    "consistencyGroup" : 1,
    "elementCount" : """ + str(self.edgesWritten * self.nEdgeAttributes) + """,
    "idCounter" : 1,
    "name" : "edgeAttributes",
    "properties" : [ ],
    "version" : "1.0"
  } ]
}""")

    def addNode(self, name):
        """Adds a node to the network. Updates the node count."""
        if name not in self.nodesdb:
            self.nodecnt += 1
            self.nodesdb[name] = self.nodecnt
            self.visibleNodes[name] = True         # All nodes are initially visible
        return self.nodesdb[name]
    
    def countEdge(self, name):
        """Increment the edge counter for node `name'."""
        if name in self.nedges:
            self.nedges[name] += 1
        else:
            self.nedges[name] = 1

    def getMaxMinDegree(self):
        """Return the highest and lowest number of edges for any node in this network."""
        maxd = 0
        mind = 100000000
        for e in self.nedges.values():
            if e > maxd:
                maxd = e
            if e < mind:
                mind = e
        return (maxd, mind)

    def collectNodes(self):
        self.edgecnt = 0
        with open(self.adjfile, "r") as f:
            for line in f:
                if line[0] != ">":
                    parsed = line.rstrip("\n").split("\t")
                    hub = parsed[0]
                    self.addNode(hub)
                    for i in range(1, len(parsed), 2):
                        self.edgecnt += 1
                        self.addNode(parsed[i])
                        self.countEdge(hub)
                        self.countEdge(parsed[i])
        sys.stderr.write("{} nodes and {} edges read from ADJ file.\n".format(len(self.nodesdb), self.edgecnt))
        sys.stderr.write("Highest / lowest degree: {}.\n".format(self.getMaxMinDegree()))

    def filterNodesByDegree(self):
        """Set to invisible all nodes with a degree smaller than self.mindegree."""
        nvisible = 0
        for (node, deg) in self.nedges.iteritems():
            if deg >= self.mindegree:
                self.visibleNodes[node] = True
                nvisible += 1
            else:
                self.visibleNodes[node] = False
        sys.stderr.write("{} nodes with degree >= {}.\n".format(nvisible, self.mindegree))

    def filterHighestDegree(self):
        """Set as visible the first self.maxnodes nodes in order of decreasing degree, and hide the rest."""
        lowdegree = 0
        degs = [(deg, node) for (node, deg) in self.nedges.iteritems() ]
        degs.sort(key=lambda d: d[0], reverse=True)
        good = degs[:self.maxnodes]
        bad = degs[self.maxnodes:]
        for g in good:
            self.visibleNodes[g[1]] = True
            lowdegree = g[0]
        for b in bad:
            self.visibleNodes[b[1]] = False
        sys.stderr.write("Top {} nodes visible, lowest degree={}.\n".format(self.maxnodes, lowdegree))

    def writeNodes(self, stream, attributes=None):
        """If `attributes' is supplied, it should be an AttrReader object."""
        self.nodesWritten = 0
        comma = False
        stream.write('{  "nodes" : [ ')
        for node, idx in self.nodesdb.iteritems():
            if self.visibleNodes[node]:
                if comma:
                    stream.write(", ")
                stream.write('{\n    "@id" : ' + str(idx) + ',\n    "n" : "' + node + '"\n  }')
                comma = True
                self.nodesWritten += 1
        stream.write(' ]\n}')

        if attributes:
            stream.write(', {  "nodeAttributes" : [ ')
            comma = False
            for node, idx in self.nodesdb.iteritems():
                if self.visibleNodes[node]:
                    attrs = attributes.attributesFor(node)
                    if attrs:
                        attrlist = [ '"{}:{}"'.format(apair[0], apair[1]) for apair in attrs ]
                        if comma:
                            stream.write(", ")
                        stream.write('{\n    "po" : ' + str(idx) + ',\n    "n" : "alias",\n    "v" : [ ' + ", ".join(attrlist) + '],\n    "d" : "list_of_string"\n}')
                        comma = True

                        # for apair in attrs:
                        #     if comma:
                        #         stream.write(", ")
                        #     #stream.write('{\n    "po" : ' + str(idx) + ',\n    "n" : "' + apair[0] + '",\n    "v" : "' + apair[1] + '"\n  }')
                        #     comma = True
                    self.nodesWritten += 1
            stream.write(' ]\n}')

        sys.stderr.write("{} nodes written to CX stream.\n".format(self.nodesWritten))

    def writeEdges(self, stream):
        self.edgesWritten = 0
        comma = False
        mis = []
        stream.write('{  "edges" : [ ')
        with open(self.adjfile, "r") as f:
            for line in f:
                if line[0] != ">":
                    parsed = line.rstrip("\n").split("\t")
                    hub = parsed[0]
                    hubi = self.nodesdb[hub]
                    for i in range(1, len(parsed), 2):
                        gene = parsed[i]
                        if self.visibleNodes[hub] and self.visibleNodes[gene]:
                            genei = self.nodesdb[gene]
                            mi = parsed[i+1]
                            if comma:
                                stream.write(", ")
                            self.edgecnt += 1
                            mis.append((self.edgecnt, mi))
                            stream.write('{\n    "@id" : ' + str(self.edgecnt) + ',\n    "s" : ' + str(hubi) + ',\n    "t" : ' + str(genei) + '\n  }')
                            comma = True
                            self.edgesWritten += 1
        stream.write('  ]\n}')

        # Now write edge attributes
        stream.write(', {  "edgeAttributes" : [ ')
        comma = False
        for m in mis:
            if comma:
                stream.write(", ")
            stream.write('{\n    "po" : ' + str(m[0]) + ',\n    "n" : "MI",\n    "v" : "' + m[1] + '"\n  }')
            comma = True
        stream.write('  ]\n}')
        sys.stderr.write("{} edges written to CX stream.\n".format(self.edgesWritten))

### Visual properties

    def scaleDegree(self, deg, mind, maxd, minr, maxr):
        return minr + (float(deg - mind) / (maxd - mind)) * (maxr - minr)

    def getFPRrange(self):
        minf = 1.0
        maxf = 0.0
        for m in self.mis:
            fpr = float(m[1][self.FPRcol])
            if fpr == 0:
                fpr = 1e6
            else:
                fpr = 1.0/fpr
            if fpr < minf:
                minf = fpr
            if fpr > maxf:
                maxf = fpr
        return (minf, maxf)

# edge defaults cause the ndex viewer to hang...
#  }, {
#    "properties_of" : "edges:default",
#    "properties" : {
#      "EDGE_WIDTH" : "3.0"
#    }

    def writeVisualProperties(self, stream):
        stream.write("""{ "cyVisualProperties" : [ {
    "properties_of" : "network",
    "properties" : {
      "NETWORK_BACKGROUND_PAINT" : "#000000"
    }
  }, {
    "properties_of" : "nodes:default",
    "properties" : {
      "NODE_BORDER_PAINT" : "#000000",
      "NODE_BORDER_STROKE" : "SOLID",
      "NODE_BORDER_WIDTH" : "1.0",
      "NODE_FILL_COLOR" : "#428BCA",
      "NODE_LABEL_COLOR" : "#DDDDDD",
      "NODE_LABEL_FONT_FACE" : "Dialog,plain,12",
      "NODE_TRANSPARENCY" : "255"
    }
  }, """)

        # Write node properties
        (maxd, mind) = self.getMaxMinDegree()
        comma = False
        for node, idx in self.nodesdb.iteritems():
            if self.visibleNodes[node]:
                if comma:
                    stream.write(", ")
                comma = True
                deg = self.nedges[node]
                size = self.scaleDegree(deg, maxd, mind, self.minNodeSize, self.maxNodeSize)
                font = self.scaleDegree(deg, maxd, mind, self.minFontSize, self.maxFontSize)
                stream.write('{\n  "properties_of" : "nodes",\n  "applies_to" : "' + str(idx) + '",\n')
                stream.write('  "properties" : {\n    "NODE_LABEL" : "' + node + '",\n    "NODE_WIDTH" : "' + str(size) + '",\n')
                stream.write('    "NODE_LABEL_FONT_SIZE" : "' + str(font) + '",\n')
                stream.write('    "NODE_HEIGHT" : "' + str(size) + '"\n    }\n  }')


        # Write edge properties
        if self.FPRcol:
            (minf, maxf) = self.getFPRrange()
            for m in self.mis:
                edgeidx = m[0]
                fpr     = float(m[1][self.FPRcol])
                if fpr == 0:
                    fpr = 1e6
                else:
                    fpr = 1.0/fpr
                size    = self.scaleDegree(fpr, minf, maxf, self.minEdgeSize, self.maxEdgeSize)
                stream.write(', {\n  "properties_of" : "edges",\n  "applies_to" : "' + str(edgeidx) + '",\n')
                stream.write('  "properties" : {\n    "EDGE_WIDTH" : "3.0",\n    "EDGE_TRANSPARENCY" : "' + str(size) + '"\n    }\n  }')
        stream.write('  ]\n}')

    def writeContext(self, stream):
        stream.write(""", { "@context" : [ {
  "AMAZE" : "http://uri.ndexbio.org/ns/15a017bb-6196-11e5-8ac5-06603eb7f303/AMAZE/",
  "BIOGRID" : "http://identifiers.org/biogrid/",
  "CAS" : "http://identifiers.org/cas/",
  "CHEMBL" : "http://identifiers.org/chembl.compound/",
  "CHEMICAL COMPONENT DICTIONARY" : "http://uri.ndexbio.org/ns/15a017bb-6196-11e5-8ac5-06603eb7f303/CHEMICAL_COMPONENT_DICTIONARY/",
  "CHEMICAL NAME" : "http://uri.ndexbio.org/ns/15a017bb-6196-11e5-8ac5-06603eb7f303/CHEMICAL_NAME/",
  "ChEBI" : "http://identifiers.org/chebi/",
  "Cosmic symbol" : "http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=",
  "DrugBank" : "http://www.drugbank.ca/drugs/",
  "ENSEMBL" : "http://identifiers.org/ensembl/",
  "Ensembl gene ID" : "http://identifiers.org/ensembl/",
  "GENATLAS" : "http://identifiers.org/genatlas/",
  "GENBANK GENE DATABASE" : "http://www.ncbi.nlm.nih.gov/nuccore/",
  "GENBANK INDENTIFIER" : "http://uri.ndexbio.org/ns/15a017bb-6196-11e5-8ac5-06603eb7f303/GENBANK_INDENTIFIER/",
  "GENBANK PROTEIN DATABASE" : "http://www.ncbi.nlm.nih.gov/protein/",
  "GENE ONTOLOGY" : "http://identifiers.org/go/",
  "GENECARDS" : "http://identifiers.org/genecards/",
  "HETEROGEN" : "http://uri.ndexbio.org/ns/15a017bb-6196-11e5-8ac5-06603eb7f303/HETEROGEN/",
  "HGNC SYMBOL" : "http://identifiers.org/hgnc.symbol/",
  "HGNC" : "http://identifiers.org/hgnc/",
  "HPRD" : "http://identifiers.org/hprd/",
  "INTACT" : "http://identifiers.org/intact/",
  "INTERPRO" : "http://identifiers.org/interpro/",
  "InChIKey" : "http://identifiers.org/inchikey/",
  "KEGG COMPOUND" : "http://identifiers.org/kegg.compound/",
  "KEGG REACTION" : "http://identifiers.org/kegg.reaction/",
  "NCBI GENE" : "http://identifiers.org/ncbigene/",
  "NCBI gene ID" : "http://identifiers.org/ncbigene/",
  "OMIM" : "http://identifiers.org/omim/",
  "PANTHER FAMILY" : "http://identifiers.org/panther.family/",
  "PHARMGKB GENE" : "http://uri.ndexbio.org/ns/15a017bb-6196-11e5-8ac5-06603eb7f303/PHARMGKB_GENE/",
  "PROTEIN DATA BANK" : "http://identifiers.org/pdb/",
  "PUBCHEM-COMPOUND" : "http://identifiers.org/pubchem.compound/",
  "PUBCHEM-SUBSTANCE" : "http://uri.ndexbio.org/ns/15a017bb-6196-11e5-8ac5-06603eb7f303/PUBCHEM-SUBSTANCE/",
  "RAT GENOME DATABASE" : "http://identifiers.org/rgd/",
  "REACTOME" : "http://identifiers.org/reactome/",
  "REFSEQ" : "http://identifiers.org/refseq/",
  "RefSeq accession" : "http://identifiers.org/refseq/",
  "UCSC gene ID" : "http://genome.ucsc.edu/cgi-bin/hgGene?org=human&db=hg38&hgg_gene=",
  "UNIPROT ACCESSION" : "http://uri.ndexbio.org/ns/15a017bb-6196-11e5-8ac5-06603eb7f303/UNIPROT_ACCESSION/",
  "UniProt accession" : "http://identifiers.org/uniprot/",
  "UniProt Isoform" : "http://identifiers.org/uniprot.isoform/",
  "UniProt Knowledgebase" : "http://identifiers.org/uniprot/",
  "VEGA" : "http://vega.archive.ensembl.org/Homo_sapiens/Gene/Summary?g=",
  "Vega gene ID" : "http://vega.archive.ensembl.org/Homo_sapiens/Gene/Summary?g="
}]}""")

    def writeClosing(self, stream):
        stream.write(""",{  "status" : [ {
    "error" : "",
    "success" : true
  } ]
}""")
        stream.write(' ]\n')

    def writeNetwork(self, stream=sys.stdout, attributes={}, nodeAttributes=None):
        self.writePreamble(stream)
        self.writeNetworkAttributes(stream, attributes)
        stream.write(", ")
        self.collectNodes()
        if self.mindegree:
            self.filterNodesByDegree()
        elif self.maxnodes:
            self.filterHighestDegree()
        self.writeNodes(stream, nodeAttributes)
        stream.write(", ")
        self.writeEdges(stream)
        stream.write(", ")
        self.writeVisualProperties(stream)
        self.writeContext(stream)
        self.writeMetadata(stream)
        self.writeClosing(stream)

## CXwriter using cytoscape file as input        

class CXwriterCyto(CXwriter):

    def collectNodes(self):
        self.edgecnt = 0
        with open(self.adjfile, "r") as f:
            for line in f:
                if line[0] != "#":
                    parsed = line.rstrip("\n").split("\t")
                    hub = parsed[0]
                    other = parsed[1]
                    self.addNode(hub)
                    self.addNode(other)
                    self.countEdge(hub)
                    self.countEdge(other)
                    self.edgecnt += 1

        sys.stderr.write("{} nodes and {} edges read from ADJ file.\n".format(len(self.nodesdb), self.edgecnt))
        sys.stderr.write("Highest / lowest degree: {}.\n".format(self.getMaxMinDegree()))

    def writeEdges(self, stream):
        self.edgesWritten = 0
        comma = False
        self.mis = []
        stream.write('{  "edges" : [ ')
        with open(self.adjfile, "r") as f:
            hdr = f.readline()
            edgeAttrs = hdr.rstrip("\n").split("\t")[2:]
            if 'FPR' in edgeAttrs:
                self.FPRcol = edgeAttrs.index('FPR')
            nattrs = len(edgeAttrs)

            for line in f:
                if line[0] != "#":
                    parsed = line.rstrip("\n").split("\t")
                    hub = parsed[0]
                    hubi = self.nodesdb[hub]
                    gene = parsed[1]
                    genei = self.nodesdb[gene]

                    if self.visibleNodes[hub] and self.visibleNodes[gene]:
                        data = parsed[2:]
                        if comma:
                            stream.write(", ")
                            
                        self.edgecnt += 1
                        self.mis.append((self.edgecnt, data))
                        stream.write('{\n    "@id" : ' + str(self.edgecnt) + ',\n    "s" : ' + str(hubi) + ',\n    "t" : ' + str(genei) + '\n  }')
                        comma = True
                        self.edgesWritten += 1
        stream.write('  ]\n}')

        # Now write edge attributes
        stream.write(', {  "edgeAttributes" : [ ')
        comma = False
        for m in self.mis:
            for i in range(nattrs):
                if comma:
                    stream.write(", ")
                data = m[1]
                stream.write('{\n    "po" : ' + str(m[0]) + ',\n    "n" : "' + edgeAttrs[i] + '",\n    "v" : "' + data[i] + '"\n  }')
                comma = True
        stream.write('  ]\n}')
        self.nEdgeAttributes = nattrs
        sys.stderr.write("Edges have {} attributes.\n".format(self.nEdgeAttributes))
        sys.stderr.write("{} edges written to CX stream.\n".format(self.edgesWritten))

### Top-level

if __name__ == "__main__":
    C = CXwriter(sys.argv[1])
    C.maxnodes = 4
    C.writeNetwork(stream=sys.stdout, attributes={'name': 'TestNet', 'description': 'A test network converted from ADJ'})
