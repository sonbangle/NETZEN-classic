# apple.py
A command-line tool that facilitates working with the [ARACNE](http://califano.c2b2.columbia.edu/aracne) program to reconstruct genetic networks from gene-level or
transcript-level expression data. 

Basic usage:

```
apple.py command arguments...
```

The following commands are available:

|Command|Description|
|-------|-----------|
|[bootstrap](#bootstrap)|Bootstrap a file into *rounds* new files, each containing *samplesize* columns.|
|[consensus](#consensus)|Generate a consensus network from multiple .adj files.|
|[convert](#convert)|Convert *infile* to a different format according to operator *op* and write the results to *outfile*.|
|[extract](#extract)|Extract edges for the genes in file genesfile from the input adj file and write them in tab-delimited format.|
|[filter](#filter)|Filter an adj file keeping only edges with MI over the threshold.|
|[histogram](#histogram)|Generate histogram of MI values from adj files.|
|[random](#random)|Generate random expression data for the genes in *genesfile* on *nsamples* samples using a negative binomial distribution.|
|[stats](#stats)|Print statistics on all supplied filenames (in adj format).|
|[translate](#translate)|Translate identifiers in *infile* writing them to *outfile*.|

Use 

```
apple.py command
```

with no additional arguments to get a description of each command and its options.

## Command descriptions
Detail usage of each apple.py command is provided below.

### Bootstrap

Usage: 

```
apple.py bootstrap [-z samplesize] filename rounds
```

This command takes as input a file containing gene expression values, and generates *rounds* new files through a bootstrap procedure.

The input file is assumed to have genes in the rows and samples in the columns. The first two columns are reserved for gene identifiers. All remaining columns contain data for different samples.

Each output file will have the same number of columns as the input file (unless a different number is specified with the -z argument), chosen at random from the input file, with replacement. Therefore a column from the input file may appear more than once (or not at all) in the output file.

### Consensus

Usage: 

```
apple.py consensus [options] outfile infiles ...
```

This command generates a consensus network from one or more .adj files (usually generated through a bootstrap procedure). The following options are available:

```
  [-c countsfile] - write a tab-delimited file with two columns: support, number of occurrences of support.
  [-d datafile]   - write a tab-delimited file with five columns: hub, gene, support, sum of MI, P-value.
  [-p pval]       - Use the specified P-value to filter edges in output (with Bonferroni correction).
  [-nb]           - If specified, disables Bonferroni correction (used with -p).
  [-s support]    - Only output edges found in at least `support' bootstrap files.
  [-f fraction]   - Like -s, but determines the support in order to have the specified fraction of edges in output.
```

### Convert

Usage: 

```
apple.py convert [options] op infile outfile
```

This command converts between different file formats, according to the specified operator *op*, which can be one of:

```
  na - convert from networkData format to adj
  nc - convert from networkData format to cytoscape
  ca - convert from cytoscape format to adj
  co - convert from cytoscape format to connections
  ac - convert from adj to cytoscape
```

The following table describes the details of each format known to apple.py. In general, all these file formats list all edges connecting pairs of genes, and may provide a measure of the strength of the relationship between the two genes (e.g., mutual information). Some formats (e.g. adj, connection) are hub-oriented: for each hub gene, they list all genes connected to it in the same entry. 

|Format|Details|
|------|-------|
|adj|ARACNE's default output format. The file begines with an optional header that contains information about the ARACNE run parameters (header lines begin with the > character). After the header, lines are tab-delimited and each line refers to a hub gene. The first entry in each line contains the id of the hub gene, while the rest of the line consists of pairs of entries: the id of the connected gene and the MI associated with this edge.|
|connections|A hub-oriented tab-delimited format with three or more columns: hub, number of connected genes, connected genes.|
|cytocscape|A tab-delimited format with three columns: gene1, gene2, mi.|
|networkData|A tab-delimited format with five columns: gene1, gene2, support (number of times this edge was observed), average MI of all observations of this edge, P-value|

### Extract

Usage: 

```
apple.py extract [-a] [-o outfile] adj genesfile [genesfile2]
```
- extract edges for the genes in file genesfile from the input adj file and write them in tab-delimited format.

This command extracts the genes specified in *genesfile* from the .adj file in input and writes their edges to *outfile*.
File *genesfile* should have a single column containing gene identifiers (one per line).

The output (sent to standard output, or to a file specified with the -o option) is tab-delimited with three columns:
hub gene, target gene, MI. The hub gene is always one of the genes specified in *genesfile*, while MI is the mutual
information of the edge connecting it to the target gene.

If the -a option is specified, both the hub gene and the target gene are required to be in *genesfile*.

This command is useful to extract sub-networks in which one gene in each edge (or both) belong to a set of interest,
e.g. a pathway or a functional class.

### Filter

Usage: 

```
apple.py filter [options] infile threshold
```

This command writes a new .adj file containing only the edges in the input .adj file *infile* with an MI value over
the specified *threshold*. The following options are available:

```
  [-o outfile] - write output to *outfile* instead of standard output.
  [-t]         - apply the threshold to the sum of all MI values for a hub, and write
                 the whole line if successful.
```

### Histogram

Usage: 

```
apple.py histogram [options] infile
```

This command computes the histogram of MI values for the edges in the specified .adj file. The following options are available:

```
  [-o outfile] - write output to *outfile* instead of standard output.
  [-n nbins]   - Specifiy number of bins to use (100 by default).
  [-r min max] - Only consider values between *min* and *max* (by default, the whole range of MIs is used).
  [-v]         - If specified, values higher than *max* are added to the last bin.
  [-s]         - If specified, the histogram is computed on the sum of the MIs of each row.
  [-m mifile]  - Write all distinct MI values to *mifile*.
```

### Random

Usage:

```
apple.py random [-o outfile] [-nb nsamples] genesfile
```

This command generates a simulated gene expression dataset, using one of two different methods:

* If -nb is specified, the program will generate *nsamples* values for each gene listed in the first
  column of file *genesfile*, using a negative binomial distribution. The output file will have a
  number of columns equal to *nsamples*+2, with the first two columns containing the gene name (for
  compatibility with ARACNE).

* Otherwise, the expression values in *genesfile* (all values in the row except for the first two)
  will be shuffled.

Output will be written to standard output or to the file specified with the -o option.

### Stats

Usage:

```
apple.py stats [-o outfile] filenames ...
```

This command prints statistics on all the .adj files supplied as arguments. For each file, the command
prints: the number of hub genes, the total number of edges, and the average number of edges per hub.

Output is in tab-delimited format, and is written to standard output or to the file specified with the
-o option.

### Translate

Usage: 

```
apple.py translate table infile outfile ... 
```

This command converts the gene identifiers in *infile* according to a supplied translation table.
Gene identifiers are looked for in the first two columns of each input file.
File *table* should e tab-delimited and contain two columns: original gene name, converted name.

Multiple pairs of input and output files may be specified on the command line. E.g.:

```
apple.py translate table.txt in1.csv out1.csv in2.csv out2.csv in3.csv out3.csv
```
