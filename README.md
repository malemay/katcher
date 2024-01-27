# Overview

The katcher program aims to extract reads containing any of a set of *k*-mers
from a BAM file. The program is written in the C language and makes use of
multithreading to improve execution speed. It uses the
[htslib](https://github.com/samtools/htslib) library for reading/writing BAM
files and its implementation of
[khash](https://github.com/samtools/htslib/blob/develop/htslib/khash.h) for
efficient hash tables.

# Usage

## katcher

The main program provided in this repository is `katcher`. You can get an overview
of command-line parameters by running `katcher --help`:

	Usage: katcher -i <input.bam> -k <kmer_list.txt>
		-i, --input:     Input file in bam format (required)
		-k, --kmers:     List of k-mers to look for in the reads (required)
		-o, --output:    Name of the output file in bam format (default: stdout)
		-t, --threads:   Number of analysis threads used in addition to decompression thread (default: 1)
		-b, --bufsize:   Number of records stored in the buffers (default: 100,000)
		-m, --maxlength: Maximum read length allowed, used for memory allocation (default: 300)
		-h, --help:      Print this help message

The list of *k*-mers should be a text file with one *k*-mer per line. You do not need
to include the reverse-complemented sequence of the *k*-mers that you are looking for.
`katcher` will automatically add them to the list of queried *k*-mers.

`katcher` will append a `KM` tag to every read it extracts indicating which of the query
*k*-mers were found. If there is more than one query *k*-mer in the read, the value
of the `KM` tag will be a list of comma-separated *k*-mers.

Because of the way multithreading is implemented at the moment, the `--bufsize`
parameter must be a multiple of the number of threads specified by `--threads`. 
`katcher` will complain accordingly if this requirement is not met.

## Other software

The `katcher` program was developed to process the output of the
[*k*-mer GWAS](https://github.com/voichek/kmersGWAS) program developed by
[Voichek and Weigel (2020)](https://doi.org/10.1038/s41588-020-0612-7), but it can be used
for any application requiring the efficient querying of *k*-mers within mapped
reads.  However, three other binaries were developed as part of the work on `katcher`
with a more specific focus on the analysis of the output of the *k*-mer GWAS
approach.

### extract_qname

The `extract_qname` program included in this repository can be used to retrieve
all the reads in a bam file that have the QNAME of any read in a subset of interest.
This can be used, for example, to retrieve all read pairs that have at least one
of the two reads in the pair containing a *k*-mer of interest. This was developed
with the idea of using those reads for local assembly. The usage for this program is:

	Usage: extract_qname <subset.bam> <complete.bam> > <output.bam>

With the output being written to stdout.

### add_pvalues

The `add_pvalues` program can be used to add p-values corresponding to significant
*k*-mers to a BAM file that has already been processed with `katcher`. The usage
for this program is:

	Usage: add_pvalues <input.bam> <pvalues.txt> > <output.bam>

With the output being written to stdout. The p-value file must be a table
of *k*-mers and their p-values as output by the *k*-mer GWAS approach, such
as the `test/kmer_pvalues.txt` file included in this repository. the p-values
will be added as a `PV` tag of comma-separated values corresponding to their
matching `KM` tag.

### list_kmers

Finally, the `list_kmers` program can be used to obtain a list of *k*-mers to
use as input to `katcher` from a presence/absence table of *k*-mers and the
name of a sample of interest. This can be useful to avoid looking for *k*-mers
that are known beforehand not to be present in a given sample. The usage
for this program is:

	Usage: list_kmers <pav_table.txt> <sample_name> > <kmer_list.txt>

An example of a PAV table file (test/significant_pav_table.txt) is included in
the repository. Such a PAV table can be obtained using the `filter_kmers`
program included in the *k*-mer GWAS suite of tools.

# Installation

## Prerequisites

You will need a Linux system with the
[htslib](https://github.com/samtools/htslib) and pthread libraries installed in
order to compile and run `katcher` and associated binaries. Users will find
that `katcher` runs faster if htslib is compiled with libdeflate support,
although this is not necessary. pthread should be installed by default on most
Linux systems.

## Download and compilation

The code can be downloaded by running `git clone https://github.com/malemay/katcher`.
Once this is done, the following commands should be sufficient to compile the programs:

	cd katcher
	make

The four compiled binaries will be located in the bin/ directory. You can add them to
your `$PATH` or call them directly from that location.

IMPORTANT: the *k*-mer length is a constant that is set at compile time. By default, it is
31. If you want to change it, you need to change the value of the `KMER_LENGTH` constant
in `include/functions.h` and re-compile the binaries.

## Testing

You can type `make test` to test that the compiled binaries work properly. This
will use sample data files in the `test` directory. The commands that show up
when running `make test` can be used as examples of how to run these binaries.

# Implementation

The development of `katcher` was motivated by the necessity to query large (>20 Gb)
BAM files for tens of thousands of *k*-mers and hundreds of samples. In order to tackle
this problem, our solution is to build a hash-table of *k*-mers to look for. We then go
along each read and look for sub-sequences of length *k* in the hash table. For this reason,
the *k*-mers that are queried must all be of a constant length that is determined at
compile time.

`katcher` implements multithreading by using a decompression thread that continuously reads
the input BAM file in blocks while analysis threads (specified by the `--threads` option)
look for *k*-mers of interest within the reads. Because the buffer filled by the decompression
thread is split evenly between threads, the `--bufsize` option must be a multiple of the
`--threads` option. This implementation is not optimal, but should be sufficient for most use
cases. Please file a bug in the issues section of this repository if you believe that the
multithreading implementation should be improved.

# Citation

If you use this software, plase cite our publication:

Lemay, M.-A., de Ronne, M., Bélanger, R., & Belzile, F. (2023). k-mer-based GWAS enhances the discovery of causal variants and candidate genes in soybean. *The Plant Genome*, 16, e20374. [doi:10.1002/tpg2.20374](https://doi.org/10.1002/tpg2.20374)

# References

Voichek, Y. and Weigel, D., 2020. Identifying genetic variants underlying
phenotypic variation in plants without complete genomes. *Nature Genetics*,
52(5), pp.534-540. <https://doi.org/10.1038/s41588-020-0612-7>

