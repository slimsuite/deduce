# Deduce v1.0.0

## Introduction

*DedUCE* is a command-line tool which efficiently finds ultra-conserved elements in a set of genomes. 

An ultra-conserved element, or UCE, is a long DNA sequence substring that appears in multiple highly-divergent genomes. For the purposes of this tool, we define a UCE as a *k*mer which is present in *p*% or more of the genomes.

### Uses

DedUCE is an excellent tool to assist researchers with the study of UCEs. DedUCE has two primary functions. The
first is to efficiently find UCEs across multiple genomes. This can be used to identify new UCEs for
further research or compare UCE ordering across different genomes. The other primary function is aligning known UCEs
with a genome. This functionality can be utilised for the validation of a whole genome assembly, as well as comparing UCEs across genome
mutations. Furthermore, it allows a user to check for known UCEs in alternate genomes.

## Getting started

### Prerequisites

Before DedUCE can function correctly, these dependencies need to be installed:

* [Jellyfish2](https://github.com/gmarcais/Jellyfish#installation)
* [Minimap2](https://github.com/lh3/minimap2#installation)
* [Bowtie2](https://github.com/BenLangmead/bowtie2#obtaining-bowtie2)
* [Python (3.6 minimum)](https://www.python.org/downloads/)

### Installing

To acquire DedUCE via command line, use the `pip` package manager:

```
$ pip install deduce_uces
$ pip3 install deduce_uces # On Ubuntu
```

## Usage

DedUCE provides a number of subcommands. Each subcommand is configured independently via command-line flags, and the overall operation is also controlled by flags common across all subcommands. 

### Accessing the built-in help

To view the available flags and their defaults, use:

```
deduce --help
deduce <subcommand> --help
```

### `info` command

`info` provides the DedUCE version and checks that all necessary dependencies are installed.

If Deduce is set up correctly, your output should look like:

```
[INFO] dedUCE version 1.0.0
[INFO] Checking dependencies...
[INFO] Found jellyfish at /usr/bin/jellyfish
[INFO] Found bowtie2 at /usr/bin/bowtie2
[INFO] Found minimap2 at /usr/bin/minimap2
```

`info` will show a `[WARNING]` for missing dependencies.

### `find` command

`find` takes a list of genomes as input and produces a list of UCEs.

**Input:** one or more FASTA files

**Output:** a single FASTA file containing a sequence per UCE (by default, this file is named `deduce_uces.fa`)

#### Implementation

`find` uses Jellyfish to generate and manipulate hashmaps containing kmer counts. Briefly, it follows these steps:

1. Generate a hash for each input genome
2. Filter the hashes for kmers with a count of `1` (to ensure uniqueness)
3. Merge the hashes, summing kmer counts
4. Filtering the combined hash to find kmers with a count greater than the required threshold

#### Options

Subcommand-specific options:

* `--dump_hashes DIRECTORY`: instead of writing Jellyfish hashes to a temporary directory, output to the given directory
* `--use_dumped`: instead of generating new Jellyfish hashes, read existing hashes from the working directory. They should be named the same as the genomes, e.g. for `human.fa` the corresponding hash should be named `human.jf`.
* `--output OUTPUT`: override the default output filename
* `--hash_size SIZE`: the size of the hash to generate with Jellyfish (default `500M`). This number specifies how many entries Jellyfish should use to store kmer counts. A larger hash uses more memory, but is faster. You should experiment with this value and set it to the largest size your system can support; for a machine with 32GB of memory this is around `2G`. For more, see the [Jellyfish documentation](https://github.com/gmarcais/Jellyfish/tree/master/doc#counting-k-mers-in-sequencing-reads). 

Useful global options:

* `--min_length`: the size of the kmers to count (default 100)
* `--min_genomes`: the threshold for inclusion in the final UCE list. A UCE must be present in at least `min_genomes` out of `total_genomes` to be included (defaults to the number of genomes supplied)

#### Example usage

To compare two genomes using the default parameters:

```
$ deduce find MT-orang.fa MT-human.fa
$ head deduce_uces.fa
>uce0001
GAACAAGAAAATCAACTTAGAAAAGAAATGTGTATTGGAAAATATCTATTTGCC
>uce0002
ATTACAGCTGTATAATTGAATGGGGAAACCATTGATCTGTAGTAGAAAGGAAAT
>uce0003
AATCATTATTTAGAACATGATTAAATATGTCTGCACCTATTGACTCTGCCTATG
>uce0004
TGATCACCTTTGTCAGTGACAATAGTGTCTTACACAGGCTGATTAGTGGAGTAA
>uce0005
CTAGTGACTAACCTTTTTTGTACTTCTGTTGAATGTTGTTCATATAACTATATC
```

To compare three genomes, requiring UCEs of length >= 500 to appear in at least two of the genomes:

```
$ deduce --min_genomes 2 --min_length 500 find a.fa b.fa c.fa \
         --hash_size 2G --dump_hashes .
```

In the above example, the hashes `a.jf`, `b.jf`, and `c.jf` will be output to the working directory.

### `align` command

`align` takes a file containing UCEs and returns the position of these UCEs (if they exist) within the supplied genomes.

**Input:** one FASTA file containing a set of UCEs; one or more genomes in FASTA format

**Output:** an alignment file for each supplied genome, where the file name is the genome name plus an extension (depending on chosen output format). In each output genome file, regardless of the output format you select, the file will include the aligned positions of each UCE.

#### Implementation

`align` will attempt to align the UCEs from the input UCEs file against the supplied genomes, using the specified algorithm (default is `bowtie`).

The current available algorithms are:

* `bowtie`: This method uses `bowtie2` configured to find exact matches. It requires constant memory of around 3.5GB, regardless of the size of input genomes, but is the slowest method.
* `minimap`: This method uses mappy (an interface to `minimap2`). It requires more memory than `bowtie2` but is faster.
* `inmemory`: This is a naive method that will scan along the genomes until UCE matches are found, and therefore, it is not recommended for long genomes. It reads the entire genome into memory.

#### Options

Subcommand-specific options:

* `--method`: The alignment method (default `bowtie`).
* `--indices`: If the `bowtie2` method is used, you can specify prebuilt indices using this flag which will be used instead of building an index for each supplied genome. This makes the alignment much faster.
* `--output`: The output format for each genome's alignment (default `gff`).

#### Example usage

Using the results of `deduce find`, by default named `deduce_uces.fa`, find the location of the UCEs using:

```
$ deduce align deduce_uces.fa MT-human.fa MT-orang.fa
```

The alignments will be saved to `MT-human.gff` and `MT-orang.gff`.

```
$ cat MT-human.gff
##gff-version 3
MT-human	.	conserved_region	181	279	100	+	ID=uce23
MT-human	.	conserved_region	182	280	100	+	ID=uce11
MT-human	.	conserved_region	183	281	100	+	ID=uce13
MT-human	.	conserved_region	184	282	100	+	ID=uce51
MT-human	.	conserved_region	185	283	100	+	ID=uce10
MT-human	.	conserved_region	186	284	100	+	ID=uce24
MT-human	.	conserved_region	187	285	100	+	ID=uce31
MT-human	.	conserved_region	188	286	100	+	ID=uce29
MT-human	.	conserved_region	189	287	100	+	ID=uce14
...

$ head MT-orang.gff
##gff-version 3
MT-orang	.	conserved_region	361	459	100	+	ID=uce23
MT-orang	.	conserved_region	362	460	100	+	ID=uce11
MT-orang	.	conserved_region	363	461	100	+	ID=uce13
MT-orang	.	conserved_region	364	462	100	+	ID=uce51
MT-orang	.	conserved_region	365	463	100	+	ID=uce10
MT-orang	.	conserved_region	366	464	100	+	ID=uce24
MT-orang	.	conserved_region	367	465	100	+	ID=uce31
MT-orang	.	conserved_region	368	466	100	+	ID=uce29
MT-orang	.	conserved_region	369	467	100	+	ID=uce14
...
```

To instead output in TSV format using the naive in-memory matching method:

```
$ deduce align deduce_uces.fa MT-human.fa MT-orang.fa \
         --output_format tsv --method inmemory
```


### `minimise` command
`minimise` attempts to merge a set of adjacent UCEs found in each genome into a single UCE.

**Input:** At least one file in GFF format for each input genome, containing information about the name of the UCEs and their start and end position in the genome.

**Output** One file for each genome, with the same name as the input and an appended `.min` extension. The files are in GFF format by default, and contain both merged UCEs (where merging is possible) and the remaining unmerged UCEs. 

#### Implementation

DedUCE uses a graph-based method for UCE merging.

For each unique UCE, we add a node. The nodes will store a unique UCE ID, and the start and end locations in each genome in which the UCE is found. For each pair of nodes corresponding to adjacent UCEs, a directed edge connects them (from the previous UCE to the next UCE). All cliques of the graph (representing a group of adjacent UCEs) are traversed from the start node to the end node. As the graph is traversed, the two nodes that are part of the edge are merged into a single node, combining the UCE start positions stored in the node with the outgoing edge and the UCE end positions stored in the node with the incoming edge. This is repeated until the entire clique is traversed and the infomation contained in the final node is returned.

#### Options

* `--output`: The output format for each genome's minimised alignment (default `gff`).

#### Example usage

These examples use the `.gff` files generated from the `align` example.

To merge the UCEs in `MT-human.gff` and `MT-orang.gff`, use:

```
$ deduce minimise MT-human.gff MT-orang.gff
$ cat MT-human.gff.min
##gff-version 3
MT-human	.	conserved_region	181	330	100	+	ID=uce153

$ cat MT-orang.gff.min 
##gff-version 3
MT-orang	.	conserved_region	361	510	100	+	ID=uce153
```

If you prefer TSV output:

```
$ deduce minimise MT-human.gff MT-orang.gff --output tsv
$ cat MT-orang.gff.min
start	end	% identity
360	510	100

$ cat MT-human.gff.min
start	end	% identity
180	330	100
```

Regardless of which output settings you choose, the results will be stored with the extension `.min`.

### `extend` command

`extend` allows the user to search a genome for non-identical UCE matches.

**Input:** A FASTA file containing one or more UCEs, one or more GFF files containing aligned UCEs and the corresponding FASTA file with the genome.

**Output:** A file for every genome supplied containing its aligned UCEs (by default, these files are in GFF format)

#### Implementation

`extend` uses the Python extension of minimap2, `mappy`, to search for instances of each UCE in the supplied FASTA file for each genome and records any instances in a corresponding GFF file.

1. Extracts all UCEs from the FASTA file
2. Read in all UCEs from the aligned UCE files for each genome
3. Build indices for each genome
4. For each UCE map the sequence onto the genome if there does not already exist a match

#### Options

Subcommand-specific options:

* `--min_similarity` : the minimum similarity the alignment must have to be classed as a UCE (default 100)
* `--output_format` : the format of the output file

Useful global options:

* `--min_length`: the size of the kmers to count (default 100)

#### Example usage

To find instances of known UCEs in a genomes using default parameters:

```
$ deduce extend deduce_uces.fa MT-orang.gff MT-orang.fa
$ head MT-orang.gff.extended
##gff version3
MT-orang  . conserved_region  361 459 100  +  .  ID=uce01                     
MT-orang  . conserved_region  362 460 100  +  .  ID=uce27                     
...
```

Finding non-identical occurrences of UCEs in multiple genomes:

```
$ deduce extend deduce_uces.fa MT-orang.gff MT-human.gff MT-orang.fa MT-human.fa \ 
         --min_similarity 90
$ head MT-orang.gff.extended
##gff version3 
MT-orang  . conserved_region  5    104  100  +  .  ID=uce01                     
MT-orang  . conserved_region  290  388  100  +  .  ID=uce27                     
MT-orang  . conserved_region  363  461  97   +  .  ID=uce16                     
```

## Summary of output formats

### `find` command

**Without any options**, the subcommand `find` ouputs a list of UCEs into a file named `deduce_uces.fa` in [FASTA format](https://zhanglab.ccmb.med.umich.edu/FASTA/).

If used with the option `--output OUTPUT`, **OUTPUT** will override the default file name.

### `align, minimise, extend` commands

**Without any options**, the following subcommands produce this output:

* `align`

    Generates file(s) with names corresponding to the genome(s) used to align the UCE set against in [GFF format](https://asia.ensembl.org/info/website/upload/gff.html).

* `minimise`

    Generates file(s) in GFF format with filenames corresponding to the alignment input file(s) with the extension `.min` appended.

* `extend`

    Generates file(s) in GFF format with filenames corresponding to the UCE input file(s) with the extension `.extended` appended.

**With options**, the following subcommands can produce output in these formats:

* GFF format
* [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md)
* [TSV format](https://en.wikipedia.org/wiki/Tab-separated_values)
