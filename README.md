# dedUCE

DedUCE is a tool for efficiently finding ultra-conserved elements across multiple genomes.

## Installation

Deduce is configured using [Flit](https://flit.readthedocs.io/en/latest/). You will need to install Flit with:

```bash
python3 -m pip install flit
```

Then you can install dedUCE locally by:

```bash
git clone git@github.com:slimsuite/deduce.git
make install
```

The tool will then be available as `deduce`. I recommend doing the installation in a virtual environment.

Use python 3.8.X (any 3.8)

To install in a virtual environment:

```
# Set up virtual environment
mkdir /home/user/.venvs/

module load python/3.8.X

python3 -m venv /home/user/.venvs/deduce

source /home/user/.venvs/deduce/bin/activate

python3 -m pip install flit

mkdir ~/tmp

export TMPDIR=~/tmp

# Install deduce
cd /home/user/.venvs/deduce

source /home/user/.venvs/deduce/bin/activate

module load python/3.8.X

git clone git@github.com:slimsuite/deduce.git

# Move to makefile and make install
cd /home/user/.venvs/deduce/deduce

make install

# Copy C++ code and give your user permissions
cp /home/user/.venvs/deduce/deduce/cpp/jf_bit_counts $HOME/bin

chmod ugo+rwx /home/user/bin/jf_bit_counts

chmod ugo+rwx /home/user/.venvs/deduce/deduce/cpp/jf_bit_counts
```


**Important:** the tool relies on some C++ code located at `cpp/jf_bit_counts`. Install in your home directory with `cp /path/to/deduce/cpp/jf_bit_counts $HOME/bin`.

### External dependencies

DedUCE requires the following external tools to be installed and available in your `$PATH`:

- `jellyfish` **>=2.3.0** (https://github.com/gmarcais/Jellyfish)
- samtools **>= 1.13**
- Either `minimap2` **2.17** (https://github.com/lh3/minimap2) or `bowtie2` **2.3.5.1** (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

Include the following lines in your job script to load modules if available:

```bash
module load bowtie/2.3.5.1
module load minimap2/2.17
module load jellyfish/2.3.0
module load samtools/1.13
```

#### Manual compilation

Compile samtools and jellyfish. The following commands will install them in `$HOME/bin`:

```bash
# SAMtools
wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2
bunzip2 samtools-1.13.tar.bz2
tar -xvf samtools-1.13.tar
cd samtools-1.13/
./configure --prefix=$HOME
make -j 4
make install

# Jellyfish
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz
tar -xzvf jellyfish-2.3.0.tar.gz
cd jellyfish-2.3.0/
./configure --prefix=$HOME
make -j 4
make install

# Minimap
wget https://github.com/lh3/minimap2/archive/refs/tags/v2.17.tar.gz
tar -xzvf v2.17.tar.gz.2
cd minimap2-2.17/
make -j 4
cp minimap2 $HOME/bin

# Bowtie2
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.1/bowtie2-2.5.1-linux-x86_64.zip
unzip bowtie2-2.5.1-linux-x86_64.zip
cd bowtie2-2.5.1-linux-x86_64/
cp bowtie2 $HOME/bin
```

To make them visible in your job scripts, add the line:

```bash
export PATH=$HOME/bin:$PATH
```

## Usage

To use dedUCE, point it at a directory containing your input genomes in FASTA format. Hopefully the defaults should be sensible; here are a couple of example invocations.


### 100% homology

To find all UCEs >= 200bp with 100% identity across all input genomes:

```bash
deduce find \
  # Size of hash used by Jellyfish, the program will fail if this is too low. If left out, dedUCE will try to calculate the best hash size based on your available virtual memory and genome sizes.
  --jf-hash-size 4G \ 
  
  # Number of parallel threads in use
  --threads 16 \
  
  # UCE DEFINITION
  # The max. number of times a UCE can appear in a genome
  --uce-max-occurrences 100 \
  
  # The minimum length of a UCE, in bp
  --uce-min-length 200 \
  
  # The minimum %identity for UCEs between genomes
  --uce-core-min-homology 1 \
 
  # OUTPUT SETTINGS
  # Prefix output with a label for this run
  --output deduce_1.0.0_uces_bej04 \
 
  # ALGORITHM SETTINGS 
  # % of genomes which unique core kmers have to appear in to be included as candidates
  --core-kmer-threshold 1 \

  # Size of core kmers used to construct UCEs
  --core-kmer-size 50 \
  
  # Mapping algorithm - either `minimap` or `bowtie`
  --mapper minimap \
  
  # Print debug information, might be quite verbose
  --debug \
  
  # Input directory
  <PATH_TO_GENOME_DIRECTORY>
```

### Inexact homology

Hints:

* DedUCE finds inexact UCEs in two stages:

    1. Core kmer matching: this stage respects the core kmer threshold and allows searching for UCEs that appear in X/Y genomes. The key parameter to set here is `core-kmer-size`: the minimum homology supported in this stage will decrease as core size decreases, at the expense of speed. If you use the `minimap` method, you will need to increase `core-kmer-non-unique-buffer` as the core size decreases, since minimap2 starts to become less accurate. For lower homologies, try the `bowtie` method.
    2. Extension: once a set of UCEs is found in X/Y genomes, they can be mapped to all other input genomes with any homology threshold. There are no parameters to tune for this stage.


To find all UCEs >= 200bp with 98% identity across at least 50% of input genomes:

```bash
deduce find \
  --jf-hash-size 4G \ 
  --threads 16 \
  --uce-max-occurrences 100 \
  --uce-min-length 200 \
  --uce-core-min-homology 0.98 \
  --output deduce_1.0.0_uces_bej04 \
  --core-kmer-threshold 0.5 \
  --core-kmer-size 75 \
  --mapper minimap \
  <PATH_TO_GENOME_DIRECTORY>
```

#### Output formats

DedUCE supports output in FASTA or BED.

The FASTA output will produce a single file ending in `.deduce.fa`, which contains all UCE sequences identified by the tool. Note that if you have a UCE identity threshold of < 100%, the sequence of the UCE might be different between genomes. You can choose which genome to print the sequences from using the `--output-reference` argument:

```bash
--output my_uces
--output-format fasta
--output-reference mouse
```

The BED output will produce a BED file for each input genome, containing the name, position, and % identity of all UCEs:

```bash
--output my_uces
--output-format bed
```

If you're interested in the regions surrounding identified UCEs, you can ask dedUCE to output the X flanking bases on each side of the UCEs using the `--output-flanks X` argument.

#### A note on mappers

I've tested the `minimap` mapping method quite extensively and it works pretty well. However, with smaller core kmers (<100bp) it only achieves 98% mapping accuracy, so there's a chance dedUCE will miss some UCEs. The alternative method, `bowtie`, achieves 100% accuracy but will make the mapping phase take about 10 times longer due to increased indexing time. However, dedUCE will output the indices generated by both mappers during operation. Subsequent runs take about the same amount of time. The `bowtie` method is less tested than `minimap` currently, so there may be bugs.


