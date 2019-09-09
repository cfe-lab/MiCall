# MiCall Lite

MiCall Lite is a fork of [MiCall](http://github.com/cfe-lab/MiCall), which is a 
bioinformatic pipeline for mapping FASTQ data to a set of reference sequences to 
generate consensus sequences, variant calls and coverage maps.  The purpose of 
MiCall Lite is to provide the core functionality of MiCall in a portable, 
lightweight version that is easy to install and use.


## Installation

To run Micall-Lite, you need the following:
* a [gcc](https://gcc.gnu.org/) build environment for compiling C source code
* Python 3.x
* [bowtie2](https://github.com/BenLangmead/bowtie2)
* the Python module [Levenshtein](https://pypi.org/project/python-Levenshtein/)
* (optional) [sierra-local](https://github.com/PoonLab/sierra-local), for HIV-1 drug resistance prediction

With these prerequisites in place, you should be able to compile the sources 
(including some C code) and install MiCall-Lite into your default Python 
directory with `sudo python3 setup.py install`. 
For more detailed instructions, please refer to the [INSTALL.md](INSTALL.md) 
Markdown document.  

> If you run into any issues trying to install MiCall-Lite, please 
[post an new issue](https://github.com/PoonLab/MiCall-Lite/issues) on our 
issue tracker so we can help!


## Usage

The most convenient way to run the MiCall-Lite pipeline is through the `bin/micall` script, which should have been installed at `/usr/local/bin`:
```
art@Kestrel:~/git/MiCall-Lite$ micall Example_S1_L001_R1_001.fastq.gz Example_S1_L001_R2_001.fastq.gz 
MiCall-Lite running sample Example_S1_L001_R1_001...
  Preliminary map
  Iterative remap
  Generating alignment file
  Generating count files
```
These paired-end [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files contain roughly 15,000 (2 x 251 nt) reads.  The pipeline required about 30 seconds to process this sample on an Intel i7-8650U processor, with a default 4 threads for bowtie2 processing.

### Compressed and uncompressed FASTQs
Note that we ran the pipeline on [gzip](https://en.wikipedia.org/wiki/Gzip) compressed files, as indicated by the conventional `.gz` file extension.  The FASTQ files compress down to roughly one-eighth of their original size.  Being able to process these files without expanding the data on your storage media is useful for conserving space.  However, you might want to process the uncompressed files instead.  You can do this by running MiCall-Lite with the `-u` (uncompressed) flag:
```
art@Kestrel:~/git/MiCall-Lite$ micall -u Example_S1_L001_R1_001.fastq Example_S1_L001_R2_001.fastq 
```

### Unpaired reads
If your sample was processed using single-read (unpaired) sequencing, then you can run the pipeline with just the one positional argument instead of two:
```
art@Kestrel:~/git/MiCall-Lite$ micall Example_S1_L001_R1_001.fastq.gz
MiCall-Lite running sample Example_S1_L001_R1_001...
```

### Filtering bad tile-cycles
The Illumina platforms produce a set of [InterOp](http://illumina.github.io/interop/index.html) files with each run.  One of these files, `ErrorMetricsOut.bin`, contains the empirical sequencing error rates associated with the [phiX174 control library](https://www.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/phix-control-v3.html) that is usually added ("spiked") into the run.  This file contains useful information because the error rates are broken down by tile and cycle, where each cycle corresponds to a specific nucleotide position in every read.  In previous work, we have observed that specific tile-cycle combinations can exhibit disproportionately high error rates --- as a result, we have developed an optional preliminary step in the MiCall pipeline for parsing this InterOp file and using the error rates to censor base-calls in the FASTQ files affected by bad tile-cycle combinations.

To use this censoring step, you need to have the `ErrorMetricsOut.bin` file from the PC linked to your Illumina instrument or from your sequencing provider, and then feed this file to the pipeline with the `-i` flag:
```
art@Kestrel:~/git/MiCall-Lite$ micall -i ErrorMetricsOut.bin Example_S1_L001_R1_001.fastq.gz Example_S1_L001_R2_001.fastq.gz
MiCall-Lite running sample Example_S1_L001_R1_001...
  Censoring bad tile-cycle combos in FASTQ
  Preliminary map
  Iterative remap
  Generating alignment file
  Generating count files
```


### HIV drug resistance prediction (optional)

The HIV-1 *pol* gene is a frequent target for genetic sequencing because it encodes 
the primary targets of antiretroviral drug treatment.
Drug resistance in HIV-1 is very well characterized and there are several algorithms
for predicting this phenotype from a *pol* sequence. 
We provide a local implementation of the 
[Stanford HIVdb Drug Resistance Database](https://hivdb.stanford.edu/) 
algorithm, which is another Python module that can be installed via `pip`:

```console
Elzar:micall artpoon$ pip3 search sierralocal
sierralocal (0.2.1)  - Local execution of HIVdb algorithm
Elzar:micall artpoon$ pip3 install sierralocal
```

To demonstrate the use of *sierra-local* to predict drug resistance levels
from next-generation sequence data processed by *MiCall-Lite*, we'll use a published 
NGS data set from a 
[study of HIV-1 drug resistance in Uganda](https://www.liebertpub.com/doi/abs/10.1089/AID.2017.0205).
(To obtain this data set, you'll need to install NCBI's 
[sratools](https://ncbi.github.io/sra-tools/), which provides the *fasterq-dump* 
file transfer program).


```console
Elzar:examples artpoon$ fasterq-dump SRS5100454
spots read      : 21,633
reads read      : 43,266
reads written   : 43,266
```

Next, we process these data through the default *MiCall-Lite* pipeline:
```console
Elzar:examples artpoon$ micall SRS5100454_1.fastq SRS5100454_2.fastq -u
Using /usr/local/bin/bowtie2 version 2.3.5.1
MiCall-Lite running sample SRS5100454_1...
  Preliminary map
  Iterative remap
  Generating alignment file
  Generating count files
```



## How it works

MiCall maps short read data to reference sequences using bowtie2.  Since it was designed for rapidly evolving RNA viruses, this mapping occurs through iterative, adaptive stages: 
1. The preliminary mapping determines the most appropriate "seed" reference sequence for a given sample.  
2. Next, Micall uses an iterative "remapping" method that updates the reference sequence using short reads that successfully mapped to the previous reference, essentially evolving the reference sequence towards the ideal reference for that sample.  This iterative stage proceeds until there are no additional reads from the original FASTQ files mapped to the current reference.
3. MiCall merges the mapped paired-end reads, using a rules-based algorithm to reconcile discordant base-calls in overlapping regions, inserts gaps where there are deletions relative to the reference sequence, and gathers identical sequences in a compact output format.
4. Finally, it determines which coordinate reference sequences are partially or fully covered by the mapped reads and outputs the nucleotide and amino acid sequences for the sample using this coordinate system.  Any insertions relative to the coordinate reference are written into a separate file by default.  

There are several optional output files that can be generated by MiCall-Lite.  To make full use of these outputs, you can manually run each Python script or adapt the `run-sample.py` driver script.

## Outputs

* `*.amino.csv` contains the amino acid frequencies at each site relative to the coordinate references that are covered by the sample.
* `*.nuc.csv` contains the nucleotide frequencies at each site relative to the coordinate references that are covered by the sample.
* `*.insert.csv` contains sequence insertions relative to the coordinate references.
* `*.conseq.csv` contains the consensus sequence of the sample with any insertions relative to the coordinate reference removed.  The `MAX` sequence is the [plurality consensus](https://www.ncbi.nlm.nih.gov/pubmed/1515745) sequence.  The additional sequences in this file contain [mixtures](https://en.wikipedia.org/wiki/Nucleic_acid_notation) (ambiguous base calls) at sites where a polymorphism exceeds the corresponding minimum frequency.  For example, the 10% consensus sequence contains mixtures wherever a polymorphism is present a frequency exceeding 10%.



### Customizing projects

[MiCall](https://github.com/cfe-lab/MiCall) was primarily developed to process Illumina FASTQ files based on samples derived from HIV-1, hepatitis C virus and human leukocyte antigen (HLA) genetic material.  These targets are determined by the *seed* and *coordinate* reference sequences, where the seed is used for the preliminary mapping of short reads, and the coordinate is used to interpret the end results of mapping.  Both types of information are stored in a [JSON](https://en.wikipedia.org/wiki/JSON) file called `projects.json`.  This JSON file is normally generated from the [CFE](http://cfenet.ubc.ca) laboratory database, but the PoonLab developers are currently working on a Python script that makes it easier to modify this JSON file, as well as a JavaScript-based graphical user interface.

