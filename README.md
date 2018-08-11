# MiCall Lite

MiCall Lite is a fork of [MiCall](http://github.com/cfe-lab/MiCall), which is a bioinformatic pipeline for mapping FASTQ data to a set of reference sequences to generate consensus sequences, variant calls and coverage maps.  The purpose of MiCall Lite is to provide the core functionality of MiCall in a portable, lightweight version that is easy to install and use.

## Installation

To run Micall-Lite, you need Python 3.x, [bowtie2](https://github.com/BenLangmead/bowtie2), and the Python module [Levenshtein](https://pypi.org/project/python-Levenshtein/).  With these prerequisites in place, you should be able to compile the sources (including some C code) and install MiCall-Lite into your default Python directory with `sudo python3 setup.py install`.  For more detailed instructions, please refer to the [INSTALL.md](INSTALL.md) Markdown document.

## Usage

The most convenient way to run the MiCall-Lite pipeline is through the `run-sample.py` script:
```
art@Kestrel:~/git/MiCall-Lite$ python3 run-sample.py Example_S1_L001_R1_001.fastq.gz Example_S2_L001_R1_001.fastq.gz 
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
art@Kestrel:~/git/MiCall-Lite$ python3 run-sample.py -u Example_S1_L001_R1_001.fastq Example_S2_L001_R1_001.fastq 
```

### Unpaired reads
If your sample was processed using single-read (unpaired) sequencing, then you can run the pipeline with just the one positional argument instead of two:
```
art@Kestrel:~/git/MiCall-Lite$ python3 run-sample.py Example_S1_L001_R1_001.fastq.gz
MiCall-Lite running sample Example_S1_L001_R1_001...
```


