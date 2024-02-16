
The purpose of this script is to merge multiple FASTQ files into a single comprehensive file while simultaneously trimming the sequences based on provided quality scores. Trimming is performed to remove reads with low quality scores, which could skew the analysis outcomes.

## Usage

To use the script, users must provide two mandatory CSV files: a trim plan and a merge plan. Optionally, users can also provide a CSV file listing compressed files or files that require compression. Detailed descriptions for each input file and its structure are provided below.

### Input Files

1. **Trim Plan (trimplan.csv):**
   - This file lists the FASTQ files that require trimming before merging.
   - It has the following column headers:
     - r1: FASTQ file for read 1
     - r2: FASTQ file for read 2
     - ...
     - rn: FASTQ file for read n
     - bad_cycles (optional): A CSV file that lists bad cycles with high error rates to be omitted during the trimming process.
   - Example
     ```csv
     r1,r2,bad_cycles
     sampleA_R1.fastq,sampleA_R2.fastq,bad_cycles_sampleA.csv
     sampleB_R1.fastq,sampleB_R2.fastq,
     ```

2. **Merge Plan (mergeplan.csv):**
   - This file maps input FASTQ files to their desired output files after merging.
   - It consists of the following columns:
     - input: Path to the input FASTQ file to be merged.
     - output: Path to the output FASTQ file after merging.
   - Example:
     ```csv
     input,output
     sampleA_R1.fastq,merged_r1.fastq.gz
     sampleB_R1.fastq,merged_r1.fastq.gz
     sampleA_R2.fastq,merged_R2.fastq
     sampleB_R2.fastq,merged_R2.fastq
     ```

3. **Zip File (zipfile.csv, optional):**
   - Lists the FASTQ files that are gzipped and need uncompressing before processing or the output files that need to be compressed after merging.
   - It has a single column:
     - file: Path to the compressed file.
   - Example:
     ```csv
     file
     sampleA_R1.fastq.gz
     sampleB_R2.fastq.gz
     merged_r1.fastq.gz
     ```

### Command-Line Arguments

The script requires the following command-line arguments to be passed:

- **trimplan:** A CSV file containing the lists of files to be trimmed.
- **mergeplan:** A CSV file containing the merge plan.
- **--zipfile (optional):** A CSV file containing a list of files that are compressed or need to be compressed.

```bash
python "micall/core/merge_fastqs.py" "trimplan.csv" "mergeplan.csv" --zipfile "zipfile.csv"
```

### Input FASTQ files

For merging two FASTQ files from the `A` and `B` sample, the input FASTQ files might look like this:

- **sampleA_R1.fastq:**
  ```
  @SEQ_ID_A1
  GGTAAC...
  +
  !!!BBB...
  ```
  
- **sampleB_R1.fastq:**
  ```
  @SEQ_ID_B1
  CCTTGG...
  +
  !!!BBB...
  ```

After running the script, `merged_r1.fastq.gz` would contain:

```
@SEQ_ID_A1
GGTAAC...
+
!!!BBB...
@SEQ_ID_B1
CCTTGG...
+
!!!BBB...
```

Resulting file may have some bases trimmed (removed from its read sequence).

## Workflow Descriptions

The script implements a series of steps to perform the merging:

1. **File Parsing:** It starts by parsing input CSV files detailing which FASTQ files need to be merged or trimmed and which files are compressed or should be post-process compression.

2. **Quality Trimming:** Reads or cycles marked for poor quality are meticulously trimmed from the FASTQ files to ensure the accuracy of the subsequent analysis steps.

3. **File Concatenation:** The script executes a concatenating routine aligning with the merge plan, producing a single, cohesive FASTQ file from multiple input sequences.

4. **Compression Handling:** When instructed, the script efficiently compresses or decompresses FASTQ files to adequately manage disk space and meet file format requirements.

5. **Finalization:** The cleanup step happens last, removing any temporary files generated throughout the process.
