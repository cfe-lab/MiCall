---
title: MiCall at BC CfE
---

**MiCall** is a free software bioinformatics pipeline for processing Illumina MiSeq FASTQ data. At BC Centre for Excellence in HIV/AIDS, it performs deep sequencing analysis for HIV, HCV, and SARS-CoV-2, producing consensus sequences, coverage metrics, and resistance interpretations.

MiCall differs from the lab's Sanger-based **ReCall** pipeline by capturing minority variants and heterogeneous viral populations through thousands of short reads per sample. At BC CfE, MiCall is used primarily for research (proviral analysis) and selectively for clinical V3 tropism testing and HCV resistance interpretation.

---

## Storage & Source-of-Truth

MiCall data exists in **three distinct stores**, each authoritative for different aspects:

### 1. RAW_DATA Filesystem

**Purpose**: Shared network storage for run folders and materialized pipeline results.

**Mount point**: `/media/raw_data` (Samba: `//raw-data.bccfe.ca/raw-data`)

**What it stores**:
* MiSeq run folders (inputs: FASTQs, sample sheets, InterOp files)
* MiCall results (outputs: CSVs, coverage maps, versioned subdirectories)
* Trigger flags (`needsprocessing`)

**Source of truth for**: File-based artifacts consumed by downstream scripts.

**Run folder naming convention**:
```
YYMMDD_<MACHINEID>_<FOLDERID>_000000000-<CELL>
Example: 230224_M04401_0258_000000000-KKRP5
```

**Results path pattern**:
```
/media/raw_data/MiSeq/runs/<run_name>/Results/<version>/
```
Multiple versions can coexist; each corresponds to a MiCall release (e.g., `version_7.17`).

**Immutability**: Do not modify files manually; managed by pipeline software and automated scripts.

### 2. Kive Provenance Store

**Purpose**: Archival storage for reproducible job execution and provenance tracking.

**What it stores**:
* All job inputs (FASTQs, sample sheets, references) uploaded by MiCall watcher
* All job outputs (CSVs, logs) produced by containerized pipeline
* Exact software versions (Singularity images), parameters, timestamps
* Job metadata (who submitted, when, success/failure status)

**Source of truth for**: Historical reproducibility and version comparison.

**Interaction model**:
* **Upload**: MiCall watcher uploads inputs via Python API (`kiveapi`)
* **Execute**: Kive runs containerized MiCall (Singularity/Docker)
* **Download**: Watcher retrieves outputs, writes to RAW_DATA

**Key property**: Results archived in Kive are **immutable**; RAW_DATA copies are **materialized views** for downstream consumption.

### 3. QAI Database (Oracle)

**Purpose**: Operational database for sample tracking, run metadata, and selected metrics accessible via web UI.

**What it stores**:
* Consensus sequences (from `conseq.csv`)
* Coverage metrics (min/max/avg coverage, positions)
* Read counts (raw, mapped)
* QC scores and acceptance flags
* Sample-to-run mappings

**Source of truth for**: Lab review interface and operational queries.

**Interaction model**:
* **Input**: QAI produces sample sheets, places in run folders
* **Output**: MiCall uploads metrics via REST API after results download
* **Consumption**: Lab staff review via QAI web interface; CFE-scripts poll for completed runs

---

## Dataflow: Run Lifecycle

MiCall processing follows a multi-stage artifact pipeline:

1. **Transport**: MiSeq instrument → RAW_DATA
   * Illumina daemon (Windows) copies completed runs to `/media/raw_data/MiSeq/runs/<run_name>/`
   * Includes FASTQs, sample sheets, InterOp files

2. **Discovery**: ScriptBunny → `needsprocessing` flag
   * ScriptBunny polls RAW_DATA every 15 minutes
   * Checks file stability (completion markers, stable sizes)
   * Creates `needsprocessing` flag in run folder

3. **Processing**: MiCall watcher → Kive job
   * MiCall watcher polls every 10 minutes, detects flag
   * Uploads inputs to Kive via Python API
   * Submits containerized job (Singularity image)
   * Removes `needsprocessing` flag

4. **Materialization**: Kive → RAW_DATA results
   * Watcher polls Kive for completion
   * Downloads all outputs from Kive archive
   * Writes to `/media/raw_data/MiSeq/runs/<run_name>/Results/<version>/`

5. **Publication**: MiCall → QAI database
   * Final pipeline step uploads metrics via REST API
   * Writes consensus, coverage, read counts to Oracle DB
   * Enables lab review via QAI web interface

6. **Consumption**: Downstream scripts → RAW_DATA
   * `miseq_gen_results.rb` polls `resistance.csv` files
   * Proviral tools read `proviral/` subdirectory
   * Report generators consume CSVs

**Parallel stream**: MiseqQC daemon also detects `needsprocessing`, parses InterOp files, uploads QC metrics independently.

---

## Input Artifacts Contract

| Artifact | Produced by | Location(s) | Consumed by | Notes |
|----------|------------|-------------|-------------|-------|
| **Run folder** | MiSeq Illumina daemon | `/media/raw_data/MiSeq/runs/<YYMMDD>_<MACHINE>_<FOLDER>_000000000-<CELL>/` | ScriptBunny, watchers | Run naming enables programmatic parsing |
| **Sample sheet** | QAI (via GUI export) | `<run_folder>/` (format: (Multi-)CSV) | MiCall, MiseqQC | Defines samples, barcodes, projects, references |
| **FASTQs** | MiSeq basecaller | Primary: `<run>/Data/Intensities/BaseCalls/`<br>Alt: `<run>/Alignment_<N>/<proj>/Fastq/`<br>Fallback: `<run>/` | MiCall | Pattern: `<sample>_S<N>_L001_R{1,2}_001.fastq.gz`<br>Also: `Undetermined_S0_*` (failed demux), `PhiX_*` |
| **needsprocessing flag** | ScriptBunny | `<run_folder>/needsprocessing` | MiCall watcher, MiseqQC daemon | Created when files stable; consumed (removed) by watcher |
| **InterOp files** | MiSeq instrument | `<run>/InterOp/*.bin` | MiseqQC | Binary format; contains tile/cycle QC metrics |
| **Metadata XMLs** | MiSeq instrument | `<run>/RunInfo.xml`, `RunParameters.xml` | MiseqQC, diagnostics | Run-level configuration and instrument state |

---

## Output Artifacts Contract

| Artifact | Produced by (MiCall stage/mode) | Location(s) | Consumers | Notes |
|----------|-------------------------------|-------------|-----------|-------|
| **conseq.csv** | Main consensus builder (both modes) | `<run>/Results/<version>/conseq.csv` | QAI upload, downstream scripts | Final consensus sequences per sample/region |
| **conseq_all.csv** | Consensus builder | `<run>/Results/<version>/conseq_all.csv` | Research analysis | Includes intermediate/all consensus variants |
| **amino.csv** | Translation stage | `<run>/Results/<version>/amino.csv` | Resistance interpretation | Translated amino acid sequences by region |
| **nuc.csv** | Nucleotide extraction | `<run>/Results/<version>/nuc.csv` | Analysis scripts | Nucleotide sequences by genomic region |
| **coverage_scores.csv** | Coverage calculator | `<run>/Results/<version>/coverage_scores.csv` | QAI upload, QC review | Min/max/avg coverage per sample |
| **genome_coverage.csv** | Coverage by region | `<run>/Results/<version>/genome_coverage.csv` | Research, visualization | Per-region coverage depth |
| **coverage_maps/** | Visualization stage | `<run>/Results/<version>/coverage_maps/` | Manual review | Graphical plots of read depth across genome |
| **concordance.csv** | Alignment QC | `<run>/Results/<version>/concordance.csv` | QC metrics | Agreement between consensus and reference |
| **cascade.csv** | Pipeline stages | `<run>/Results/<version>/cascade.csv` | Debugging | Sequence evolution through filtering/assembly |
| **mutations.csv** | Mutation caller | `<run>/Results/<version>/mutations.csv` | Resistance, research | All detected nucleotide mutations |
| **resistance.csv** | Resistance interpreter | `<run>/Results/<version>/resistance.csv` | CFE-scripts (`miseq_gen_results.rb`) | Drug resistance scores (HIVdb, HCV databases) |
| **denovo/ subdir** | De novo assembly mode | `<run>/Results/<version>/denovo/*.csv` | Research workflows | Separate outputs when de novo mode used |
| **proviral/ subdir** | Proviral pipeline (cfeproviral) | `<run>/Results/<version>/proviral/*.csv` | Research, intactness analysis | Defect classification, gene extraction, primers |

**Storage duality**: All outputs exist in both Kive archive (immutable, with provenance) and RAW_DATA (materialized view for file-based consumption).

---

## MiCall Modes

MiCall supports two assembly strategies, selectable via sample sheet configuration:

## MiCall Modes

MiCall supports two assembly strategies, selectable via sample sheet configuration:

### Remapping Mode
* **Strategy**: Align reads to reference sequences (e.g., HXB2), iteratively refine, build consensus
* **Process**: Map → update reference → remap → converge → consensus + coverage
* **Outputs**: Standard result files in `Results/<version>/`
* **Performance**: Faster; suitable when sample similar to references
* **Config requirement**: Reference sequences specified in sample sheet

### De Novo Assembly Mode
* **Strategy**: Assemble sequences from scratch using k-mer hashing, then align contigs to identify
* **Process**: K-mer index → overlap detection → contig assembly → stitching → alignment identification
* **Outputs**: Standard files + additional artifacts in `denovo/` subdirectory
* **Performance**: Slower; better for divergent samples, large indels, novel variants
* **Config requirement**: De novo flag in sample sheet

### Common Elements (Both Modes)
* **Preprocessing**: Phred-based quality filtering, adapter stripping, primer removal
* **QC outputs**: `coverage_scores.csv`, `concordance.csv`, `cascade.csv`
* **Consensus formats**: `conseq.csv` (final), `conseq_all.csv` (intermediate)
* **Proviral extension**: `proviral/` subdir produced when proviral targets specified

**Mode selection impact**: Runtime (de novo ~2-3× longer), output subdirectories, sensitivity to novel variants.

---

## System Integrations

### QAI (LIMS Database)

**Technology**: Oracle database + REST API

**MiCall → QAI**:
* **Trigger**: Final step in pipeline after results download
* **Uploaded fields** (from CSVs):
  * Consensus sequences (`conseq.csv`)
  * Coverage metrics: min/max/avg, positions
  * Read counts: raw, mapped
* **Failure surface**: Network errors, API unavailable, malformed CSV, DB lock
* **Debugging**: Check Kive job logs for REST API errors

**QAI → MiCall**:
* **Sample sheets**: Placed in run folder, define samples/projects/references
* **Format**: (Multi-)CSV

**QAI Review UI**:
* Lab staff query by run ID, view metrics, mark acceptance
* Not consumed programmatically by MiCall

### Kive (Pipeline Orchestration)

**Technology**: Python API (`kiveapi`), Singularity containers

**MiCall watcher → Kive**:
* **Upload**: FASTQs, sample sheet via API
* **Submit**: Specify pipeline version (e.g., MiCall 7.17), parameters
* **Poll**: Check job status (queued/running/success/failed)
* **Download**: Retrieve all outputs on success

**Kive archive properties**:
* Immutable storage of all inputs/outputs
* Software version (container hash) recorded
* Enables historical reprocessing: re-run old data with new version

**Failure surfaces**:
* Network/authentication errors
* Cluster resource exhaustion
* Container runtime errors (check Kive job logs)

### CFE-Scripts (Downstream Report Generation)

**Primary consumer**: `miseq_gen_results.rb` (Ruby)

**Polling behavior**:
* Scans `/media/raw_data/MiSeq/runs/*/Results/*/resistance.csv`
* Checks for new/unprocessed files (timestamp or DB flag)
* Extracts resistance scores, generates PDF reports

**Not MiCall's responsibility**:
* Report formatting, clinical thresholds, delivery to clinicians
* Separation enables independent versioning/validation

### MiseqQC (Parallel QC Pipeline)

**Trigger**: Same `needsprocessing` flag (independent consumption)

**Process**:
* Perl daemon parses InterOp binary files (tile QC, Phred distributions)
* Uploads metrics to Oracle (`MISEQQC_*` tables)
* Python tool generates HTML/PDF reports (Levey-Jennings, Westgard rules)

**Not consumed by MiCall**: QC metrics inform manual run acceptance, not pipeline logic.

---

## Troubleshooting Entry Points

## Troubleshooting Entry Points

| Symptom | System/Stage | Check |
|---------|-------------|-------|
| **FASTQs missing** | RAW_DATA transport | 1. Verify run folder copied from MiSeq<br>2. Check hierarchical paths (see Input Artifacts table)<br>3. Confirm sample sheet names match FASTQ prefixes<br>4. Check file permissions |
| **`needsprocessing` missing** | Discovery | 1. Verify ScriptBunny daemon running (15-min cycle)<br>2. Check if files still copying (unstable sizes)<br>3. Look for flag already consumed (check Kive for submitted job)<br>4. Verify ScriptBunny write permissions |
| **Results folder missing** | Kive job | 1. Check Kive job status via UI (queued/running/failed)<br>2. Review Kive job logs for errors<br>3. Verify watcher connectivity to Kive API<br>4. Check cluster resource availability |
| **QAI not updated** | QAI upload| 1. Verify `conseq.csv` and `coverage_scores.csv` exist in Results<br>2. Check Kive job logs for REST API errors<br>3. Confirm QAI endpoint reachable<br>4. Check Oracle DB status (locks, space) |
| **Reports missing** | Downstream scripts | 1. Verify `resistance.csv` exists in Results<br>2. Check `miseq_gen_results.rb` daemon status<br>3. Review CFE-scripts logs<br>4. Not a MiCall issue—contact report team |

---

## Further Reading

* **cfeproviral Documentation**: [https://cfe-lab.github.io/proviral/](https://cfe-lab.github.io/proviral/) — proviral intactness analysis
* **Kive Platform**: [https://cfe-lab.github.io/kive/](https://cfe-lab.github.io/kive/) — job management and reproducibility platform
* **QAI System**: Contact lab staff for training
