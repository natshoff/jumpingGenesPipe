# jumpingGenesPipe

A bioinformatics pipeline for estimating transposable element (TE) composition from low-coverage shotgun sequencing data.

## Overview

This pipeline processes Illumina paired-end reads (150 bp) to generate transposable element estimates for *Sorghum halepense*. The workflow includes quality control, decontamination, and TE annotation using Transposome with bootstrapping for robust estimates.

**Study Goal**: Characterizing intra- and inter-population genomic variation through TE composition analysis.

## Pipeline Workflow

The pipeline consists of four main stages:

### Stage 1: Sample Job Submission (`1_callTEpipe.pl`)
- Reads a sample list file
- Locates paired-end FASTQ files for each sample
- Submits batch jobs to SLURM scheduler

### Stage 2: Quality Control & Decontamination (`2_TEpipe_fastp.sh`)
Processes each sample through four steps:

**Step 1 - Quality Filtering & Adapter Trimming (fastp)**
- Removes adapters and low-quality bases
- Trims first 5 bases from forward reads
- Removes PCR duplicates
- Removes polyG tails (common in NovaSeq data)

**Step 2 - Read Decontamination (Bowtie2)**
- **2a**: Removes organellar reads (chloroplast/mitochondria)
- **2c**: Removes bacterial contamination
- Retains only nuclear genomic reads

**Step 3 - Read Interleaving**
- Converts paired-end FASTQ to interleaved FASTA format
- Prepares reads for Transposome analysis

**Step 4 - TE Analysis Setup**
- Initiates Transposome analysis via `3_multiTransposome.pl`

### Stage 3: Bootstrapped TE Estimation (`3_multiTransposome.pl`)
- Generates multiple random subsamples from cleaned reads
- **Default**: 100 iterations of 200,000 read subsamples
- Creates configuration files for each iteration
- Submits individual Transposome jobs

### Stage 4: TE Annotation (`4_transposome.srun`)
- Runs Transposome on each subsample
- Performs BLAST against grass TE database
- Clusters similar sequences
- Annotates repeat families

### Stage 5: Results Summary (`5_processedDataPath.pl`, `6_callSummary.sh`, `7_pullAnnotation.pl`)

**Step 5 - Generate Sample List (`5_processedDataPath.pl`)**
- Scans processed data directory for valid samples
- Filters samples by expected number of transposome runs (e.g., 100)
- Outputs valid samples and excluded samples lists

**Step 6 - Call Summary (`6_callSummary.sh`)**
- SBATCH wrapper script for HPC job submission
- Calls the annotation pulling script

**Step 7 - Pull Annotations (`7_pullAnnotation.pl`)**
- Aggregates TE annotations across all bootstrap runs
- Calculates average and standard deviation of GenomeFraction
- Supports three classification levels: Order, Superfamily, Family
- Outputs tab-separated results files

## Requirements

### Software Dependencies

| Software | Version | Purpose |
|----------|---------|---------|
| Perl | 5.x | Pipeline scripts |
| fastp | [version] | Quality control & adapter trimming |
| Bowtie2 | 2.3 | Contamination removal |
| Transposome | 0.12.1 | TE annotation |
| SLURM | - | Job scheduling |

### Reference Databases

1. **Sorghum organellar genomes** (Bowtie2 index)
   - Location: `/mrm/Sorghum_halepense_Utilities/Organellar/Sorghum_organellar_genomes`

2. **NCBI Bacterial genomes** (Bowtie2 index)
   - Location: `/mrm/bin/Bacterial_Genome_DB/NCBI_Bacterial_bowtie`

3. **Grass TE database** (RepBase 21.10)
   - Location: `/mrm/RepBase21.10.fasta/grasrep.ref`
   - Required for TE annotation in `GRStransposome_config.yml`

### System Requirements

- **Memory**: 100GB per sample (main pipeline)
- **CPUs**: 4 cores per sample
- **Storage**: Sufficient scratch space for intermediate files

## Directory Structure

```
jumpingGenesPipe/
├── scripts/
│   ├── 1_callTEpipe_mapped.pl        # Job submission script (with location map)
│   ├── 2_TEpipe_fastp.sh             # Main processing pipeline
│   ├── 3_multiTransposome.pl         # Transposome bootstrapping
│   ├── 4_transposome.srun            # Transposome execution
│   ├── 5_processedDataPath.pl        # Generate valid sample list
│   ├── 6_callSummary.sh              # SBATCH wrapper for results summary
│   ├── 7_pullAnnotation.pl           # Aggregate TE annotations
│   ├── fastq2fasta_interleaved.pl    # Format conversion
│   ├── GRStransposome_config.yml     # Transposome configuration
│   ├── helper/                       # Utility scripts
│   └── deprecated/                   # Old/replaced scripts
├── samples/
│   ├── samples.txt                   # List of sample IDs
│   └── outputRuns/                   # Sample validation output
│       ├── *_valid.txt               # Samples with correct run count
│       └── *_excluded.txt            # Samples with wrong run count
├── data/
│   ├── raw/                          # Raw FASTQ files
│   ├── processed/                    # Processed sample directories
│   │   └── SAMPLE_NAME/
│   │       ├── 1_TrimmedReads/
│   │       ├── 2_Bowtie/
│   │       │   ├── 2a_Organellar/
│   │       │   └── 2c_Bacterial/
│   │       ├── 3_Interleaved/
│   │       └── 4_Transposome/
│   └── output/                       # Final summarized results
│       ├── *_order_TE_annotation_results_average.txt
│       ├── *_order_TE_annotation_results_stdev.txt
│       ├── *_superfamily_TE_annotation_results_average.txt
│       ├── *_superfamily_TE_annotation_results_stdev.txt
│       ├── *_family_TE_annotation_results_average.txt
│       └── *_family_TE_annotation_results_stdev.txt
└── logs/
    └── slurm-*.out                   # SLURM log files
```

## Usage

### 1. Setup

Clone the repository to your HPC scratch directory:

```bash
cd /scratch/YOUR_USERNAME
git clone git@github.com:natshoff/jumpingGenesPipe.git
cd jumpingGenesPipe
```

### 2. Prepare Sample List

Create a text file with sample identifiers (one per line):

```bash
# samples/samples.txt
JG_CO_011_C
JG_CO_099_A
JG_CO_150_B
```

### 3. Generate Sample Location Map

Since samples are spread across multiple directories, generate a map file that locates each sample:

```bash
# Run as SLURM job (recommended for large sample lists)
sbatch scripts/generate_sample_map.sh samples/samples.txt samples/sample_locations.txt

# Or run directly
chmod +x scripts/generate_sample_map.sh
./scripts/generate_sample_map.sh samples/samples.txt samples/sample_locations.txt
```

This creates a tab-delimited file mapping each sample to its directory:
```
/grps2/mrmckain/Sequence_Vault/NovaSeq.022122	JG_CO_011_C
/scratch/nphofford/jumpingGenesPipe/data/raw/NovaSeq.042622	JG_CO_099_A
/grps2/mrmckain/Sequence_Vault/NovaSeq.110321	JG_CO_150_B
```

**Note**: The script searches the vault directories, but automatically remaps NovaSeq.042622 samples to the scratch location (where converted .gz files are stored).

### 4. Run Pipeline

Submit jobs for all samples using the location map:

```bash
# Run full pipeline (all steps)
perl scripts/1_callTEpipe_mapped.pl samples/sample_locations.txt scripts/2_TEpipe_fastp.sh

# Or start from a specific step (e.g., Step 3)
perl scripts/1_callTEpipe_mapped.pl samples/sample_locations.txt scripts/2_TEpipe_fastp.sh 3
```

**Available steps:**
- Step 1: Quality Filtering & Adapter Trimming
- Step 2: Contamination Removal  
- Step 3: Read Interleaving
- Step 4: Transposome Analysis

This will submit one SLURM job per sample.

**Alternative**: If all your samples are in one directory, use the original script:
```bash
perl scripts/1_callTEpipe.pl samples/samples.txt scripts/2_TEpipe_fastp.sh
```
(You'll need to update the search path in line 23 of that script)

### 5. Monitor Progress

Check job status:

```bash
# View active jobs
squeue -u $USER

# Monitor a specific job log
tail -f logs/slurm-<jobid>.out
```

### 6. Restart Failed Jobs (Optional)

If a job fails at a specific step, you can restart from that step without reprocessing earlier steps:

```bash
# Option 1: Create a location map for just the failed samples
echo "JG_CO_011_C" > samples/failed_samples.txt
./scripts/generate_sample_map.sh samples/failed_samples.txt samples/failed_locations.txt

# Restart from Step 3 (skipping Steps 1-2)
perl scripts/1_callTEpipe_mapped.pl samples/failed_locations.txt scripts/2_TEpipe_fastp.sh 3

# Option 2: Or manually edit sample_locations.txt to include only failed samples
# Then run with the starting step
perl scripts/1_callTEpipe_mapped.pl samples/sample_locations_failed.txt scripts/2_TEpipe_fastp.sh 3
```

The pipeline will skip completed steps and use existing intermediate files.

### 7. Summarize Results

After all Transposome jobs complete, aggregate the annotations:

#### Step 5: Generate Valid Sample List

Filter samples that have exactly 100 transposome output files (or your expected count):

```bash
perl scripts/5_processedDataPath.pl <output_name> <expected_count>

# Example: Find all samples with exactly 100 runs
perl scripts/5_processedDataPath.pl feb1samples 100
```

**Output** (in `/scratch/nphofford/jumpingGenesPipe/samples/outputRuns/`):
- `feb1samples_valid.txt` - Samples with exactly 100 annotation files
- `feb1samples_excluded.txt` - Samples with different counts (for debugging)

#### Step 6: Run Annotation Summary

Submit a batch job to aggregate TE annotations across all bootstrap runs:

```bash
sbatch scripts/6_callSummary.sh <valid_samples_file> <level> <output_name>

# Example: Summarize at superfamily level
sbatch scripts/6_callSummary.sh \
  /scratch/nphofford/jumpingGenesPipe/samples/outputRuns/feb1samples_valid.txt \
  superfamily \
  feb1Superfamily
```

**Classification Levels:**

| Level | Description | Example Values |
|-------|-------------|----------------|
| `order` | Highest level | ltr_retrotransposon, dna_transposon, satellite |
| `superfamily` | Mid level | Gypsy, Copia, hAT |
| `family` | Most specific | ATHILA-1_SBi-I, Copia-141_SB-I |

**Output** (in `/scratch/nphofford/jumpingGenesPipe/data/output/`):
- `<output_name>_<level>_TE_annotation_results_average.txt` - Mean GenomeFraction per TE type
- `<output_name>_<level>_TE_annotation_results_stdev.txt` - Standard deviation across runs

#### Run All Three Levels

To get complete results at all classification levels:

```bash
# Generate sample list once
perl scripts/5_processedDataPath.pl feb1samples 100

# Submit all three levels
sbatch scripts/6_callSummary.sh /scratch/nphofford/jumpingGenesPipe/samples/outputRuns/feb1samples_valid.txt order feb1Order
sbatch scripts/6_callSummary.sh /scratch/nphofford/jumpingGenesPipe/samples/outputRuns/feb1samples_valid.txt superfamily feb1Superfamily
sbatch scripts/6_callSummary.sh /scratch/nphofford/jumpingGenesPipe/samples/outputRuns/feb1samples_valid.txt family feb1Family
```

## Configuration

### Adjusting Transposome Sampling Parameters

In `scripts/3_multiTransposome.pl`, you can modify:

```perl
# Line 14: Number of bootstrap iterations (default: 100)
for (my $i=0; $i<=99; $i++){

# Line 19 & 22: Number of reads per subsample (default: 200,000)
$start_rand = int(rand($total_pairs-200000));
`tail -n $start_rand $ARGV[0] | head -n 200000 > ...`;
```

### Modifying Transposome Parameters

Edit `scripts/GRStransposome_config.yml`:

```yaml
clustering_options:
  - percent_identity:   90      # Clustering threshold
  - fraction_coverage:  0.55    # Minimum alignment coverage

annotation_options:
  - cluster_size:       100     # Minimum cluster size for annotation
```

See the [Transposome documentation](https://github.com/sestaton/Transposome/wiki) for more details.


## Output Files

### Per-Sample Intermediate Files

| Location | Description |
|----------|-------------|
| `data/processed/SAMPLE/1_TrimmedReads/` | Quality-filtered FASTQ files |
| `data/processed/SAMPLE/2_Bowtie/2a_Organellar/` | Reads after organellar removal |
| `data/processed/SAMPLE/2_Bowtie/2c_Bacterial/` | Final cleaned reads (nuclear only) |
| `data/processed/SAMPLE/3_Interleaved/` | Interleaved FASTA for Transposome |
| `data/processed/SAMPLE/4_Transposome/` | 100 transposome output directories |

### Per-Run Transposome Output

Each bootstrap run creates a directory like `SAMPLE_transposome_JOBID_out/` containing:

| File | Description |
|------|-------------|
| `*_cluster_report.txt` | Clustering results |
| `*_cluster_report_annotations.tsv` | Annotated clusters |
| `*_cluster_report_annotations_summary.tsv` | **Main output**: TE type, superfamily, family, GenomeFraction |
| `*_cluster_report_singletons_annotations.tsv` | Singleton read annotations |

### Final Summary Output

Located in `data/output/`:

| File | Description |
|------|-------------|
| `*_order_TE_annotation_results_average.txt` | Mean GenomeFraction by Order |
| `*_order_TE_annotation_results_stdev.txt` | Std dev by Order |
| `*_superfamily_TE_annotation_results_average.txt` | Mean GenomeFraction by Superfamily |
| `*_superfamily_TE_annotation_results_stdev.txt` | Std dev by Superfamily |
| `*_family_TE_annotation_results_average.txt` | Mean GenomeFraction by Family |
| `*_family_TE_annotation_results_stdev.txt` | Std dev by Family |

**Output Format**: Tab-separated with samples as rows and TE types as columns. Values represent the proportion of the genome (GenomeFraction) occupied by each TE type.

## Troubleshooting

### Log Files

- **Pipeline logs**: `logs/slurm-<jobid>.out`
- **Transposome logs**: `data/processed/SAMPLE/4_Transposome/*_log.txt`

## Citation

If you use this pipeline, please cite:

**Transposome:**
- Staton SE, Burke JM (2015). Transposome: A toolkit for annotation of transposable element families from unassembled sequence reads. *Bioinformatics* 31(11): 1827-1829.

**fastp:**
- Chen S, Zhou Y, Chen Y, Gu J (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34(17): i884-i890.

**Bowtie2:**
- Langmead B, Salzberg SL (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods* 9(4): 357-359.

## Author

NP Hofford  
University of Alabama

## License

[Add license information]

## Acknowledgments

- McKain Lab for computational resources
- RepBase for TE reference database

