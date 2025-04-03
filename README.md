# JointCalledSet

A Python toolkit for merging and analyzing concordance between genetic variant datasets from different sequencing technologies.

## Overview

JointCalledSet provides a comprehensive set of tools for analyzing genetic variant data from multiple sources, with a specific focus on comparing and merging exome sequencing data with imputation data. The toolkit enables researchers to assess concordance between different sequencing technologies, identify discrepancies, and create high-quality merged datasets.

Key features include:
- Loading and merging variant call files from different sources (exome and imputation)
- Detailed concordance analysis between datasets
- Stratification of concordance by allele count (AC) and minor allele frequency (MAF)
- Visualization of concordance metrics
- Multi-allelic variant handling
- Flexible export options for downstream analysis

## Dependencies

The tool requires the following Python packages:
- [Hail](https://hail.is/)
- Pandas
- NumPy
- Bokeh
- PySpark

## Installation

```bash
pip install hail pandas numpy bokeh

git clone https://github.com/gnchau/JointCalledSet.git
cd JointCalledSet
```

## Quick Start

```python
import hail as hl
from JointCalledSet import JointCalledSet

# Initialize Hail
hl.init()

# Create a joint set from exome and imputation data
joint_set = JointCalledSet(
    exo_path='path/to/exome.vcf.gz',
    imp_path='path/to/imputation.vcf.gz',
    output_directory='output/'
)

# Analyze the data
joint_set.concordance()

# Export the merged set
joint_set.export_table('vcf', file_name_prefix='merged_variants')

# For detailed concordance analysis
results = joint_set.analyze_concordance(joint_set.exo, joint_set.imp)
```

## Core Components

### JointCalledSet Class

The primary class for loading, merging, and analyzing datasets:

```python
joint_set = JointCalledSet(
    exo_path='path/to/exome.vcf.gz',
    imp_path='path/to/imputation.vcf.gz',
    output_directory='output/',
    priority='exome',
    keep_all_fields=False,
    impute_nonsimilar_samples=False
)
```

#### Parameters:

- `exo_path`: Path to the exome variant call file in VCF format
- `imp_path`: Path to the imputation variant call file in VCF format
- `output_directory`: Directory where the output files will be saved
- `log_outname`: Name of the log file (default: generated based on date)
- `priority`: The priority dataset ('exome' or 'imputation') for merging
- `keep_all_fields`: Whether to retain all fields from the original datasets
- `impute_nonsimilar_samples`: Whether to impute genotypes for samples not present in both datasets

### Key Methods

#### Loading and Merging

- `load_sets(exo_path, imp_path)`: Loads VCF files into Hail MatrixTables
- `merge()`: Merges datasets giving priority to specified dataset
- `set_priority(new)`: Changes merge priority ('exome' or 'imputation')

#### Analysis

- `concordance(output_name=None)`: Calculates concordance between datasets
- `calc_overlap(output_name=None)`: Calculates overlapping variants and samples
- `plot_concordance(col_conc_fname)`: Visualizes concordance metrics
- `split_multi()`: Splits multi-allelic variants for downstream analysis
- `remove_discordant(which, thresh=None)`: Filters out discordant variants or samples

#### Export

- `export_table(out_type, file_name_prefix='merged', overwrite=False)`: Exports merged dataset
- `export_samples(sample_file_prefix='jc_samples')`: Exports sample IDs
- `save_log()`: Saves operations log to file

### Advanced Concordance Analysis

The `analyze_concordance` function provides detailed concordance metrics stratified by allele count and minor allele frequency:

```python
results = analyze_concordance(exome_mt, imputation_mt, filtered=False)
display_result_tables(results)
```

This function:
1. Filters datasets to bi-allelic variants
2. Annotates variants with AC and MAF bins
3. Computes concordance matrices for each bin
4. Calculates multiple concordance metrics including:
   - Overall concordance rate
   - Non-reference concordance
   - Heterozygote concordance
5. Produces summary tables for interpretation

## Output Formats

The toolkit supports various output formats:

- **VCF**: Standard format for variant calls
- **PLINK**: Binary format for genetic analysis
- **Hail MatrixTable**: Native format for continued analysis with Hail
- **TSV/CSV**: Concordance and sample information in tabular format
- **JSON**: Detailed concordance metrics and results

## Concordance Metrics

The toolkit calculates several key metrics:

- **Overall concordance**: Percentage of genotypes that match between datasets
- **Non-reference concordance**: Concordance considering only non-reference genotypes
- **Heterozygote concordance**: How well heterozygous calls are preserved

These metrics are stratified by:
- **Allele Count (AC) bins**: "1", "2", "3", "4", "5", "6-10", "10+"
- **Minor Allele Frequency (MAF) bins**: "1-2%", "2-5%", "5%+"

## Example Workflow

```python
import hail as hl
from JointCalledSet import JointCalledSet

# Initialize Hail
hl.init()

# Load and merge datasets
joint_set = JointCalledSet(
    exo_path='exome_data.vcf.gz',
    imp_path='imputation_data.vcf.gz',
    output_directory='results/',
    priority='exome'
)

# Check overlap between datasets
joint_set.calc_overlap(output_name='overlapping_variants')

# Basic concordance analysis
joint_set.concordance(output_name='basic_concordance')

# Remove highly discordant variants
joint_set.remove_discordant('variant', thresh=10)

# Export merged dataset
joint_set.export_table('vcf', file_name_prefix='high_quality_merged')

# Detailed concordance analysis
results = joint_set.analyze_concordance(joint_set.exo, joint_set.imp)
display_result_tables(results)
```

## Some Notes
