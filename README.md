# zUMIs_dev

**zUMIs_dev** is a comprehensive and flexible pipeline designed for processing high-sensitivity full-length transcriptome data. It is an enhanced version of the zUMIs pipeline for MGI, automated to handle workflows from raw FASTQ data to gene expression quantification, specifically optimized for MGI sequencing chemistry.

## Features

- **End-to-End Pipeline**: Handles filtering, mapping, counting, and statistical analysis.
- **Flexible Configuration**: Supports manual, auto, and custom sample types via YAML configuration.
- **Multi-Species Support**: Optimized for Human and Mouse genomes.
- **Efficient Processing**: Distinct stages for Filtering, Mapping, Counting, and Summarising.
- **Detailed Reporting**: Generates comprehensive statistics and analysis results.

## Prerequisites

Ensure you have the following installed on your system:

- **Python 3.6+**
- **PyYAML**: Install via pip:
  ```bash
  pip install PyYAML
  ```
- **STAR**: For genome alignment.
- **Samtools**: For BAM file manipulation.
- **Pigz**: For parallel gzip compression.
- **Seqkit**: For FASTQ manipulation.
- **FeatureCounts**: For read counting.

## Installation

1.  **Clone the Repository**
    ```bash
    git clone https://github.com/Yzh25/Mhsflt_toolkit.git
    cd Mhsflt_toolkit
    ```

2.  **Environment Setup**
    The toolkit relies on specific directory structures and environment settings. Run the release script to initialize the environment:
    ```bash
    sh release_env.sh
    ```

## Configuration

1.  **Build Reference Index**
    Download the genome sequence (FASTA) and gene annotation (GTF) files from [GENCODE](https://www.gencodegenes.org). Build the genome index using STAR.

2.  **Create Config File**
    Copy the example configuration and modify it for your project:
    ```bash
    cp yaml/example_config.yaml my_config.yaml
    ```
    
    Edit `my_config.yaml` to set your paths, sample IDs, and analysis parameters. Key sections include:
    - `project`: Project name (used for output files).
    - `sample`: Define species (`human` or `mouse`) and type (`auto`, `manual`, `custom`).
    - `sequence_files`: Paths to your FASTQ files and read structure definitions.
    - `reference`: Paths to your STAR index and GTF file.
    - `out_dir`: Directory for analysis outputs.

## Usage

Run the analysis pipeline by pointing to your configuration file:

```bash
python3 run_analysis_pipeline.py -y my_config.yaml
```

The pipeline will execute the stages defined in your config (Filtering -> Mapping -> Counting -> Summarising).

## Output Structure

Results will be saved in the `results` subdirectory within your specified `out_dir`.

- `results/`: Final analysis reports and gene expression matrices.
- `analysis/`: Intermediate files including BAMs and stats.
- `config/`: Configuration files used for the run.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.