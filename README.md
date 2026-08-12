# Plasmid Conservation Analysis in *Cupriavidus*

This repository contains the custom computational pipeline and scripts used to assess the structural conservation of the replicons carrying symbiotic genes across the *Burkholderiaceae* family. 

This work is part of the manuscript:
**"Insights into the evolution of symbiosis in the *Cupriavidus* genus: a highly conserved symbiotic island and a patchy phylogenetic distribution"** 
*(Magallanes Alba et al., 2026 - Environmental Microbiology, Under Review).*

## Overview

To trace the evolutionary origins of the cohesive symbiotic machinery in *Cupriavidus*, this pipeline evaluates whether the beta-rhizobial replicon carrying symbiotic genes acts as a single evolutionary unit or if it is the result of a chimeric assembly (mosaicism). 

The script processes BLASTn alignments and intersects them with predefined genomic coordinates to strictly separate the coverage of the "symbiotic island" from the non-symbiotic "plasmid backbone".

## Calculated Metrics

For each subject genome analyzed, the pipeline calculates and outputs the following metrics (which correspond to **Table S3** of the manuscript):

* **`nt_total`**: Total nucleotides aligned to the query plasmid.
* **`nt_in_island`**: Nucleotides aligned strictly within the coordinates of the symbiotic island.
* **`nt_out_island`**: Nucleotides aligned outside the island (representing the plasmid backbone).
* **`pct_total`**: Total percentage of coverage of the query plasmid.
* **`pct_isla_plasmid`**: Coverage of the plasmid explained exclusively by the symbiotic island.
* **`pct_out_plasmid`**: Coverage of the plasmid outside the island (backbone coverage).
* **`pct_isla_relative`**: Fraction of the symbiotic island covered by the subject sequence.

## Prerequisites

To run this pipeline, the following tools must be installed in your environment:
* [BLAST+ suite](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [Bedtools](https://bedtools.readthedocs.io/)
* [Python 3.x](https://www.python.org/) / Bash (depending on execution)

## Usage

*Note: Update this section with the exact command needed to run your script.*

```bash
bash plasmid_coverage_postblast.sh \
   -b resultados_uymma_blast.tsv \
   -i isla_uymma02A.tsv \
   -o resultados_mosaico_uymma

## Arguments:

    -b: BLASTn output file in tabular format (outfmt 6/7) containing the alignments of the query plasmid against the subject database. The query length (qlen) must be included in column 11.

    -i: Symbiotic island coordinates file (1-based TSV) defining the genomic boundaries of the island within the query plasmid.

    -o: Name of the output directory where the final coverage matrices and intermediate BED files will be saved.
