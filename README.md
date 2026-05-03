# CS823 Bioinformatics Project

## Regulatory Variant Analysis with AlphaGenome and Evo2

This repository contains the code used for my CS823 Introduction to Bioinformatics course project. The project investigates how single-nucleotide variants in regulatory DNA near T-cell-related transcription factors affect predicted regulatory activity and genomic sequence likelihood.

The analysis focuses on three genes:

- `FOXP3`
- `BACH2`
- `SATB1`

The workflow combines two genomic foundation models:

- **AlphaGenome**: used to predict regulatory genomic tracks from DNA sequence and estimate how mutations affect predicted regulatory activity.
- **Evo2**: used to score reference and alternate DNA sequences and estimate whether a mutation makes the sequence less likely under a genomic sequence model.

The goal is not to treat Evo2 as a validation of AlphaGenome. Instead, the project compares two different model outputs: predicted regulatory impact from AlphaGenome and sequence-likelihood disruption from Evo2.

---

## Repository Contents

```text
compute_gene_correlations.py
regulatory_region_discovery_alphagenome_scan.ipynb
variant_scoring_ui.ipynb
score_gene_evo2.py
README.md
```

### `regulatory_region_discovery_alphagenome_scan.ipynb`

Main AlphaGenome workflow notebook.

This notebook performs the region-discovery and AlphaGenome mutation-scanning portion of the analysis. For each gene, it:

1. Retrieves gene coordinates.
2. Defines a local search region around the gene.
3. Retrieves regulatory annotations from genome annotation sources.
4. Builds candidate regulatory scan windows.
5. Generates single-nucleotide variant candidates.
6. Runs AlphaGenome scoring on candidate variants.
7. Produces ranked variant summaries and heatmaps.

Output files:

```text
<GENE>_gene_coordinates.json
<GENE>_search_locus.json
<GENE>_candidate_annotations.csv
<GENE>_candidate_scan_windows.csv
<GENE>_candidate_variants.csv
<GENE>_alphagenome_variant_summary_bio_ranked.csv
<GENE>_alphagenome_track_variance_bio_ranked.csv
heatmaps/
```

The key output used by the Evo2 stage is:

```text
<GENE>_alphagenome_variant_summary_bio_ranked.csv
```

---

### `variant_scoring_ui.ipynb`

AlphaGenome single-variant visualization notebook.

This notebook is used after AlphaGenome mutation scanning to visualize a selected variant in genomic context. It compares reference and alternate allele predictions and plots regulatory tracks such as:

- ATAC
- DNase
- CAGE
- histone ChIP-seq tracks

This notebook was used to generate REF vs ALT plots for top variants in `FOXP3`, `BACH2`, and `SATB1`.

Example variants used in the project included:

```text
FOXP3  chrX:49269631:T>C
BACH2  chr6:90296863:G>C
SATB1  chr3:18444811:G>C
```

---

### `score_gene_evo2.py`

Command-line Evo2 scoring script.

This script scores variants using Evo2 by comparing reference and alternate DNA sequence likelihoods.

For each variant, the script:

1. Reads variants from an AlphaGenome-ranked CSV or accepts a single variant.
2. Fetches an hg38 reference sequence window from the UCSC sequence API.
3. Confirms the reference base matches the input variant.
4. Creates the alternate sequence.
5. Scores the reference and alternate sequences with Evo2.
6. Computes:

```text
evo2_delta_score = evo2_alt_score - evo2_ref_score
```

Interpretation:

- `evo2_delta_score < 0`: alternate sequence is less likely than the reference sequence.
- `evo2_delta_score > 0`: alternate sequence is more likely than the reference sequence.
- `evo2_delta_score ≈ 0`: mutation has little effect on Evo2 sequence likelihood.

Example:

```bash
python score_gene_evo2.py \
  --gene-symbol FOXP3 \
  --alphagenome-csv FOXP3_alphagenome_variant_summary_bio_ranked.csv \
  --top-n 100 \
  --model-name evo2_7b \
  --output-dir FOXP3_evo2_outputs_top100
```

Outputs:

```text
<GENE>_evo2_outputs_top100/
  input_variants_for_evo2.csv
  sequence_metadata.csv
  <GENE>_evo2_scores.csv
  <GENE>_alphagenome_vs_evo2.png
```

---

### `compute_gene_correlations.py`

Cross-model comparison script.

This script merges AlphaGenome and Evo2 results by `variant_id_custom` and computes Pearson and Spearman correlations.

It compares:

1. AlphaGenome absolute regulatory score vs signed Evo2 delta.
2. AlphaGenome absolute regulatory score vs absolute Evo2 delta.

Example:

```bash
python compute_gene_correlations.py \
  --gene-symbol FOXP3 \
  --alphagenome-csv FOXP3_alphagenome_variant_summary_bio_ranked.csv \
  --evo2-csv FOXP3_evo2_outputs_top100/FOXP3_evo2_scores.csv \
  --output-dir FOXP3_final_correlation \
  --top-n-label top100
```

Outputs:

```text
<GENE>_final_correlation/
  <GENE>_alphagenome_evo2_merged_top100.csv
  <GENE>_correlation_summary_top100.csv
  <GENE>_descriptive_summary_top100.csv
  <GENE>_top10_alphagenome_with_evo2.csv
  <GENE>_top10_evo2_disruptive.csv
  <GENE>_alphagenome_vs_evo2_signed_delta_top100.png
  <GENE>_alphagenome_vs_evo2_abs_delta_top100.png
```

---

## Workflow

### Step 1: Run AlphaGenome region discovery and mutation scanning

Run:

```text
regulatory_region_discovery_alphagenome_scan.ipynb
```

Set the target gene symbol in the notebook:

```python
gene_symbol = "FOXP3"
```

Repeat for:

```python
gene_symbol = "BACH2"
gene_symbol = "SATB1"
```

This produces candidate regulatory regions, mutation sensitivity heatmaps, and ranked AlphaGenome variant summaries.

---

### Step 2: Visualize top variants with AlphaGenome

Run:

```text
variant_scoring_ui.ipynb
```

Use top variants from:

```text
<GENE>_alphagenome_variant_summary_bio_ranked.csv
```

Recommended plot settings used in this project:

```python
ontology_terms = ["EFO:0001187"]

plot_atac = True
plot_dnase = True
plot_cage = True
plot_chip_histone = True

plot_rna_seq = False
plot_chip_tf = False
plot_splice_sites = False
plot_splice_site_usage = False
plot_splice_junctions = False
plot_contact_maps = False

transcription_factors = None
filter_to_positive_strand = False
filter_to_negative_strand = False
```

---

### Step 3: Run Evo2 scoring

Copy the AlphaGenome ranked CSV to the system where Evo2 is installed. Then run:

```bash
python score_gene_evo2.py \
  --gene-symbol <GENE> \
  --alphagenome-csv <GENE>_alphagenome_variant_summary_bio_ranked.csv \
  --top-n 100 \
  --model-name evo2_7b \
  --output-dir <GENE>_evo2_outputs_top100
```

Example:

```bash
python score_gene_evo2.py \
  --gene-symbol BACH2 \
  --alphagenome-csv BACH2_alphagenome_variant_summary_bio_ranked.csv \
  --top-n 100 \
  --model-name evo2_7b \
  --output-dir BACH2_evo2_outputs_top100
```

---

### Step 4: Compute AlphaGenome vs Evo2 correlations

Run:

```bash
python compute_gene_correlations.py \
  --gene-symbol <GENE> \
  --alphagenome-csv <GENE>_alphagenome_variant_summary_bio_ranked.csv \
  --evo2-csv <GENE>_evo2_outputs_top100/<GENE>_evo2_scores.csv \
  --output-dir <GENE>_final_correlation \
  --top-n-label top100
```

Example:

```bash
python compute_gene_correlations.py \
  --gene-symbol SATB1 \
  --alphagenome-csv SATB1_alphagenome_variant_summary_bio_ranked.csv \
  --evo2-csv SATB1_evo2_outputs_top100/SATB1_evo2_scores.csv \
  --output-dir SATB1_final_correlation \
  --top-n-label top100
```

---

## Data and Model Inputs
The analysis uses:

- Human reference genome sequence from hg38/GRCh38.
- Gene and regulatory annotations from genome annotation resources.
- AlphaGenome's pretrained regulatory predictions, including predicted regulatory tracks such as ATAC, DNase, CAGE, and histone marks.
- Evo2 sequence likelihood scoring on reference and alternate DNA sequence windows.

AlphaGenome provides predicted regulatory outputs. Evo2 provides sequence-likelihood scores. The project compares those two model-derived signals.

---

## Environment

### AlphaGenome

The AlphaGenome notebooks were run in Google Colab/Jupyter using the AlphaGenome tutorial as a starting template.

### Evo2

Evo2 was run using NVIDIA H100 GPU nodes on the Waterfield HPC cluster using a dedicated conda environment.

---

## Interpretation of Model Outputs

### AlphaGenome

AlphaGenome estimates how a variant changes predicted regulatory signals. A high AlphaGenome score indicates that the mutation is predicted to alter one or more regulatory tracks.

### Evo2

Evo2 scores the reference and alternate sequence windows. The delta score compares the alternate sequence to the reference sequence:

```text
evo2_delta_score = evo2_alt_score - evo2_ref_score
```

A negative value means the alternate allele is less likely than the reference under Evo2. This does **not** mean the mutation is unlikely to occur in a population. It means the alternate sequence is less consistent with the sequence patterns learned by the model.

### Cross-model comparison

The correlation analysis asks whether variants with high AlphaGenome regulatory impact also show high Evo2 sequence-likelihood disruption. In this project, the models often prioritized different variants, suggesting that regulatory impact and sequence likelihood capture related but distinct properties of noncoding DNA.

---

## Reproducibility Notes

To reproduce the full analysis for a gene:

1. Run `regulatory_region_discovery_alphagenome_scan.ipynb`.
2. Use `variant_scoring_ui.ipynb` to visualize top variants.
3. Run `score_gene_evo2.py` on the ranked AlphaGenome CSV.
4. Run `compute_gene_correlations.py` on the AlphaGenome and Evo2 outputs.

For consistency across genes, the project used the top 100 AlphaGenome-prioritized variants for Evo2 scoring and correlation analysis.

---

## Project Summary

This project demonstrates a reproducible workflow for comparing regulatory-effect predictions from AlphaGenome with sequence-likelihood scores from Evo2 across regulatory regions near `FOXP3`, `BACH2`, and `SATB1`.

The main interpretation is that the two models do not measure the same biological quantity:

- AlphaGenome asks whether a mutation changes predicted regulatory behavior.
- Evo2 asks whether the alternate DNA sequence is less likely than the reference sequence.

Together, they provide complementary views of noncoding regulatory variants.
