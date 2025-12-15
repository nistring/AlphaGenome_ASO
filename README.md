# AlphaGenome ASO Workflow

This repository provides a reproducible workflow for antisense oligonucleotide (ASO) masking experiments using [AlphaGenome](https://github.com/google-deepmind/alphagenome). The analysis proceeds through:

1. Visualize the gene and context-specific model outputs
2. Generate ASO-masked sequence variants across a target window
3. Predict effects and compute ASO impact scores
4. Visualize ASO scores as sequence logos around the target exon

The `alphagenome` library is used as-is via its public API.

---
## 1. Installation

Get an AlphaGenome API key: https://deepmind.google.com/science/alphagenome/

Create and activate a Python environment (Python ≥3.10 recommended):

```bash
conda create -y -n alphagenome_env python=3.11
conda activate alphagenome_env
pip install -e alphagenome biopython pandas numpy matplotlib
```

---
## 2. Reference Data (GRCh38 / GENCODE v46)

Download primary assembly FASTA:

```bash
cd data
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
cd ..
```

Note: AlphaGenome currently supports GENCODE up to v46.

---
## 3. Configuration

Edit parameters (e.g., `config_private/SMN2.json`) before running `aso.ipynb`:

- **gene_symbol**: Target gene (e.g., `SMN2`)
- **dna_api_key**: AlphaGenome API key
- **ontology_terms**: Context identifiers (UBERON/CL ontology). See `results/metadata*.csv`
- **requested_outputs**: Prediction types (e.g., `RNA_SEQ`, `SPLICE_SITE_USAGE`)
- **exon_intervals**: Genomic start/end of target exon
- **flank**: Bases around exon to include in ASO window (<250 recommended)
- **ASO_length**: Length of the ASO masking window
- **strand**: Track strand filter (`+`, `-`, `stranded`, `unstranded`, `all`)
- **track_filter**: Substring to filter tracks (optional)
- **SNV**: Single nucleotide variant positions (0-based index; currently SNV is the only supported variant type).

Tips:
- Ensure `results_dir` exists or let the notebook create it; ASO outputs are saved under `results_dir/ASO/`.
- For SNVs, set `position` relative to the resized interval start (0-based); the notebook converts to 1-based for AlphaGenome.

---
## 4. How it works — ASO

ASO experiments mask a short window (length = `ASO_length`) with `N` bases while sliding across a target region around the exon. For each masked variant, AlphaGenome predicts selected outputs (e.g., RNA-seq coverage, splice-site usage). Differences relative to the reference are aggregated into ASO impact scores and visualized as sequence logos over the window.

### Outputs
- `${results_dir}/ASO/{config}_ASO_scores.csv`: Per-ASO scores (one row per mask position) for each requested output type.
- `${results_dir}/ASO/{config}_ASO_{OUTPUT}.bed`: Top and bottom ASOs colored by effect size for genome browser tracks.
- `${results_dir}/ASO/{config}_ASO_{OUTPUT}_full.bed`: Full set of ASO windows (all positions) for the output type.
- SeqLogo plots and overlaid tracks are rendered inline in `aso.ipynb` for quick inspection.