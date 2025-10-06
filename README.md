# EGFR Allosteric Drug Discovery Workflow

A complete computational pipeline for unbiased, site-agnostic ranking of EGFR allosteric inhibitor candidates from Boltz-2 structure predictions.

## Overview

This workflow processes molecular docking results from two compound libraries (cmolgpt and enamine_real), discovers binding sites without structural bias, selects optimal poses, extracts protein-ligand interactions, and ranks molecules site-agnostically based on potency, interactions, and confidence.

### Key Features

- **Unbiased site discovery**: DBSCAN clustering on ligand centroids without prior knowledge
- **Deterministic pose selection**: Cascade-based selection prioritizing affinity probability, pIC50, and confidence
- **Site-agnostic ranking**: No preference for specific binding sites
- **Diversity enforcement**: Chemical similarity clustering to avoid redundant hits
- **Reproducible**: Fixed algorithms ensure same input → same output

## Directory Structure

```
egfr_allosteric/
├── data/                     # Raw Boltz-2 outputs (mmCIF + JSON)
│   ├── cmolgpt/
│   │   └── EGFR_ALLO_XXX/
│   │       ├── boltz2_output.json
│   │       └── structures/
│   │           ├── structure_1.mmcif
│   │           ├── structure_2.mmcif
│   │           └── structure_3.mmcif
│   └── enamine_real/
│       └── ZXXXXXXXXXX/
│           ├── boltz2_output.json
│           └── structures/
├── data_std/                 # Standardized PDB files + metrics
│   ├── cmolgpt/
│   └── enamine_real/
├── outputs/                  # Analysis results
│   ├── samples_long.csv
│   ├── site_samples.csv
│   ├── site_clusters.csv
│   ├── master_per_molecule.csv
│   ├── master_with_interactions.csv
│   └── ranked_shortlist.csv
├── scripts/                  # Python pipeline scripts (numbered 00-05)
│   ├── 00_prepare_inputs_from_mmcif.py    # mmCIF conversion + PDB fixing
│   ├── 01_build_samples_long_ic50.py      # IC50 normalization
│   ├── 02_detect_sites_unsupervised.py    # Binding site discovery
│   ├── 03_select_best_pose_unbiased.py    # Best pose selection
│   ├── 04_plip_contacts.py                # Interaction extraction
│   ├── 05_score_and_cluster_unbiased.py   # Final ranking
│   ├── utils.py                           # Shared utilities
│   └── visualize_results.py               # Result visualization
├── requirements.txt
└── README.md
```

## Installation

### Requirements

- Python 3.10 or higher
- PLIP (Protein-Ligand Interaction Profiler) - optional but recommended

### Setup

```bash
# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Optional: Install PLIP for detailed interaction analysis
pip install plip
```

## Usage

### Quick Start

Run the complete pipeline sequentially:

```bash
# 1. Prepare inputs (convert mmCIF to PDB, extract metrics)
python scripts/00_prepare_inputs_from_mmcif.py \
    --cmolgpt data/cmolgpt \
    --enamine data/enamine_real \
    --out data_std

# 2. Build samples table with IC50 normalization
python scripts/01_build_samples_long_ic50.py \
    --data_std data_std \
    --out outputs/samples_long.csv

# 3. Discover binding sites using DBSCAN
python scripts/02_detect_sites_unsupervised.py \
    --samples outputs/samples_long.csv \
    --out_samples outputs/site_samples.csv \
    --out_clusters outputs/site_clusters.csv \
    --eps 6.0 \
    --min_samples 8

# 4. Select best pose per molecule
python scripts/03_select_best_pose_unbiased.py \
    --samples outputs/site_samples.csv \
    --clusters outputs/site_clusters.csv \
    --out outputs/master_per_molecule.csv

# 5. Extract protein-ligand interactions (using PLIP)
python scripts/04_plip_contacts.py \
    --master outputs/master_per_molecule.csv \
    --out outputs/master_with_interactions.csv \
    --use_plip  # Omit if PLIP not installed

# 6. Score and rank molecules site-agnostically
python scripts/05_score_and_cluster_unbiased.py \
    --master outputs/master_with_interactions.csv \
    --out_csv outputs/ranked_shortlist.csv \
    --similarity_threshold 0.6
```

### Pipeline Scripts

#### Script 1: `00_prepare_inputs_from_mmcif.py`

Converts mmCIF files to PDB format and extracts prediction metrics from Boltz-2 JSON outputs.

**Inputs**:
- Raw Boltz-2 structure files (mmCIF)
- `boltz2_output.json` with affinity predictions

**Outputs**:
- `data_std/*/sample_*/complex.pdb` - PDB format structures
- `data_std/*/sample_*/metrics.json` - Extracted metrics (IC50, confidence, pIC50)

**Key metrics extracted**:
- `affinity_pred_value`: Predicted IC50 value
- `affinity_probability`: Binding probability
- `structure_confidence`: Overall structure confidence
- `protein_confidence`: Protein-specific confidence

#### Script 2: `01_build_samples_long_ic50.py`

Creates long-format table of all samples with normalized IC50 values and pIC50 calculations.

**Normalization**:
- Converts all IC50 values to Molarity (M)
- Unit conversions: nM→1e-9, μM→1e-6, mM→1e-3
- Calculates pIC50 = -log10(IC50_M)

**Output**: `outputs/samples_long.csv`

Columns: `molecule_id`, `sample_id`, `source`, `complex_pdb`, `affinity_probability`, `ic50_raw`, `ic50_unit`, `ic50_M`, `pIC50`, `structure_confidence`, `protein_confidence`

#### Script 3: `02_detect_sites_unsupervised.py`

Discovers binding sites using DBSCAN clustering on ligand centroids.

**Algorithm**:
1. Compute ligand center of mass for each sample
2. Identify contact residues (within 5Å)
3. Cluster centroids using DBSCAN (default: eps=6.0Å, min_samples=8)
4. Generate site statistics

**Parameters**:
- `--eps`: DBSCAN epsilon (distance threshold in Angstroms)
- `--min_samples`: Minimum samples to form a cluster

**Outputs**:
- `outputs/site_samples.csv` - Samples with site assignments
- `outputs/site_clusters.csv` - Site statistics (median pIC50, sample counts)

#### Script 4: `03_select_best_pose_unbiased.py`

Selects one best pose per molecule using deterministic cascade.

**Selection cascade** (strict ordering):
1. **Maximize** `affinity_probability` (binding confidence)
2. **If tied** (Δp ≤ 0.05): Maximize `pIC50` (potency)
3. **If tied** (ΔpIC50 ≤ 0.30): Maximize `structure_confidence`
4. **If still tied**: Choose sample from site with highest `median_pIC50`

**Output**: `outputs/master_per_molecule.csv` - One row per molecule

#### Script 5: `04_plip_contacts.py`

Extracts protein-ligand interactions using PLIP or fallback method.

**Interaction types counted**:
- Hydrogen bonds (`hbonds`)
- Hydrophobic contacts (`hydrophobics`)
- Salt bridges (`saltbridges`)
- Pi-stacking (`pi_stacks`)
- Pi-cation (`pi_cations`)
- Halogen bonds (`halogens`)
- Water bridges (`waterbridges`)
- Metal coordination (`metals`)

**Output**: `outputs/master_with_interactions.csv`

**Note**: Requires PLIP installation for accurate results. Use `--use_plip` flag if available.

#### Script 6: `05_score_and_cluster_unbiased.py`

Final site-agnostic ranking with diversity enforcement.

**Scoring formula**:
```
composite_score = 0.40*S_pot + 0.35*S_int + 0.25*S_conf - 0.20*(P_pains + P_noise)
```

Where:
- **S_pot**: Normalized pIC50 (potency)
- **S_int**: Normalized interaction count (hbonds + hydrophobics + pi_stacks)
- **S_conf**: Normalized confidence (0.6×affinity_probability + 0.4×structure_confidence)
- **P_pains**: PAINS penalty (1 if flagged, 0 otherwise)
- **P_noise**: Noise site penalty (1 if site_cluster_id=-1, 0 otherwise)

**Diversity clustering**:
- ECFP4 fingerprint similarity (Tanimoto > 0.6)
- Provides both global and within-cluster rankings

**Output**: `outputs/ranked_shortlist.csv`

Columns include: `global_rank`, `cluster_rank`, `composite_score`, `S_pot`, `S_int`, `S_conf`, `diversity_cluster`, plus all previous columns.

## Output Files

### Final Outputs

1. **`ranked_shortlist.csv`** - Primary deliverable
   - Ranked list of all molecules
   - Composite scores and normalized subscores
   - Site assignments and diversity clusters
   - Ready for experimental validation prioritization

2. **`site_clusters.csv`** - Site characterization
   - Statistics for each discovered binding site
   - Median/mean pIC50 per site
   - Sample and molecule counts
   - Centroid coordinates

### Intermediate Files

- `samples_long.csv`: All structure predictions with pIC50
- `site_samples.csv`: Samples with site assignments
- `master_per_molecule.csv`: Best pose per molecule
- `master_with_interactions.csv`: With PLIP interaction counts

## Key Concepts

### Site-Agnostic Ranking

Unlike traditional workflows that prioritize specific binding sites (e.g., ATP pocket), this pipeline:
- Discovers sites purely from data
- Ranks molecules independently of binding location
- Applies minimal penalty to "noise" (uncharacterized) sites

### Deterministic Pose Selection

The cascade ensures:
- Reproducibility (same input → same output)
- Interpretable criteria
- Balanced consideration of affinity, potency, and structure quality

### Diversity Enforcement

Chemical clustering prevents:
- Over-representation of similar scaffolds
- Redundant experimental validation
- Bias toward prolific chemical series

## Mathematical Details

### IC50 Normalization
```
IC50_M = IC50_value × unit_factor
pIC50 = -log10(IC50_M)
```

### Min-Max Normalization
```
x̃ = (x - min(x)) / (max(x) - min(x))
```

### DBSCAN Parameters
- **eps**: Maximum distance between points in a cluster (Angstroms)
- **min_samples**: Minimum points to form a dense region

### Tanimoto Similarity
```
T(A,B) = |A ∩ B| / |A ∪ B|
```

## Customization

### Adjust Scoring Weights

Edit `05_score_and_cluster_unbiased.py`:

```python
df = calculate_composite_score(
    df,
    w_pot=0.40,   # Potency weight
    w_int=0.35,   # Interaction weight
    w_conf=0.25,  # Confidence weight
    p_pains=0.20, # PAINS penalty
    p_noise=0.20  # Noise site penalty
)
```

### Modify Site Discovery

Adjust DBSCAN parameters:

```bash
python scripts/02_detect_sites_unsupervised.py \
    --eps 8.0 \            # Larger clusters
    --min_samples 5        # Fewer samples required
```

### Change Diversity Threshold

```bash
python scripts/05_score_and_cluster_unbiased.py \
    --similarity_threshold 0.7  # More strict clustering
```

## Troubleshooting

### PLIP Installation Issues

If PLIP fails to install:
```bash
# Use conda
conda install -c conda-forge plip

# Or run without PLIP
python scripts/04_plip_contacts.py --master ... --out ...
# (omit --use_plip flag for fallback method)
```

### Memory Issues with Large Datasets

Process in batches by splitting input directories or reduce parallelism in MDAnalysis operations.

### Missing mmCIF Files

Ensure Boltz-2 outputs are complete:
```bash
# Check for missing structures
find data/cmolgpt -name "structure_*.mmcif" | wc -l
# Should be 3× number of molecules
```

## Citation

If you use this pipeline, please cite:

```
[Your paper citation here]
Boltz-2: [Boltz-2 citation]
PLIP: Salentin et al., Nucleic Acids Res, 2015
```

## License

 MIT, Apache 2.0

## Contact

For questions or issues, please contact [your email] or open an issue on GitHub.

## Acknowledgments

- Boltz-2 for structure predictions
- PLIP for interaction profiling
- RDKit for cheminformatics utilities
