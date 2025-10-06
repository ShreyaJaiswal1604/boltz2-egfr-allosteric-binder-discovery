#!/usr/bin/env python3
"""
Script 2: Build samples table with IC50 normalization
Traverses standardized data and creates long-format table with pIC50
"""

import argparse
import json
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Optional

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Unit conversion factors to Molarity (M)
UNIT_CONVERSIONS = {
    'M': 1.0,
    'mM': 1e-3,
    'μM': 1e-6,
    'uM': 1e-6,
    'nM': 1e-9,
    'pM': 1e-12,
}


def normalize_ic50(ic50_value: float, ic50_unit: str) -> Optional[float]:
    """
    Normalize IC50 to Molarity.

    Args:
        ic50_value: IC50 value
        ic50_unit: Unit string

    Returns:
        IC50 in Molarity (M) or None if invalid
    """
    if ic50_value is None or ic50_value <= 0:
        return None

    factor = UNIT_CONVERSIONS.get(ic50_unit)
    if factor is None:
        # Try smart inference or default to nM
        logger.warning(f"Unknown unit '{ic50_unit}', defaulting to nM")
        factor = 1e-9

    return ic50_value * factor


def calculate_pic50(ic50_M: Optional[float]) -> Optional[float]:
    """
    Calculate pIC50 = -log10(IC50_M).

    Args:
        ic50_M: IC50 in Molarity

    Returns:
        pIC50 value or None
    """
    if ic50_M is None or ic50_M <= 0:
        return None

    return -np.log10(ic50_M)


def extract_sample_data(
    sample_dir: Path,
    molecule_id: str,
    source: str
) -> Optional[Dict]:
    """
    Extract data for a single sample.

    Args:
        sample_dir: Path to sample directory
        molecule_id: Molecule identifier
        source: Data source (cmolgpt or enamine_real)

    Returns:
        Dictionary with sample data or None
    """
    metrics_file = sample_dir / 'metrics.json'
    pdb_file = sample_dir / 'complex.pdb'

    if not metrics_file.exists() or not pdb_file.exists():
        return None

    try:
        with open(metrics_file, 'r') as f:
            metrics = json.load(f)

        sample_id = sample_dir.name  # e.g., sample_1

        # Extract fields
        ic50_raw = metrics.get('affinity_pred_value')
        ic50_unit = metrics.get('affinity_pred_unit', 'nM')
        affinity_prob = metrics.get('affinity_probability')
        structure_conf = metrics.get('structure_confidence')
        protein_conf = metrics.get('protein_confidence')

        # Normalize IC50
        ic50_M = normalize_ic50(ic50_raw, ic50_unit)
        pic50 = calculate_pic50(ic50_M)

        return {
            'molecule_id': molecule_id,
            'sample_id': f"{molecule_id}_{sample_id}",
            'source': source,
            'complex_pdb': str(pdb_file.resolve()),
            'affinity_probability': affinity_prob,
            'ic50_raw': ic50_raw,
            'ic50_unit': ic50_unit,
            'ic50_M': ic50_M,
            'pIC50': pic50,
            'structure_confidence': structure_conf,
            'protein_confidence': protein_conf,
        }

    except Exception as e:
        logger.error(f"Failed to process {sample_dir}: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description='Build samples table with IC50 normalization'
    )
    parser.add_argument(
        '--data_std',
        type=Path,
        default=Path('data_std'),
        help='Standardized data directory'
    )
    parser.add_argument(
        '--out',
        type=Path,
        default=Path('outputs/samples_long.csv'),
        help='Output CSV file'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Collect all samples
    samples = []

    for source in ['cmolgpt', 'enamine_real']:
        source_dir = args.data_std / source

        if not source_dir.exists():
            logger.warning(f"Source directory not found: {source_dir}")
            continue

        logger.info(f"Processing {source}...")

        mol_dirs = sorted([d for d in source_dir.iterdir() if d.is_dir()])

        for mol_dir in mol_dirs:
            molecule_id = mol_dir.name

            # Process each sample
            for sample_num in [1, 2, 3]:
                sample_dir = mol_dir / f'sample_{sample_num}'

                if not sample_dir.exists():
                    continue

                sample_data = extract_sample_data(sample_dir, molecule_id, source)

                if sample_data:
                    samples.append(sample_data)

    # Create DataFrame
    df = pd.DataFrame(samples)

    if len(df) == 0:
        logger.error("No samples found!")
        return

    # Sort by source and molecule_id
    df = df.sort_values(['source', 'molecule_id', 'sample_id'])

    # Save to CSV
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    logger.info(f"Processed {len(df)} samples")
    logger.info(f"Output saved to: {args.out}")

    # Print summary statistics
    logger.info("\nSummary Statistics:")
    logger.info(f"  cmolgpt samples: {len(df[df['source'] == 'cmolgpt'])}")
    logger.info(f"  enamine_real samples: {len(df[df['source'] == 'enamine_real'])}")
    logger.info(f"  pIC50 range: {df['pIC50'].min():.2f} - {df['pIC50'].max():.2f}")
    logger.info(f"  Mean affinity probability: {df['affinity_probability'].mean():.3f}")


if __name__ == '__main__':
    main()
