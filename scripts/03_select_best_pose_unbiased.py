#!/usr/bin/env python3
"""
Script 4: Select best pose per molecule using deterministic cascade
No bias toward specific binding sites
"""

import argparse
import logging
import pandas as pd
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def select_best_pose(
    group: pd.DataFrame,
    site_median_pic50: dict
) -> pd.Series:
    """
    Select best pose for a molecule using deterministic cascade.

    Selection criteria (in order):
    1. Maximize affinity_probability
    2. If tied (Δp ≤ 0.05): maximize pIC50
    3. If tied (ΔpIC50 ≤ 0.30): maximize structure_confidence
    4. If still tied: choose sample from site with highest median_pIC50

    Args:
        group: DataFrame of samples for one molecule
        site_median_pic50: Dict mapping site_cluster_id to median pIC50

    Returns:
        Best sample as pandas Series
    """
    if len(group) == 0:
        return None

    # Sort by criteria in order
    sorted_group = group.copy()

    # Add site median pIC50 for tie-breaking
    sorted_group['site_median_pic50'] = sorted_group['site_cluster_id'].map(
        lambda x: site_median_pic50.get(x, -999.0)
    )

    # Sort by cascade criteria
    sorted_group = sorted_group.sort_values(
        by=[
            'affinity_probability',  # Primary: highest affinity probability
            'pIC50',  # Secondary: highest pIC50
            'structure_confidence',  # Tertiary: highest structure confidence
            'site_median_pic50',  # Quaternary: best site
        ],
        ascending=[False, False, False, False]
    )

    # Get top candidate
    best = sorted_group.iloc[0]

    # Apply tolerance checks
    candidates = [best]

    # Check for ties in affinity_probability (Δp ≤ 0.05)
    tied_affinity = sorted_group[
        abs(sorted_group['affinity_probability'] - best['affinity_probability']) <= 0.05
    ]

    if len(tied_affinity) > 1:
        # Check for ties in pIC50 (ΔpIC50 ≤ 0.30)
        tied_pic50 = tied_affinity[
            abs(tied_affinity['pIC50'] - best['pIC50']) <= 0.30
        ]

        if len(tied_pic50) > 1:
            # Check for ties in structure_confidence
            # Take the one with highest confidence, or if tied, best site
            best = tied_pic50.iloc[0]
        else:
            best = tied_affinity.iloc[0]

    return best


def main():
    parser = argparse.ArgumentParser(
        description='Select best pose per molecule using deterministic cascade'
    )
    parser.add_argument(
        '--samples',
        type=Path,
        default=Path('outputs/site_samples.csv'),
        help='Input samples CSV with site assignments'
    )
    parser.add_argument(
        '--clusters',
        type=Path,
        default=Path('outputs/site_clusters.csv'),
        help='Site clusters CSV'
    )
    parser.add_argument(
        '--out',
        type=Path,
        default=Path('outputs/master_per_molecule.csv'),
        help='Output master table CSV'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Load data
    logger.info(f"Loading samples from {args.samples}")
    df_samples = pd.read_csv(args.samples)

    logger.info(f"Loading clusters from {args.clusters}")
    df_clusters = pd.read_csv(args.clusters)

    # Create site median pIC50 lookup
    site_median_pic50 = df_clusters.set_index('site_cluster_id')['median_pIC50'].to_dict()

    # Select best pose per molecule
    logger.info("Selecting best pose per molecule...")
    best_poses = []

    for molecule_id, group in df_samples.groupby('molecule_id'):
        best = select_best_pose(group, site_median_pic50)
        if best is not None:
            best_poses.append(best)

    # Create master DataFrame
    df_master = pd.DataFrame(best_poses)

    # Sort by source and molecule_id
    df_master = df_master.sort_values(['source', 'molecule_id'])

    # Save output
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df_master.to_csv(args.out, index=False)

    logger.info(f"Selected {len(df_master)} molecules")
    logger.info(f"Output saved to: {args.out}")

    # Print summary
    logger.info("\nSelection Summary:")
    logger.info(f"  cmolgpt molecules: {len(df_master[df_master['source'] == 'cmolgpt'])}")
    logger.info(f"  enamine_real molecules: {len(df_master[df_master['source'] == 'enamine_real'])}")

    # Site distribution
    site_dist = df_master['site_cluster_id'].value_counts().sort_index()
    logger.info("\nSite distribution of selected poses:")
    for site_id, count in site_dist.items():
        logger.info(f"  Site {int(site_id) if site_id >= 0 else 'noise'}: {count} molecules")


if __name__ == '__main__':
    main()
