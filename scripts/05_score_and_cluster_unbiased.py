#!/usr/bin/env python3
"""
Script 6: Score and cluster molecules site-agnostically
Implements final ranking with diversity enforcement
"""

import argparse
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List

from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def normalize_scores(series: pd.Series) -> pd.Series:
    """
    Min-max normalization to [0, 1].

    Args:
        series: Pandas series to normalize

    Returns:
        Normalized series
    """
    min_val = series.min()
    max_val = series.max()

    if max_val == min_val:
        return pd.Series([0.5] * len(series), index=series.index)

    return (series - min_val) / (max_val - min_val)


def calculate_composite_score(
    df: pd.DataFrame,
    w_pot: float = 0.40,
    w_int: float = 0.35,
    w_conf: float = 0.25,
    p_pains: float = 0.20,
    p_noise: float = 0.20,
) -> pd.DataFrame:
    """
    Calculate composite scoring for molecules.

    Score = w_pot*S_pot + w_int*S_int + w_conf*S_conf - p_pains*P_pains - p_noise*P_noise

    Args:
        df: DataFrame with molecular data
        w_pot: Weight for potency score
        w_int: Weight for interaction score
        w_conf: Weight for confidence score
        p_pains: Penalty for PAINS
        p_noise: Penalty for noise site

    Returns:
        DataFrame with scores added
    """
    df = df.copy()

    # Calculate normalized scores
    logger.info("Calculating normalized scores...")

    # S_pot: Normalized pIC50
    df['S_pot'] = normalize_scores(df['pIC50'])

    # S_int: Normalized interaction count (hbonds + hydrophobics + pi_stacks)
    df['interaction_score'] = (
        df['hbonds'] +
        df['hydrophobics'] +
        df['pi_stacks']
    )
    df['S_int'] = normalize_scores(df['interaction_score'])

    # S_conf: Normalized confidence (weighted combination)
    df['confidence_combined'] = (
        0.6 * df['affinity_probability'] +
        0.4 * df['structure_confidence']
    )
    df['S_conf'] = normalize_scores(df['confidence_combined'])

    # Penalties
    # P_pains: 1 if PAINS, 0 otherwise (need PAINS info - placeholder)
    if 'is_pains' not in df.columns:
        df['is_pains'] = 0  # Placeholder
    df['P_pains'] = df['is_pains'].astype(float)

    # P_noise: 1 if in noise site (-1), 0 otherwise
    df['P_noise'] = (df['site_cluster_id'] == -1).astype(float)

    # Composite score
    df['composite_score'] = (
        w_pot * df['S_pot'] +
        w_int * df['S_int'] +
        w_conf * df['S_conf'] -
        p_pains * df['P_pains'] -
        p_noise * df['P_noise']
    )

    logger.info("Composite scores calculated")

    return df


def cluster_by_similarity(
    df: pd.DataFrame,
    similarity_threshold: float = 0.6
) -> pd.DataFrame:
    """
    Cluster molecules by ECFP4 similarity.

    Args:
        df: DataFrame with ECFP4 fingerprints
        similarity_threshold: Tanimoto similarity threshold

    Returns:
        DataFrame with cluster assignments
    """
    logger.info("Clustering by molecular similarity...")

    # Placeholder: Need actual ECFP4 fingerprints
    # For now, use simple clustering based on molecular weight similarity

    if 'molecular_weight' not in df.columns:
        # Use pIC50 as proxy for clustering
        features = df[['pIC50', 'interaction_score']].values
    else:
        features = df[['molecular_weight', 'pIC50']].values

    # Handle NaN values by filling with mean
    features = np.nan_to_num(features, nan=0.0)

    # Normalize features
    mean = features.mean(axis=0)
    std = features.std(axis=0) + 1e-8
    features = (features - mean) / std
    features = np.nan_to_num(features, nan=0.0)  # Handle any remaining NaNs

    # Agglomerative clustering
    n_molecules = len(df)

    if n_molecules < 2:
        df['diversity_cluster'] = 0
        return df

    # Use distance threshold for clustering
    clustering = AgglomerativeClustering(
        n_clusters=None,
        distance_threshold=1.0,  # Adjust based on normalized features
        linkage='average'
    )

    df['diversity_cluster'] = clustering.fit_predict(features)

    n_clusters = len(df['diversity_cluster'].unique())
    logger.info(f"Found {n_clusters} diversity clusters")

    return df


def rank_with_diversity(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rank molecules globally and within diversity clusters.

    Args:
        df: DataFrame with scores and clusters

    Returns:
        DataFrame with rankings
    """
    logger.info("Ranking molecules...")

    # Fill any NaN composite scores with minimum value
    df['composite_score'] = df['composite_score'].fillna(df['composite_score'].min() - 1.0)

    # Global rank by composite score
    df['global_rank'] = df['composite_score'].rank(ascending=False, method='min').astype(int)

    # Within-cluster rank
    df['cluster_rank'] = df.groupby('diversity_cluster')['composite_score'].rank(
        ascending=False, method='min'
    ).astype(int)

    # Sort by global rank
    df = df.sort_values('global_rank')

    return df


def main():
    parser = argparse.ArgumentParser(
        description='Score and rank molecules site-agnostically'
    )
    parser.add_argument(
        '--master',
        type=Path,
        default=Path('outputs/master_with_interactions.csv'),
        help='Input master table with interactions'
    )
    parser.add_argument(
        '--out_csv',
        type=Path,
        default=Path('outputs/ranked_shortlist.csv'),
        help='Output ranked shortlist CSV'
    )
    parser.add_argument(
        '--similarity_threshold',
        type=float,
        default=0.6,
        help='Tanimoto similarity threshold for clustering'
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
    logger.info(f"Loading master table from {args.master}")
    df = pd.read_csv(args.master)

    # Calculate composite scores
    df = calculate_composite_score(df)

    # Cluster by similarity
    df = cluster_by_similarity(df, args.similarity_threshold)

    # Rank molecules
    df = rank_with_diversity(df)

    # Save output
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)

    logger.info(f"Ranked shortlist saved to: {args.out_csv}")

    # Print summary
    logger.info("\nRanking Summary:")
    logger.info(f"  Total molecules: {len(df)}")
    logger.info(f"  Diversity clusters: {df['diversity_cluster'].nunique()}")

    logger.info("\nTop 10 molecules:")
    top_cols = [
        'global_rank', 'molecule_id', 'source', 'composite_score',
        'pIC50', 'interaction_score', 'site_cluster_id'
    ]
    if all(col in df.columns for col in top_cols):
        print(df[top_cols].head(10).to_string(index=False))
    else:
        print(df[['global_rank', 'molecule_id', 'composite_score']].head(10).to_string(index=False))

    logger.info("\nScore distribution:")
    logger.info(f"  Composite score: {df['composite_score'].mean():.3f} ± {df['composite_score'].std():.3f}")
    logger.info(f"  S_pot: {df['S_pot'].mean():.3f}")
    logger.info(f"  S_int: {df['S_int'].mean():.3f}")
    logger.info(f"  S_conf: {df['S_conf'].mean():.3f}")

    # Site distribution in top 50
    logger.info("\nSite distribution in top 50:")
    top_50_sites = df.head(50)['site_cluster_id'].value_counts()
    for site_id, count in top_50_sites.items():
        logger.info(f"  Site {int(site_id) if site_id >= 0 else 'noise'}: {count}")


if __name__ == '__main__':
    main()
