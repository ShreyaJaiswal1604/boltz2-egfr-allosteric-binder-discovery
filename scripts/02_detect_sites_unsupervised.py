#!/usr/bin/env python3
"""
Script 3: Detect binding sites using unsupervised clustering
Uses DBSCAN on ligand centroids without structural priors
"""

import argparse
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Optional

import MDAnalysis as mda
from sklearn.cluster import DBSCAN

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def compute_ligand_centroid(pdb_file: Path) -> Optional[np.ndarray]:
    """
    Compute centroid of ligand atoms.

    Args:
        pdb_file: Path to PDB file

    Returns:
        3D centroid coordinates or None if failed
    """
    try:
        # Use format='PDB' explicitly and ignore warnings
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(str(pdb_file), format='PDB', guess_bonds=False)

        # Select ligand residues - specifically target chains B and C (renamed ATP and L)
        # Chain A is protein, chains B and C are ligands
        ligand_sel = u.select_atoms("segid B or segid C or chainID B or chainID C")

        if len(ligand_sel) == 0:
            # Fallback: try selecting by resname
            ligand_sel = u.select_atoms("resname LIG1 or resname LIG2 or resname LIG or resname L or resname ATP")

        if len(ligand_sel) == 0:
            # Last resort: select non-protein residues
            ligand_sel = u.select_atoms("not protein and not resname HOH and not resname WAT and not resname TIP3")

        if len(ligand_sel) == 0:
            logger.warning(f"No ligand found in {pdb_file}")
            return None

        centroid = ligand_sel.center_of_geometry()  # Use geometry instead of mass
        return centroid

    except Exception as e:
        logger.error(f"Failed to compute centroid for {pdb_file}: {e}")
        return None


def get_contact_residues(pdb_file: Path, cutoff: float = 5.0) -> List[int]:
    """
    Get residue numbers within cutoff distance of ligand.

    Args:
        pdb_file: Path to PDB file
        cutoff: Distance cutoff in Angstroms

    Returns:
        List of residue numbers in contact with ligand
    """
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(str(pdb_file), format='PDB', guess_bonds=False)

        # Select ligand - try chain-based selection first
        ligand_sel = u.select_atoms("segid B or segid C or chainID B or chainID C")

        if len(ligand_sel) == 0:
            ligand_sel = u.select_atoms("resname LIG1 or resname LIG2 or resname LIG or resname L or resname ATP")

        if len(ligand_sel) == 0:
            ligand_sel = u.select_atoms("not protein and not resname HOH and not resname WAT and not resname TIP3")

        if len(ligand_sel) == 0:
            return []

        # Select protein atoms near ligand
        protein_sel = u.select_atoms("protein")
        contact_residues = set()

        for atom in protein_sel:
            try:
                distances = mda.lib.distances.distance_array(
                    atom.position[np.newaxis, :],
                    ligand_sel.positions,
                    box=None  # Don't use periodic boundary conditions
                )

                if distances.min() <= cutoff:
                    contact_residues.add(atom.resid)
            except:
                continue

        return sorted(list(contact_residues))

    except Exception as e:
        logger.error(f"Failed to get contact residues for {pdb_file}: {e}")
        return []


def cluster_sites(
    centroids: np.ndarray,
    eps: float = 6.0,
    min_samples: int = 8
) -> np.ndarray:
    """
    Cluster binding sites using DBSCAN.

    Args:
        centroids: Nx3 array of ligand centroids
        eps: DBSCAN epsilon parameter (Angstroms)
        min_samples: Minimum samples per cluster

    Returns:
        Cluster labels (-1 for noise)
    """
    if len(centroids) < min_samples:
        logger.warning(f"Too few centroids ({len(centroids)}) for clustering")
        return np.array([-1] * len(centroids))

    dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric='euclidean')
    labels = dbscan.fit_predict(centroids)

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = list(labels).count(-1)

    logger.info(f"DBSCAN results: {n_clusters} clusters, {n_noise} noise points")

    return labels


def main():
    parser = argparse.ArgumentParser(
        description='Detect binding sites using unsupervised clustering'
    )
    parser.add_argument(
        '--samples',
        type=Path,
        default=Path('outputs/samples_long.csv'),
        help='Input samples CSV file'
    )
    parser.add_argument(
        '--out_samples',
        type=Path,
        default=Path('outputs/site_samples.csv'),
        help='Output samples CSV with site assignments'
    )
    parser.add_argument(
        '--out_clusters',
        type=Path,
        default=Path('outputs/site_clusters.csv'),
        help='Output site clusters CSV'
    )
    parser.add_argument(
        '--eps',
        type=float,
        default=6.0,
        help='DBSCAN epsilon parameter (Angstroms)'
    )
    parser.add_argument(
        '--min_samples',
        type=int,
        default=8,
        help='DBSCAN minimum samples per cluster'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Load samples
    logger.info(f"Loading samples from {args.samples}")
    df = pd.read_csv(args.samples)

    # Compute centroids and contact residues
    logger.info("Computing ligand centroids and contact residues...")
    centroids = []
    contact_res_list = []
    valid_indices = []

    for idx, row in df.iterrows():
        pdb_path = Path(row['complex_pdb'])

        centroid = compute_ligand_centroid(pdb_path)
        contact_res = get_contact_residues(pdb_path, cutoff=5.0)

        if centroid is not None:
            centroids.append(centroid)
            contact_res_list.append(contact_res)
            valid_indices.append(idx)

    logger.info(f"Valid centroids: {len(centroids)}/{len(df)}")

    if len(centroids) == 0:
        logger.error("No valid centroids found!")
        return

    # Cluster sites
    logger.info("Clustering binding sites...")
    centroids_array = np.array(centroids)
    labels = cluster_sites(centroids_array, eps=args.eps, min_samples=args.min_samples)

    # Add cluster assignments to DataFrame
    df['site_cluster_id'] = -1  # Initialize all as noise
    df.loc[valid_indices, 'site_cluster_id'] = labels
    df['centroid_x'] = np.nan
    df['centroid_y'] = np.nan
    df['centroid_z'] = np.nan
    df.loc[valid_indices, 'centroid_x'] = centroids_array[:, 0]
    df.loc[valid_indices, 'centroid_y'] = centroids_array[:, 1]
    df.loc[valid_indices, 'centroid_z'] = centroids_array[:, 2]

    # Add contact residues as string
    df['contact_residues'] = ''
    df.loc[valid_indices, 'contact_residues'] = [
        ','.join(map(str, res)) for res in contact_res_list
    ]

    # Save samples with site assignments
    args.out_samples.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_samples, index=False)
    logger.info(f"Saved site samples to: {args.out_samples}")

    # Generate site cluster statistics
    cluster_stats = []
    unique_clusters = sorted(set(labels))

    for cluster_id in unique_clusters:
        if cluster_id == -1:
            continue  # Skip noise

        cluster_mask = df['site_cluster_id'] == cluster_id
        cluster_samples = df[cluster_mask]

        stats = {
            'site_cluster_id': cluster_id,
            'n_samples': len(cluster_samples),
            'n_molecules': cluster_samples['molecule_id'].nunique(),
            'median_pIC50': cluster_samples['pIC50'].median(),
            'mean_pIC50': cluster_samples['pIC50'].mean(),
            'std_pIC50': cluster_samples['pIC50'].std(),
            'mean_affinity_probability': cluster_samples['affinity_probability'].mean(),
            'mean_structure_confidence': cluster_samples['structure_confidence'].mean(),
            'centroid_x': centroids_array[labels == cluster_id, 0].mean(),
            'centroid_y': centroids_array[labels == cluster_id, 1].mean(),
            'centroid_z': centroids_array[labels == cluster_id, 2].mean(),
        }
        cluster_stats.append(stats)

    # Create cluster DataFrame
    df_clusters = pd.DataFrame(cluster_stats)
    df_clusters = df_clusters.sort_values('median_pIC50', ascending=False)

    df_clusters.to_csv(args.out_clusters, index=False)
    logger.info(f"Saved cluster statistics to: {args.out_clusters}")

    # Print summary
    logger.info("\nSite Cluster Summary:")
    logger.info(f"  Total clusters: {len(df_clusters)}")
    logger.info(f"  Noise samples: {len(df[df['site_cluster_id'] == -1])}")

    if len(df_clusters) > 0:
        logger.info("\nTop 5 sites by median pIC50:")
        for _, row in df_clusters.head(5).iterrows():
            logger.info(
                f"  Site {int(row['site_cluster_id'])}: "
                f"n={int(row['n_samples'])}, "
                f"median_pIC50={row['median_pIC50']:.2f}"
            )


if __name__ == '__main__':
    from typing import Optional
    main()
