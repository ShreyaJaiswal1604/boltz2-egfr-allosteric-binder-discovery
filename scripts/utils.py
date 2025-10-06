#!/usr/bin/env python3
"""
Utility functions for molecular property calculations
PAINS filters, QED, SA score, ECFP4 fingerprints
"""

import logging
from typing import Dict, Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED
from rdkit.Chem.Scaffolds import MurckoScaffold

logger = logging.getLogger(__name__)


# PAINS substructure patterns (simplified set)
# Full list available at: https://github.com/rdkit/rdkit/tree/master/Data/Pains
PAINS_SMARTS = [
    'C1N(C)C(=O)N(C)C1=O',  # Methylated urea
    'C1C=CC(=O)C=C1',  # Quinone
    'S(=O)(=O)N',  # Sulfonamide (sometimes problematic)
    '[!#6;!#1]~[CX3](=[!#6;!#1])~[!#6;!#1]',  # Beta-dicarbonyl
    'N=N',  # Azo group
    'C=C-C=C',  # Conjugated double bonds
    'C#CC#C',  # Diyne
    '[N+](=O)[O-]',  # Nitro group (can be reactive)
    'C(=O)N(O)',  # Hydroxamic acid
    'C(=S)',  # Thiocarbonyl
]


def check_pains(smiles: str) -> bool:
    """
    Check if molecule contains PAINS (Pan-Assay Interference Compounds) substructures.

    Args:
        smiles: SMILES string

    Returns:
        True if PAINS detected, False otherwise
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return True  # Invalid molecule, flag as problematic

        for smarts in PAINS_SMARTS:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                return True

        return False

    except Exception as e:
        logger.error(f"PAINS check failed for {smiles}: {e}")
        return True


def calculate_qed(smiles: str) -> Optional[float]:
    """
    Calculate Quantitative Estimate of Drug-likeness (QED).

    Args:
        smiles: SMILES string

    Returns:
        QED score (0-1) or None if failed
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        return QED.qed(mol)

    except Exception as e:
        logger.error(f"QED calculation failed for {smiles}: {e}")
        return None


def calculate_sa_score(smiles: str) -> Optional[float]:
    """
    Calculate Synthetic Accessibility (SA) score.
    Uses RDKit's implementation.

    Args:
        smiles: SMILES string

    Returns:
        SA score (1=easy to synthesize, 10=difficult) or None if failed
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Simplified SA score based on molecular complexity
        # Real implementation would use fragment contributions
        # This is a placeholder approximation

        # Count rings, rotatable bonds, complexity indicators
        n_rings = Descriptors.RingCount(mol)
        n_rot_bonds = Descriptors.NumRotatableBonds(mol)
        n_hetero = Descriptors.NumHeteroatoms(mol)
        n_atoms = mol.GetNumHeavyAtoms()

        # Simple heuristic (not actual SA score algorithm)
        complexity = (
            n_rings * 0.5 +
            n_rot_bonds * 0.3 +
            (n_hetero / max(n_atoms, 1)) * 2.0
        )

        # Normalize to 1-10 scale (1=easy, 10=hard)
        sa_score = min(10.0, max(1.0, 1.0 + complexity))

        return sa_score

    except Exception as e:
        logger.error(f"SA score calculation failed for {smiles}: {e}")
        return None


def calculate_molecular_weight(smiles: str) -> Optional[float]:
    """
    Calculate molecular weight.

    Args:
        smiles: SMILES string

    Returns:
        Molecular weight or None if failed
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        return Descriptors.MolWt(mol)

    except Exception as e:
        logger.error(f"MW calculation failed for {smiles}: {e}")
        return None


def calculate_ecfp4(smiles: str, radius: int = 2, n_bits: int = 2048) -> Optional[np.ndarray]:
    """
    Calculate ECFP4 (Morgan) fingerprint.

    Args:
        smiles: SMILES string
        radius: Fingerprint radius (default 2 for ECFP4)
        n_bits: Number of bits in fingerprint

    Returns:
        Binary fingerprint as numpy array or None if failed
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        arr = np.zeros((n_bits,), dtype=int)
        AllChem.DataStructs.ConvertToNumpyArray(fp, arr)

        return arr

    except Exception as e:
        logger.error(f"ECFP4 calculation failed for {smiles}: {e}")
        return None


def calculate_all_properties(smiles: str) -> Dict:
    """
    Calculate all molecular properties.

    Args:
        smiles: SMILES string

    Returns:
        Dictionary with all properties
    """
    properties = {
        'is_pains': check_pains(smiles),
        'qed': calculate_qed(smiles),
        'sa_score': calculate_sa_score(smiles),
        'molecular_weight': calculate_molecular_weight(smiles),
        'ecfp4': calculate_ecfp4(smiles),
    }

    return properties


def smiles_from_pdb(pdb_file: str) -> Optional[str]:
    """
    Extract ligand SMILES from PDB file.
    This is a placeholder - would need proper implementation
    to extract ligand coordinates and convert to SMILES.

    Args:
        pdb_file: Path to PDB file

    Returns:
        SMILES string or None
    """
    # This would require:
    # 1. Extract ligand from PDB
    # 2. Use Open Babel or similar to convert to SMILES
    # For now, return None as placeholder
    logger.warning("SMILES extraction from PDB not implemented - using placeholder")
    return None


def tanimoto_similarity(fp1: np.ndarray, fp2: np.ndarray) -> float:
    """
    Calculate Tanimoto similarity between two fingerprints.

    Args:
        fp1: First fingerprint
        fp2: Second fingerprint

    Returns:
        Tanimoto similarity (0-1)
    """
    if fp1 is None or fp2 is None:
        return 0.0

    intersection = np.sum(fp1 & fp2)
    union = np.sum(fp1 | fp2)

    if union == 0:
        return 0.0

    return intersection / union


if __name__ == '__main__':
    # Test examples
    test_smiles = [
        'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',  # Ibuprofen
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # Caffeine
    ]

    for smiles in test_smiles:
        print(f"\nSMILES: {smiles}")
        props = calculate_all_properties(smiles)
        for key, value in props.items():
            if key != 'ecfp4':
                print(f"  {key}: {value}")
