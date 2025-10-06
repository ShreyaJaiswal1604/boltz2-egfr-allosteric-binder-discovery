#!/usr/bin/env python3
"""
Script 1: Prepare inputs from mmCIF files
Converts .mmcif to .pdb and extracts metrics from boltz2_output.json
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

from Bio import PDB

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def fix_pdb_file(input_pdb: Path, output_pdb: Path) -> bool:
    """
    Fix PDB file formatting issues to be compatible with MDAnalysis.

    Args:
        input_pdb: Input PDB file path
        output_pdb: Output PDB file path

    Returns:
        True if successful
    """
    try:
        with open(input_pdb, 'r') as f:
            lines = f.readlines()

        fixed_lines = []
        for line in lines:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                # PDB format specification:
                # Extract fields
                record = line[0:6]
                serial = line[6:11]
                atom_name = line[12:16]
                alt_loc = line[16:17] if len(line) > 16 else ' '
                res_name = line[17:20] if len(line) > 19 else '   '
                chain_id = line[21:22] if len(line) > 21 else ' '
                res_seq = line[22:26] if len(line) > 25 else '    '
                i_code = line[26:27] if len(line) > 26 else ' '

                # Try to parse coordinates
                try:
                    x_str = line[30:38].strip() if len(line) > 37 else '0.000'
                    y_str = line[38:46].strip() if len(line) > 45 else '0.000'
                    z_str = line[46:54].strip() if len(line) > 53 else '0.000'

                    x = float(x_str)
                    y = float(y_str)
                    z = float(z_str)

                    occupancy = line[54:60].strip() if len(line) > 59 else '1.00'
                    temp_factor = line[60:66].strip() if len(line) > 65 else '0.00'
                    element = line[76:78].strip() if len(line) > 77 else '  '

                    # Reformat the line correctly
                    new_line = f"{record:6s}{serial:>5s} {atom_name:4s}{alt_loc:1s}{res_name:3s} {chain_id:1s}{res_seq:>4s}{i_code:1s}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:>6s}{temp_factor:>6s}          {element:>2s}  \n"
                    fixed_lines.append(new_line)
                except ValueError:
                    # If parsing fails, keep original
                    fixed_lines.append(line)
            else:
                fixed_lines.append(line)

        with open(output_pdb, 'w') as f:
            f.writelines(fixed_lines)

        return True

    except Exception as e:
        logger.error(f"Error fixing {input_pdb}: {e}")
        return False


def convert_mmcif_to_pdb_fixed(mmcif_path: Path, pdb_path: Path) -> bool:
    """
    Convert mmCIF to PDB with chain ID renaming and residue name truncation.
    Renames chains: ATP->B, L->C to fit PDB 1-char chain ID limit.
    Truncates residue names to 3 characters to fit PDB format.

    Args:
        mmcif_path: Path to input mmCIF file
        pdb_path: Path to output PDB file

    Returns:
        True if successful, False otherwise
    """
    try:
        parser = PDB.MMCIFParser(QUIET=True)
        structure = parser.get_structure('complex', str(mmcif_path))

        # Rename chain IDs and fix residue names
        chain_mapping = {'ATP': 'B', 'L': 'C'}
        for model in structure:
            for chain in model:
                # Rename chains
                if chain.id in chain_mapping:
                    chain.id = chain_mapping[chain.id]

                # Fix residue names - truncate to 3 characters for PDB format
                for residue in chain.get_residues():
                    resname = residue.resname
                    if len(resname) > 3:
                        # Truncate to 3 characters (e.g., LIG1 -> LIG)
                        residue.resname = resname[:3]

        # Save as PDB with proper formatting
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(str(pdb_path))

        logger.debug(f"Converted {mmcif_path.name} to PDB (renamed chains)")
        return True
    except Exception as e:
        logger.error(f"Failed to convert {mmcif_path}: {e}")
        return False


def extract_metrics(boltz_output: Path, sample_num: int) -> Optional[Dict]:
    """
    Extract metrics from boltz2_output.json for a specific sample.

    Args:
        boltz_output: Path to boltz2_output.json
        sample_num: Sample number (1-indexed)

    Returns:
        Dictionary with metrics or None if failed
    """
    try:
        with open(boltz_output, 'r') as f:
            data = json.load(f)

        # Extract affinity data for ligand 'L'
        affinities = data.get('affinities', {})
        ligand_affinity = affinities.get('L', {})

        # Get confidence scores (0-indexed for sample_num)
        idx = sample_num - 1
        confidence_scores = data.get('confidence_scores', [])
        protein_conf = data.get('protein_iptm_scores', [])

        # Extract affinity values
        # Using ensemble affinity (averaged across models)
        affinity_pred = ligand_affinity.get('affinity_pred_value', [None])[0]
        affinity_prob = ligand_affinity.get('affinity_probability_binary', [None])[0]
        pic50 = ligand_affinity.get('affinity_pic50', [None])[0]

        # Structure confidence for this sample
        structure_conf = confidence_scores[idx] if idx < len(confidence_scores) else None
        protein_conf_val = protein_conf[idx] if idx < len(protein_conf) else None

        metrics = {
            'affinity_pred_value': affinity_pred,  # IC50 value
            'affinity_pred_unit': 'nM',  # Boltz-2 uses nM by default
            'affinity_probability': affinity_prob,
            'structure_confidence': structure_conf,
            'protein_confidence': protein_conf_val,
            'pic50': pic50
        }

        logger.debug(f"Extracted metrics for sample {sample_num}: {metrics}")
        return metrics

    except Exception as e:
        logger.error(f"Failed to extract metrics from {boltz_output}: {e}")
        return None


def process_molecule(mol_dir: Path, output_base: Path, source: str) -> int:
    """
    Process a single molecule directory.

    Args:
        mol_dir: Path to molecule directory
        output_base: Base output directory
        source: Source type (cmolgpt or enamine_real)

    Returns:
        Number of successfully processed samples
    """
    mol_id = mol_dir.name
    structures_dir = mol_dir / 'structures'
    boltz_output = mol_dir / 'boltz2_output.json'

    if not structures_dir.exists():
        logger.warning(f"No structures directory for {mol_id}")
        return 0

    if not boltz_output.exists():
        logger.warning(f"No boltz2_output.json for {mol_id}")
        return 0

    # Output directory for this molecule
    mol_output_dir = output_base / source / mol_id

    processed_count = 0

    # Process each structure sample
    for sample_num in [1, 2, 3]:
        mmcif_file = structures_dir / f'structure_{sample_num}.mmcif'

        if not mmcif_file.exists():
            logger.warning(f"Missing {mmcif_file}")
            continue

        # Create sample output directory
        sample_dir = mol_output_dir / f'sample_{sample_num}'
        sample_dir.mkdir(parents=True, exist_ok=True)

        # Convert mmCIF to PDB with chain renaming
        pdb_file_temp = sample_dir / 'complex_temp.pdb'
        pdb_file = sample_dir / 'complex.pdb'

        if not convert_mmcif_to_pdb_fixed(mmcif_file, pdb_file_temp):
            continue

        # Fix PDB formatting for MDAnalysis compatibility
        if not fix_pdb_file(pdb_file_temp, pdb_file):
            logger.warning(f"Failed to fix PDB format for {pdb_file_temp}, using original")
            pdb_file_temp.rename(pdb_file)
        else:
            # Remove temporary file
            pdb_file_temp.unlink()

        # Extract metrics
        metrics = extract_metrics(boltz_output, sample_num)
        if metrics is None:
            continue

        # Save metrics
        metrics_file = sample_dir / 'metrics.json'
        with open(metrics_file, 'w') as f:
            json.dump(metrics, f, indent=2)

        processed_count += 1

    if processed_count > 0:
        logger.info(f"Processed {mol_id}: {processed_count}/3 samples")

    return processed_count


def main():
    parser = argparse.ArgumentParser(
        description='Prepare inputs from mmCIF files'
    )
    parser.add_argument(
        '--cmolgpt',
        type=Path,
        required=True,
        help='Path to cmolgpt data directory'
    )
    parser.add_argument(
        '--enamine',
        type=Path,
        required=True,
        help='Path to enamine_real data directory'
    )
    parser.add_argument(
        '--out',
        type=Path,
        required=True,
        help='Output directory for standardized data'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Create output directory
    args.out.mkdir(parents=True, exist_ok=True)

    total_processed = 0

    # Process cmolgpt molecules
    logger.info("Processing cmolgpt molecules...")
    if args.cmolgpt.exists():
        mol_dirs = sorted([d for d in args.cmolgpt.iterdir() if d.is_dir()])
        for mol_dir in mol_dirs:
            total_processed += process_molecule(mol_dir, args.out, 'cmolgpt')
    else:
        logger.error(f"cmolgpt directory not found: {args.cmolgpt}")

    # Process enamine_real molecules
    logger.info("Processing enamine_real molecules...")
    if args.enamine.exists():
        mol_dirs = sorted([d for d in args.enamine.iterdir() if d.is_dir()])
        for mol_dir in mol_dirs:
            total_processed += process_molecule(mol_dir, args.out, 'enamine_real')
    else:
        logger.error(f"enamine_real directory not found: {args.enamine}")

    logger.info(f"Total samples processed: {total_processed}")
    logger.info(f"Output written to: {args.out}")


if __name__ == '__main__':
    main()
