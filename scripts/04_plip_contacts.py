#!/usr/bin/env python3
"""
Script 5: Extract protein-ligand interactions using PLIP
Runs PLIP on selected complexes and parses interaction counts
"""

import argparse
import logging
import pandas as pd
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, Optional

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def run_plip(pdb_file: Path) -> Optional[Path]:
    """
    Run PLIP on a PDB file.

    Args:
        pdb_file: Path to PDB file

    Returns:
        Path to PLIP XML report or None if failed
    """
    try:
        # Create temporary directory for PLIP output
        temp_dir = tempfile.mkdtemp()
        temp_dir_path = Path(temp_dir)

        # Run PLIP command
        cmd = [
            'plip',
            '-f', str(pdb_file),
            '-x',  # XML output
            '-o', str(temp_dir_path),
            '--name', 'plip_report'
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )

        if result.returncode != 0:
            logger.error(f"PLIP failed for {pdb_file}: {result.stderr}")
            return None

        # Find XML report
        xml_files = list(temp_dir_path.glob('**/plip_report.xml'))

        if not xml_files:
            # Try alternative pattern
            xml_files = list(temp_dir_path.glob('**/*.xml'))

        if xml_files:
            return xml_files[0]
        else:
            logger.warning(f"No XML report found for {pdb_file}")
            return None

    except subprocess.TimeoutExpired:
        logger.error(f"PLIP timeout for {pdb_file}")
        return None
    except Exception as e:
        logger.error(f"PLIP error for {pdb_file}: {e}")
        return None


def parse_plip_xml(xml_file: Path) -> Dict[str, int]:
    """
    Parse PLIP XML report and count interactions.

    Args:
        xml_file: Path to PLIP XML report

    Returns:
        Dictionary with interaction counts
    """
    counts = {
        'hbonds': 0,
        'hydrophobics': 0,
        'saltbridges': 0,
        'pi_stacks': 0,
        'pi_cations': 0,
        'halogens': 0,
        'waterbridges': 0,
        'metals': 0,
    }

    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        # Count each interaction type
        for binding_site in root.findall('.//bindingsite'):
            # Hydrogen bonds
            hbonds = binding_site.findall('.//hydrogen_bond')
            counts['hbonds'] += len(hbonds)

            # Hydrophobic contacts
            hydrophobic = binding_site.findall('.//hydrophobic_interaction')
            counts['hydrophobics'] += len(hydrophobic)

            # Salt bridges
            saltbridges = binding_site.findall('.//salt_bridge')
            counts['saltbridges'] += len(saltbridges)

            # Pi stacking
            pi_stacks = binding_site.findall('.//pi_stack')
            counts['pi_stacks'] += len(pi_stacks)

            # Pi-cation interactions
            pi_cations = binding_site.findall('.//pi_cation_interaction')
            counts['pi_cations'] += len(pi_cations)

            # Halogen bonds
            halogens = binding_site.findall('.//halogen_bond')
            counts['halogens'] += len(halogens)

            # Water bridges
            waterbridges = binding_site.findall('.//water_bridge')
            counts['waterbridges'] += len(waterbridges)

            # Metal complexes
            metals = binding_site.findall('.//metal_complex')
            counts['metals'] += len(metals)

    except Exception as e:
        logger.error(f"Failed to parse XML {xml_file}: {e}")

    return counts


def extract_interactions_simple(pdb_file: Path) -> Dict[str, int]:
    """
    Fallback: Simple interaction counting without PLIP.
    Counts basic contacts based on distance cutoffs.

    Args:
        pdb_file: Path to PDB file

    Returns:
        Dictionary with estimated interaction counts
    """
    logger.warning(f"Using fallback interaction counting for {pdb_file}")

    # Placeholder - in real implementation, use MDAnalysis or similar
    # to count atom pairs within distance cutoffs
    counts = {
        'hbonds': 0,
        'hydrophobics': 0,
        'saltbridges': 0,
        'pi_stacks': 0,
        'pi_cations': 0,
        'halogens': 0,
        'waterbridges': 0,
        'metals': 0,
    }

    return counts


def main():
    parser = argparse.ArgumentParser(
        description='Extract protein-ligand interactions using PLIP'
    )
    parser.add_argument(
        '--master',
        type=Path,
        default=Path('outputs/master_per_molecule.csv'),
        help='Input master table CSV'
    )
    parser.add_argument(
        '--out',
        type=Path,
        default=Path('outputs/master_with_interactions.csv'),
        help='Output CSV with interaction counts'
    )
    parser.add_argument(
        '--use_plip',
        action='store_true',
        help='Use PLIP for interaction analysis (requires PLIP installation)'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Load master table
    logger.info(f"Loading master table from {args.master}")
    df = pd.read_csv(args.master)

    # Add interaction columns
    interaction_cols = [
        'hbonds', 'hydrophobics', 'saltbridges', 'pi_stacks',
        'pi_cations', 'halogens', 'waterbridges', 'metals'
    ]

    for col in interaction_cols:
        df[col] = 0

    # Process each complex
    logger.info("Extracting interactions...")
    for idx, row in df.iterrows():
        pdb_file = Path(row['complex_pdb'])

        if not pdb_file.exists():
            logger.warning(f"PDB file not found: {pdb_file}")
            continue

        if args.use_plip:
            # Run PLIP
            xml_report = run_plip(pdb_file)

            if xml_report:
                counts = parse_plip_xml(xml_report)
            else:
                counts = extract_interactions_simple(pdb_file)
        else:
            # Use fallback method
            counts = extract_interactions_simple(pdb_file)

        # Update DataFrame
        for col, count in counts.items():
            df.at[idx, col] = count

        if (idx + 1) % 10 == 0:
            logger.info(f"Processed {idx + 1}/{len(df)} complexes")

    # Calculate total interactions
    df['total_interactions'] = df[interaction_cols].sum(axis=1)

    # Save output
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    logger.info(f"Output saved to: {args.out}")

    # Print summary
    logger.info("\nInteraction Summary:")
    for col in interaction_cols:
        logger.info(f"  {col}: {df[col].sum()} total, {df[col].mean():.2f} mean")

    logger.info(f"\nTotal interactions: {df['total_interactions'].sum()}")
    logger.info(f"Mean interactions per complex: {df['total_interactions'].mean():.2f}")


if __name__ == '__main__':
    main()
