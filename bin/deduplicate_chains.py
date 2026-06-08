#!/usr/bin/env python3
"""
Deduplicate protein/RNA chains across all input FASTA files.

Reads all FASTA files from a samplesheet, extracts protein and RNA sequences,
deduplicates by sequence content, and outputs:
  - Per-chain FASTA files for each unique protein/RNA sequence
  - A chain_map.json recording which chains belong to which sample

DNA, CCD ligands, SMILES ligands, and ions are NOT deduplicated
(they don't need MSA search) — they are recorded in the map as-is.
"""

import argparse
import hashlib
import json
import os
import re
import string
from Bio import SeqIO


VALID_CHAIN_IDS = list(string.ascii_uppercase) + list(string.ascii_lowercase) + [str(x) for x in range(0, 10)]


def parse_args(args=None):
    Description = "Deduplicate protein/RNA chains across input FASTA files for AlphaFold3."
    Epilog = "Example: deduplicate_chains.py --input fastas.tsv --output_dir unique_chains/ --chain_map_out chain_map.json"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument(
        "--input",
        required=True,
        help="TSV file with columns: sample_id\\tfasta_path"
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to write unique chain FASTA files"
    )
    parser.add_argument(
        "--chain_map_out",
        required=True,
        help="Path to write the chain_map.json file"
    )

    return parser.parse_args(args)


def infer_entity_type(header, sequence):
    """
    Infer the entity type from the FASTA header and sequence content.

    Returns one of: "protein", "rna", "dna", "ccd", "smiles", "unknown"
    """
    ENTITY_TYPES = ["protein", "ccd", "smiles", "dna", "rna"]

    header_lower = header.lower()
    for entity in ENTITY_TYPES:
        if entity in header_lower:
            return entity

    seq = sequence.strip()
    seq_set = set(seq)

    # RNA: only A,C,U,G,N
    if len(seq_set - set("ACUGN")) == 0:
        return "rna"
    # DNA: only A,C,T,G,N
    if len(seq_set - set("ACTGN")) == 0:
        return "dna"
    # Protein: only 20 AA, not just A,C,T,G,U,N
    protein_letters = set("ACDEFGHIKLMNPQRSTVWY")
    if len(seq_set - protein_letters) == 0 and not (seq_set <= set("ACUGTN")):
        return "protein"
    # SMILES: fallback
    if re.fullmatch(r"[A-Za-z0-9@+\\-\\[\\]\\(\\)=#\\\$%]+", seq):
        return "smiles"
    return "unknown"


def make_chain_hash(entity_type, sequence):
    """
    Create a deterministic hash key for a chain.

    For protein and RNA, the hash is based on entity type + sequence.
    DNA, CCD, and SMILES chains are not deduplicated, so they get
    unique hashes that include a counter to prevent collisions.
    """
    if entity_type in ("protein", "rna"):
        return f"{entity_type}_{sequence}"
    else:
        # Non-deduplicatable types get a unique key via content hash
        # (these will never collide with each other in unique_chains
        #  because we don't deduplicate them)
        return None


def read_input_tsv(input_path):
    """
    Read the input TSV file and return a list of (sample_id, fasta_path) tuples.
    """
    samples = []
    with open(input_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                raise ValueError(f"Invalid input line (expected 2 tab-separated columns): {line}")
            samples.append((parts[0], parts[1]))
    return samples


def parse_fasta_chains(fasta_path):
    """
    Parse a FASTA file and return a list of chain descriptors.

    Each descriptor is a dict with:
      - entity_type: str
      - sequence: str
      - header: str (original FASTA header)
      - chain_id: str (assigned from VALID_CHAIN_IDS based on position)
    """
    chains = []
    for i, record in enumerate(SeqIO.parse(fasta_path, "fasta")):
        sequence = str(record.seq)
        header = record.description
        entity_type = infer_entity_type(header, sequence)
        chains.append({
            "entity_type": entity_type,
            "sequence": sequence,
            "header": header,
            "chain_id": VALID_CHAIN_IDS[i] if i < len(VALID_CHAIN_IDS) else f"chain_{i}",
        })
    return chains


def main(args=None):
    args = parse_args(args)

    os.makedirs(args.output_dir, exist_ok=True)

    # Read input samples
    samples = read_input_tsv(args.input)

    # Parse all FASTA files and collect chains per sample
    unique_chains = {}  # hash -> {type, sequence, fasta_path}
    sample_chains = {}  # sample_id -> list of chain descriptors for the map

    for sample_id, fasta_path in samples:
        chains = parse_fasta_chains(fasta_path)
        sample_chain_entries = []

        for chain in chains:
            entity_type = chain["entity_type"]
            sequence = chain["sequence"]
            chain_id = chain["chain_id"]

            chain_hash = make_chain_hash(entity_type, sequence)

            if chain_hash is not None and chain_hash not in unique_chains:
                # New unique protein/RNA chain — write a FASTA for it
                fasta_filename = f"unique_{entity_type}_{len(unique_chains)}.fasta"
                fasta_out_path = os.path.join(args.output_dir, fasta_filename)
                with open(fasta_out_path, 'w') as f:
                    f.write(f">{entity_type}\n{sequence}\n")
                unique_chains[chain_hash] = {
                    "type": entity_type,
                    "sequence": sequence,
                    "fasta_path": fasta_filename,
                    "data_json_path": None,  # filled after data pipeline runs
                }

            # Build the chain entry for this sample's map
            if chain_hash is not None:
                entry = {
                    "entity_type": entity_type,
                    "hash": chain_hash,
                    "chain_id": chain_id,
                }
            elif entity_type == "dna":
                entry = {
                    "entity_type": "dna",
                    "sequence": sequence,
                    "chain_id": chain_id,
                }
            elif entity_type == "ccd":
                entry = {
                    "entity_type": "ccd",
                    "ccdCodes": [sequence],
                    "chain_id": chain_id,
                }
            elif entity_type == "smiles":
                entry = {
                    "entity_type": "smiles",
                    "smiles": sequence,
                    "chain_id": chain_id,
                }
            else:
                entry = {
                    "entity_type": entity_type,
                    "sequence": sequence,
                    "chain_id": chain_id,
                }

            sample_chain_entries.append(entry)

        sample_chains[sample_id] = {
            "model_seeds": [11],  # default, overridden by pipeline params
            "chains": sample_chain_entries,
        }

    # Write chain_map.json
    chain_map = {
        "unique_chains": unique_chains,
        "samples": sample_chains,
    }

    with open(args.chain_map_out, 'w') as f:
        json.dump(chain_map, f, indent=2)

    print(f"Deduplicated {sum(len(v['chains']) for v in sample_chains.values())} total chains "
          f"across {len(sample_chains)} samples into {len(unique_chains)} unique protein/RNA chains")


if __name__ == "__main__":
    main()
