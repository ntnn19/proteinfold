#!/usr/bin/env python3
"""
Assemble multimeric AlphaFold3 JSON files from pre-computed MSAs.

Reads chain_map.json and a directory of *_data.json files produced by the
AlphaFold3 data pipeline (--run_inference=false). For each sample in the
chain_map, assembles a complete AF3 input JSON with pre-filled MSAs and
templates for protein/RNA chains, and direct entries for DNA/ligand chains.

The output JSONs can be used as input for inference-only runs
(--run_data_pipeline=false).
"""

import argparse
import json
import os
import string


def parse_args(args=None):
    Description = "Assemble multimeric AlphaFold3 JSONs from pre-computed MSAs."
    Epilog = "Example: assemble_multimer_json.py --chain_map chain_map.json --data_json_dir data_outputs/ --output_dir assembled/"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument(
        "--chain_map",
        required=True,
        help="Path to chain_map.json from deduplicate_chains.py"
    )
    parser.add_argument(
        "--data_json_dir",
        required=True,
        help="Directory containing *_data.json files from the data pipeline"
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to write assembled multimer JSON files"
    )
    parser.add_argument(
        "--model_seeds",
        type=int,
        nargs='+',
        default=[11],
        help="Model seeds to use in assembled JSONs (default: [11])"
    )

    return parser.parse_args(args)


def sanitised_name(name):
    """
    Sanitize a name for use as an AF3 job name / filename.
    Matches the AlphaFold3 sanitisation logic.
    """
    lower_spaceless_name = name.lower().replace(' ', '_')
    allowed_chars = set(string.ascii_lowercase + string.digits + '_-.')
    return ''.join(l for l in lower_spaceless_name if l in allowed_chars)


def find_data_json(data_json_dir, chain_hash, unique_chains):
    """
    Find the *_data.json file for a given chain hash.

    The data pipeline writes to <output_dir>/<sanitised_name>/<sanitised_name>_data.json
    where sanitised_name is derived from the JSON name field, which in turn comes
    from the FASTA filename (without extension).

    We look for the data JSON by scanning the directory for subdirectories matching
    the sanitised name of the chain.
    """
    chain_info = unique_chains[chain_hash]
    fasta_filename = chain_info["fasta_path"]
    # The FASTA filename is like "unique_protein_0.fasta"
    # The JSON name will be the same stem
    chain_name = os.path.splitext(fasta_filename)[0]
    s_name = sanitised_name(chain_name)

    # Look for the data JSON in the expected location
    data_json_path = os.path.join(data_json_dir, s_name, f"{s_name}_data.json")
    if os.path.exists(data_json_path):
        return data_json_path

    # Fallback: scan all subdirectories for matching *_data.json
    for entry in os.listdir(data_json_dir):
        candidate = os.path.join(data_json_dir, entry, f"{entry}_data.json")
        if os.path.exists(candidate):
            # Check if this data JSON contains a chain with the same sequence
            try:
                with open(candidate, 'r') as f:
                    data = json.load(f)
                for seq_entry in data.get("sequences", []):
                    for chain_type in ("protein", "rna"):
                        if chain_type in seq_entry:
                            if seq_entry[chain_type].get("sequence") == chain_info["sequence"]:
                                return candidate
            except (json.JSONDecodeError, KeyError):
                continue

    return None


def load_msa_and_templates(data_json_path):
    """
    Load unpairedMsa, pairedMsa, and templates from a *_data.json file.

    Returns a dict with keys: unpairedMsa, pairedMsa, templates
    """
    with open(data_json_path, 'r') as f:
        data = json.load(f)

    # The data JSON has the same structure as an AF3 input JSON,
    # but with MSA and template fields filled in.
    # It contains a single chain (monomer), so we take the first sequence entry.
    sequences = data.get("sequences", [])
    if not sequences:
        return {"unpairedMsa": None, "pairedMsa": None, "templates": None}

    first_seq = sequences[0]

    if "protein" in first_seq:
        chain_data = first_seq["protein"]
        return {
            "unpairedMsa": chain_data.get("unpairedMsa"),
            "pairedMsa": chain_data.get("pairedMsa"),
            "templates": chain_data.get("templates"),
        }
    elif "rna" in first_seq:
        chain_data = first_seq["rna"]
        return {
            "unpairedMsa": chain_data.get("unpairedMsa"),
            "pairedMsa": None,  # RNA doesn't have pairedMsa in AF3
            "templates": None,  # RNA doesn't have templates in AF3
        }

    return {"unpairedMsa": None, "pairedMsa": None, "templates": None}


def build_protein_entry(chain_id, sequence, msa_data):
    """
    Build a protein sequence entry for the assembled JSON.

    Uses the current nf-core format: {"protein": {"id": ..., "sequence": ...}}
    with MSA/template fields added when available.
    """
    entry = {
        "id": chain_id,
        "sequence": sequence,
    }

    # Add MSA fields if available from pre-computed data
    if msa_data.get("unpairedMsa") is not None:
        entry["unpairedMsa"] = msa_data["unpairedMsa"]
    if msa_data.get("pairedMsa") is not None:
        entry["pairedMsa"] = msa_data["pairedMsa"]
    if msa_data.get("templates") is not None:
        entry["templates"] = msa_data["templates"]

    return {"protein": entry}


def build_rna_entry(chain_id, sequence, msa_data):
    """
    Build an RNA sequence entry for the assembled JSON.

    Uses the current nf-core format: {"rna": {"id": ..., "sequence": ...}}
    with unpairedMsa added when available.
    """
    entry = {
        "id": chain_id,
        "sequence": sequence,
    }

    if msa_data.get("unpairedMsa") is not None:
        entry["unpairedMsa"] = msa_data["unpairedMsa"]

    return {"rna": entry}


def build_dna_entry(chain_id, sequence):
    """
    Build a DNA sequence entry for the assembled JSON.
    """
    return {"dna": {"id": chain_id, "sequence": sequence}}


def build_ligand_ccd_entry(chain_id, ccd_codes):
    """
    Build a ligand entry using CCD codes.

    Uses the AF3-spec-compliant "ligand" key with "ccdCodes" sub-field.
    """
    return {"ligand": {"id": chain_id, "ccdCodes": ccd_codes}}


def build_ligand_smiles_entry(chain_id, smiles):
    """
    Build a ligand entry using a SMILES string.

    Uses the AF3-spec-compliant "ligand" key with "smiles" sub-field.
    """
    return {"ligand": {"id": chain_id, "smiles": smiles}}


def main(args=None):
    args = parse_args(args)

    os.makedirs(args.output_dir, exist_ok=True)

    # Load chain map
    with open(args.chain_map, 'r') as f:
        chain_map = json.load(f)

    unique_chains = chain_map["unique_chains"]
    samples = chain_map["samples"]

    # Pre-load data JSONs for all unique chains
    chain_msa_cache = {}
    for chain_hash, chain_info in unique_chains.items():
        data_json_path = find_data_json(args.data_json_dir, chain_hash, unique_chains)
        if data_json_path is not None:
            chain_msa_cache[chain_hash] = load_msa_and_templates(data_json_path)
            # Update the chain_map with the data_json_path
            chain_info["data_json_path"] = data_json_path
        else:
            chain_msa_cache[chain_hash] = {
                "unpairedMsa": None,
                "pairedMsa": None,
                "templates": None,
            }

    # Assemble JSON for each sample
    for sample_id, sample_info in samples.items():
        sequences = []

        for chain in sample_info["chains"]:
            entity_type = chain["entity_type"]
            chain_id = chain["chain_id"]

            if entity_type == "protein":
                chain_hash = chain["hash"]
                chain_info = unique_chains[chain_hash]
                msa_data = chain_msa_cache[chain_hash]
                sequences.append(
                    build_protein_entry(chain_id, chain_info["sequence"], msa_data)
                )

            elif entity_type == "rna":
                chain_hash = chain["hash"]
                chain_info = unique_chains[chain_hash]
                msa_data = chain_msa_cache[chain_hash]
                sequences.append(
                    build_rna_entry(chain_id, chain_info["sequence"], msa_data)
                )

            elif entity_type == "dna":
                sequences.append(
                    build_dna_entry(chain_id, chain["sequence"])
                )

            elif entity_type == "ccd":
                sequences.append(
                    build_ligand_ccd_entry(chain_id, chain["ccdCodes"])
                )

            elif entity_type == "smiles":
                sequences.append(
                    build_ligand_smiles_entry(chain_id, chain["smiles"])
                )

            else:
                # Unknown type — skip with warning
                print(f"WARNING: Skipping unknown entity type '{entity_type}' "
                      f"in sample '{sample_id}', chain '{chain_id}'")

        # Build the complete AF3 JSON
        assembled_json = {
            "name": sample_id,
            "sequences": sequences,
            "modelSeeds": args.model_seeds,
            "dialect": "alphafold3",
            "version": 1,
        }

        # Write to output file
        out_path = os.path.join(args.output_dir, f"{sanitised_name(sample_id)}.json")
        with open(out_path, 'w') as f:
            json.dump(assembled_json, f, indent=4)

        print(f"Assembled JSON for sample '{sample_id}' -> {out_path} "
              f"({len(sequences)} chains)")

    # Update chain_map.json with data_json_path values
    with open(args.chain_map, 'w') as f:
        json.dump(chain_map, f, indent=2)


if __name__ == "__main__":
    main()
