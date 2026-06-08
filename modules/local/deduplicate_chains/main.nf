/*
 * Deduplicate protein/RNA chains across all input FASTA files.
 * Takes a TSV mapping sample_id -> fasta_path, deduplicates
 * protein/RNA chains by sequence, and outputs per-chain FASTAs
 * plus a chain_map.json for downstream assembly.
 */
process DEDUPLICATE_CHAINS {
    tag "deduplicate_chains"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7dabd3f132a613fb11ee27c66e9517eb7649eee64f4e4f63747841105883b40/data' :
        'community.wave.seqera.io/library/biopython_python:06582b7b722f3db3' }"

    input:
    path "fastas/*"
    path input_tsv

    output:
    path "chain_map.json"              , emit: chain_map
    path "unique_chains/*.fasta"       , emit: unique_fasta
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    deduplicate_chains.py \\
        --input ${input_tsv} \\
        --output_dir unique_chains/ \\
        --chain_map_out chain_map.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch chain_map.json
    mkdir -p unique_chains
    touch unique_chains/stub.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>/dev/null | sed 's/Python //g' || echo "unknown")
    END_VERSIONS
    """
}
