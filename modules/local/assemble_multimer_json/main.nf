/*
 * Assemble multimeric AlphaFold3 JSON files from pre-computed MSAs.
 * Reads chain_map.json and per-chain *_data.json files, then
 * assembles per-sample multimer JSONs with MSAs filled in.
 *
 * Outputs all assembled JSONs at once. The workflow splits them
 * into per-sample channels for downstream inference.
 */
process ASSEMBLE_MULTIMER_JSON {
    tag "assemble_multimer_json"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7dabd3f132a613fb11ee27c66e9517eb7649eee64f4e4f63747841105883b40/data' :
        'community.wave.seqera.io/library/biopython_python:06582b7b722f3db3' }"

    input:
    path chain_map
    path "data_jsons/*"

    output:
    path "assembled/*.json"  , emit: json
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p assembled

    assemble_multimer_json.py \\
        --chain_map ${chain_map} \\
        --data_json_dir data_jsons/ \\
        --output_dir assembled/ \\
        --model_seeds 11

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p assembled
    touch assembled/stub.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>/dev/null | sed 's/Python //g' || echo "unknown")
    END_VERSIONS
    """
}
