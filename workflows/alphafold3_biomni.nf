/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { FASTA_TO_ALPHAFOLD3_JSON          } from '../modules/local/fasta_to_alphafold3_json'
include { RUN_ALPHAFOLD3                    } from '../modules/local/run_alphafold3'
include { MMCIF2PDB as MMCIF2PDB_TOP_RANKED } from '../modules/local/mmcif2pdb/main.nf'
include { MMCIF2PDB as MMCIF2PDB_MODELS     } from '../modules/local/mmcif2pdb/main.nf'

//
// MODULE: Deduplicated MSA workflow modules (opt-in via --alphafold3_deduplicate_msa)
//
include { DEDUPLICATE_CHAINS                } from '../modules/local/deduplicate_chains'
include { ASSEMBLE_MULTIMER_JSON            } from '../modules/local/assemble_multimer_json'
include { RUN_ALPHAFOLD3_DATA_PIPELINE      } from '../modules/local/run_alphafold3_data_pipeline'
include { RUN_ALPHAFOLD3_INFERENCE          } from '../modules/local/run_alphafold3_inference'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ALPHAFOLD3 {

    take:
    input 
    ch_samplesheet       // channel: samplesheet read in from --input
    ch_versions          // channel: [ path(versions.yml) ]
    ch_alphafold3_params // channel: path(alphafold2_params)
    ch_small_bfd         // channel: path(small_bfd)
    ch_mgnify            // channel: path(mgnify)
    ch_mmcif_files       // channel: path(mmcif_files)
    ch_uniref90          // channel: path(uniref90)
    ch_pdb_seqres        // channel: path(pdb_seqres)
    ch_uniprot           // channel: path(uniprot)

    main:
    ch_pdb_final      = channel.empty()
    ch_top_ranked_pdb = channel.empty()
    ch_msa_final      = channel.empty()
    ch_multiqc_report = channel.empty()

    // In DSL2, a queue channel can only be consumed once.
    // collect() produces a VALUE channel (multi-readable) that both
    // conditional branches can derive what they need from.
    ch_ss_collected = ch_samplesheet
    ch_ss_collected
        .view { println "BIOMNI: $it" }

    if (params.alphafold3_deduplicate_msa) {
        //
        // DEDUPLICATED MSA WORKFLOW
        //

        // Step 1: Derive input TSV and FASTA file list from the collected value channel
        ch_input_tsv = ch_ss_collected
            .map { meta, fasta ->
                "${meta.id}\t${fasta}"
            }
            .collect()
            .map { lines ->
                file("input.tsv") << lines.join('\n')
            }

        ch_fasta_collect = ch_input_tsv
            .flatMap { entries -> entries.collect { meta, fasta -> fasta } }
            .collect()

        DEDUPLICATE_CHAINS(ch_fasta_collect, ch_input_tsv)
        ch_versions = ch_versions.mix(DEDUPLICATE_CHAINS.out.versions)

        // Step 2: Convert each unique chain FASTA to monomer JSON
        ch_unique_chain_meta = DEDUPLICATE_CHAINS.out.unique_fasta
            .map { fasta ->
                def name = fasta.baseName
                def meta = [:]
                meta.id = name
                [ meta, fasta ]
            }

        FASTA_TO_ALPHAFOLD3_JSON(ch_unique_chain_meta)
        ch_versions = ch_versions.mix(FASTA_TO_ALPHAFOLD3_JSON.out.versions)

        // Step 3: Run data pipeline (MSA + templates) for each unique chain
        RUN_ALPHAFOLD3_DATA_PIPELINE(
            FASTA_TO_ALPHAFOLD3_JSON.out.json,
            ch_alphafold3_params,
            ch_small_bfd,
            ch_mgnify,
            ch_mmcif_files,
            ch_uniref90,
            ch_pdb_seqres,
            ch_uniprot
        )
        ch_versions = ch_versions.mix(RUN_ALPHAFOLD3_DATA_PIPELINE.out.versions)

        // Step 4: Collect all data JSONs and assemble multimeric JSONs
        ch_data_jsons_collect = RUN_ALPHAFOLD3_DATA_PIPELINE.out.data_json
            .map { meta, data_json -> data_json }
            .collect()

        ASSEMBLE_MULTIMER_JSON(
            DEDUPLICATE_CHAINS.out.chain_map,
            ch_data_jsons_collect
        )
        ch_versions = ch_versions.mix(ASSEMBLE_MULTIMER_JSON.out.versions)

        // Split assembled JSONs into per-sample channels
        ch_assembled_json = ASSEMBLE_MULTIMER_JSON.out.json
            .flatten()
            .map { json_file ->
                def name = json_file.baseName
                def meta = [:]
                meta.id = name
                [ meta, json_file ]
            }

        // Step 5: Run inference-only on assembled multimer JSONs
        RUN_ALPHAFOLD3_INFERENCE(
            ch_assembled_json,
            ch_alphafold3_params,
            ch_small_bfd,
            ch_mgnify,
            ch_mmcif_files,
            ch_uniref90,
            ch_pdb_seqres,
            ch_uniprot
        )
        ch_versions = ch_versions.mix(RUN_ALPHAFOLD3_INFERENCE.out.versions)

        // Post-processing: same as standard flow
        MMCIF2PDB_MODELS (
            RUN_ALPHAFOLD3_INFERENCE
                .out
                .cif
                .groupTuple()
                .map {
                    meta, files ->
                    [ meta, files.flatten() ]
                }
        )
        ch_versions = ch_versions.mix(MMCIF2PDB_MODELS.out.versions)

        MMCIF2PDB_MODELS
            .out
            .pdb
            .map { it ->
                def meta   = it[0].clone();
                meta.model = "alphafold3";
                def files = (it[1] instanceof List) ? it[1] : [ it[1] ]
                [ meta, files ]
            }
            .set { ch_pdb_final }

        MMCIF2PDB_TOP_RANKED (
            RUN_ALPHAFOLD3_INFERENCE
                .out
                .top_ranked_cif
        )
        ch_versions = ch_versions.mix(MMCIF2PDB_TOP_RANKED.out.versions)

        MMCIF2PDB_TOP_RANKED
            .out
            .pdb
            .map { it ->
                def meta = it[0].clone();
                meta.model = "alphafold3";
                [ meta, it[1] ]
            }
            .set { ch_top_ranked_pdb }

        RUN_ALPHAFOLD3_INFERENCE
            .out
            .msa
            .map { it ->
                def meta = it[0].clone();
                meta.model = "alphafold3";
                [ meta, it[1] ]
            }
            .set { ch_msa_final }

        RUN_ALPHAFOLD3_INFERENCE
            .out
            .multiqc
            .map { it -> it[1] }
            .toSortedList()
            .map { it ->
                [ [ "model": "alphafold3" ], it.flatten() ]
            }
            .set { ch_multiqc_report }

        RUN_ALPHAFOLD3_INFERENCE
            .out
            .pae
            .map { it ->
                def meta = it[0].clone();
                meta.model = "alphafold3";
                [ meta, it[1] ]
            }
            .set { ch_pae_final }

    } else {
        //
        // STANDARD WORKFLOW (unchanged)
        //

        // Split the collected value channel back into individual [meta, fasta] tuples
        ch_ss_split = ch_ss_collected.flatMap { it }

        FASTA_TO_ALPHAFOLD3_JSON(ch_ss_split)
        ch_versions       = ch_versions.mix(FASTA_TO_ALPHAFOLD3_JSON.out.versions)

        RUN_ALPHAFOLD3 (
            FASTA_TO_ALPHAFOLD3_JSON.out.json,
            ch_alphafold3_params,
            ch_small_bfd,
            ch_mgnify,
            ch_mmcif_files,
            ch_uniref90,
            ch_pdb_seqres,
            ch_uniprot
        )
        ch_versions = ch_versions.mix(RUN_ALPHAFOLD3.out.versions)

        RUN_ALPHAFOLD3
                .out
                .cif
                .groupTuple()
                .map {
                    meta, files ->
                    [ meta, files.flatten() ]
                }

        MMCIF2PDB_MODELS (
            RUN_ALPHAFOLD3
                .out
                .cif
                .groupTuple()
                .map {
                    meta, files ->
                    [ meta, files.flatten() ]
                }
        )
        ch_versions = ch_versions.mix(MMCIF2PDB_MODELS.out.versions)

        MMCIF2PDB_MODELS
            .out
            .pdb
            .map { it ->
                def meta   = it[0].clone();
                meta.model = "alphafold3";
                def files = (it[1] instanceof List) ? it[1] : [ it[1] ]
                [ meta, files ]
            }
            .set { ch_pdb_final }

        MMCIF2PDB_TOP_RANKED (
            RUN_ALPHAFOLD3
                .out
                .top_ranked_cif
        )
        ch_versions = ch_versions.mix(MMCIF2PDB_TOP_RANKED.out.versions)

        MMCIF2PDB_TOP_RANKED
            .out
            .pdb
            .map { it ->
                def meta = it[0].clone();
                meta.model = "alphafold3";
                [ meta, it[1] ]
            }
            .set { ch_top_ranked_pdb }

        RUN_ALPHAFOLD3
            .out
            .msa
            .map { it ->
                def meta = it[0].clone();
                meta.model = "alphafold3";
                [ meta, it[1] ]
            }
            .set { ch_msa_final }

        RUN_ALPHAFOLD3
            .out
            .multiqc
            .map { it -> it[1] }
            .toSortedList()
            .map { it ->
                [ [ "model": "alphafold3" ], it.flatten() ]
            }
            .set { ch_multiqc_report }

        RUN_ALPHAFOLD3
            .out
            .pae
            .map { it ->
                def meta = it[0].clone();
                meta.model = "alphafold3";
                [ meta, it[1] ]
            }
            .set { ch_pae_final }

    } // end if/else deduplicate_msa

    emit:
    top_ranked_pdb = ch_top_ranked_pdb
    pdb            = ch_pdb_final
    msa            = ch_msa_final
    pae            = ch_pae_final
    multiqc_report = ch_multiqc_report
    versions       = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
