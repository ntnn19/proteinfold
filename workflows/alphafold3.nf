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

    if (params.alphafold3_deduplicate_msa) {
        //
        // DEDUPLICATED MSA WORKFLOW
        // Run data pipeline only for unique protein/RNA chains,
        // then assemble multimeric JSONs with pre-computed MSAs,
        // then run inference-only.
        //

        // Step 1: Build input TSV and collect all FASTAs for deduplication
        // The TSV maps sample_id -> fasta_path so that chain_map.json
        // uses the original sample IDs from the samplesheet.
        // We fork the samplesheet channel to consume it twice.
        ch_samplesheet_dup = ch_samplesheet.dup()

        ch_input_tsv = ch_samplesheet_dup
            .map { meta, fasta -> "${meta.id}\t${fasta}" }
            .collectFile(name: 'input.tsv', storeDir: '.', sort: true)

        ch_fasta_collect = ch_samplesheet_dup
            .map { meta, fasta -> fasta }
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
        // Each JSON filename is <sanitised_sample_id>.json
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
        // Convert models mmcifs to pdbs
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

        // Convert top ranked mmcif to pdb
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

        // Prepare msa input
        RUN_ALPHAFOLD3_INFERENCE
            .out
            .msa
            .map { it ->
                def meta = it[0].clone();
                meta.model = "alphafold3";
                [ meta, it[1] ]
            }
            .set { ch_msa_final }

        // Prepare report input
        RUN_ALPHAFOLD3_INFERENCE
            .out
            .multiqc
            .map { it -> it[1] }
            .toSortedList()
            .map { it ->
                [ [ "model": "alphafold3" ], it.flatten() ]
            }
            .set { ch_multiqc_report }

        // Prepare pae input
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

        FASTA_TO_ALPHAFOLD3_JSON(ch_samplesheet)
        ch_versions       = ch_versions.mix(FASTA_TO_ALPHAFOLD3_JSON.out.versions)

        //
        // SUBWORKFLOW: Run AlphaFold3
        //
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

        // Convert mmcif to pdbs
        RUN_ALPHAFOLD3
                .out
                .cif
                .groupTuple()
                .map {
                    meta, files ->
                    [ meta, files.flatten() ]
                }

        // Convert models mmcifs to pdbs
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

        // Convert top ranked mmcif to pdb
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

        // Prepare msa input
        RUN_ALPHAFOLD3
            .out
            .msa
            .map { it ->
                def meta = it[0].clone();
                meta.model = "alphafold3";
                [ meta, it[1] ]
            }
            .set { ch_msa_final }

        // Prepare report input
        RUN_ALPHAFOLD3
            .out
            .multiqc
            .map { it -> it[1] }
            .toSortedList()
            .map { it ->
                [ [ "model": "alphafold3" ], it.flatten() ]
            }
            .set { ch_multiqc_report }

        // Prepare dummy pae input
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
    top_ranked_pdb = ch_top_ranked_pdb // channel: [ id, /path/to/*.pdb ]
    pdb            = ch_pdb_final      // channel: [ meta, /path/to/*.pdb, ...,/path/to/*.pdb ]
    msa            = ch_msa_final      // channel: [ meta, /path/to/*.pdb, /path/to/*_coverage.png ]
    pae            = ch_pae_final      // channel: [ meta, path/to/*_pae.tsv ]
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
