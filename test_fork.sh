#!/bin/bash
# Create test FASTA files
cat > test_sample1.fasta << 'EOF'
>protein
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAVKNT
>ccd
ATP
EOF

cat > test_sample2.fasta << 'EOF'
>protein
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAVKNT
>protein
MAKEDNIKLQQTFSVLISDGRTVGFSGRQLVAGKQVVNIPYKLENLGLDHIVKQTLKNS
EOF

# Create samplesheet pointing to the FASTAs
cat > test_samplesheet.csv << 'EOF'
id,fasta
sample1,test_sample1.fasta
sample2,test_sample2.fasta
EOF

# Run from your fork (replace <branch> with your actual branch name)
nextflow run main.nf \
    -profile singularity \
    --input test_samplesheet.csv \
    --outdir results_test_dedup \
    --mode alphafold3 \
    --alphafold3_db /scratch1/common/databases/af3_db_structured \
    --alphafold3_deduplicate_msa \
    --use_gpu
