set -ev

methplotlib -h

# bam with MM/ML tags, bedMethyl of modkit (read with tabix) and bedgraph,
# synthetic data made for methplotlib-rs (scripts/make_test_data.py there)
methplotlib -m examples/synthetic_reads.bam examples/synthetic_modkit.bed.gz examples/synthetic_modkit.bedgraph \
            -n reads frequencies values \
            -w chrS:5,000-25,000 \
            --split

# the same, restricted to 5mC
methplotlib -m examples/synthetic_reads.bam examples/synthetic_modkit.bed.gz \
            -n reads frequencies \
            -w chrS:5,000-25,000 \
            --mods C+m

# nanopolish output, a legacy format, with gtf and bed annotation
methplotlib -m examples/NA19240-methylation_ACTB_calls.tsv.gz examples/NA19240-methylation_ACTB_frequency.tsv.gz \
            -n calls frequencies \
            -w chr7:5,525,542-5,543,028 \
            -g examples/GRCh38-ACTB-locus.gtf.gz \
            --simplify \
            -b examples/DNase_cluster_ACTB.bed.gz

# the same annotation without --simplify, which takes a different code path
# through parse_annotation and used to fail on gtf records without a transcript_id
methplotlib -m examples/NA19240-methylation_ACTB_frequency.tsv.gz \
            -n frequencies \
            -w chr7:5,525,542-5,543,028 \
            -g examples/GRCh38-ACTB-locus.gtf.gz \
            -b examples/DNase_cluster_ACTB.bed.gz
