set -ev

methplotlib -h
methplotlib -m examples/NA19240-methylation_ACTB_calls.tsv.gz examples/NA19240-methylation_ACTB_frequency.tsv.gz \
            -n calls frequencies \
            -w chr7:5,525,542-5,543,028 \
            -g examples/GRCh38-ACTB-locus.gtf.gz \
            --simplify \
            -b examples/DNase_cluster_ACTB.bed.gz
