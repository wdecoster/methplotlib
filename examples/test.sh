set -ev

methplotlib -h
methplotlib -m test/NA19240-methylation_ACTB_calls.tsv.gz test/NA19240-methylation_ACTB_frequency.tsv.gz \
            -n calls frequencies \
            -w chr7:5,525,542-5,543,028 \
            -g test/GRCh38-ACTB-locus.gtf.gz \
            --simplify \
            -b test/DNase_cluster_ACTB.bed.gz
