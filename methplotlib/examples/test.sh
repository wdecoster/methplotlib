set -ev

methplotlib -h
methplotlib -m examples/ACTB_calls.tsv.gz examples/meth_freq.tsv.gz \
            -n calls frequencies \
            -w chr7:5,525,542-5,543,028 \
            -g examples/g38_locus.gtf.gz \
            --simplify \
            -b examples/DNAse_cluster.bed.gz
