from argparse import ArgumentParser
import bgzip
import pandas as pd


def main():
    args = get_args()

    if args.naive:
        phased_calls = bgzip.open(args.phased_methylation, 'rt')

        phase1 = bgzip.open(args.prefix + "_phase1.tsv.gz", 'wt')
        phase2 = bgzip.open(args.prefix + "_phase2.tsv.gz", 'wt')
        uphase = bgzip.open(args.prefix + "_unphased.tsv.gz", 'wt')

        header = next(phased_calls)
        for output in [phase1, phase2, uphase]:
            output.write(header)

        for line in phased_calls:
            if line.endswith('1.0\n') or line.endswith('1\n'):
                phase1.write(line)
            elif line.endswith('2.0\n') or line.endswith('2\n'):
                phase2.write(line)
            else:
                uphase.write(line)
    else:
        df = pd.read_csv(args.phased_methylation, sep="\t")
        df[df["HP"].isna()].to_csv(args.prefix + "_calls_unphased.tsv.gz", sep="\t", index=False)
        df = df[df["HP"].notna()]
        ps_counts = df.loc[:, ["PS", "HP"]].drop_duplicates()["PS"].value_counts()
        hom_blocks = ps_counts[ps_counts == 1].index
        df[df.loc[:, "PS"].isin(hom_blocks)] \
            .to_csv(args.prefix + "_calls_homozygous.tsv.gz", sep="\t", index=False)
        df = df[~df.loc[:, "PS"].isin(hom_blocks)]
        df[df["HP"] == 1.0].to_csv(args.prefix + "_calls_phase1.tsv.gz", sep="\t", index=False)
        df[df["HP"] == 2.0].to_csv(args.prefix + "_calls_phase2.tsv.gz", sep="\t", index=False)


def get_args():
    parser = ArgumentParser(description="Split file with phased calls by phase.")
    parser.add_argument("phased_methylation", help="File created by annotate_calls_by_phase.py")
    parser.add_argument("-p", "--prefix", help="Prefix for output files", required=True)
    parser.add_argument("--naive",
                        help="Naively split reads, not taking homozygous regions into account.")
    return parser.parse_args()


if __name__ == '__main__':
    main()
