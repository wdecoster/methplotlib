import pandas as pd
import pyranges as pr
import numpy as np
import sys
import logging
from methplotlib.utils import file_sniffer, flatten
from itertools import repeat


class Modification(object):
    def __init__(self, table, data_type, name, called_sites, start_end_table=None):
        self.table = table
        self.data_type = data_type
        self.name = name
        self.called_sites = called_sites
        self.start_end_table = start_end_table


def get_data(args, window):
    """
    Import methylation data from all files in the list methylation_files

    Data can in various formats
    - nanopolish raw, phased, frequency
    - cram
    - bam
    - nanocompore
    - bedgraph
    - bedmethyl in the extended form (as obtained from modbam2bed --extended)

    data is extracted within the window args.window
    Frequencies are smoothened using a sliding window
    """
    return flatten([read_mods(f, n, window, args) for f, n in zip(args.methylation, args.names)])


def read_mods(filename, name, window, args):
    """
    converts a file from nanopolish to a pandas dataframe
    input can be from calculate_methylation_frequency
    which will return a dataframe with 'chromosome', 'pos', 'methylated_frequency'
    smoothening the result by a rolling average

    input can also be raw data per read, optionally phased
    which will return a dataframe with 'read', 'chromosome', 'pos', 'log_lik_ratio', 'strand'
    """
    file_type = file_sniffer(filename)
    logging.info(f"File {filename} is of type {file_type}")
    try:
        if file_type.startswith("nanopolish"):
            return parse_nanopolish(filename, file_type, name, window, smoothen=args.smoothen)
        elif file_type == "nanocompore":
            return [parse_nanocompore(filename, name, window)]
        elif file_type in ["cram", "bam"]:
            return parse_cram(filename, file_type, name, window, args.mods)
        elif file_type == "bedgraph":
            return [parse_bedgraph(filename, name, window)]
        elif file_type == "bedmethyl_extended":
            return parse_bedmethyl(filename, name, window, smoothen=args.smooth, flavor="modbam2bed", mods_of_interest=args.mods)
        elif file_type == "bedmethyl":
            return parse_bedmethyl(filename, name, window, smoothen=args.smooth, flavor="modkit", mods_of_interest=args.mods)
    except Exception as e:
        logging.error(f"Error processing {filename}.")
        logging.error(e, exc_info=True)
        sys.stderr.write(f"\n\n\nError processing {filename}!\n")
        sys.stderr.write("\n\n\nDetailed error:\n")
        raise


def parse_nanopolish(filename, file_type, name, window, smoothen=5):
    if window:
        from pathlib import Path

        if Path(filename + ".tbi").is_file():
            import subprocess
            import gzip

            logging.info(f"Reading {filename} using a tabix stream.")
            region = f"{window.chromosome}:{window.begin}-{window.end}"
            try:
                tabix_stream = subprocess.Popen(
                    ["tabix", filename, region],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
            except FileNotFoundError as e:
                logging.error("Error when opening a tabix stream.")
                logging.error(e, exc_info=True)
                sys.stderr.write("\n\nERROR when opening a tabix stream.\n")
                sys.stderr.write("Is tabix installed and on the PATH?.")
                raise
            header = gzip.open(filename, "rt").readline().rstrip().split("\t")
            table = pd.read_csv(tabix_stream.stdout, sep="\t", header=None, names=header)
        else:
            logging.info(f"Reading {filename} slowly by splitting the file in chunks.")
            sys.stderr.write(f"\nReading {filename} would be faster with bgzip and tabix.\n")
            if file_type in ["nanopolish_call", "nanopolish_phased"]:
                sys.stderr.write("Please index with 'tabix -S1 -s1 -b3 -e4'.\n")
            else:
                sys.stderr.write("Please index with 'tabix -S1 -s1 -b2 -e3'.\n")
            iter_csv = pd.read_csv(filename, sep="\t", iterator=True, chunksize=1e6)
            table = pd.concat(
                [chunk[chunk["chromosome"] == window.chromosome] for chunk in iter_csv]
            )
    else:
        table = pd.read_csv(filename, sep="\t")
    gr = pr.PyRanges(
        table.rename(
            columns={
                "start": "Start",
                "chromosome": "Chromosome",
                "end": "End",
                "Strand": "strand",
            }
        )
    )
    logging.info("Read the file in a dataframe.")

    if window:
        gr = gr[str(window.chromosome), window.begin : window.end]
    try:
        gr.pos = np.floor(gr.drop().df[["Start", "End"]].mean(axis=1))
    except KeyError:
        sys.stderr.write(f"\n\n\nProblem parsing nanopolish file {filename}!\n")
        sys.stderr.write(
            "Could it be that there are no calls in your selected window "
            f"{window.chromosome}:{window.begin}-{window.end}?\n"
        )
        sys.stderr.write("\n\n\nDetailed error:\n")
        raise

    table = gr.df

    if file_type in ["nanopolish_call", "nanopolish_phased"]:
        table = table.drop(
            columns=[
                "Start",
                "End",
                "log_lik_methylated",
                "log_lik_unmethylated",
                "num_calling_strands",
                "num_motifs",
                "sequence",
            ]
        )
        if "motif" in table:
            return [
                Modification(
                    table=sub_df,
                    data_type=file_type,
                    name=f"{name}_{mod}",
                    called_sites=len(sub_df),
                )
                for mod, sub_df in table.groupby("motif")
            ]
        else:
            return [
                Modification(
                    table=table.sort_values(["read_name", "pos"]),
                    data_type=file_type,
                    name=name,
                    called_sites=len(table),
                )
            ]
    if file_type == "nanopolish_freq":
        called_sites = table.called_sites
        chromosome = table.Chromosome.values[0]

        table = table.drop(
            columns=[
                "Chromosome",
                "Start",
                "End",
                "num_motifs_in_group",
                "called_sites",
                "called_sites_methylated",
                "group_sequence",
            ]
        )
        return [
            Modification(
                table=table.sort_values("pos")
                .groupby("pos")
                .mean()
                .rolling(window=smoothen, center=True)
                .mean()
                .assign(Chromosome=chromosome),
                data_type=file_type,
                name=name,
                called_sites=called_sites.sum(),
            )
        ]


def parse_nanocompore(filename, name, window):
    def nanocompore_columns_of_interest(column):
        if column in ["pos", "ref_id"]:
            return True
        elif column.endswith("pvalue_context_2") or column.endswith("pvalue"):
            return True
        else:
            return False

    table = pd.read_csv(filename, sep="\t", usecols=nanocompore_columns_of_interest)
    if window:
        table = table[table["ref_id"] == window.chromosome]
    return Modification(
        table=table.sort_values("pos")
        .append({"pos": window.end}, ignore_index=True)
        .drop(columns="ref_id")
        .fillna(1.0),
        data_type="nanocompore",
        name=name,
        called_sites=len(table),
    )


def parse_bedgraph(filename, name, window):
    if window:
        from pathlib import Path

        if Path(filename + ".tbi").is_file():
            import subprocess

            logging.info(f"Reading {filename} using a tabix stream.")
            region = f"{window.chromosome}:{window.begin}-{window.end}"
            try:
                tabix_stream = subprocess.Popen(
                    ["tabix", filename, region],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
            except FileNotFoundError as e:
                logging.error("Error when opening a tabix stream.")
                logging.error(e, exc_info=True)
                sys.stderr.write("\n\nERROR when opening a tabix stream.\n")
                sys.stderr.write("Is tabix installed and on the PATH?.")
                raise
            table = pd.read_csv(
                tabix_stream.stdout,
                sep="\t",
                header=None,
                names=["Chromosome", "Start", "End", "Value"],
            )
        else:
            logging.info(f"Reading {filename} slowly by splitting the file in chunks.")
            sys.stderr.write(
                f"\nReading {filename} would be faster with bgzip and 'tabix -p bed'.\n"
            )
            iter_csv = pd.read_csv(
                filename,
                sep="\t",
                iterator=True,
                chunksize=1e6,
                header=None,
                names=["Chromosome", "Start", "End", "Value"],
            )
            table = pd.concat(
                [chunk[chunk["Chromosome"] == window.chromosome] for chunk in iter_csv]
            )
    else:
        table = pd.read_csv(
            filename,
            sep="\t",
            header=None,
            names=["Chromosome", "Start", "End", "Value"],
        )
    gr = pr.PyRanges(table)
    logging.info("Read the file in a dataframe.")
    if window:
        gr = gr[window.chromosome, window.begin : window.end]
        if len(gr.df) == 0:
            sys.exit(f"No records for {filename} in {window.string}!\n")
    return Modification(
        table=gr.df.sort_values("Start"),
        data_type="bedgraph",
        name=name,
        called_sites=len(table),
    )


def parse_bedmethyl(filename, name, window, smoothen=5, flavor="modkit", mods_of_interest=None):
    if not flavor in ["modkit", "modbam2bed"]:
        sys.exit(f"ERROR: flavor {flavor} not supported for bedmethyl.")
    if flavor == "modkit":
        colnames = {
            0: "Chromosome",
            1: "Start",
            2: "End",
            3: "Modification",
            10: "Frequency",
        }
    else:
        colnames = {
            0: "Chromosome",
            1: "Start",
            2: "End",
            3: "Modification",
            11: "canonical",
            12: "modified",
        }
    usecols = colnames.keys()
    if window:
        from pathlib import Path

        if Path(filename + ".tbi").is_file():
            import subprocess

            logging.info(f"Reading {filename} using a tabix stream.")
            region = f"{window.chromosome}:{window.begin}-{window.end}"
            try:
                tabix_stream = subprocess.Popen(
                    ["tabix", filename, region],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
            except FileNotFoundError as e:
                logging.error("Error when opening a tabix stream.")
                logging.error(e, exc_info=True)
                sys.stderr.write("\n\nERROR when opening a tabix stream.\n")
                sys.stderr.write("Is tabix installed and on the PATH?.")
                raise
            table = pd.read_csv(
                tabix_stream.stdout,
                sep="\t",
                header=None,
                usecols=usecols,
            ).rename(columns=colnames)
        else:
            logging.info(f"Reading {filename} slowly by splitting the file in chunks.")
            sys.stderr.write(
                f"\nReading {filename} would be faster with bgzip and 'tabix -p bed'.\n"
            )
            iter_csv = pd.read_csv(
                filename,
                sep="\t",
                iterator=True,
                chunksize=1e6,
                header=None,
                usecols=usecols,
            )
            table = pd.concat([chunk[chunk[0] == window.chromosome] for chunk in iter_csv]).rename(
                columns=colnames
            )
    else:
        table = pd.read_csv(
            filename,
            sep="\t",
            header=None,
            usecols=usecols,
        ).rename(columns=colnames)
    gr = pr.PyRanges(table)
    logging.info("Read the file in a dataframe.")
    if window:
        gr = gr[window.chromosome, window.begin : window.end]
        if len(gr.df) == 0:
            sys.exit(f"No records for {filename} in {window.string}!\n")
    table = gr.df.sort_values("Start")
    if flavor == "modkit":
        table["modified_frequency"] = table["Frequency"] / 100
        table.drop(columns=["Frequency"], inplace=True)
    else:
        table["modified_frequency"] = table["modified"] / (table["canonical"] + table["modified"])
    table.drop(
        columns=["Chromosome", "End", "canonical", "modified"], inplace=True, errors="ignore"
    )
    if mods_of_interest:
        # mods of interest is in the form of C+m or C+h (as in the cram file)
        # however, the bedmethyl file has the modifications in the form of m or h
        # so we need to split the mods in mods_of_interest on the '+' and select the second part
        table = table[table["Modification"].isin(m.split('+')[1] for m in mods_of_interest.split(","))]

    return [
        Modification(
            table=sub_df.drop(columns=["Modification"]).rolling(window=smoothen, center=True).mean(),
            data_type="bedmethyl_extended",
            name=f"{name}_{mod}",
            called_sites=len(sub_df),
        )
        for mod, sub_df in table.groupby("Modification") if len(sub_df) > 0
    ]


def parse_cram(filename, filetype, name, window, mods_of_interest=None):
    """
    Extracts modified positions from a CRAM file

    :param filename: str, path to the CRAM file
    :param filetype: str, either 'cram' or 'bam'
    :param name: str, name for the trace/sample
    :param window: Region object to extract data for from the file
    :param mods_of_interest: list of str, optional, list of modifications to extract
    """
    import pysam

    mode = "rc" if filetype == "cram" else "rb"
    cram = pysam.AlignmentFile(filename, mode)
    data = []
    start_stops = []
    for read in cram.fetch(reference=str(window.chromosome), start=window.begin, end=window.end):
        if not read.is_supplementary and not read.is_secondary:
            start_stops.append((read.query_name, read.reference_start, read.reference_end))
            mod_positions = get_modified_reference_positions(read)
            if mod_positions:
                data.extend(mod_positions)
    df = (
        pd.DataFrame(data, columns=["read_name", "strand", "pos", "quality", "mod"])
        .astype(dtype={"mod": "category", "quality": "float"})
        .sort_values(["read_name", "pos"])
    )


    if mods_of_interest:
        mods_found = df["mod"].unique().categories.values
        df = df[df["mod"].isin(mods_of_interest.split(","))]
        if len(df) == 0:
            sys.exit(f"No more records after selecting --mods!\nDetected modifications: {mods_found}\n")
    return [
        Modification(
            table=sub_df,
            data_type="ont-cram",
            name=f"{name}_{mod}",
            called_sites=len(sub_df),
            start_end_table=pd.DataFrame(
                start_stops, columns=["read_name", "posmin", "posmax"]
            ).set_index("read_name"),
        )
        for mod, sub_df in df.groupby("mod") if len(sub_df) > 0
    ]


def get_modified_reference_positions(read):
    mod_positions = []
    if read.has_tag("Mm") or read.has_tag("MM"):
        modtag = "Mm" if read.has_tag("Mm") else "MM"
        qualtag = "Ml" if read.has_tag("Ml") else "ML"
        offset = 0
        for context in read.get_tag(modtag).split(";"):
            if context:
                basemod = context.split(",", 1)[0].replace('?', '')
                if "-" in basemod:
                    sys.exit(
                        "ERROR: modifications on negative strand currently unsupported.\n"
                        "Please contact me if this would be of interest for you."
                    )
                base, mod = basemod.split("+")
                # code below does not work with for ambiguous N and will search for a literal N
                if base == "N":
                    sys.exit(
                        "ERROR: modifications of N nucleotides currently unsupported.\n"
                        "Please contact me if this would be of interest for you."
                    )
                # The positions are encoded by specifying the number of non-modified occurences
                # of that specific bases to skip in the read sequence
                deltas = [int(i) for i in context.split(",")[1:]]
                if not deltas:
                    continue
                # Make an array with all positions for which a modification is reported
                # In coordinates relative to occurences of the nucleotide in the read
                locations = np.cumsum(deltas) + np.arange(len(deltas))
                # Nake an array with all read positions for which that nucleotide exists in the read
                # Oddly, get_forward_sequence() returns the sequence as it came from the sequencer
                # This is the same direction as the Mm/Ml tags are specified in
                base_index = np.array(
                    [i for i, letter in enumerate(read.get_forward_sequence()) if letter == base]
                )
                # Based on both arrays, find those read coordinates that have a modification
                modified_bases = base_index[locations]
                # Convert the array of read coordinates to reference coordinates
                # full_length=True will make sure the list is of the same length as the read
                # by adding 'None' for softclipped nucleotides
                refpos = np.array(read.get_reference_positions(full_length=True))
                # read.get_reference_positions() returns incrementing coordinates
                # regardless of strand, so reversing the order for reverse_stranded reads
                if read.is_reverse:
                    refpos = refpos[::-1]
                # The likelihoods are in an array of length of all Mm/MM deltas,
                # and are not separated by context/modified nucleotide type
                if read.has_tag(qualtag):
                    likelihoods = read.get_tag(qualtag).tolist()[offset : offset + len(deltas)]
                else:
                    likelihoods = [255] * len(deltas)
                offset += len(deltas)
                mod_positions.extend(
                    zip(
                        repeat(read.query_name),  # read_name
                        repeat("-" if read.is_reverse else "+"),  # strand
                        refpos[modified_bases],  # pos
                        likelihoods,  # quality
                        repeat(basemod),  # mod
                    )
                )
        return mod_positions
    else:
        return None


def errs_tab(n):
    """Generate list of error rates for qualities less than equal than n."""
    return [10 ** (q / -10) for q in range(n + 1)]


def phred_to_probability(quals, tab=errs_tab(128)):
    return [tab[ord(q) - 33] for q in quals]
