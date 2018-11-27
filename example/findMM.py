#!/usr/bin/env python

# ./findMM.py -p 4 -i ./index/example_transcriptome -k 30 -ref example_transcriptome.fasta -out out_test

# could add 'save_kmers'. Action.... In kmer_function(args, save_kmers):


import argparse
import sys
import os
import subprocess
import pandas as pd
from string import Template
from Bio import SeqIO

# %%
# ==============================================================================
# // k-mer scripts - in progress
# ==============================================================================


def gen_reads(seq, k):
    """
    gen_reads generates list of k base "reads" sliding across transcripts.
    """
    k2 = k - 1
    reads = [""] * (len(seq) - k2)
    for i in range(len(seq) - k2):
        reads[i] = seq[i : i + k]
    return reads


def read_to_fasta(reads, transcript_name, k):
    """
    read_to_fasta converts reads to fasta format for each transcript
    """
    fasta = [""] * (len(reads) * 2)
    j = 0
    for n, s in enumerate(reads):
        # info listed before each read.
        # position of alignment in SAM file is 1-based. So position is n+1
        fasta[j] = ">{}.{}".format(transcript_name, str(n + 1))
        fasta[j + 1] = s
        j = j + 2
    return fasta


def kmer_gen(k, reference, output_filename):
    """
    inputs - k (int), reference transcriptome file, output k-mer dataset filename
    kmer_gene generates the k-mer fasta formatted dataset and saves the reads as output_filename
    """
    # import and parse fasta file with transcripts. returns seqs: list of sequences
    k = int(k)
    seqs = []
    ids = []
    with open(reference) as handle:
        fasta_sequences = SeqIO.parse(handle, "fasta")
        for fasta in fasta_sequences:
            seqs.append(str(fasta.seq))
            ids.append(str(fasta.id))

    # save reads to new file
    with open(output_filename, "w+") as output_handle:
        for sequence, name in zip(seqs, ids):
            fasta_reads = read_to_fasta(gen_reads(sequence, k), name, k)
            for i in fasta_reads:
                output_handle.write("{}\n".format(i))


# %%
# ==============================================================================
# // main scripts
# ==============================================================================


def generate_filenames(args, message=True):
    """
    generate filenames for rest of script from command line args
    return dictionary of filenames/arguments
    """
    # parse input fasta file name to get output filename base
    temp = os.path.basename(args.reference)
    basenoext = os.path.splitext(temp)[0]

    # could put extension check here
    # turn command line arguments into params dictionary
    params = dict(args._get_kwargs())

    new_filename = "{}-{}-mers{}".format(basenoext, str(args.k), ".fa")
    k_mer_file = os.path.join(args.out, new_filename)
    params["k_mer_file"] = k_mer_file

    # use k_mer_file as base for other filenames
    k_mer_noext, ext = os.path.splitext(k_mer_file)

    params["sam_path"] = "{}.sam".format(k_mer_noext)
    params["max_output"] = "{}-multi-mapping_kmers.fa".format(k_mer_noext)
    params["MM_sam_path"] = "{}_MMers.sam".format(k_mer_noext)
    params["MM_sam_path_noH"] = "{}_MMers_noH.sam".format(k_mer_noext)
    params["MM_sorted_bam_path"] = "{}_sorted.bam".format(k_mer_noext)
    params["depth"] = "{}_depth.txt".format(k_mer_noext)
    params["idxstats"] = "{}_idx.txt".format(k_mer_noext)
    params["MM_network_file"] = "{}-MM_network_all_connections.csv".format(k_mer_noext)
    params[
        "MM_network_file_duplicates_removed"
    ] = "{}-MM_network_unique_connections.csv".format(k_mer_noext)
    params["table_file"] = "{}-multi-mapping_transcripts_table.csv".format(k_mer_noext)

    if message:
        parameter_statement = Template(
            """OUTPUT FILES:
- k-mer filename:\n    $k_mer_file
- multi-mapping (MMing) k-mers:\n    $max_output
- network csv file:\n    $MM_network_file
- network csv file with duplicate connections removed:\n    $MM_network_file_duplicates_removed
- mulit-mapping transcripts table:\n    $table_file

temperary files:
- k-mer alignment sam (removed):\n    $sam_path
- MMing k_mer alignment sam (removed):\n    $MM_sam_path
- MMing k_mer alignment sam, no header (removed):\n    $MM_sam_path_noH
- MMing k_mer alignment BAM, sorted (removed):\n    $MM_sorted_bam_path
- samtools depth (removed):\n    $depth
- samtools idxstats (removed):\n    $idxstats
        """
        )
        print(parameter_statement.substitute(params))
    return params


def bowtie_alignment1(parameter_dictionary):
    """
    args - argparse parsed command line arguments
    k_mer_file - path to fasta formatted k-mer reads
    """
    # TODO - change to log file and display

    # bowtie alignment 1
    bowtie_command1 = Template(
        "$b -S -p $threads --norc -v $mismatches -f -m 1 --max $max_output $index $k_mer_file $sam_path"
    )
    print(
        "running bowtie alignment with k-mers:\n"
        + bowtie_command1.substitute(parameter_dictionary)
    )
    subprocess.call(bowtie_command1.substitute(parameter_dictionary), shell=True)
    # subprocess.call("rm {}".format(parameter_dictionary["sam_path"]), shell=True)

    # bowtie alignment 2
    bowtie_command2 = Template(
        "$b -S -p $threads --norc -v $mismatches -f -a $index $max_output $MM_sam_path"
    )  # realign just multi-mapping k-mers
    print(
        "running bowtie alignment with multi-mapping k-mers: \n"
        + bowtie_command2.substitute(parameter_dictionary)
    )
    subprocess.call(bowtie_command2.substitute(parameter_dictionary), shell=True)

    # need sorted bam file for samtools depth
    print("converting to BAM, sorting, and indexing")
    subprocess.call(
        "samtools view -u {} | samtools sort -o {}".format(
            parameter_dictionary["MM_sam_path"],
            parameter_dictionary["MM_sorted_bam_path"],
        ),
        shell=True,
    )
    subprocess.call(
        "samtools index {}".format(parameter_dictionary["MM_sorted_bam_path"]),
        shell=True,
    )

    # get coverage depth
    subprocess.call(
        "samtools depth {} > {}".format(
            parameter_dictionary["MM_sorted_bam_path"], parameter_dictionary["depth"]
        ),
        shell=True,
    )

    # calculate idxstats
    subprocess.call(
        "samtools idxstats {} > {}".format(
            parameter_dictionary["MM_sorted_bam_path"], parameter_dictionary["idxstats"]
        ),
        shell=True,
    )

    # remove header from sam file
    subprocess.call(
        "samtools view -o {} {}".format(
            parameter_dictionary["MM_sam_path_noH"], parameter_dictionary["MM_sam_path"]
        ),
        shell=True,
    )


def multi_map_network(alignment_filename):
    df = pd.read_table(alignment_filename, sep="\t", header=None)
    df = df[[0, 2, 3]]
    df.columns = ["read", "gene", "position (1-based leftmost mapping POSition)"]
    read_split = df["read"].str.split(".")
    df["read_origin"] = read_split.str[0]
    df["read_origin_position"] = read_split.str[1].astype(int)
    df["pos_diff"] = (
        df["position (1-based leftmost mapping POSition)"] - df["read_origin_position"]
    )
    df["map_pos=read_pos?"] = df["pos_diff"] == 0

    # self_maps = transcript of origin and mapped transcript are equal.
    # Includes correct maps and internal multi-maps
    self_maps = df[df["gene"] == df["read_origin"]]
    # true selfMM = internal multi-maps
    true_selfMM = self_maps[self_maps["map_pos=read_pos?"] == False]
    true_selfMM2 = true_selfMM[["gene", "read_origin"]].copy()
    otherMM = df[df["gene"] != df["read_origin"]]
    MMnetwork = pd.concat(
        [otherMM[["gene", "read_origin"]], true_selfMM2], ignore_index=True
    )
    MMnetwork["alignments"] = "blank"
    network2 = MMnetwork.groupby(["gene", "read_origin"]).count().reset_index()
    network2.columns = ["gene", "read_origin", "MM-kmer-alignments"]
    return network2


def rm_duplicates(network):
    # remove duplicate connections
    test = []
    for i in network.index:
        j = network.loc[i, "gene"]
        k = network.loc[i, "read_origin"]
        if [k, j, network.loc[i, "MM-kmer-alignments"]] not in test:
            test.append(list(network.loc[i, :]))
    network3 = pd.DataFrame(test, columns=["gene", "read_origin", "MM-kmer-alignments"])
    return network3


def importidx(filename):
    df = pd.read_table(
        filename,
        sep="\t",
        header=None,
        names=["transcript", "length", "alignments", "unmapped"],
    )
    df = df[df["transcript"] != "*"]
    df = df.drop(["unmapped", "alignments"], axis=1)
    df = df.reset_index(drop=True)
    return df


def importdepth(filename):
    depth = pd.read_table(
        filename, header=None, names=["transcript", "position", "coverage"]
    )
    depth = depth.drop("position", axis=1)
    return depth


def num_of_MM(df):
    df2 = df[df["coverage"] > 0]
    counts = df2.groupby("transcript").count()
    counts = counts.reset_index()
    counts.columns = ["transcript", "MM_bps"]
    return counts


def fractionMM(depth_file, idxstats_file):
    depth = importdepth(depth_file)
    depth_count = num_of_MM(depth)
    idx = importidx(idxstats_file)
    dep_idx = pd.merge(depth_count, idx, on="transcript")
    dep_idx["percentMM"] = (dep_idx["MM_bps"] / dep_idx["length"]) * 100
    return dep_idx


# TODO - print fraction of transcriptome that multi-maps
# TODO - external multi-maps -> int
# TODO - turn multimaps to multi-maps
def build_table(depth_file, idxstats_file, MM_network_file):
    # create name of multi-mapper table csv file

    # build network using samtools depth, samtools idxstats, and network file
    net = pd.read_csv(MM_network_file)
    net = net.drop("MM-kmer-alignments", axis=1)
    otherMM = net[net["gene"] != net["read_origin"]]
    internalMM = net[net["gene"] == net["read_origin"]]
    table = otherMM.groupby("gene").count()
    table = table.reset_index()
    internalMM = internalMM.drop(["read_origin"], axis=1)
    internalMM["internal multi-maps?"] = "yes"
    table = pd.merge(table, internalMM, how="outer")
    table["internal multi-maps?"] = table["internal multi-maps?"].fillna("no")
    table["read_origin"] = table["read_origin"].fillna(0)
    table.columns = ["transcript", "external multimaps", "internal multi-maps?"]
    frac = fractionMM(depth_file, idxstats_file)
    table = pd.merge(frac, table, on="transcript")
    table = table.drop("MM_bps", axis=1)
    table.columns = [
        "transcript",
        "length (bp)",
        "percent of transcript that multimaps",
        "external multimaps",
        "internal multi-maps?",
    ]
    table = table.round({"percent of transcript that multimaps": 2})
    table = table.sort_values("percent of transcript that multimaps", ascending=False)
    return table


def clean_files(file_list):
    for file in file_list:
        subprocess.call("rm {}".format(file), shell=True)


def driver(params):
    # make output directory if it doesn't exist
    subprocess.call("mkdir {}".format(params["out"]), shell=True)
    print("generating transcriptome k-mers...")
    kmer_gen(params["k"], params["reference"], params["k_mer_file"])
    bowtie_alignment1(params)

    print("creating multi-mapping network...")
    network = multi_map_network(params["MM_sam_path_noH"])
    network.to_csv(params["MM_network_file"], index=False)

    network3 = rm_duplicates(network)
    network3.to_csv(params["MM_network_file_duplicates_removed"], index=False)

    print("building multi-mapper table...")
    table = build_table(params["depth"], params["idxstats"], params["MM_network_file"])
    table.to_csv(params["table_file"], index=False)

    print("cleaning up temporary files...")
    erase_list = [
        params[i]
        for i in [
            "MM_sam_path_noH",
            "MM_sorted_bam_path",
            "depth",
            "idxstats",
            "MM_sam_path",
            "sam_path",
        ]
    ]
    erase_list.append(params["MM_sorted_bam_path"] + ".bai")
    clean_files(erase_list)
    print("done!")


def main():
    """
    parse command line arguments
    execute functions with inputs from command line
    """
    parser = argparse.ArgumentParser(
        description="""identify multi-mapping transcripts and their connectivity in a reference transcriptome fasta file"""
    )
    parser.add_argument(
        "-b",
        metavar="<bowtie path>",
        default="bowtie",
        help="""the path to the bowtie executable. Default: bowtie in PATH""",
    )
    parser.add_argument(
        "-v",
        "--mismatches",
        metavar="<INT>",
        default=2,
        type=int,
        help="""maximum number of mismatches allowed in alignments (given directly to bowtie parameter -v). Default: 2""",
    )
    parser.add_argument(
        "-p",
        "--threads",
        metavar="<INT>",
        default=1,
        type=int,
        help="""Number of alignment threads for bowtie to use (given directly to bowtie parameter -p). Default: 1""",
    )
    parser.add_argument(
        "-out",
        metavar="<directory>",
        default="./",
        help="""Directory where script outputs are saved. Default is current directory""",
    )
    parser.add_argument(
        "-i",
        "--index",
        required=True,
        metavar="<ebwt_base>",
        help="""location of bowtie index made from input reference transcriptome.
        Use bowtie-build <reference_in> <ebwt_outfile_base>
        example:
        bowtie-build ./example_transcriptome.fasta ./index/example_transcriptome
        future:
        make argument optional and build index in ./index by default""",
    )
    parser.add_argument(
        "-k",
        metavar="<INT>",
        required=True,
        type=int,
        help="""k-mer read lengths (i.e. k). Recommend k be equal to the smallest read length in an experimental dataset""",
    )
    parser.add_argument(
        "-ref",
        "--reference",
        metavar="<file>",
        required=True,
        help="""input reference transcriptome file in fasta format""",
    )
    args = parser.parse_args()
    print(
        """
Running findMM with the following parameters:
     - bowtie path: '{}'
     - number of allowed mismatches (-v): '{}'
     - threads used in alignment: '{}'
     - output directory: '{}'
     - bowtie index: '{}'
     - k: '{}'
     - reference transcriptome: '{}'
    """.format(
            args.b,
            args.mismatches,
            args.threads,
            args.out,
            args.index,
            args.k,
            args.reference,
        )
    )
    params = generate_filenames(args)
    driver(params)


if __name__ == "__main__":
    main()
