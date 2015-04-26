"""
calculate coverage across a list of regions
"""
import os

import six
from argparse import ArgumentParser
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# import matplotlib
# import seaborn as sns
from ichwrapper import cluster, arguments
import pandas as pd
from collections import Counter, defaultdict

import pybedtools

from bcbio.utils import rbind, file_exists
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction

from ecov.total import _calc_total_exome_coverage
from ecov.bias import calculate_bias_over_multiple_regions
from ecov.variants import calc_variants_stats

def calculate_genes_per_vcf(args):
    """
    count number of genes with variants
    """

def calculate_cg_depth_coverage(args):
    resources = {'name': 'vcf_stats', 'mem': 1, 'cores': 1}
    cluster.send_job(calc_variants_stats, args.bams, args, resources)

def calculate_bam(args):
    """
    samtools flagstat output
    """
    stats = defaultdict(list)
    samples = []
    for bam in args.bams:
        # sample = os.path.splitext(bam)[0].split("-")[0]
        sample = os.path.basename(bam).split("-")[0]
        out_sample = sample + ".flagstat"
        if not file_exists(out_sample):
            with file_transaction(out_sample) as tx_out:
                cmd = ("samtools flagstat {bam} > {tx_out}")
                do.run(cmd.format(**locals()), "bam stats for %s" % bam)
        with open(out_sample) as in_handle:
            for line in in_handle:
                if line.find("mapQ") == -1:
                    stats[line.strip().split(" + 0 ")[1].split("(")[0].strip()].append(line.strip().split(" + ")[0])
                # print stats[sample]
        samples.append(sample)
    with open(args.out, 'w') as out_handle:
        out_handle.write("\t".join(['measure'] + samples) + '\n')
        for feature in stats:
            out_handle.write("\t".join([feature] + stats[feature]) + "\n")

def calculate_tstv(args):
    """
    get tstv from bcftools stat for all, known and new variants
    """
    tstv = defaultdict(list)
    for in_vcf in args.bams:
        out_file = os.path.splitext(in_vcf)[0] + ".stats"
        known_file = os.path.splitext(in_vcf)[0] + ".known.stats"
        new_file = os.path.splitext(in_vcf)[0] + ".new.stats"
        sample = os.path.basename(in_vcf).split("-")[0]
        if not file_exists(out_file):
            with file_transaction(out_file) as tx_out:
                cmd = ("bcftools stats {in_vcf} > {tx_out}")
                do.run(cmd.format(**locals()), "ts/tv ratio for %s" % in_vcf)
        if not file_exists(new_file):
            with file_transaction(new_file) as tx_new:
                cmd = ("bcftools filter -i DB=0 {in_vcf} | bcftools stats /dev/stdin > {tx_new}")
                do.run(cmd.format(**locals()), "ts/tv ratio for %s" % in_vcf)
        if not file_exists(known_file):
            with file_transaction(known_file) as tx_known:
                cmd = ("bcftools filter -i DB=1 {in_vcf} | bcftools stats /dev/stdin > {tx_known}")
                do.run(cmd.format(**locals()), "ts/tv ratio for %s" % in_vcf)
        for fn, name in zip([out_file, known_file, new_file], ['all', 'known', 'new']):
            with open(fn) as in_handle:
                for line in in_handle:
                    if line.startswith("TSTV"):
                        tstv[sample].append(line.split()[4])
                        break
    df = pd.DataFrame(tstv, index=['all', 'known', 'new'])
    df.to_csv(args.out)

def bias_exome_coverage(args):
    resources = {'name': 'bias', 'mem': 1, 'cores': 1}
    cluster.send_job(calculate_bias_over_multiple_regions, args.bams, args, resources)

def average_exome_coverage(args):
    out_file = args.out
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tx_out_file:
        # dfs = [_calc_total_exome_coverage(bam, bed_file) for bam in in_bams]
        resources = {'name': 'bedtools', 'mem': 12, 'cores': 1}
        cluster.send_job(_calc_total_exome_coverage, args.bams, args, resources)
        # df = rbind(dfs)
        # df.to_csv(tx_out_file, mode='a', index=False, header=["r10", "r25", "r50", "region", "size", "sample"])
    return out_file

def _calc_regional_coverage(in_bam, chrom, start, end, samplename, work_dir):
    """
    given a BAM and a region, calculate the coverage for each base in that
    region. returns a pandas dataframe of the format:

    chrom position coverage name

    where the samplename column is the coverage at chrom:position
    """
    region_bt = pybedtools.BedTool("%s\t%s\t%s\n" % (chrom, start, end), from_string=True).saveas()
    region_file = region_bt.fn
    coords = "%s:%s-%s" % (chrom, start, end)
    tx_tmp_file = os.path.join(work_dir, "coverage-%s-%s.txt" % (samplename, coords.replace(":", "_")))
    cmd = ("samtools view -b {in_bam} {coords} | "
           "bedtools coverage -abam - -b {region_file} -d > {tx_tmp_file}")
    do.run(cmd.format(**locals()), "Plotting coverage for %s %s" % (samplename, coords))
    names = ["chom", "start", "end", "offset", "coverage"]
    df = pd.io.parsers.read_table(tx_tmp_file, sep="\t", header=None,
                                  names=names)
    os.remove(tx_tmp_file)
    df["sample"] = samplename
    df["chrom"] = chrom
    df["position"] = df["start"] + df["offset"] - 1
    return df[["chrom", "position", "coverage", "sample"]]

def _combine_regional_coverage(in_bams, samplenames, chrom, start, end, work_dir):
    """
    given a list of bam files, sample names and a region, calculate the
    coverage in the region for each of the samples and return a tidy pandas
    dataframe of the format:

    chrom position coverage name
    """
    dfs = [_calc_regional_coverage(bam, chrom, start, end, sample, work_dir) for bam, sample
           in zip(in_bams, samplenames)]
    return rbind(dfs)

def save_multiple_regions_coverage(samples, out_file, region_bed=None):
    """
    given a list of bcbio samples and a bed file or BedTool of regions,
    makes a plot of the coverage in the regions for the set of samples

    if given a bed file or BedTool of locations in stem_bed with a label,
    plots lollipops at those locations
    Adapted form Rory Kirchner (@roryk)
    """
    PAD = 100
    if file_exists(out_file):
        return out_file
    in_bams = samples
    samplenames = [os.path.splitext(os.path.basename(fn))[0] for fn in samples]
    if isinstance(region_bed, six.string_types):
        region_bed = pybedtools.BedTool(region_bed)
    with file_transaction(out_file) as tx_out_file:
        for line in region_bed:
            chrom = line.chrom
            start = max(line.start - PAD, 0)
            end = line.end + PAD
            df = _combine_regional_coverage(in_bams, samplenames, chrom,
                                            start, end, os.path.dirname(tx_out_file))
            df.to_csv(tx_out_file, mode='a', index=False, header=None)
    return out_file

if __name__ == "__main__":
    parser = ArgumentParser(description="Create file with coverage of a region")
    parser = arguments.myargs(parser)
    parser.add_argument("--region", help="bed file with regions.")
    parser.add_argument("--reference", help="genome fasta file.")
    parser.add_argument("--out", required=True, help="output file.")
    parser.add_argument("bams", nargs="*", help="Bam files.")
    parser.add_argument("--basic-bam", action="store_true", help="Calculate bam stats")
    parser.add_argument("--basic-vcf", action="store_true", help="Calculate vcf stats")
    parser.add_argument("--stats", action="store_true", help="Calculate stats for all regions.")
    parser.add_argument("--tstv", action="store_true", help="Calculate ts/tv ratio.")
    parser.add_argument("--bias", action="store_true", help="Calculate bias for all regions.")
    parser.add_argument("--plot", action="store_true", help="Calculate nt coverage for given regions.")
    parser.add_argument("--n_sample", default=1000, help="sample bed files with this number of lines")
    parser.add_argument("--seed", help="replication of sampling")
    args = parser.parse_args()

    if os.path.exists(args.out):
        os.remove(args.out)

    if args.stats:
        average_exome_coverage(args)
    elif args.bias:
        bias_exome_coverage(args)
    elif args.tstv:
        calculate_tstv(args)
    elif args.basic_bam:
        calculate_bam(args)
    elif args.basic_vcf:
        calculate_cg_depth_coverage(args)
    elif args.plot:
        save_multiple_regions_coverage(args.bams, args.out, args.region)
