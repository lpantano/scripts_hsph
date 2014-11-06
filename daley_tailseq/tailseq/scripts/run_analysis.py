from argparse import ArgumentParser
import os
from collections import OrderedDict
from itertools import izip
from tailseq import align
from tailseq import cluster
from tailseq import detect


def get_sample(line, sample_map_filename):
    keys = ["sample_id", "r1_path", "r2_path"]
    sample_id, r1_filename, r2_filename = line.split(",")
    return dict(zip(keys, [sample_id, r1_filename, r2_filename]))


def get_samples_to_process(sample_file):
    with open(sample_file) as in_handle:
        return [get_sample(x, sample_file) for x in in_handle]


def get_r2_prepped_outfile(sample, alignment_dir):
    return os.path.join(alignment_dir,
                        ".".join([sample["sample_id"], sample["subsample_id"]]))


def get_star_prefix(fastq_file):
    base, _ = os.path.splitext(fastq_file)
    return base


def get_cleaned_outfile(align_file):
    base, ext = os.path.splitext(align_file)
    return base + ".cleaned" + ext


if __name__ == "__main__":
    parser = ArgumentParser(description="Run a single cell analysis.")
    parser.add_argument("--multimappers", action="store_true",
                        default=False, help="Keep multimappers")
    parser.add_argument("--sample-map", required=True, help="Sample map file.")
    parser.add_argument("--aligner-index", help="Path to aligner index.")
    parser.add_argument("--alignment-dir", help="Output directory")
    parser.add_argument("--gtf-file", required=True, help="GTF file")
    parser.add_argument("--num-jobs", type=int,
                        default=1, help="Number of concurrent jobs to process.")
    parser.add_argument("--cores-per-job", type=int,
                        default=1, help="Number of cores to use.")
    parser.add_argument("--memory-per-job", default=2, help="Memory in GB to reserve per job.")
    parser.add_argument("--timeout", default=15, help="Time to wait before giving up starting.")
    parser.add_argument("--scheduler", default=None, help="Type of scheduler to use.",
                        choices=["lsf", "slurm", "torque", "sge"])
    parser.add_argument("--resources", default=None, help="Extra scheduler resource flags.")
    parser.add_argument("--queue", default=None, help="Queue to submit jobs to.")
    parser.add_argument("--parallel", choices = ["local", "ipyton"], default="local",
                        help="Run in parallel on a local machine.")

    args = parser.parse_args()

    samples = get_samples_to_process(args.sample_map)
    prepped = []

    print "Beginning alignment."
    for sample in samples:
        data[sample['sample_id']] = [sample['r1_path'], args.aligner_index, get_star_prefix(sample['r1_path']), args.cores_per_job]
    aligned = cluster.send_job(align.star_align, data, args)
    print "Finished alignment."

    print "Beginning QC."
    for sample, sam_file in izip(samples,aligned):
        data[sample['sample_id']] = [sam_file]
    qc.append(cluster.send_job(align.qc, data, args))
    print "Finished QC."

    print "Begin cleaning of poorly mapped reads."
    for sample, sam_file in izip(samples, aligned):
        data[sample['sample_id']] = [sam_file, get_cleaned_outfile(sam_file)]
    cleaned.append(cluster.send_job(align.clean_align, data, args)
    print "Finished cleaning."

    print "Detecting polyA and modifications."
    for sample in samples:
        data[sample['sample_id']] = [sample['r2_path'], sample['sample_id']]
    polyA.append(cluster.send_job(detect.detect, data, args))
