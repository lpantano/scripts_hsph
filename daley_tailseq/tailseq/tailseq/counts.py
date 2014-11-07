from utils import file_transaction
import os
import gzip
from collections import defaultdict, Counter
from tailseq import do
from tailseq import logger
from tailseq.detect import tune


def counts(data, args):
    logger.my_logger.info("Counting sample %s:" % data['sample_id'])
    in_file = data['clean']
    prefix = data['sample_id'] + "_"
    gtf_file = args.gtf_file
    cores = args.cores_per_job
    out_file = prefix + "counts"
    data['counts'] = _cmd_counts(in_file, out_file, gtf_file, cores)
    data['assign'] = _assign_gene(data['counts'], prefix)
    _summarize(data['detect'], data['assign'], prefix + "summary.dat")
    return data


def _cmd_counts(in_file, out_file, gtf_file, cores):
    if not os.path.exists(out_file):
        cmd = "featureCounts -R -T {cores} --primary -a {gtf_file} -o {out_file} {in_file}"
        do.run(cmd.format(**locals()))
    return in_file + ".featureCounts"


def _assign_gene(in_file, prefix):
    """read featureCounts output and assign each read a gene"""
    out_file = prefix + "assign.dat"
    if not os.path.exists(out_file):
        with open(in_file) as handle, file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, 'w') as out:
                for line in handle:
                    cols = line.strip().split("\t")
                    if cols[1] == "Assigned":
                        out.write("%s\t%s\n" % (cols[0], cols[2]))
    return out_file


def _summarize(in_file, count_file, out_file):
    read_gene, counts_gene = _get_first_read(count_file)
    stats = defaultdict(Counter)
    with gzip.open(in_file) as handle_polya:
        for line in handle_polya:
            cols = line.strip().split("\t")
            find = tune(cols[3], cols[4])
            read = cols[0].split(" ")[0].replace("@", "")
            print "%s %s -----> %s " % (cols[3], cols[4], find)
            if find:
                if read in read_gene:
                    gene = read_gene[read]
                    stats[gene]["polyA"] += 1
                    if find[0] != "":
                        stats[gene][find[0]] += 1
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, 'w') as out:
            for gene in counts_gene:
                out.write("%s %s\n" % (gene, counts_gene[gene]))
                if gene in stats:
                    for mod, c in stats[gene].iteritems():
                        out.write("%s %s\n" % (mod, c))


def _get_first_read(in_file):
    gene = {}
    stats = Counter()
    with open(in_file) as counts:
        for line in counts:
            cols = line.strip().split("\t")
            gene[cols[0]] = cols[1]
            stats[cols[1]] += 1
    return gene, stats
