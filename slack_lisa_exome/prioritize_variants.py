import os
import sys
import vcf
import gzip
import pandas as pd
from collections import defaultdict

from argparse import ArgumentParser

def _moderate(record):
    if record.find("MODERATE") > -1 or record.find("HIGH") > -1:
        return True

def _gtp(info):
    return info.split(":")[0]

def _info_sample(name, info, tag):
    k = name.split(":")
    v = info.split(":")
    d = dict(zip(k, v))
    return d[tag]

def _is_in_any_control(samples, t, n):
    c = set()
    for s in n:
       c.add(_gtp(samples[s]))
    for s in t:
        if _gtp(samples[s]) in c:
            return True
    return False

def _get_info(info):
    cols = info.split(";")
    effect = []
    genes = []
    for field in cols:
        if field.find("=") == -1:
            continue
        k, v = field.split("=")
        if k == "STATUS":
            status = v
        elif k == "EFF":
            for e in v.split(","):
                type_e = e.split("(")[0]
                gene_e = e.split("(")[1].split("|")[5]
                level_e = e.split("(")[1].split("|")[0]
                effect.append("_".join([type_e, level_e]))
                genes.append(gene_e)
    return [status, effect, genes]

def _frmt(info, t, n, genes=None):
    """Get information for the variant"""
    in_gene = "No"
    cols = info.split("\t")
    samples = cols[9:]
    chrom, pos, rs = cols[:3]
    change = "%s:%s" % (cols[3], cols[4])
    status, eff, gene = _get_info(cols[7])
    effect = ";".join(list(set(["%s=%s" % (e, g) for e, g in zip(eff, gene)])))
    if any([g in genes for g in gene]):
        in_gene = "Yes"
    gt_t = ":".join([_gtp(samples[s]) for s in t])
    gt_n = ":".join([_gtp(samples[s]) for s in n])
    af_t = ":".join([_info_sample(cols[8], samples[s], "AF") for s in t])
    af_n = ":".join([_info_sample(cols[8], samples[s], "AF") for s in n])
    line = ("{chrom}\t{pos}\t{rs}\t{change}\t{gt_t}|{gt_n}\t{af_t}|{af_n}\t{status}\t{effect}\t{in_gene}")
    return line.format(**locals())

def _read_genes(fn):
    if fn:
        df = pd.read_csv(fn)
        return set(df["Gene"])

def _add_to_genes(record, t, n, d):
    cols = record.strip().split("\t")
    samples = cols[9:]
    status, eff, gene = _get_info(cols[7])
    for s in [0, 1, 2]:
        if not _gtp(samples[t[s]]) == "./." and not _gtp(samples[n[s]]) == "./.":
            if _gtp(samples[t[s]]) != _gtp(samples[n[s]]):
                for g in gene:
                    if g not in d:
                        d[g].update({s: []})
                    if s not in d[g]:
                        d[g].update({s: []})
                    d[g][s].append(record)
    return d

def _clean_variants(fn, genes):
    """
    Read variants and keep the ones in controls
    """
    t = [0, 2, 4]
    n = [1, 3, 5]
    genes = _read_genes(genes)
    print "\t".join(["chrom", "pos", "rs", "change", "gt:tumor|normal",
                     "af:tumor|normal", "status", "effects",
                     "in_gene_list", "num_samples_seen", "filter"])
    var_by_genes = defaultdict(dict)
    with gzip.open(fn[0], 'rb') as vcf_reader:
        for record in vcf_reader:
            skip = 0
            if not record.startswith("#"):
                samples = record.strip().split("\t")[9:]
                for s in [0, 1, 2]:
                    skip += _gtp(samples[t[s]]) == "./." or _gtp(samples[n[s]]) == "./."
                if skip > 1:
                    if _moderate(record):
                        var_by_genes = _add_to_genes(record, t, n, var_by_genes)
                    continue
                # if _is_in_any_control(samples, t, n):
                #    continue
                frmt = _frmt(record, t, n, genes)
                print "%s\t%s\t%s" % (frmt, 3 - skip, "variant")
    var_seen = set()
    for g in var_by_genes:
        if len(var_by_genes[g]) > 2:
            for s in var_by_genes[g]:
                for var in var_by_genes[g][s]:
                    var_frmt = _frmt(var, t, n, genes)
                    if var_frmt not in var_seen:
                        print "%s\t%s\t%s" % (var_frmt, 1, "gene")
                        var_seen.add(var_frmt)

if __name__ == "__main__":
    parser = ArgumentParser(description="Detect enhancer activity.")
    parser.add_argument("--genes", help="List of genes to annotate.")
    parser.add_argument("files", nargs="*", help="vcf files.")
    args = parser.parse_args()

    clean_vcf = _clean_variants(args.files, args.genes)
