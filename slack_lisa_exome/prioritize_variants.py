import os
import sys
import vcf
import gzip

from argparse import ArgumentParser


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

def _frmt(info, t, n):
    """Get information for the variant"""
    cols = info.split("\t")
    samples = cols[9:]
    chrom, pos, rs = cols[:3]
    change = "%s:%s" % (cols[3], cols[4])
    status, eff, gene = _get_info(cols[7])
    effect = ";".join(list(set(["%s=%s" % (e, g) for e, g in zip(eff, gene)])))
    gt_t = ":".join([_gtp(samples[s]) for s in t])
    gt_n = ":".join([_gtp(samples[s]) for s in n])
    af_t = ":".join([_info_sample(cols[8], samples[s], "AF") for s in t])
    af_n = ":".join([_info_sample(cols[8], samples[s], "AF") for s in n])
    line = ("{chrom}\t{pos}\t{rs}\t{change}\t{gt_t}|{gt_n}\t{af_t}|{af_n}\t{status}\t{effect}")
    return line.format(**locals())

def _clean_variants(fn):
    """
    Read variants and keep the ones in controls
    """
    t = [0, 2, 4]
    n = [1, 3, 5]
    print "\t".join(["chrom", "pos", "rs", "change", "gt:tumor|normal",
                     "af:tumor|normal", "status", "effects", "num_samples_seen"])
    with gzip.open(fn[0], 'rb') as vcf_reader:
        for record in vcf_reader:
            skip = 0
            if not record.startswith("#"):
                samples = record.strip().split("\t")[9:]
                for s in [0, 1, 2]:
                    skip += _gtp(samples[t[s]]) == "./." or _gtp(samples[n[s]]) == "./."
                if skip > 1:
                    continue
                # if _is_in_any_control(samples, t, n):
                #    continue
                frmt = _frmt(record, t, n)
                print "%s\t%s" % (frmt, 3 - skip)

if __name__ == "__main__":
    parser = ArgumentParser(description="Detect enhancer activity.")
    parser.add_argument("files", nargs="*", help="vcf files.")
    args = parser.parse_args()

    clean_vcf = _clean_variants(args.files)
