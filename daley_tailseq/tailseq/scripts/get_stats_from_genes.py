from argparse import ArgumentParser
from collections import defaultdict

class Gene:
    def __init__(self):
        self.seq = []
        self.qual = []
        self.counts = 0
    def add(self, seq, qual):
        self.seq.append(seq)
        self.qual.append(qual)
        self.counts += 1
    def get_seq(self):
        seq = [0] * 150
        for s in self.seq:
            for i, letter in enumerate(s):
                if letter == "T":
                    seq[i] += 1
        return " ".join([str(1.0 * s / self.counts) for s in seq])
    def get_qual(self):
        seq = [0] * 150
        for q in self.qual:
            for i, phred in enumerate(q):
                seq[i] += ord(phred)-33
        return " ".join([str(1.0 * s / self.counts) for s in seq])


def read_genes(fn_file):
    data = []
    with open(fn_file) as in_handle:
        for line in in_handle:
            data.append(line.strip())
    return data


if __name__ == "__main__":
    parser = ArgumentParser(description="Run a single cell analysis.")
    parser.add_argument("--quality-log-file", required=True, help="Sample map file.")
    parser.add_argument("--genes", required=True, help="Sample map file.")
    parser.add_argument("--out-file", required=True, help="Sample map file.")

    args = parser.parse_args()

    genes = read_genes(args.genes)
    data = defaultdict(Gene)
    with open(args.quality_log_file) as in_handle:
        for line in in_handle:
            cols = line.strip().split(" ")
            if cols[0] in genes:
                data[cols[0]].add(cols[3], cols[4])
    with open(args.out_file, "w") as out_handle:
        for g in data:
            out_handle.write("seq %s\n" % " ".join([g, data[g].get_seq()]))
            out_handle.write("qual %s\n" % " ".join([g, data[g].get_qual()]))
