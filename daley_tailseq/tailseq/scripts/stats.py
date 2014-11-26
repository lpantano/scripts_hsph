from tailseq.stats import summarize_stats
from argparse import ArgumentParser


if __name__ == "__main__":
    parser = ArgumentParser(description="Run a single cell analysis.")
    parser.add_argument("--count-file", required=True, help="Sample map file.")
    parser.add_argument("--summary-log-file", required=True, help="Sample map file.")
    parser.add_argument("--polya-log-file", required=True, help="Sample map file.")
    parser.add_argument("--out-file", required=True, help="Sample map file.")

    args = parser.parse_args()
    print args
    summarize_stats(args.count_file, args.summary_log_file,
                           args.polya_log_file, args.out_file)

