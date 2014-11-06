from tailseq import cluster
from argparse import ArgumentParser


if __name__ == "__main__":
    parser = ArgumentParser(description="Run a single cell analysis.")
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
    parser.add_argument("--parallel", choices = ["local", "ipython"], default="ipython",
                        help="Run in parallel on a local machine.")
    parser.add_argument("--local", action="store_true",
                        default=True, help="Run parallel locally")

    args = parser.parse_args()

    res = []
    with cluster.get_cluster_view(args) as view:
        for sample in [1, 2, 3]:
            res.append(view.apply_async(sum, [1, 2, 3]))
    res = cluster.wait_until_complete(res)
    print res
