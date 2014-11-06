from cluster_helper import cluster as ipc

def get_cluster_view(args):
    return ipc.cluster_view(args.scheduler, args.queue,
                          args.num_jobs, args.cores_per_job,
                          start_wait=args.timeout,
                          extra_params={"resources": args.resources,
                                        "mem": args.memory_per_job,
                                        "tag": "tailseq",
                                        "run_local": args.local})


def wait_until_complete(jobs):
    return [j.get() for j in jobs]


def send_jobs(fn, data, args):
    """decide if send jobs with ipython or run locally"""
    res = []
    if args.parallel:
        with get_cluster_view(args) as view:
            for sample in data:
                res.append(view.apply_async(fn, *data[sample]))
        res = wait_until_complete(res)
    else:
        for sample in data:
            res.append(fn(*data[sample]))
    return res
