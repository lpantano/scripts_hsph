import os
import logger
from cluster_helper import cluster as ipc


resources = {"align": [45, 12],
             "qc": [8, 1],
             "cleaning": [4, 1],
             "polyA": [4, 1]}


def get_cluster_view(args):
    return ipc.cluster_view(args.scheduler, args.queue,
                          args.num_jobs, args.cores_per_job,
                          start_wait=args.timeout,
                          extra_params={"resources": args.resources,
                                        "mem": args.memory_per_job,
                                        "tag": "ts",
                                        "run_local": args.local})


def wait_until_complete(jobs):
    return [j.get() for j in jobs]


def is_done(step):
    if os.path.exists(step):
        return True
    return False


def flag_done(step):
    with open(step,"w") as handle:
        handle.write("done")


def send_job(fn, data, args, step):
    """decide if send jobs with ipython or run locally"""
    res = []
    logger.my_logger.debug("doing %s" % step)
    if step not in resources:
        raise ValueError("step not in resources %s" % step)
    else:
        args.memory_per_job = resources[step][0]
        args.cores_per_job = resources[step][1]
    if args.parallel == "ipython":
        if not is_done(step):
            with get_cluster_view(args) as view:
                for sample in data:
                    res.append(view.apply_async(fn, *data[sample]))
            res = wait_until_complete(res)
            flag_done(step)
            return res
    for sample in data:
        res.append(fn(*data[sample]))
    return res
