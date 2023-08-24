# !/home/keppen/miniconda3/bin/python3
import sys
import time

from rama_new import Visualize

truncate = None
options = None
step = None
no_plot = True
to_csv = True
log_file = None
structure = None
cluster = None

start_time = time.time()


def get_log(log_file, **kargs):
    with open(log_file, 'r') as file:
        print("{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}\
        {:<8}{:<8}{:<8}{:<8}{:<8}".format("Data id",
                                          "Structure",
                                          "Cluster",
                                          "Truncate",
                                          "Step",
                                          "Plot type",
                                          "plot name",
                                          "Plot labels",
                                          "Plot limits",
                                          "Plot res",
                                          "Data npz")
              )
        logs = file.readlines()
        for log in logs:
            log = log.split(';')
            msg = ""
            msg += f"{log[0]:<8}"
            msg += f"{log[1]:<8}"
            msg += f"{log[2]:<8}"
            msg += f"{log[3]:<8}"
            msg += f"{log[4]:<8}"
            msg += f"{log[5]:<8}"
            msg += f"{log[6]:<8}"
            msg += f"{log[7]:<8}"
            msg += f"{log[8]:<8}"
            msg += f"{log[9]:<8}"
            msg += f"{log[10]:<8}"
            msg += f"{log[11]:<8}"
            print(msg)


def main(*arguments):
    """
    Docstring
    """
    option_actions = {
        "-h": lambda: print(Visualize.__doc__) or exit(0),
        "--help": lambda: print(Visualize.__doc__) or exit(0),
        "readlog": lambda: read_log(arguments),
        "plot": lambda: run(arguments)
    }

    for arg in arguments:
        print(arg)
        if arg in option_actions:
            option_actions[arg]()


def read_log(arguments):
    for index, opt in enumerate(arguments):
        if opt in ["-l", "--log-file"]:
            log_file = arguments[index + 1]
        if opt in ["-s", "--strucutre"]:
            structure = arguments[index + 1]
        if opt in ["-c", "--cluster"]:
            cluster = arguments[index + 1]
        if opt in ["-rd", "--raw-data"]:
            plot_type = "rawdata"
        if opt in ["-pm", "--plot-name"]:
            plot_name = arguments[index + 1]
        if opt in ["-lb", "--plot-labels"]:
            plot_labels = arguments[index + 1]
        if opt in ["-lm", "--plot-limits"]:
            plot_limits = arguments[index + 1]
        if opt in ["-rs", "--resolution"]:
            plot_resolution = arguments[index + 1]
        if opt in ["-pt", "--plot-type"]:
            plot_type = arguments[index + 1]
        if opt in ["-id", "--data-id"]:
            data_id = arguments[index + 1]

    get_log(
        log_file,
        structure=structure or None,
        cluster=cluster or None,
        plot_name=plot_name or None,
        plot_labels=plot_labels or None,
        plot_limits=plot_limits or None,
        plot_resolution=plot_resolution or None,
        plot_type=plot_type or None,
        data_id=data_id or None,
    )

    exit(0)


def run(arguments):

    for index, opt in enumerate(arguments):
        if opt in ["-t", "--truncate"]:
            truncate = arguments[index + 1]
        if opt in ["-o", "--options"]:
            options = arguments[index + 1]
        if opt in ["--step"]:
            step = arguments[index + 1]
        if opt in ["--no-plot"]:
            no_plot = True
        if opt in ["-l", "--log-file"]:
            log_file = arguments[index + 1]
        if opt in ["-s", "--structure"]:
            structure = arguments[index + 1]
        if opt in ["-c", "--cluster"]:
            cluster = arguments[index + 1]
        if opt in ["--no-csv"]:
            to_csv = False

    R = Visualize(
        structure,
        cluster,
        truncate=truncate,
        step=step,
        options=options,
        no_plot=no_plot,
        to_csv=to_csv,
        log_file=log_file
    )
    R._init_run()

    exit(0)


if __name__ == "__main__":
    main(*sys.argv)
