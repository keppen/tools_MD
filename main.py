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


def main(*arguments):
    """
    Docstring
    """
    option_actions = {
        "-h": lambda: print(main.__doc__) or exit(0),
        "--help": lambda: print(main.__doc__) or exit(0),
        "readlog": lambda: read_log(arguments),
        "run": lambda: run(arguments),
        "rerun": lambda: rerun(arguments)
    }

    if arguments[1] in option_actions:
        option_actions[arguments[1]]()


def read_log(arguments):

    structure = None
    cluster = None
    plot_name = None
    plot_labels = None
    plot_limits = None
    plot_resolution = None
    plot_type = None
    data_id = None

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
        if opt in ["-h", "--help"]:
            print(read_log.__doc__)
            exit(0)

    print_log(
        log_file,
        structure=structure,
        cluster=cluster,
        plot_name=plot_name,
        plot_labels=plot_labels,
        plot_limits=plot_limits,
        plot_resolution=plot_resolution,
        plot_type=plot_type,
        data_id=data_id,
    )

    exit(0)


def run(arguments):

    truncate = None
    step = None
    options = None
    no_plot = None
    to_csv = None
    log_file = None

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
        if opt in ["-h", "--help"]:
            print(Visualize.__doc__)
            exit(0)

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


def rerun(arguments):
    for index, opt in enumerate(arguments):
        if opt in ["-l", "--log-file"]:
            log_file = arguments[index + 1]
        if opt in ["-id", "--data-id"]:
            data_id = arguments[index + 1]
        if opt in ["-h", "--help"]:
            print(rerun.__doc__)
            exit(0)

    logs = open(log_file)
    logs_lines = logs.readlines()
    logs.close()

    for entry in logs_lines:
        fileds = entry.split(';')
        if int(data_id) == str(fileds[0]):
            run_arguments = fileds
            break


def print_log(log_file, **kargs):
    with open(log_file, 'r') as file:
        titles = [
                "Data id",
                "Structure",
                "Cluster",
                "Truncate",
                "Step",
                "Plot type",
                "plot name",
                "Plot labels",
                "Plot limits",
                "Plot res",
                "Data npz"
                    ]
        widths = [8, 16, 25, 9, 9, 12, 18, 18, 25, 19, 18]

        msg = [f"{t:<{s}}" for t, s in zip(titles, widths)]
        print("".join(msg))
        print("-" * sum(widths))

        logs = file.readlines()
        msg = ""
        for log in logs:
            fields = [field.strip() for field in log.split(';')]
            formatted_fields = [f"{field:<{width}}" for field, width in zip(fields, widths)]
            msg += "".join(formatted_fields) + '\n'
        print(msg)


if __name__ == "__main__":
    main(*sys.argv)
