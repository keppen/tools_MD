# !/home/keppen/miniconda3/bin/python3
import sys
import time
import pandas as pd
import numpy as np
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
Entry point function for the program.

Parameters:
    run             | Run the analysis.
    rerun           | Rerun from a run from log entry.
    readlog         | Print formatted logfile.
    -h, --help      | Display this help message and exit.
    """

    option_actions = {
        "-h": lambda: print(main.__doc__),
        "--help": lambda: print(main.__doc__),
        "readlog": lambda: read_log(arguments),
        "run": lambda: run(arguments),
        "rerun": lambda: rerun(arguments)
    }

    if arguments[1] in option_actions:
        option_actions[arguments[1]]()


def read_log(arguments):
    """
Function to read and process log file.

Parameters:
    -l, --log-file      | Log file name.
    -s, --structure     | Structure pdb file name.
    -c, --cluster       | Cluster pdb file name.
    -rd, --raw-data     | Use raw data for plotting. Flag
    -pm, --plot-name    | Filter by name of the plot.
    -lb, --plot-labels  | Filter by labels for the plot.
    -lm, --plot-limits  | Filter by limits for the plot.
    -rs, --resolution   | Filter by resolution for the plot.
    -pt, --plot-type    | Filter by type of the plot.
    -id, --data-id      | Filter by data ID.
    -h, --help          | Display this help message and exit.
s
    Filters does not work now!
    """

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
    """
Function to run analysis and visualization.

Parameters:
    -t, --truncate   | Truncate option.
    -o, --options    | Options for the analysis.
    --step           | Step option.
    --no-plot        | Disable plotting.
    -l, --log-file   | Path to the log file.
    -s, --structure  | Path to the structure file.
    -c, --cluster    | Cluster data.
    --no-csv         | Disable CSV output.
    -h, --help       | Display this help message and exit.
    """

    truncate = None
    step = None
    options = None
    no_plot = False
    to_csv = True
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
    """
Function to rerun analysis from log entry.

Parameters:
    -l, --log-file     | Path to the log file.
    -id, --data-id     | Data ID for the log entry.
    -dim, --dimension  | Plot dimension.
    -n, --plot_name    | Plot name.
    -h, --help         | Display this help message and exit.
    """

    plot_dim = None
    plot_name = None

    for index, opt in enumerate(arguments):
        if opt in ["-l", "--log-file"]:
            log_file = arguments[index + 1]
        if opt in ["-id", "--data-id"]:
            data_id = arguments[index + 1]
        if opt in ["-dim", "--dimmention"]:
            plot_dim = arguments[index + 1]
        if opt in ["-n", "--plot_name"]:
            plot_name = arguments[index + 1]
        if opt in ["-h", "--help"]:
            print(rerun.__doc__)
            exit(0)

    logs = open(log_file)
    logs_lines = logs.readlines()
    logs.close()

    for entry in logs_lines:
        fileds = entry.split(';')
        if str(data_id) == str(fileds[0]):
            run_opts = fileds
            break

    if "rawdata" in run_opts:
        from_raw_data(run_opts, plot_name=plot_name,
                      plot_dim=plot_dim, log_file=log_file)
    else:
        from_processed_data(run_opts, plot_name=plot_name,
                            plot_dim=plot_dim, log_file=log_file)


def from_processed_data(args, **kargs):
    from libplot import mavi_contour, plt_ramachandran, plt_distribution

    plotting_functions = {
        '1': lambda density, coords, limits, name, labels: plt_distribution(
            density, coords, limits=limits, name=name, labels=plot_labels
            ),
        '2': lambda density, coords, limits, name, labels: plt_ramachandran(
            density, coords, limits=limits, name=name, labels=plot_labels
            ),
        '3': lambda density, coords, limits, name, labels: mavi_contour(
            density, coords, limits=limits, name=name, labels=plot_labels
            ),
    }

    plot_type = args[5]
    plot_dim = kargs["plot_dim"] or plot_type[-2]
    plot_dim = int(plot_dim)

    plot_name = args[6]
    plot_labels = args[7].split()
    # print(plot_labels)

    plot_limits = args[8].split()
    plot_limits = np.array([int(x) for x in plot_limits])
    print(plot_limits)
    if plot_dim in [1, 3]:
        plot_limits = plot_limits.reshape(3, 2)
    else:
        plot_limits = plot_limits.reshape(plot_dim, 2)

    # print(plot_limits)

    if plot_dim == 1:
        loaded_data = np.load(args[-1].strip())

        density = {key: loaded_data[key] for key in loaded_data.files if "dens" in key}
        coords = {key: loaded_data[key] for key in loaded_data.files if "coords" in key}

    else:
        loaded_data = np.load(args[-1].strip())
        density = loaded_data["data"]
        coords = loaded_data["coords"]

    for plot in plotting_functions:
        if int(plot) == plot_dim:
            plotting_functions[plot](
                density,
                coords,
                limits=plot_limits,
                name=plot_name,
                labels=plot_labels,
            )
            exit(0)


def from_raw_data(args, **kargs):
    plot_axes_titles = {
        "dihedral": ["model", "residue index", "Phi", "Psi", "Xi", "Chi"],
        "geometry": ["model", "residue index", "Alpha", "Theta1", "Theta2"],
        "hbonds": ["model", "residue acceptor", "residue donor", "angle", "lenght"],
        "axis": ["model", "index", "distance", "resid"],
        "contact": ["models", "index 1", "index 2", "distance"],
    }

    plotting_functions = {
        'd1': lambda data: R._do_rama_1d(data=data),
        'd2': lambda data: R._do_rama_2d(data=data),
        'd3': lambda data: R._do_rama_3d(data=data),
        'g1': lambda data: R._do_geom_1d(data=data),
        'g2': lambda data: R._do_geom_2d(data=data),
        'g3': lambda data: R._do_geom_3d(data=data),
        'h': lambda data: R._do_hbond(data=data),
        'c': lambda data: R._do_contact(data=data),
        'a': lambda data: R._do_axis(data=data),
    }
    structure = args[1]
    cluster = args[2]

    # plot_type = args[5]
    plot_name = args[6]
    plot_dim = kargs["plot_dim"]

    DF_data = read_npz(args[-1].strip(), plot_axes_titles[plot_name])

    options = f"{plot_name[0]}{plot_dim or ''}"
    print(args[3])
    if args[3] == "None":
        truncate = None
    else:
        truncate = args[3]

    R = Visualize(
        structure,
        cluster,
        log_file=kargs["log_file"],
        no_plot=False,
        truncate=truncate,
        step=args[4]
    )
    R._init_log()

    for plot in plotting_functions:
        print(plot)
        if options == plot:
            plotting_functions[plot](DF_data)
            exit(0)


def read_npz(data_npz, columns=None):

    data_raw = np.load(data_npz)
    data_raw = data_raw['rawdata']
    return pd.DataFrame(data_raw, columns=columns)


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
        widths = [8, 16, 25, 9, 9, 12, 30, 25, 30, 19, 18]

        msg = [f"{t:<{s}}" for t, s in zip(titles, widths)]
        print("".join(msg))
        print("-" * sum(widths))

        logs = file.readlines()
        msg = ""
        for log in logs:
            fields = [field.strip() for field in log.split(';')]
            formatted_fields = [
                f"{field:<{width}}" for field, width in zip(fields, widths)]
            msg += "".join(formatted_fields) + '\n'
        print(msg)


if __name__ == "__main__":
    main(*sys.argv)
