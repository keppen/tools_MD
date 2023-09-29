import os
import numpy as np

directory_path = "."  # Change this to the directory path you want to search

chirality_dict = {
        "ou1": 'sssss',
        "ou2": 'srsss',
        "ou3": 'ssrss',
        "ou4": 'sssrs',
        "ou5": 'srrss',
        "ou6": 'srsrs',
        "ou7": 'ssrrs',
        }


new_dict = {}


def read_npz(loaded_data):
    # loaded_data = np.load(npz_file.strip())

    density = {key: loaded_data[key] for key in loaded_data.files if "dens" in key}
    coords = {key: loaded_data[key] for key in loaded_data.files if "coords" in key}
    return density, coords


def get_run(filename):
    run_id = [char for char in filename if char.isdigit()]
    run_id = "".join(run_id)
    log_file = open(directory_path + "/log_file.log", 'r')
    logs = log_file.readlines()
    log_file.close()

    run_data = None

    for log in logs:
        log = log.strip().split(";")
        if log[0] == run_id:
            run_data = log
            break

    if not run_data:
        return None

    if run_data[5] == "rama_1d":
        return run_data[6]
    return None


def plt_distribution(data_dict: dict, limits=None, name="max_prob.png",):
    # import necessary libraries and modules
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import matplotlib.font_manager as font_manager
    import matplotlib
    from matplotlib.gridspec import GridSpec
    from utils import latin_2_greek

    data_keys = data_dict.keys()
    data_keys = [key.split() for key in data_keys]
    data_keys = zip(*data_keys)
    data_keys = [set(keys) for keys in data_keys]
    var1_names, var2_names = data_keys
    var1_names = sorted(var1_names)
    var2_names = sorted(var2_names)

    print(var1_names, var2_names)

    vnum_plots = len(var1_names)
    hnum_plots = len(var2_names)

    # Configure font and font size for the plot
    font_path = "/home/mszatko/.local/share/fonts/ARIALUNI.TTF"
    custom_font = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = custom_font.get_name()
    plt.rcParams['font.size'] = 8

    # Create a figure and subplots
    fig = plt.figure(figsize=(3.45, 3.45/2), layout='constrained')
    GS = GridSpec(1, vnum_plots, figure=fig)

    # fig, axes = plt.subplots(
    #         1, num_plots,
    #         figsize=(3.45, 3.45),
    #         sharey=True,
    #         layout="constrained"
    #         )

    # Change lables to greek
    labels = latin_2_greek(var1_names)
    print(labels, limits)

    # Set colors
    # colormap = matplotlib.colormaps["viridis"]
    # bins = [x/(len(var2_names) - 1) for x in range(len(var2_names))]
    # colors = [colormap(x) for x in bins]
    colors = ["green", "blue"]
    markers = ["x", "+"]

    for i, var1_name in enumerate(var1_names):

        ax = fig.add_subplot(GS[0, i])
        # ax = axes[i]

        # Plot the data
        for j, var2_name in enumerate(var2_names):
            y = data_dict[f"{var1_name} {var2_name}"]
            x = np.zeros(len(y))
            if j == 0:
                x += 0.1
            else:
                x -= 0.1
            print(x, y)
            ax.scatter(
                    x,
                    y,
                    # coordinates_dict[f"{var1_name} {var2_name} coords"][0],
                    color=colors[j],
                    linewidth=0.75,
                    # label=f"{var2_name.upper()}",
                    s=5,
                    marker=markers[j],
                    alpha=0.3,
                    lw=0.5,
                    )

        # Set labels and title for each subplot
        ax.set_xlabel(labels[i])
        if i == 0:
            ax.set_ylabel("Torional angle [\u00b0]")

        # Extract limits from variable
        xg = limits[i]
        xmin, xmax = xg

        # Set the limits of the x and y axes to (-180, 180)
        plt.ylim(xmin, xmax)
        plt.xlim(-1.5, 1.5)

        # Customize axes properties
        ax.spines['bottom'].set_linewidth(0.5)
        ax.spines['left'].set_linewidth(0.5)
        ax.spines['top'].set_linewidth(0.5)
        ax.spines['right'].set_linewidth(0.5)

        # Customize tick parameters for major ticks
        ax.set_xticks([])

        num_ticks = (xmax - xmin) / 6
        major_locator = ticker.MultipleLocator(
            base=num_ticks)  # Set major tick every 30 degrees
        ax.yaxis.set_major_locator(major_locator)
        ax.yaxis.set_tick_params(width=0.5, length=3,
                                        direction='out', which='both', right=False)

        # Customize tick parameters for minor ticks
        minor_locator = ticker.MultipleLocator(
            base=num_ticks/4)  # Set minor tick every 10 degrees
        ax.yaxis.set_minor_locator(minor_locator)
        ax.yaxis.set_tick_params(width=0.3, length=1.5,
                                        direction='out', which='minor', right=False)

        # Set the position of the subplot
        ax.set_position(GS[i].get_position(fig))

    # Adjust layout to avoid overlapping labels
    # plt.tight_layout(rect=[0, 0, 1, 1])

    # Create a separate scatter plot with opaque markers for the legend
    legend_marker_1 = plt.scatter(
            [], [],
            label=f"{var2_names[0].upper()}",
            s=5,
            marker=markers[0],
            lw=.5,
            alpha=1.0,
            color=colors[0]
            )

    legend_marker_2 = plt.scatter(
            [], [],
            label=f"{var2_names[1].upper()}",
            s=9,
            marker=markers[1],
            lw=.5,
            alpha=1.0,
            color=colors[1]
            )


    # Add a legend above the plot figure in the upper right corner
    plt.legend(
            loc='upper right',
            bbox_to_anchor=(1.05, 1.24),
            ncol=len(labels),
            frameon=False,
            columnspacing=0.7,
            )


    # Increase spacing between subplots and the surrounding area
    plt.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.15, wspace=0.6, hspace=0.5)

    # Show the plot
    # plt.show()

    # Save the plot as an image
    plt.savefig(
        name + ".png",
        dpi=1200
    )
    # Close the plot
    plt.close()


def get_chirality(filename):
    for code in chirality_dict:
        if code in run_name:
            return chirality_dict[code]


def get_maximum(data: dict, coords: dict, chirality):
    print(chirality)
    for key, item in data.items():
        angle, res_id, _ = key.split()
        new_key = f"{angle} {chirality[int(res_id) - 1]}"
        key_to_coords = f"{angle} {res_id} coords"
        max_item_index = np.argmax(item)
        if new_key not in new_dict:
            new_dict[new_key] = []
        new_dict[new_key].append(coords[key_to_coords][0][max_item_index])
    # print(new_dict)


for filename in os.listdir(directory_path):
    if filename.endswith(".npz"):
        file_path = os.path.join(directory_path, filename)
        run_name = get_run(filename)
        if not run_name:
            continue
        chirality = get_chirality(filename)
        data = np.load(file_path)
        density, coords = read_npz(data)
        get_maximum(density, coords, chirality)
        
print(new_dict)
        # Process the loaded data as needed

plt_distribution(data_dict=new_dict, limits=3*[[-180, 180]])
