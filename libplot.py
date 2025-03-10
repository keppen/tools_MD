from numpy import float16


def mavi_contour(
    data,
    coordinates,
    limits,
    name,
    labels,
    downsample_factor=1,
    alpha_power=0.8,
    threshold=1e-8,
):
    """
    3D density visualization using an optimized scatter plot.

    - Colors are scaled manually from the minimum (nonzero) to maximum value.
    - Points with values close to 0 (below 'threshold') are omitted.
    - The alpha channel is modulated based on the original data values.
    """

    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import matplotlib.font_manager as font_manager
    import matplotlib
    from matplotlib.gridspec import GridSpec
    from utils import latin_2_greek

    import matplotlib.pyplot as plt
    import numpy as np

    # font_path = "/home/mszatko/.local/share/fonts/ARIALUNI.TTF"
    # custom_font = font_manager.FontProperties(fname=font_path)
    # plt.rcParams["font.family"] = custom_font.get_name()
    plt.rcParams["font.size"] = 8

    labels = latin_2_greek(labels)

    fig = plt.figure(figsize=(6.9, 6.9))
    ax = fig.add_subplot(111, projection="3d")

    # Downsample data and convert to float32 for sufficient precision
    ds_data = data[::downsample_factor, ::downsample_factor, ::downsample_factor]
    print("Downsampled shape:", ds_data.shape)
    ds_data = ds_data.astype(np.float32)
    coordinates = coordinates.astype(np.float32)

    # Add a small epsilon to avoid zeros (will be filtered later)
    epsilon = 1e-8
    ds_data = ds_data + epsilon

    # Unpack coordinate grids (assumed to be provided as a tuple: (xx, yy, zz))
    xx, yy, zz = coordinates

    # Flatten arrays for plotting
    x_flat = xx.ravel()
    y_flat = yy.ravel()
    z_flat = zz.ravel()
    d_flat = ds_data.ravel()

    # Filter out points that are close to 0 (below threshold)
    mask = d_flat > threshold
    x_flat = x_flat[mask]
    y_flat = y_flat[mask]
    z_flat = z_flat[mask]
    d_flat = d_flat[mask]

    # Manual scaling: scale d_flat values to [0,1] based on the filtered data
    d_min = d_flat.min()
    d_max = d_flat.max()
    scaled = (d_flat - d_min) / (d_max - d_min)

    # Map the scaled values to colors using the colormap
    cmap = plt.cm.viridis
    colors = cmap(scaled)

    # Compute custom alpha values based on the original d_flat values.
    # Here, each alpha value is computed as (d / d_max)^alpha_power.
    alphas = np.power(d_flat / d_max, alpha_power)
    colors[:, -1] = alphas  # Replace the alpha channel

    print("Final number of points:", x_flat.shape[0])
    print("Data range used for colors:", d_min, d_max)

    # Create scatter plot using the computed RGBA colors
    sc = ax.scatter(
        x_flat,
        y_flat,
        z_flat,
        c=colors,
        s=1,  # Adjust point size as needed
        edgecolor="none",
    )
    # Identify the maximum density point and its coordinates
    max_index = np.argmax(d_flat)
    max_x, max_y, max_z = x_flat[max_index], y_flat[max_index], z_flat[max_index]
    print("Maximum density point at:", (max_x, max_y, max_z))

    # Draw dashed lines from the maximum point to the lower bound of each axis.
    # Assuming limits[i][0] represents the lower bound for the respective axis.
    ax.plot(
        [max_x, max_x],
        [max_y, max_y],
        [max_z, limits[2][0]],
        color="red",
        linestyle="--",
    )
    ax.plot(
        [max_x, max_x],
        [max_y, limits[1][0]],
        [max_z, max_z],
        color="red",
        linestyle="--",
    )
    ax.plot(
        [max_x, limits[0][0]],
        [max_y, max_y],
        [max_z, max_z],
        color="red",
        linestyle="--",
    )

    # Configure axes
    ax.set_xlim(limits[0])
    ax.set_ylim(limits[1])
    ax.set_zlim(limits[2])
    ax.set_xlabel(labels[0], fontsize=12, labelpad=15)
    ax.set_ylabel(labels[1], fontsize=12, labelpad=15)
    ax.set_zlabel(labels[2], fontsize=12, labelpad=15)

    # Customize axes properties

    plt.gca().spines["bottom"].set_linewidth(0.5)
    plt.gca().spines["left"].set_linewidth(0.5)
    plt.gca().spines["top"].set_linewidth(0.5)
    plt.gca().spines["right"].set_linewidth(0.5)

    # Customize tick parameters for major ticks

    xmax = limits[0][1]
    xmin = limits[0][0]

    ymax = limits[1][1]
    ymin = limits[1][0]

    zmax = limits[2][1]
    zmin = limits[2][0]

    num_ticks_x = (xmax - xmin) / 6
    num_ticks_y = (ymax - ymin) / 6
    num_ticks_z = (zmax - zmin) / 6

    major_locator_x = ticker.MultipleLocator(base=num_ticks_x)
    # Set major tick every 30 degrees
    major_locator_y = ticker.MultipleLocator(base=num_ticks_y)
    # Set major tick every 30 degrees

    major_locator_z = ticker.MultipleLocator(
        base=num_ticks_z
    )  # Set major tick every 30 degrees
    plt.gca().xaxis.set_major_locator(major_locator_x)
    plt.gca().yaxis.set_major_locator(major_locator_y)
    plt.gca().zaxis.set_major_locator(major_locator_z)
    plt.gca().xaxis.set_tick_params(
        width=0.5, length=3, direction="out", which="both", top=True
    )
    plt.gca().yaxis.set_tick_params(
        width=0.5, length=3, direction="out", which="both", right=True
    )

    # Customize tick parameters for minor ticks
    minor_locator_x = ticker.MultipleLocator(
        base=num_ticks_x / 4
    )  # Set minor tick every 10 degrees

    minor_locator_y = ticker.MultipleLocator(
        base=num_ticks_y / 4
    )  # Set minor tick every 10 degrees

    minor_locator_z = ticker.MultipleLocator(
        base=num_ticks_z / 4
    )  # Set minor tick every 10 degrees

    plt.gca().xaxis.set_minor_locator(minor_locator_x)
    plt.gca().yaxis.set_minor_locator(minor_locator_y)
    plt.gca().zaxis.set_minor_locator(minor_locator_z)

    plt.gca().xaxis.set_tick_params(
        width=0.3, length=1.5, direction="out", which="minor", top=True
    )
    plt.gca().yaxis.set_tick_params(
        width=0.3, length=1.5, direction="out", which="minor", right=True
    )

    plt.gca().zaxis.set_tick_params(
        width=0.3, length=1.5, direction="out", which="minor", right=True
    )

    # Add colorbar (Note: the colorbar now reflects the manual scaling)
    # cbar = fig.colorbar(sc, ax=ax, shrink=0.6)
    # cbar.set_label("Density", rotation=270, labelpad=20)

    # Set view angle
    ax.view_init(elev=25, azim=45)
    # plt.show()

    plt.savefig(name + ".png", dpi=1200)


def debug_geometry(normal1, normal2, point1, point2):
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D

    Axes3D.text()

    # Normalize the normal vectors
    normal1 /= np.linalg.norm(normal1)
    normal2 /= np.linalg.norm(normal2)

    # Create a grid of points in 3D space
    x, y = np.meshgrid(np.linspace(-5, 5, 10), np.linspace(-5, 5, 10))
    z1 = (-normal1[0] * x - normal1[1] * y + np.dot(normal1, point1)) * 1.0 / normal1[2]
    z2 = (-normal2[0] * x - normal2[1] * y + np.dot(normal2, point2)) * 1.0 / normal2[2]

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Plot the planes
    ax.plot_surface(x, y, z1, color="blue", alpha=0.5)
    ax.plot_surface(x, y, z2, color="green", alpha=0.5)

    # Plot the normal vectors as arrows with bigger size
    # Plot the normal vectors as lines
    ax.plot(
        [point1[0], point1[0] + normal1[0]],
        [point1[1], point1[1] + normal1[1]],
        [point1[2], point1[2] + normal1[2]],
        color="red",
    )

    ax.plot(
        [point2[0], point2[0] + normal2[0]],
        [point2[1], point2[1] + normal2[1]],
        [point2[2], point2[2] + normal2[2]],
        color="red",
    )

    arrow_scale = 0.5
    ax.quiver(
        point1[0] + normal1[0],
        point1[1] + normal1[1],
        point1[2] + normal1[2],
        normal1[0],
        normal1[1],
        normal1[2],
        color="red",
        arrow_length_ratio=arrow_scale,
    )
    ax.quiver(
        point2[0] + normal2[0],
        point2[1] + normal2[1],
        point2[2] + normal2[2],
        normal2[0],
        normal2[1],
        normal2[2],
        color="red",
        arrow_length_ratio=arrow_scale,
    )

    normal_text1 = f"""\
    Normal 1: ({normal1[0]:6.2f}, \
    {normal1[1]:6.2f}, \
    {normal1[2]:6.2f})\
    """
    normal_text2 = f"""\
    Normal 2: ({normal2[0]:6.2f}, \
    {normal2[1]:6.2f}, \
    {normal2[2]:6.2f})\
    """

    # Plot the points corresponding to the origins of the normal vectors
    ax.scatter(
        [point1[0], point2[0]],
        [point1[1], point2[1]],
        [point1[2], point2[2]],
        color="black",
        s=50,
    )

    # Add labels to display point and normal vector coordinates

    # ax.text(point1[0],
    #         point1[1],
    #         point1[2],
    #         f'Point 1: ({point1[0]:6.2f}, {point1[1]:6.2f}, {point1[2]:6.2f})',
    #         color='black')
    #
    # ax.text(point2[0],
    #         point2[1],
    #         point2[2],
    #         f'Point 2: ({point2[0]:6.2f}, {point2[1]:6.2f}, {point2[2]:6.2f})',
    #         color='black')

    ax.text(
        point1[0] + normal1[0],
        point1[1] + normal1[1],
        point1[2] + normal1[2],
        normal_text1,
        color="black",
    )
    ax.text(
        point2[0] + normal2[0],
        point2[1] + normal2[1],
        point2[2] + normal2[2],
        normal_text2,
        color="black",
    )

    # Set labels and title
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Plot of Two Planes with Normal Vectors")

    # Set the scale of the plot (adjust these values as needed)
    ax.set_xlim([-5, 5])
    ax.set_ylim([-5, 5])
    ax.set_zlim([-5, 5])

    # Show the plot
    # plt.show()
    plt.close()


def plt_ramachandran(data, coordinates, limits=None, name=None, labels=None):
    """
    Generate a Ramachandran plot using kernel density estimation (KDE).

    Parameters:
    data (numpy.ndarray): The data to generate the plot from.
    labels (list): Labels for the X and Y axes.
    name (str): The name of the output plot file (without extension).

    Returns:
    None

    Requires:
    - numpy
    - matplotlib.pyplot
    - matplotlib.ticker
    """
    # import necessary libraries and modules
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import matplotlib.font_manager as font_manager
    from utils import latin_2_greek

    # Define the grid for the contour plot

    if labels[0].startswith("Theta"):
        limits[0] = [0, 90]
    if labels[1].startswith("Theta"):
        limits[1] = [0, 90]
    xg, yg = limits
    xmin, xmax = xg
    ymin, ymax = yg

    # Configure font and font size for the plot
    font_path = "/home/mszatko/.local/share/fonts/ARIALUNI.TTF"
    custom_font = font_manager.FontProperties(fname=font_path)
    plt.rcParams["font.family"] = custom_font.get_name()
    plt.rcParams["font.size"] = 8

    # Create a new figure
    plt.figure(figsize=(3.45, 3.45))

    # Create a filled contour plot
    plt.contourf(
        coordinates[0],
        coordinates[1],
        data,
        cmap="Greys",
        levels=100,
        vmin=0,
    )

    # Customize labels and title
    labels = latin_2_greek(labels)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])

    # Set the limits of the x and y axes to (-180, 180)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    # Customize axes properties
    plt.gca().spines["bottom"].set_linewidth(0.5)
    plt.gca().spines["left"].set_linewidth(0.5)
    plt.gca().spines["top"].set_linewidth(0.5)
    plt.gca().spines["right"].set_linewidth(0.5)

    # Customize tick parameters for major ticks
    num_ticks_x = (xmax - xmin) / 6
    num_ticks_y = (ymax - ymin) / 6
    major_locator_x = ticker.MultipleLocator(
        base=num_ticks_x
    )  # Set major tick every 30 degrees
    major_locator_y = ticker.MultipleLocator(
        base=num_ticks_y
    )  # Set major tick every 30 degrees
    plt.gca().xaxis.set_major_locator(major_locator_x)
    plt.gca().yaxis.set_major_locator(major_locator_y)
    plt.gca().xaxis.set_tick_params(
        width=0.5, length=3, direction="out", which="both", top=True
    )
    plt.gca().yaxis.set_tick_params(
        width=0.5, length=3, direction="out", which="both", right=True
    )

    # Customize tick parameters for minor ticks
    minor_locator_x = ticker.MultipleLocator(
        base=num_ticks_x / 4
    )  # Set minor tick every 10 degrees
    minor_locator_y = ticker.MultipleLocator(
        base=num_ticks_y / 4
    )  # Set minor tick every 10 degrees
    plt.gca().xaxis.set_minor_locator(minor_locator_x)
    plt.gca().yaxis.set_minor_locator(minor_locator_y)
    plt.gca().xaxis.set_tick_params(
        width=0.3, length=1.5, direction="out", which="minor", top=True
    )
    plt.gca().yaxis.set_tick_params(
        width=0.3, length=1.5, direction="out", which="minor", right=True
    )

    # Ensure a tight layout
    plt.tight_layout()

    # Save the plot as an image
    plt.savefig(name + ".png", dpi=1200)
    # Close the plot
    plt.close()


def plt_distribution(
    data_dict: dict, coordinates_dict, limits=None, name=None, labels=None
):
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
    var1_names, var2_names, _ = data_keys
    var1_names = sorted(var1_names)
    var2_names = sorted([int(v) for v in var2_names])

    num_plots = len(var1_names)

    # Configure font and font size for the plot
    # font_path = "/home/keppen/.local/share/fonts/ARIALUNI.TTF"
    # custom_font = font_manager.FontProperties(fname=font_path)
    # plt.rcParams["font.family"] =
    plt.rcParams["font.size"] = 12

    # Create a figure and subplots
    # fig = plt.figure(figsize=(3.45, 3.45 / 2), layout="constrained")
    fig = plt.figure(figsize=(3.45, 3.45 / 2))
    GS = GridSpec(1, num_plots, figure=fig)

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
    # colormap = matplotlib.colormaps["Greens"]
    # colormap = matplotlib.colormaps["Blues"]
    colormap = matplotlib.colormaps["Reds"]
    bins = [(x + 1) / (len(var2_names) - 0) for x in range(len(var2_names))]
    print(bins)
    colors = [colormap(x) for x in bins]

    for i, var1_name in enumerate(var1_names):
        ax = fig.add_subplot(GS[0, i])
        # ax = axes[i]

        # Plot the data
        for j, var2_name in enumerate(var2_names):
            ax.plot(
                data_dict[f"{var1_name} {var2_name} dens"],
                coordinates_dict[f"{var1_name} {var2_name} coords"][0],
                color=colors[j],
                linewidth=0.75,
                label=f"{var2_name}",
            )

        # Set labels and title for each subplot
        ax.set_xlabel(labels[i])
        if i == 0:
            ax.set_ylabel("Angle distribution [\u00b0]")

        # Add a legend above the plot figure in the upper right corner
        if i == len(var2_names) - 1:
            plt.legend(
                loc="upper right",
                bbox_to_anchor=(1.2, 1.3),
                ncol=len(labels) + 1,
                frameon=False,
                columnspacing=1,
                # title="Mer: ",
            )

        # Extract limits from variable
        xg = limits[i]
        xmin, xmax = xg

        # Set the limits of the x and y axes to (-180, 180)
        plt.ylim(xmin, xmax)

        # Customize axes properties
        ax.spines["bottom"].set_linewidth(0.5)
        ax.spines["left"].set_linewidth(0.5)
        ax.spines["top"].set_linewidth(0.5)
        ax.spines["right"].set_linewidth(0.5)

        # Customize tick parameters for major ticks
        ax.set_xticks([])

        if i == 0:
            num_ticks = (xmax - xmin) / 4
            major_locator = ticker.MultipleLocator(
                base=num_ticks
            )  # Set major tick every 30 degrees
            ax.yaxis.set_major_locator(major_locator)
            ax.yaxis.set_tick_params(
                width=0.5, length=3, direction="out", which="both", right=False
            )

            # Customize tick parameters for minor ticks
            minor_locator = ticker.MultipleLocator(
                base=num_ticks / 4
            )  # Set minor tick every 10 degrees
            ax.yaxis.set_minor_locator(minor_locator)
            ax.yaxis.set_tick_params(
                width=0.3, length=1.5, direction="out", which="minor", right=False
            )
        else:
            ax.set_yticks([])

        # Set the position of the subplot
        ax.set_position(GS[i].get_position(fig))

    # Adjust layout to avoid overlapping labels
    # plt.tight_layout(rect=[0, 0, 1, 1])

    # Increase spacing between subplots and the surrounding area
    # plt.subplots_adjust(
    #     left=0.15, right=0.95, top=0.80, bottom=0.15, wspace=1, hspace=0.5
    # )
    plt.subplots_adjust(
        left=0.25, right=0.98, top=0.85, bottom=0.15, wspace=0, hspace=0
    )

    # Show the plot
    # plt.show()

    # Save the plot as an image
    plt.savefig(name + ".png", dpi=1200)
    # Close the plot
    plt.close()


def plt_heatmap(data, limits=None, name=None, labels=None, **kargs):
    # import necessary libraries and modules
    import matplotlib.pyplot as plt
    import numpy as np

    # Extract X and Y coordinates
    # x_coords, y_coords = coordinates

    # Configure font and font size for the plot
    plt.rcParams["font.sans-serif"] = "Arial"
    plt.rcParams["font.size"] = 8

    palette = None
    if "palette" in kargs:
        palette = kargs["palette"]

    # Create a new figure
    fig, ax = plt.subplots(figsize=(3.45, 3.45))

    # Create a heatmap
    ax.imshow(
        data,
        cmap=palette or "viridis",
        origin="lower",
        # extent=[x_coords.min(), x_coords.max(), y_coords.min(), y_coords.max()],
        aspect="auto",
    )

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(labels[0])), labels=labels[0])
    ax.set_yticks(np.arange(len(labels[1])), labels=labels[1])

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    if "ax_labels" in kargs:
        ax.set_xlabel(kargs["ax_labels"][0])
        ax.set_ylabel(kargs["ax_labels"][1])

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.to_numpy().shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.to_numpy().shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Set the limits of the x and y axes
    # plt.xlim(*limits[0])
    # plt.ylim(*limits[1])

    # Ensure a tight layout
    plt.tight_layout()

    # Save the plot as an image
    plt.savefig(name + ".png", dpi=1200)
    # Close the plot
    plt.close()


def plt_lineplot(data, limits=None, name=None, labels=None):
    # import necessary libraries and modules
    import matplotlib.pyplot as plt

    # Extract X and Y coordinates
    # x_coords, y_coords = coordinates

    # Configure font and font size for the plot
    plt.rcParams["font.sans-serif"] = "Arial"
    plt.rcParams["font.size"] = 8

    # Create a new figure
    plt.figure(figsize=(3.45, 3.45))

    # Create a line plot
    plt.plot(data.index, data.distance)

    # Customize labels and title
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])

    # Ensure a tight layout
    plt.tight_layout()

    # Save the plot as an image
    plt.savefig(name + ".png", dpi=1200)
    # Close the plot
    plt.close()
