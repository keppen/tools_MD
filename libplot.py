def mavi_contour(data, coordinates, limits=None, name="test", labels=None):
    from mayavi import mlab
    # import matplotlib.pyplot as plt

    # Configure Matplotlib font settings
    # plt.rcParams['font.sans-serif'] = 'Arial'
    # plt.rcParams['font.size'] = 10

    # mlab.options.offscreen = True
    # fig = mlab.figure(figure=None)
    # mlab.options.backend = 'envisage'

    xg, yg, zg = limits
    xi, yi, zi = coordinates

    # Create a white scene background
    fig = mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(963, 963))

    # Manually create a bounding box
    xmin, xmax = xg
    ymin, ymax = yg
    zmin, zmax = zg

    color_axes = (0, 0, 0)

    mlab.plot3d(
        [xmin, xmax, xmax, xmin, xmin],
        [ymin, ymin, ymax, ymax, ymin],
        [zmin, zmin, zmin, zmin, zmin],
        color=color_axes,
        tube_radius=None,  # Solid lines
        line_width=2.0
    )

    mlab.plot3d(
        [xmin, xmax, xmax, xmin, xmin],
        [ymin, ymin, ymax, ymax, ymin],
        [zmax, zmax, zmax, zmax, zmax],
        color=color_axes,
        tube_radius=None,  # Solid lines
        line_width=2.0
    )

    mlab.plot3d(
        [xmin, xmax, xmax, xmin, xmin],
        [ymax, ymax, ymax, ymax, ymax],
        [zmin, zmin, zmax, zmax, zmin],
        color=color_axes,
        tube_radius=None,  # Solid lines
        line_width=2.0
    )

    mlab.plot3d(
        [xmin, xmax, xmax, xmin, xmin],
        [ymin, ymin, ymin, ymin, ymin],
        [zmin, zmin, zmax, zmax, zmin],
        color=color_axes,
        tube_radius=None,  # Solid lines
        line_width=2.0
    )

    contour = mlab.contour3d(xi, yi, zi, data, opacity=0.5)

    # Create a Mayavi text actor for the axes labels
    ax = mlab.axes(
        contour,  # The Mayavi figure
        # Specify the ranges for each axis
        ranges=(xg[0], xg[1], yg[0], yg[1], zg[0], zg[1]),
        extent=(xg[0], xg[1], yg[0], yg[1], zg[0], zg[1]),
        # xlabel='\u03c6',  # Label for the X axis
        # ylabel='\u03bf',  # Label for the Y axis
        # zlabel='\u03ba',  # Label for the Z axis
        color=color_axes,    # Axes color
        line_width=1.0,     # Line width of the axes
        nb_labels=5,        # Number of labels on each axis
    )

    # Customize the font properties of the scale (ticks and labels)
    # for tick in ax.xaxis.get_major_ticks():
    #     tick.label1.set_fontname('arial')  # Change the font name
    #     tick.label1.set_fontsize(10)       # Change the font size
    #
    # for tick in ax.yaxis.get_major_ticks():
    #     tick.label1.set_fontname('arial')  # Change the font name
    #     tick.label1.set_fontsize(10)       # Change the font size
    #
    # for tick in ax.zaxis.get_major_ticks():
    #     tick.label1.set_fontname('arial')  # Change the font name
    #     tick.label1.set_fontsize(10)       # Change the font size

    ax.axes.font_factor = 1
    ax.axes.label_format = '    %4.0f'
    ax.axes.x_label = labels[0]
    ax.axes.y_label = labels[1]
    ax.axes.z_label = labels[2]

    ax.label_text_property.bold = False
    ax.label_text_property.italic = False
    ax.label_text_property.color = (0, 0, 0)
    ax.label_text_property.font_size = 8
    ax.label_text_property.font_family = 'arial'

    ax.title_text_property.bold = False
    ax.title_text_property.italic = False
    ax.title_text_property.font_size = 8
    ax.title_text_property.font_family = 'arial'
    ax.property.color = (0, 0, 0)

    # Adjust the layout settings for a tight layout
    fig.scene.disable_render = True  # Disable rendering during layout adjustments

    mlab.view(30, 30,                # Adjust camera distance
              distance=(xmax - xmin)*3.6
              )
    mlab.roll(0)                     # Reset camera roll

    fig.scene.disable_render = False  # Re-enable renderingv

    # mlab.savefig("test.png",
    # size=(600, 600)
    # )
    # mlab.close()
    mlab.show()


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
    z1 = (-normal1[0] * x - normal1[1] * y +
          np.dot(normal1, point1)) * 1. / normal1[2]
    z2 = (-normal2[0] * x - normal2[1] * y +
          np.dot(normal2, point2)) * 1. / normal2[2]

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the planes
    ax.plot_surface(x, y, z1, color='blue', alpha=0.5)
    ax.plot_surface(x, y, z2, color='green', alpha=0.5)

    # Plot the normal vectors as arrows with bigger size
    # Plot the normal vectors as lines
    ax.plot([point1[0], point1[0] + normal1[0]],
            [point1[1], point1[1] + normal1[1]],
            [point1[2], point1[2] + normal1[2]], color='red')

    ax.plot([point2[0], point2[0] + normal2[0]],
            [point2[1], point2[1] + normal2[1]],
            [point2[2], point2[2] + normal2[2]], color='red')

    arrow_scale = 0.5
    ax.quiver(
        point1[0] + normal1[0],
        point1[1] + normal1[1],
        point1[2] + normal1[2],
        normal1[0],
        normal1[1],
        normal1[2],
        color='red',
        arrow_length_ratio=arrow_scale
    )
    ax.quiver(
        point2[0] + normal2[0],
        point2[1] + normal2[1],
        point2[2] + normal2[2],
        normal2[0],
        normal2[1],
        normal2[2],
        color='red',
        arrow_length_ratio=arrow_scale
    )

    normal_text1 = f'''\
    Normal 1: ({normal1[0]:6.2f}, \
    {normal1[1]:6.2f}, \
    {normal1[2]:6.2f})\
    '''
    normal_text2 = f'''\
    Normal 2: ({normal2[0]:6.2f}, \
    {normal2[1]:6.2f}, \
    {normal2[2]:6.2f})\
    '''

    # Plot the points corresponding to the origins of the normal vectors
    ax.scatter([point1[0], point2[0]], [point1[1], point2[1]],
               [point1[2], point2[2]], color='black', s=50)

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

    ax.text(point1[0] + normal1[0],
            point1[1] + normal1[1],
            point1[2] + normal1[2],
            normal_text1,
            color='black')
    ax.text(point2[0] + normal2[0],
            point2[1] + normal2[1],
            point2[2] + normal2[2],
            normal_text2,
            color='black')

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Plot of Two Planes with Normal Vectors')

    # Set the scale of the plot (adjust these values as needed)
    ax.set_xlim([-5, 5])
    ax.set_ylim([-5, 5])
    ax.set_zlim([-5, 5])

    # Show the plot
    plt.show()
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

    xg, yg = limits
    xmin, xmax = xg
    ymin, ymax = yg

    # Configure font and font size for the plot
    font_path = "/home/mszatko/.local/share/fonts/ARIALUNI.TTF"
    custom_font = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = custom_font.get_name()
    plt.rcParams['font.size'] = 8

    # Create a new figure
    plt.figure(figsize=(3.45, 3.45))

    # Create a filled contour plot
    plt.contourf(
        coordinates[0],
        coordinates[1],
        data,
        cmap='Greys',
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
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)

    # Customize tick parameters for major ticks
    num_ticks = (xmax - xmin) / 6
    major_locator = ticker.MultipleLocator(
        base=num_ticks)  # Set major tick every 30 degrees
    plt.gca().xaxis.set_major_locator(major_locator)
    plt.gca().yaxis.set_major_locator(major_locator)
    plt.gca().xaxis.set_tick_params(width=0.5, length=3,
                                    direction='out', which='both', top=True)
    plt.gca().yaxis.set_tick_params(width=0.5, length=3,
                                    direction='out', which='both', right=True)

    # Customize tick parameters for minor ticks
    minor_locator = ticker.MultipleLocator(
        base=num_ticks/4)  # Set minor tick every 10 degrees
    plt.gca().xaxis.set_minor_locator(minor_locator)
    plt.gca().yaxis.set_minor_locator(minor_locator)
    plt.gca().xaxis.set_tick_params(width=0.3, length=1.5,
                                    direction='out', which='minor', top=True)
    plt.gca().yaxis.set_tick_params(width=0.3, length=1.5,
                                    direction='out', which='minor', right=True)

    # Ensure a tight layout
    plt.tight_layout()

    # Save the plot as an image
    plt.savefig(
        name+".png",
        dpi=1200
    )
    # Close the plot
    plt.close()


def plt_distribution(data_list, coordinates_list,  limits=None, name=None, labels=None):
    # import necessary libraries and modules
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import matplotlib.font_manager as font_manager
    from utils import latin_2_greek

    num_plots = len(data_list)

    # Configure font and font size for the plot
    font_path = "/home/mszatko/.local/share/fonts/ARIALUNI.TTF"
    custom_font = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = custom_font.get_name()
    plt.rcParams['font.size'] = 8

    # Create a figure and subplots
    fig, axes = plt.subplots(1, num_plots, figsize=(3.45, 3.45), sharey=True)

    # Extract limits from variable
    xg = limits[0]
    xmin, xmax = xg

    # Change lables to greek
    labels = latin_2_greek(labels)

    for i, packed in enumerate(zip(data_list, coordinates_list, labels)):

        # unpack variables
        data, coords, label = packed
        ax = axes[i]

        # Plot the data
        ax.plot(data, coords[0], color="b", linewidth=0.75)

        # Set labels and title for each subplot
        ax.set_xlabel(label)
        if i == 0:
            ax.set_ylabel("Angle distibution")

        # Set the limits of the x and y axes to (-180, 180)
        plt.ylim(xmin, xmax)

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

    # Adjust layout to avoid overlapping labels
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Show the plot
    # plt.show()

    # Save the plot as an image
    plt.savefig(
        name + ".png",
        dpi=1200
    )
    # Close the plot
    plt.close()


def plt_heatmap(data, limits=None, name=None, labels=None, **kargs):
    # import necessary libraries and modules
    import matplotlib.pyplot as plt
    import numpy as np

    # Extract X and Y coordinates
    # x_coords, y_coords = coordinates

    # Configure font and font size for the plot
    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.size'] = 8

    palette = None
    if "palette" in kargs:
        palette = kargs["palette"]

    # Create a new figure
    fig, ax = plt.subplots(figsize=(3.45, 3.45))

    # Create a heatmap
    ax.imshow(
        data,
        cmap=palette or 'viridis',
        origin='lower',
        # extent=[x_coords.min(), x_coords.max(), y_coords.min(), y_coords.max()],
        aspect='auto'
    )

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(labels[0])), labels=labels[0])
    ax.set_yticks(np.arange(len(labels[1])), labels=labels[1])

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    if "ax_labels" in kargs:
        ax.set_xlabel(kargs["ax_labels"][0])
        ax.set_ylabel(kargs["ax_labels"][1])

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.to_numpy().shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.to_numpy().shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Set the limits of the x and y axes
    # plt.xlim(*limits[0])
    # plt.ylim(*limits[1])

    # Ensure a tight layout
    plt.tight_layout()

    # Save the plot as an image
    plt.savefig(
        name + ".png",
        dpi=1200
    )
    # Close the plot
    plt.close()


def plt_lineplot(data, limits=None, name=None, labels=None):
    # import necessary libraries and modules
    import matplotlib.pyplot as plt

    # Extract X and Y coordinates
    # x_coords, y_coords = coordinates

    # Configure font and font size for the plot
    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.size'] = 8

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
    plt.savefig(
        name + ".png",
        dpi=1200
    )
    # Close the plot
    plt.close()
