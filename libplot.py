def mavi_contour(data):
    from libmath import calculate_kde
    import numpy as np
    from mayavi import mlab
    import matplotlib.pyplot as plt

    # mlab.options.offscreen = True
    # fig = mlab.figure(figure=None)
    # mlab.options.backend = 'envisage'
    grid = np.array([
        [-180, 180],
        [-180, 180],
        [-180, 180]
            ])

    xg, yg, zg = grid
    density, coords = calculate_kde(data, grid, 10)
    xi, yi, zi = coords

    # Create a white scene background
    fig = mlab.figure(bgcolor=(1, 1, 1))

    # print(coords, density)
    mlab.contour3d(xi, yi, zi, density, opacity=0.5)

    color_axes = (0, 0, 0)
    mlab.axes(
        fig,  # The Mayavi figure
        ranges=(xg[0], xg[1], yg[0], yg[1], zg[0], zg[1]),  # Specify the ranges for each axis
        extent=(xg[0], xg[1], yg[0], yg[1], zg[0], zg[1]),
        # xlabel='X Label',  # Label for the X axis
        # ylabel='Y Label',  # Label for the Y axis
        # zlabel='Z Label',  # Label for the Z axis
        color=color_axes,    # Axes color
        line_width=1.0,     # Line width of the axes
        # nb_labels=5,        # Number of labels on each axis
    )

    # # Customize Matplotlib style and add labels
    # title = "Mayavi Scene with Matplotlib Labels"
    # xlabel, ylabel, zlabel = "X Label", "Y Label", "Z Label"
    # plt.title(title, color='black', fontsize=16)  # Set title font color and size
    #
    # # Create a subplots arrangement with a single subplot (Mayavi scene)
    # fig, ax = plt.subplots()
    # ax.set_xlabel(xlabel, color='blue', fontsize=14)  # Set xlabel color and size
    # ax.set_ylabel(ylabel, color='green', fontsize=14)  # Set ylabel color and size
    # ax.set_zlabel(zlabel, color='red', fontsize=14)  # Set zlabel color and size


    # Manually create a bounding box
    xmin, xmax = xg
    ymin, ymax = yg
    zmin, zmax = zg

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

    # fig.scene.aspect_ratio = (1, 1, 1)

    mlab.axes()
    mlab.show()
    # mlab.savefig("test.png")
    # mlab.close()


def debug_geometry(normal1, normal2, point1, point2):
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    # Normalize the normal vectors
    normal1 /= np.linalg.norm(normal1)
    normal2 /= np.linalg.norm(normal2)

    # Create a grid of points in 3D space
    x, y = np.meshgrid(np.linspace(-5, 5, 10), np.linspace(-5, 5, 10))
    z1 = (-normal1[0] * x - normal1[1] * y + np.dot(normal1, point1)) * 1. / normal1[2]
    z2 = (-normal2[0] * x - normal2[1] * y + np.dot(normal2, point2)) * 1. / normal2[2]

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
    ax.scatter([point1[0], point2[0]], [point1[1], point2[1]], [point1[2], point2[2]], color='black', s=50)

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
