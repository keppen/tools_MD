import numpy as np


def compare_kde(kde1, kde2, x_values):
    diff_kde = np.abs(kde2 - kde1)
    return abs(integral(x_values, diff_kde) / 2 * 100 - 100)


def integral(x_values, y_values):
    x_values = np.array(x_values)
    y_values = np.array(y_values)

    dx = np.diff(x_values)

    integral_value = np.sum((y_values[:-1] + y_values[1:]) * dx / 2)

    return integral_value


def helix_find(coords_subject, density_subject, density_proper, density_loose, *args):
    keys_density = list(density_subject.keys())
    keys_coords = list(coords_subject.keys())

    similarities = []
    msg = [args[-5].strip()]
    for index, _ in enumerate(keys_coords):
        x_values = coords_subject[keys_coords[index]]
        kde1 = density_proper[keys_density[index]]
        kde2 = density_subject[keys_density[index]]
        normalized_kde1 = kde1 / integral(x_values, kde1)
        normalized_kde2 = kde2 / integral(x_values, kde2)

        similarity = compare_kde(normalized_kde1, normalized_kde2, x_values)
        similarities.append(similarity)
    msg.append(str(np.sum(similarities) / 9))

    similarities = []
    for index, _ in enumerate(keys_coords):
        x_values = coords_subject[keys_coords[index]]
        kde1 = density_loose[keys_density[index]]
        kde2 = density_subject[keys_density[index]]
        normalized_kde1 = kde1 / integral(x_values, kde1)
        normalized_kde2 = kde2 / integral(x_values, kde2)

        similarity = compare_kde(normalized_kde1, normalized_kde2, x_values)
        similarities.append(similarity)
    msg.extend([str((np.sum(similarities) / 9)), "\n"])

    write_to_file("helix_similarities.txt", *msg)


def write_to_file(file_name, *args):
    with open(file_name, "a") as file:
        file.writelines(" ".join(args))
