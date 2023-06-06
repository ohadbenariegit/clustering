import numpy as np
from spkmeansmodule import fit_goal, fit_kmeans
import sys

np.random.seed(0)


def validate_input(condition):
    if not condition:
        print('Invalid Input!', end='')
        exit(1)


def convert_arg(arg, convert):
    try:
        return convert(arg)
    except:
        validate_input(False)


def get_initial_centroid_indexes():
    n, m = datapoints.shape

    indexes = [np.random.choice(n)]

    for i in range(1, m):
        # for each datapoint, calculate the minimum squared euclidean distance of datapoint from existing centroid
        distances = np.apply_along_axis(
            lambda datapoint: (np.linalg.norm(datapoint - datapoints[indexes], axis=1) ** 2).min(),
            1,
            datapoints
        )

        probabilities = distances / (distances.sum() or 1)

        indexes.append(np.random.choice(n, p=probabilities))

    return indexes


def print_matrix(mat, jacobi):
    # jacobi eigenvalues are supposed to be positive
    for row in range(len(mat)):
        print(','.join(format(0 if jacobi and row == 0 and abs(val) < 0.0001 else val, '.4f') for val in mat[row]))


try:
    args = sys.argv[1:]

    validate_input(len(args) == 3)

    k, goal, filename = args

    try:
        k = int(k)
    except:
        validate_input(False)

    output = fit_goal(k, goal, filename)

    if goal == 'spk':
        datapoints = np.array(output)

        initial_centroid_indexes = get_initial_centroid_indexes()

        print(','.join(map(str, initial_centroid_indexes)))

        initial_centroids = datapoints[initial_centroid_indexes]

        output = fit_kmeans(*datapoints.shape, datapoints.tolist(), initial_centroids.tolist())

    print_matrix(output, goal == 'jacobi')
except Exception:
    print('An Error Has Occurred', end='')
    exit(1)
