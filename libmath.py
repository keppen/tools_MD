import numpy as np
import itertools


def calcVector(p0, p1) -> np.array:
    p0 = np.array(p0)
    p1 = np.array(p1)
    return p1 - p0


def calcAngle(v0, v1) -> np.array:
    len0 = np.linalg.norm(v0)
    len1 = np.linalg.norm(v1)
    return np.dot(v0, v1) / len0 / len1


def calcPlane(p0, p1, p2) -> np.array:
    v0 = calcVector(p0, p1)
    v1 = calcVector(p1, p2)
    return np.cross(v0, v1)


def calcMiddlePoint(p0, p1):
    v = calcVector(p0, p1) / 2
    return p0 + v


def calcCentroid(points):
    return sum(points) / len(points)


# def formulaPowerIter(direction, points, geomteric_centre):
#     nextd = np.zeros(3)
#     for p in points:
#         centeredpoint = p - geomteric_centre
#         nextd += np.dot(centeredpoint, direction) * centeredpoint
#     return nextd/np.linalg.norm(nextd)


def calcDistance_form_Vector(point, vector):
    cross = np.cross(vector, point)
    return np.linalg.norm(cross) / np.linalg.norm(vector)


def calcProjectPoint(point, origin, direction):
    """Produce an orthagonal vector |x - z|, where z = t * u
    A dot product is a projection of v1 onto v2!

    x, point - a point to project
    x0, origin - origination of point and direction
    u ,direction - a vector on which a point projected
    z, projection - projected point
    t, scalar - np.dot((x - x0), u)
    """

    direction = direction / np.linalg.norm(direction)

    vector = calcVector(origin, point)
    scalar = np.dot(vector, direction)
    return origin + scalar * direction


def calcLinearRegression_PowerIteration(points):
    """Linear regresion using power iteration.
    Points need to have its geometrical centre at origin.

    Power interation approximates main eigenvector of a martix.
    It does it in recurrent series.
    b_1 = Ab_0 / ||Ab_0||, but here Ab_0 = sum( dot(p, b_0) * p)

    p, point = a point to be regressed
    A = a matrix, its eigenvector is approximated by power iteration
    b_x = approximated eigenvector of A
    """
    eigenvector = np.ones(3)
    err = 1
    while err > 1e-6:
        previous_eigenvector = eigenvector
        matrix = np.zeros(3)
        for p in points:
            matrix += np.dot(p, eigenvector) * p
        eigenvector = matrix / np.linalg.norm(matrix)
        err = np.linalg.norm(eigenvector - previous_eigenvector)
    return eigenvector


def calcDihedral(p0, p1, p2, p3):
    """Calculated a dihedral provided with four points.

    Parameters:
        px (list): List of 3 element list [x, y, z]
                    converted into numpy.array object.

    Returns:
        angle (float)
    """
    p0 = np.array(p0)
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    b0 = p1 - p0
    b1 = p2 - p1
    b2 = p3 - p2
    b1xb0 = np.cross(b1, b0)
    b2xb1 = np.cross(b2, b1)
    b1 /= np.linalg.norm(b1)
    b1xb0 /= np.linalg.norm(b1xb0)
    b2xb1 /= np.linalg.norm(b2xb1)
    b1_b0_x_b2xb1 = np.cross(b1xb0, b2xb1)
    cos = np.dot(b1xb0, b2xb1)
    sin = np.dot(b1_b0_x_b2xb1, b1)
    return np.degrees(np.arctan2(sin, cos))


def deriv_1st(p0, p1, p2, p3, h):
    return (p0 - 8 * p1 + 8 * p2 - p3) / 12 / h


def deriv_2nd(p0, p1, p2, h):
    return (p0 - 2 * p1 + p2) / h / h


class NumericAnalysis:
    def __init__(
        self,
        data: np.array,
        resolution: float,
        point: np.array,
        dimension: int,
        debug=False,
    ):
        self.debug = debug
        self.data = data
        self.resolution = resolution
        self.start_point = point
        self.start_point.shape = (len(point), 1)
        self.point = self.start_point
        self.dimension = dimension
        if self.debug:
            print(
                f"""Starting steepestdescent:
point: {np.transpose(self.start_point)}^T
energy: {self.lin_inter(point)}
dimension: {self.dimension}
resolution: {self.resolution}
"""
            )

    def deriv_1D(self, point, alpha, resolution):
        return deriv_1st(
            self.projection_1D(self.point, alpha - 2 * resolution),
            self.projection_1D(self.point, alpha - 1 * resolution),
            self.projection_1D(self.point, alpha + 1 * resolution),
            self.projection_1D(self.point, alpha + 2 * resolution),
            resolution,
        )

    def gradient(self):
        grad = np.zeros((self.dimension, 1))
        for i in range(self.dimension):
            vector = np.zeros((self.dimension, 1))
            vector[i] = self.resolution
            grad[i] = deriv_1st(
                self.lin_inter(self.point - 2 * vector),
                self.lin_inter(self.point - 1 * vector),
                self.lin_inter(self.point + 1 * vector),
                self.lin_inter(self.point + 2 * vector),
                self.resolution,
            )
        return grad

    def projection_1D(self, point, alpha):
        return self.lin_inter(point + alpha * self.direction)

    def bisection(self, a, b):
        step = 0
        tolerance = 1e-6
        a_deriv = self.deriv_1D(self.point, a, self.step_size)
        trial = 0
        if a_deriv > 0:
            hi, low = a, b
        else:
            hi, low = b, a
        while step < 1000:
            now_trial = (hi + low) / 2
            trial_deriv = self.deriv_1D(self.point, now_trial, self.step_size)
            if trial_deriv > 0:
                hi = now_trial
            else:
                low = now_trial
            if abs(trial - now_trial) < tolerance:
                return now_trial
            trial = now_trial
        return False

    def descent_1D(self):
        alpha = 0
        step = 0
        self.step_size = self.resolution / 10
        point = self.point
        point_value = self.lin_inter(point)
        while True:
            step += 1
            alpha += self.step_size
            new_point_value = self.projection_1D(point, alpha)
            if new_point_value > point_value:
                a = alpha
                b = alpha - self.step_size
                return self.bisection(a, b)
            point_value = new_point_value

    def steepest_descent(self):
        step = 0
        tolerance = 1e-5
        self.direction = -1 * self.gradient() / np.linalg.norm(self.gradient())

        while step < 500:
            step_size = self.descent_1D()

            if self.debug:
                print(
                    f"""Step {step}:
point: {np.transpose(self.point)}^T
energy: {self.lin_inter(self.point)}
step size: {step_size}
"""
                )

            point_now = self.point + self.direction * step_size

            point_value = self.lin_inter(self.point)
            point_now_value = self.lin_inter(point_now)
            if abs(point_now_value - point_value) < tolerance:
                break

            self.point = point_now
            self.direction = -1 * self.gradient() / np.linalg.norm(self.gradient())
            step += 1
        return self.point

    def calculate_point(self, point, params):
        combination = [(1)]
        for i in range(len(point)):
            combination.extend(list(itertools.combinations(point, i + 1)))
        combination = [np.prod(i) for i in combination]
        return sum(np.multiply(np.transpose([combination]), params))

    def lin_inter(self, point):
        point_before = np.trunc(point)
        A_matrix = np.zeros((2**self.dimension, 2**self.dimension))
        b_matrix = np.zeros((2**self.dimension, 1))
        for index, bool in enumerate(
            itertools.product([False, True], repeat=self.dimension)
        ):
            combinations = [[1]]
            next_to = list(
                map(
                    lambda values, condition: values + 1 if condition else values,
                    point_before,
                    bool,
                )
            )
            for j in range(self.dimension):
                combinations.extend(
                    list(itertools.combinations(next_to, j + 1)))
            combinations = [np.prod(set) for set in combinations]
            A_matrix[index] = combinations
            b_matrix[index] = self.data[
                int(next_to[0]), int(next_to[1]), int(next_to[2])
            ]
        params = np.linalg.solve(A_matrix, b_matrix)
        return self.calculate_point(point, params)
