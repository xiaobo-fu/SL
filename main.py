import numpy as np
import math
from matplotlib import pyplot as plt
from scipy.spatial.transform import Rotation


# points
class Point:
    def __init__(self, x, y, z=0):
        self.x = x
        self.y = y
        self.z = z
        self.cord = (x, y)
        self.cord3d = (x, y, z)
        self.width = z


# strokes
class Strokes:
    def __init__(self, points, mus, sigmas, deltas, alphas, t0s):
        self.points = points
        self.mus = mus
        self.sigmas = sigmas
        self.deltas = deltas
        self.alphas = alphas
        self.t0s = t0s

    def get_trejectory(self):
        trejectory = []
        for t in range(300):
            point_cord = np.array([0.0, 0.0, 0.0])
            for j in range(len(self.points) - 1):
                cord = get_arc_cord(self.points[j], self.points[j+1], t/100, mu=self.mus[j], sigma=self.sigmas[j], delta=self.deltas[j], alpha=self.alphas[j], t0=self.t0s[j])
                cord = np.array(cord)
                point_cord += cord
            point_cord += np.array(self.points[j].cord3d)
            trejectory.append(point_cord)
        return trejectory


if __name__ == "__main__":
    def get_arc_cord(p_oirgin, p_end, t, mu=0, sigma=0, delta=0.001, alpha=0, t0=0):
        norm = get_norm(p_oirgin, p_end)
        theta = get_theta(p_oirgin, p_end)
        D = norm * delta / 2 / np.sin(delta/2)
        if 0 >= t0:
            w = 1 / 2 * (1 + math.erf((np.log(t - t0) - mu) / (sigma * np.sqrt(2))))
        else:
            w = 0
        # print(w)
        phi = theta - delta / 2 + delta * w
        # print(phi)
        dx, dy = get_displacement_cord(delta, D, theta, w, phi)
        dz = w * (p_end.width - p_oirgin.width)

        #rotation
        axis = (np.array(p_end.cord3d) - np.array(p_oirgin.cord3d))/norm
        v = [dx, dy, dz]
        new_v = np.dot(rotation_matrix(axis, alpha), v)

        return new_v[0], new_v[1], new_v[2]

    def get_norm(p_oirgin, p_end):
        return np.sqrt((p_oirgin.x - p_end.x) ** 2 + (p_oirgin.y - p_end.y) ** 2 + (p_oirgin.z - p_end.z) ** 2)

    def get_theta(p_oirgin, p_end):
        dx = p_end.x - p_oirgin.x
        dy = p_end.y - p_oirgin.y
        theta = np.arctan2(dy, dx)
        return theta


    def get_displacement_cord(delta, D, theta, w, phi):
        if delta != 0:
            displacement_x = D / delta * (np.sin(phi) - np.sin(theta - delta / 2))
            displacement_y = D / delta * -(np.cos(phi) - np.cos(theta - delta / 2))
        else:
            displacement_x = w * np.cos(theta)
            displacement_y = w * np.sin(theta)
        return displacement_x, displacement_y


    def rotation_matrix(axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        axis = axis / math.sqrt(np.dot(axis, axis))
        a = math.cos(theta / 2.0)
        b, c, d = -axis * math.sin(theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                         [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                         [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


    p1 = Point(500, 500, 50)
    p2 = Point(100, 100, 10)
    p3 = Point(0, 200, 10)

    points = [p1, p2, p3]
    mus = [0.1, 0.1]
    sigmas = [0.4, 0.4]
    deltas = [1, 0.1]
    alphas = [-0.2, 0.2]
    t0s = [0, 0]

    stroke = Strokes(points, mus, sigmas, deltas, alphas, t0s)
    trejectory = stroke.get_trejectory()

    plt.rcParams["figure.figsize"] = [8.00, 8.50]
    plt.rcParams["figure.autolayout"] = True
    x = [5]
    y = [5]
    plt.xlim(-200, 600)
    plt.ylim(-200, 600)
    plt.grid()
    # plt.plot(p1.x, p1.y, marker="o", markersize=p1.width, markeredgecolor="black", markerfacecolor="black")
    # plt.plot(p2.x, p2.y, marker="o", markersize=p2.width, markeredgecolor="black", markerfacecolor="black")
    for i in trejectory:
        plt.plot(i[0], i[1], marker="o", markersize=20, markeredgecolor="black", markerfacecolor="black")
    plt.show()
