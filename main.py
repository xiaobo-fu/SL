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
        for t in range(800):
            point_cord = np.array([0.0, 0.0, 0.0])
            for j in range(len(self.points) - 1):
                cord = get_arc_cord(self.points[j], self.points[j+1], t/100, mu=self.mus[j], sigma=self.sigmas[j], delta=self.deltas[j], alpha=self.alphas[j], t0=self.t0s[j])
                cord = np.array(cord)
                point_cord += cord
            point_cord += np.array(self.points[0].cord3d)
            trejectory.append(point_cord)
        return trejectory


if __name__ == "__main__":
    def get_arc_cord(p_oirgin, p_end, t, mu=0.1, sigma=0.4, delta=0.001, alpha=0, t0=0):
        norm = get_norm(p_oirgin, p_end)
        theta = get_theta(p_oirgin, p_end)
        D = norm * delta / 2 / np.sin(delta/2)
        if t >= t0:
            w = 1 / 2 * (1 + math.erf((np.log(t - t0) - mu) / (sigma * np.sqrt(2))))
        else:
            return 0, 0, 0
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


    p11 = Point(109, 500, 1)
    p12 = Point(170, 479, 70)
    p13 = Point(160, 480, 1)

    points1 = [p11, p12, p13]
    mus1 = [0.1, 0.1]
    sigmas1 = [0.4, 0.4]
    deltas1 = [-0.2, -0.6]
    alphas1 = [-0.2, -0.2]
    t0s1 = [0, 0.6]

    stroke1 = Strokes(points1, mus1, sigmas1, deltas1, alphas1, t0s1)
    trejectory1 = stroke1.get_trejectory()

    p21 = Point(72, 367, 1)
    p22 = Point(102, 326, 30)
    p23 = Point(176, 112, 50)
    p24 = Point(214, 340, 1)

    points2 = [p21, p22, p23, p24]
    mus2 = [0.1, 0.1, 0.1]
    sigmas2 = [0.4, 0.4, 0.4]
    deltas2 = [-1, 0.6, -0.6]
    alphas2 = [-0.2, -0.2, -0.2]
    t0s2 = [0, 0.4, 3]

    stroke2 = Strokes(points2, mus2, sigmas2, deltas2, alphas2, t0s2)
    trejectory2 = stroke2.get_trejectory()

    p31 = Point(235, 450, 1)
    p32 = Point(269, 400, 50)
    p33 = Point(283, 300, 30)
    p34 = Point(520, 450, 30)
    p35 = Point(376, 390, 0)

    points3 = [p31, p32, p33, p34, p35]
    mus3 = [0.1, 0.1, 0.1, 0.1]
    sigmas3 = [0.4, 0.4, 0.4, 0.4]
    deltas3 = [-2.5, 0.6, -2.5, -2.5]
    alphas3 = [-0.2, -0.2, -0.2, 0.2]
    t0s3 = [0, 0.4, 3, 4.7]

    stroke3 = Strokes(points3, mus3, sigmas3, deltas3, alphas3, t0s3)
    trejectory3 = stroke3.get_trejectory()

    p41 = Point(383, 500, 0)
    p42 = Point(409, 480, 100)
    p43 = Point(390, 380, 10)
    p44 = Point(500, 280, 50)
    p45 = Point(230, 200, 0)

    points4 = [p41, p42, p43, p44, p45]
    mus4 = [0.1, 0.1, 0.1, 0.1]
    sigmas4 = [0.4, 0.4, 0.4, 0.4]
    deltas4 = [-2, 0.6, -1, -1.8]
    alphas4 = [-0.2, -0.2, -0.2, 0.2]
    t0s4 = [0, 0.4, 2, 2.5]

    stroke4 = Strokes(points4, mus4, sigmas4, deltas4, alphas4, t0s4)
    trejectory4 = stroke4.get_trejectory()

    p51 = Point(320, 180, 1)
    p52 = Point(380, 200, 2)
    p53 = Point(540, 60, 4)

    points5 = [p51, p52, p53]
    mus5 = [0.1, 0.1]
    sigmas5 = [0.4, 0.4]
    deltas5 = [-1, -1]
    alphas5 = [1, 1]
    t0s5 = [0, 0.1]

    stroke5 = Strokes(points5, mus5, sigmas5, deltas5, alphas5, t0s5)
    trejectory5 = stroke5.get_trejectory()




    plt.rcParams["figure.figsize"] = [8.00, 8.50]
    plt.rcParams["figure.autolayout"] = True
    x = [5]
    y = [5]
    plt.xlim(-0, 600)
    plt.ylim(-0, 600)
    plt.grid()
    # plt.plot(p1.x, p1.y, marker="o", markersize=p1.width, markeredgecolor="black", markerfacecolor="black")
    # plt.plot(p2.x, p2.y, marker="o", markersize=p2.width, markeredgecolor="black", markerfacecolor="black")
    for i in trejectory1:
        plt.plot(i[0], i[1], marker="o", markersize=i[2], markeredgecolor="black", markerfacecolor="black")
    for i in trejectory2:
        plt.plot(i[0], i[1], marker="o", markersize=i[2], markeredgecolor="black", markerfacecolor="black")
    for i in trejectory3:
        plt.plot(i[0], i[1], marker="o", markersize=i[2], markeredgecolor="black", markerfacecolor="black")
    for i in trejectory4:
        plt.plot(i[0], i[1], marker="o", markersize=i[2], markeredgecolor="black", markerfacecolor="black")
    for i in trejectory5:
        plt.plot(i[0], i[1], marker="o", markersize=i[2], markeredgecolor="black", markerfacecolor="black")
    plt.show()
