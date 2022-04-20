import numpy as np
import math
from matplotlib import pyplot as plt


# points
class Point:
    def __init__(self, x, y, z=0):
        self.x = x
        self.y = y
        self.z = z
        self.cord = (x, y)
        self.width = z

# strokes
class Strokes:
    def __int__(self, points, parameters):
        self.points = points
        self.parameters = parameters

    def get_trejector(self):
        for i in range(len(self.points) - 1):
            trejectory = 0
        return trejectory



if __name__ == "__main__":
    def get_arc_cord(p_oirgin, p_end, t, mu=0, sigma=0, delta=0.001):
        norm = get_norm(p_oirgin, p_end)
        theta = get_theta(p_oirgin, p_end)
        alpha = np.arcsin((p_end.z - p_oirgin.z) / norm)
        D = norm * delta / 2 / np.sin(delta/2)
        t0 = 0
        w = 1 / 2 * (1 + math.erf((np.log(t - t0) - mu) / (sigma * np.sqrt(2))))
        # print(w)
        phi = theta - delta / 2 + delta * w
        # print(phi)
        dx, dy = get_displacement_cord(delta, D, theta, w, phi)
        dz = w * (p_end.width - p_oirgin.width)
        return p_oirgin.x + dx, p_oirgin.y + dy, p_oirgin.z + dz

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





    p1 = Point(100, 100, 0)
    p2 = Point(400, 450, 60)
    p3 = []
    for i in range(1000):
        x, y, z = get_arc_cord(p1, p2, i/100, mu=0.1, sigma=0.4, delta=1)
        p3.append(Point(x, y, z))

    plt.rcParams["figure.figsize"] = [8.00, 8.50]
    plt.rcParams["figure.autolayout"] = True
    x = [5]
    y = [5]
    plt.xlim(0, 600)
    plt.ylim(0, 600)
    plt.grid()
    plt.plot(p1.x, p1.y, marker="o", markersize=p1.width, markeredgecolor="black", markerfacecolor="black")
    plt.plot(p2.x, p2.y, marker="o", markersize=p2.width, markeredgecolor="black", markerfacecolor="black")
    for i in range(1000):
        plt.plot(p3[i].x, p3[i].y, marker="o", markersize=p3[i].width, markeredgecolor="black", markerfacecolor="black")
    plt.show()
