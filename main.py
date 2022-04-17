import numpy as np
import math
from matplotlib import pyplot as plt

if __name__ == "__main__":
    def get_arc_cord(p_oirgin, p_end, t, mu=0, sigma=0, delta=2):
        theta = get_theta(p_oirgin, p_end)
        D = 1
        t0 = 0
        w = 1 / 2 * (1 + math.erf((np.log(t - t0) - mu) / (sigma * np.sqrt(2))))
        # print(w)
        phi = theta - delta / 2 + delta * w
        # print(phi)
        dx, dy = get_displacement_cord(delta, D, theta, w, phi)
        return p_oirgin.x + dx, p_oirgin.y + dy

    def get_theta(p_oirgin, p_end):
        dx = p_end.x - p_oirgin.x
        dy = p_end.y - p_oirgin.y
        theta = np.arctan2(dy, dx)
        return theta

    def get_displacement_cord(delta, D, theta, w, phi):
        if delta != 0:
            displacement_x = D / delta * (np.sin(phi) - np.sin(theta - delta / 2))
            displacement_y = D / delta * (np.cos(phi) - np.cos(theta - delta / 2))
        else:
            displacement_x = w * np.cos(theta)
            displacement_y = w * np.sin(theta)
        return displacement_x, displacement_y

    # points
    class Point:
        def __init__(self, x, y, z=0):
            self.x = x
            self.y = y
            self.z = z
            self.cord = (x, y)
            self.width = z

    p1 = Point(100, 100, 3)
    p2 = Point(500, 200, 20)
    p3 = Point(500, 500, 3)
    x, y = get_arc_cord(p1, p2, 1.0, mu=0.1, sigma=0.4, delta=0.1)
    print(x, y)
    p3 = Point(x, y, 3)

    plt.rcParams["figure.figsize"] = [8.00, 8.50]
    plt.rcParams["figure.autolayout"] = True
    x = [5]
    y = [5]
    plt.xlim(0, 600)
    plt.ylim(0, 600)
    plt.grid()
    plt.plot(p1.x, p1.y, marker="o", markersize=p1.width, markeredgecolor="black", markerfacecolor="black")
    plt.plot(p2.x, p2.y, marker="o", markersize=p2.width, markeredgecolor="black", markerfacecolor="black")
    plt.plot(p3.x, p3.y, marker="o", markersize=p3.width, markeredgecolor="red", markerfacecolor="red")
    plt.show()

    print(get_theta(p1, p3)/np.pi*180)

