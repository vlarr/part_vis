import matplotlib.pyplot as plt
import sys
import math
import numpy as np


def read_xyz_file(path, first_line_offset=2):
    file = open(path, "r")
    slines = file.readlines()[first_line_offset:]
    return [[line.split('\t')[0] for line in slines],
            [[float(coord) for coord in line.split('\t')[1:]] for line in slines]]


def write_r_file(path, data):
    open(path, "w").writelines('\n'.join([str(rv) for rv in data]))


def r(c1, c2):
    return math.sqrt(math.pow(c1[0] - c2[0], 2) + math.pow(c1[1] - c2[1], 2) + math.pow(c1[2] - c2[2], 2))


def calc_center_of_mass(coords, masses):
    cm = [0, 0, 0]
    for i in range(len(coords)):
        cm[0] += coords[i][0] * masses[i]
        cm[1] += coords[i][1] * masses[i]
        cm[2] += coords[i][2] * masses[i]
    return [c / sum(masses) for c in cm]


# def calc_r_array(coords, num_r_point=0):
#     result = []
#     for i in range(len(coords)):
#         if i != num_r_point:
#             result.append(r(coords[i], coords[num_r_point]))
#     return sorted(result)


def calc_cm_array(coords, masses):
    result = []
    cm = calc_center_of_mass(coords, masses)
    for i in range(len(coords)):
        result.append(r(coords[i], cm))
    return sorted(result)


def gauss(x, sigma, mu):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi))


def plot_graph(points):
    offset = 0
    for parr in points:
        plt.plot(parr, [offset] * len(parr), 'bd')
        offset += 0.5

    axes = plt.gca()
    axes.set_ylim([-1, offset + 0.5])
    plt.show()


def calc_parr(parr, r_min, r_max, step, sigma):
    r_len = r_max - r_min
    npoints = int(r_len / step) + 1
    x_arr = np.linspace(r_min, r_max, npoints)
    y_arr = [0] * npoints

    for i in range(npoints):
        for pt in parr:
            y_arr[i] += gauss(x_arr[i], sigma, pt)

    return [x_arr, y_arr]

def find_r_min_max(points):
    g_min = 10000.0
    g_max = 0
    for parr in points:
        g_min = min(g_min, parr[0])
        g_max = max(g_max, parr[len(parr) - 1])
    return [g_min, g_max]

########################################################################################################################

masses_lib = {
    'Ir': 192.217,
    'Ag': 107.8682,
    'Au': 196.966
}

points = []

if __name__ == "__main__":
    for filename in sys.argv[1:]:
        read_result = read_xyz_file(filename)
        coords = read_result[1]  # координаты частиц [][x,y,z]
        masses = [masses_lib[mass_id] for mass_id in read_result[0]]  # массы частиц []
        rarr = calc_cm_array(coords, masses)  # массив расстояний от центра массы до частиц |r_cm - r_i|
        points.append(rarr)

# for filename in ["example1xyz"]:
#     read_result = read_xyz_file(filename)
#     coords = read_result[1]
#     masses = [masses_lib[mass_id] for mass_id in read_result[0]]
#     rarr = calc_cm_array(coords, masses)
#     points.append(rarr)


c_result = []

r_min_max = find_r_min_max(points)

for parr in points:
    c_result.append(calc_parr(parr, r_min_max[0], r_min_max[1], 0.01, 0.001))

for ln in c_result:
    plt.plot(ln[0], ln[1], '-')

plt.show()